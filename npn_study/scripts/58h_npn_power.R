#!/usr/bin/env Rscript
# Step 58h (NPN grass curing study; laptop): how much ground monitoring would it take to see a
# weather signal in the year-to-year departures? Built from 53's variance components (all grasses):
# a plant's departure in a year = regional year signal (year + region:year) + site-year (site:year)
# + species at the site (site:species:year) + plant (residual).
# For a design of n_sites sites in a region, n_species species per site and n_plants plants per
# species, the regional mean departure in a year measures the regional year signal with
#   noise variance = v_site:year / n_sites + v_site:species:year / (n_sites n_species)
#                    + v_plant / (n_sites n_species n_plants)
#   reliability    = v_signal / (v_signal + noise variance),   v_signal = v_year + v_region:year
# If a weather driver explains a share f of the true regional year signal (true r = sqrt(f)), the
# observed correlation is r sqrt(reliability), and the years needed to detect it (two-sided 5 %,
# power 80 %, Fisher z) are n = ((1.96 + 0.84) / atanh(r_obs))^2 + 3.
# Also for one site (the site-year signal = v_signal + v_site:year).
# Output <run_dir>/npn_study/record/power.csv and report/tables/npn_power.csv (grid of designs) and
# npn_power_current.csv (the record's typical design: from 53, plants per series-year, and series
# per site-year).
#
# Usage: Rscript npn_study/scripts/58h_npn_power.R [profile]   (after 53)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("58h_npn_power", profile)
base <- file.path(cfg$run_dir, "npn_study"); rdir <- file.path(base, "record"); dir.create(rdir, showWarnings = FALSE)
log_msg("=== 58h_npn_power | profile: ", profile)
C <- utils::read.csv(file.path(base, "variance", "components.csv"), stringsAsFactors = FALSE)
Sm <- utils::read.csv(file.path(base, "variance", "summary.csv"), stringsAsFactors = FALSE)
C <- C[C$group == "all", ]
vv <- function(tg, k) { x <- C$var[C$target == tg & C$component == k]; if (length(x)) x else 0 }
years_needed <- function(r) ifelse(r <= 0, Inf, ((stats::qnorm(0.975) + stats::qnorm(0.8)) / atanh(pmin(r, 0.999)))^2 + 3)
power_at <- function(r, n) stats::pnorm(atanh(pmin(r, 0.999)) * sqrt(n - 3) - stats::qnorm(0.975))
n_rec <- 11                                                      # years in the NPN year-to-year record

grid <- expand.grid(n_sites = c(1, 2, 3, 5, 10, 20), n_species = c(1, 2, 3), n_plants = c(1, 3, 5, 10, 20), f = c(0.25, 0.5))
OUT <- list()
for (tg in unique(C$target)) {
  sig <- vv(tg, "year") + vv(tg, "region:year"); sy <- vv(tg, "site:year"); ssy <- vv(tg, "site:species:year"); pl <- vv(tg, "plant")
  g <- grid
  g$noise_var <- sy / g$n_sites + ssy / (g$n_sites * g$n_species) + pl / (g$n_sites * g$n_species * g$n_plants)
  g$reliability <- sig / (sig + g$noise_var)
  g$r_obs <- sqrt(g$f) * sqrt(g$reliability)
  g$years_80 <- ceiling(years_needed(g$r_obs)); g$power_11y <- power_at(g$r_obs, n_rec)
  g$se_days <- sqrt(g$noise_var)
  OUT[[tg]] <- cbind(target = tg, signal_sd_days = sqrt(sig), g)
  log_msg(sprintf("%-8s regional year signal SD %.1f d; site-year %.1f, species-site-year %.1f, plant %.1f d", tg, sqrt(sig), sqrt(sy), sqrt(ssy), sqrt(pl)))
  for (ns_ in c(1, 3, 10)) { z <- g[g$n_sites == ns_ & g$n_species == 2 & g$f == 0.5, ]
    log_msg(sprintf("    %2d site(s), 2 species: reliability %s by plants per species %s; years for 80 %% power (driver explains half): %s", ns_,
                    paste(sprintf("%.2f", z$reliability), collapse = "/"), paste(z$n_plants, collapse = "/"), paste(z$years_80, collapse = "/"))) }
}
P <- do.call(rbind, OUT)
utils::write.csv(P, file.path(rdir, "power.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE); utils::write.csv(P, file.path(tdir, "npn_power.csv"), row.names = FALSE)

# the record's typical design: plants per series-year (53), about 2 species per NEON site, regions with 1-3 sites
CU <- do.call(rbind, lapply(unique(C$target), function(tg) {
  npl <- Sm$plants_per_series_year[Sm$target == tg & Sm$group == "all"][1]
  do.call(rbind, lapply(c(1, 2, 3), function(ns_) {
    sig <- vv(tg, "year") + vv(tg, "region:year"); nv <- vv(tg, "site:year") / ns_ + vv(tg, "site:species:year") / (ns_ * 2) + vv(tg, "plant") / (ns_ * 2 * npl)
    rel <- sig / (sig + nv)
    data.frame(target = tg, n_sites = ns_, n_species = 2, n_plants = npl, reliability = rel,
               years_80_half = ceiling(years_needed(sqrt(0.5 * rel))), years_80_quarter = ceiling(years_needed(sqrt(0.25 * rel))),
               power_11y_half = power_at(sqrt(0.5 * rel), n_rec), power_11y_quarter = power_at(sqrt(0.25 * rel), n_rec))
  }))
}))
utils::write.csv(CU, file.path(tdir, "npn_power_current.csv"), row.names = FALSE)
for (i in seq_len(nrow(CU))) with(CU[i, ], log_msg(sprintf("  record-like design %-8s %d site(s) x 2 species x %d plants: reliability %.2f; years for 80 %% power: %d (driver explains half), %d (a quarter); power with 11 years %.0f %% / %.0f %%",
                                                           target, n_sites, n_plants, reliability, years_80_half, years_80_quarter, 100 * power_11y_half, 100 * power_11y_quarter)))
log_msg("=== done | ", file.path(rdir, "power.csv"), " | log: ", log_file)
