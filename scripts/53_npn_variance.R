#!/usr/bin/env Rscript
# Step 53 (NPN grass curing study, Q3 place versus year): variance partition of plant dates.
# Unit: plant-year (one tagged plant, one year; >= npn_study: min_readings readings; crops out).
# Each date is unwrapped around its series' (site x species) circular median; |departure| > 90 d
# is a season mismatch and dropped. Mixed model per date (green-up = 50 % up, curing = 50 % down,
# cured), all grasses and C3 / C4 separately:
#   date ~ 1 + (1|site) + (1|species) + (1|site:species)          place and species
#            + (1|year) + (1|region:year) + (1|site:year)          year signal shared by species
#            + (1|site:species:year)                                year signal of one species
#            + residual                                             plant (and observer) within series-year
# Place = site + species + site:species. Year-to-year = the rest. "Shared weather" = year +
# region:year + site:year: the part every grass at a site moves together, the most any site
# weather record could explain. Also: expected reading-interval error (gap^2 / 12 for a crossing
# interpolated inside a gap of that many days) against the plant residual, and the reliability
# of a site-year mean.
# Inputs: npn_curing_metrics.rds (GSI_OBS_DIR or config study_dir), data/npn_species_pathway.csv,
# region from <run_dir>/eval/target_check/site_years.csv. Needs lme4. Runs anywhere (< 1 min).
# Output <run_dir>/npn_study/variance/components.csv, summary.csv.
#
# Usage: Rscript scripts/53_npn_variance.R [profile] [sensitivity]
#   sensitivity: a name under config npn_study: sensitivity (e.g. no1617); outputs then go to
#   <run_dir>/npn_study/sens/<name>/variance/

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
source(here::here("R", "npn.R"))

cfg <- load_cfg(profile)
log_file <- log_init("53_npn_variance", profile)
ns <- cfg$npn_study %||% list()
sens <- npn_sensitivity(ns, args)
odir <- file.path(cfg$run_dir, "npn_study", sens$sub, "variance"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 53_npn_variance | profile: ", profile, if (nzchar(sens$name)) sprintf(" | sensitivity %s: excluding %s", sens$name, paste(sens$exclude_years, collapse = ", ")) else "")
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs the lme4 package")

ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
cm <- readRDS(file.path(ddir, "npn_curing_metrics.rds"))
sp <- utils::read.csv(here::here("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
u <- cm[cm$usable %in% TRUE & cm$n_obs >= (ns$min_readings %||% 12) & !(cm$year %in% sens$exclude_years), ]
m <- match(u$species_id, sp$species_id)
u$pathway <- sp$pathway[m]; u <- u[!(sp$crop[m] %in% "yes"), ]
tc <- file.path(cfg$run_dir, "eval", "target_check", "site_years.csv")
u$region <- NA_integer_
if (file.exists(tc)) {
  r <- utils::read.csv(tc, stringsAsFactors = FALSE)
  r$site_id <- as.integer(sub("^npn(\\d+)\\|.*$", "\\1", r$site)); r <- r[!is.na(r$site_id) & !is.na(r$region), ]
  u$region <- r$region[match(u$site_id, r$site_id)]
}
u$region[is.na(u$region)] <- 0                              # unknown region: one shared level
log_msg(sprintf("plant-years: %d (%d plants, %d sites, %d species); region known for %.0f %%",
                nrow(u), length(unique(u$individual_id)), length(unique(u$site_id)), length(unique(u$species_id)), 100 * mean(u$region > 0)))

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
targets <- list(greenup = list(col = "g50_up_doy", gap = NA), curing = list(col = "g50_down_doy", gap = "g50_down_gap"),
                cured = list(col = "cured_doy", gap = "leaves_off_gap"))
groups <- list(all = c("C3", "C4", NA), C3 = "C3", C4 = "C4")
fml <- y ~ 1 + (1 | site) + (1 | species) + (1 | site:species) + (1 | year) + (1 | region:year) + (1 | site:year) + (1 | site:species:year)
place <- c("site", "species", "site:species"); shared <- c("year", "region:year", "site:year")
comp_rows <- list(); sum_rows <- list()
for (tg in names(targets)) {
  x <- u[!is.na(u[[targets[[tg]]$col]]), ]
  x$series <- paste(x$site_id, x$species_id)
  med <- c(tapply(x[[targets[[tg]]$col]], x$series, cmed))
  d <- cdiff(x[[targets[[tg]]$col]], med[x$series])
  x$y <- as.numeric(med[x$series]) + d
  mis <- abs(d) > 90; mis_pw <- x$pathway[mis]; x <- x[!mis, ]
  for (gn in names(groups)) {
    g <- x[x$pathway %in% groups[[gn]], ]
    if (nrow(g) < 100) next
    g$site <- factor(g$site_id); g$species <- factor(g$species_id); g$year <- factor(g$year); g$region <- factor(g$region)
    fit <- lme4::lmer(fml, data = g, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
    vc <- as.data.frame(lme4::VarCorr(fit))
    v <- stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "plant", vc$grp))
    tot <- sum(v); within <- sum(v[setdiff(names(v), place)])
    comp_rows[[length(comp_rows) + 1]] <- data.frame(target = tg, group = gn, component = names(v), var = v, sd = sqrt(v),
      share_total = v / tot, share_year_to_year = ifelse(names(v) %in% place, NA, v / within))
    # reading-interval error: crossing interpolated inside a gap -> uniform error, var = gap^2 / 12
    gv <- if (!is.na(targets[[tg]]$gap)) { gp <- g[[targets[[tg]]$gap]]; mean(gp^2 / 12, na.rm = TRUE) } else NA_real_
    # plants per series-year, and reliability of a series-year mean of the shared signal
    npl <- stats::median(table(paste(g$series, g$year)))
    sh <- sum(v[shared]); spec <- v[["site:species:year"]]; res <- v[["plant"]]
    sum_rows[[length(sum_rows) + 1]] <- data.frame(target = tg, group = gn, plant_years = nrow(g), series = length(unique(g$series)),
      sites = length(unique(g$site)), mismatches_dropped = sum(mis_pw %in% groups[[gn]]),
      sd_total = sqrt(tot), sd_year_to_year = sqrt(within), place_share = sum(v[place]) / tot,
      shared_share = sh / within, species_year_share = spec / within, plant_share = res / within,
      region_year_share = (v[["year"]] + v[["region:year"]]) / within,
      reading_gap_var = gv, reading_gap_share_of_plant = gv / res, plants_per_series_year = npl,
      shared_of_series_year_mean = sh / (sh + spec + res / npl))
    with(sum_rows[[length(sum_rows)]], log_msg(sprintf(
      "%-7s %-3s: %d plant-years, %d series | SD total %.0f d, year-to-year %.0f d | place %.0f %% of total | of year-to-year: shared %.0f %% (region-wide %.0f %%), species %.0f %%, plant %.0f %% | reading gaps %s of plant | shared in a series-year mean %.0f %%",
      target, group, plant_years, series, sd_total, sd_year_to_year, 100 * place_share, 100 * shared_share, 100 * region_year_share,
      100 * species_year_share, 100 * plant_share, ifelse(is.na(reading_gap_share_of_plant), "-", sprintf("%.0f %%", 100 * reading_gap_share_of_plant)),
      100 * shared_of_series_year_mean)))
  }
}
C <- do.call(rbind, comp_rows); rownames(C) <- NULL
utils::write.csv(C, file.path(odir, "components.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, sum_rows), file.path(odir, "summary.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
