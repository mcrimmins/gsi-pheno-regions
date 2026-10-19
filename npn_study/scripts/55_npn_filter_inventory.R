#!/usr/bin/env Rscript
# Step 55 (NPN grass curing study, Q1): what each filter removes, NEON versus other sites.
# Purpose: decide whether relaxing filters for the record, climatology and spatial questions
# (Q1, Q2, Q4) recovers enough non-NEON sites to be worth it. Changes nothing upstream.
#   funnel.csv     cumulative filters on plant-years (usable -> crops out -> min readings -> PRISM
#                  site and study years -> date present -> season mismatch -> series with >= min_years
#                  years); per date; plant-years, plants, sites, species, series; NEON / other /
#                  managed
#   readings.csv   min-readings sensitivity (1, 4, 6, 8, 10, 12, 16, 20): plant-years and sites kept,
#                  NEON versus other, per date; series with >= 1, 2, 3 years of the date
#   by_region.csv  sites and plant-years per region (EPA Level II code; none = unknown) and state: current
#                  analysis set (min readings, >= min_years) versus a relaxed Q2/Q4 set
#                  (relaxed_min_readings, >= 1 year)
#   density.csv    readings per plant-year and median reading gap, NEON versus other
# Site names come from data/npn_site_names.csv (local only, 54a) and are used only for the NEON /
# managed flags; no names are written. Runs anywhere (seconds).
# Inputs: npn_curing_metrics.rds (GSI_OBS_DIR or config study_dir), data/npn_species_pathway.csv,
# <run_dir>/npn_study/npn_site_daily.rds (50, optional), data/npn_site_ecoregions.csv (59).
# Output <run_dir>/npn_study/inventory/.
#
# Usage: Rscript npn_study/scripts/55_npn_filter_inventory.R [profile] [relaxed_min_readings, default 6]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "full")
relaxed_n <- if (length(args) >= 2) as.integer(args[2]) else 6L

npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("55_npn_filter_inventory", profile)
ns <- cfg$npn_study %||% list()
sdir <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(sdir, "inventory")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 55_npn_filter_inventory | profile: ", profile, " | relaxed min readings: ", relaxed_n)
min_n <- ns$min_readings %||% 12; min_years <- ns$min_years %||% 3
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
log_msg(sprintf("npn_curing_metrics: %d rows, columns: %s", nrow(cm), paste(names(cm), collapse = ", ")))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id)
cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"

# flags (names used only here, never written)
nf <- npn_path("data", "npn_site_names.csv")
managed_re <- "garden|campus|universit|college|school|park|yard|home|house|residen|arboret|cemeter|golf|lawn|church|office|museum|zoo|nursery|farm|ranch hq|backyard|front|street|apartment|condo|library|center"
if (file.exists(nf)) {
  sn <- utils::read.csv(nf, stringsAsFactors = FALSE)
  nm <- tolower(sn$site_name[match(cm$site_id, sn$site_id)])
  cm$neon <- grepl("\\.phenology\\.phe", nm)
  cm$managed <- !cm$neon & grepl(managed_re, nm)
} else {
  log_msg("no data/npn_site_names.csv: NEON / managed flags unavailable (run 54a_npn_site_names.R)")
  cm$neon <- NA; cm$managed <- NA
}
cm$group <- ifelse(is.na(cm$neon), "unknown", ifelse(cm$neon, "NEON", ifelse(cm$managed, "other_managed", "other")))

cm$region <- npn_site_regions(cm$site_id)                  # EPA Level II ecoregion (59); "none" = unknown
dfile <- file.path(sdir, "npn_site_daily.rds")
prism_sites <- if (file.exists(dfile)) readRDS(dfile)$sites$site_id else NULL
if (is.null(prism_sites)) log_msg("no npn_site_daily.rds: PRISM-site filter skipped")

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
targets <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
groups <- c("NEON", "other", "other_managed", "unknown")

count <- function(x, step, tg) {
  base <- data.frame(target = tg, step = step, group = c("all", groups))
  cbind(base, do.call(rbind, lapply(base$group, function(g) {
    y <- if (g == "all") x else x[x$group == g, ]
    data.frame(plant_years = nrow(y), plants = length(unique(y$individual_id)), sites = length(unique(y$site_id)),
               species = length(unique(y$species_id)), series = length(unique(paste(y$site_id, y$species_id))))
  })))
}
# series with >= k years of the date (series-year = site x species x year)
keep_years <- function(x, k) {
  ny <- tapply(x$year, paste(x$site_id, x$species_id), function(v) length(unique(v)))
  x[as.numeric(ny[paste(x$site_id, x$species_id)]) >= k, ]
}
drop_mismatch <- function(x, col) {
  if (!nrow(x)) return(x)
  s <- paste(x$site_id, x$species_id); med <- c(tapply(x[[col]], s, cmed))
  x[abs(cdiff(x[[col]], med[s])) <= 90, ]
}

# 1. funnel
F <- list()
x0 <- cm;                                  F[[1]] <- count(x0, "0 all plant-years", "any")
x1 <- x0[x0$usable %in% TRUE, ];           F[[2]] <- count(x1, "1 usable", "any")
x2 <- x1[!x1$crop, ];                      F[[3]] <- count(x2, "2 crops out", "any")
x3 <- x2[x2$n_obs >= min_n, ];             F[[4]] <- count(x3, sprintf("3 >= %d readings", min_n), "any")
x4 <- x3[x3$year %in% yrs & (is.null(prism_sites) | x3$site_id %in% prism_sites), ]
F[[5]] <- count(x4, "4 PRISM site and study years", "any")
for (tg in names(targets)) {
  col <- targets[[tg]]
  x5 <- x4[!is.na(x4[[col]]), ];           F[[length(F) + 1]] <- count(x5, "5 date present", tg)
  x6 <- drop_mismatch(x5, col);            F[[length(F) + 1]] <- count(x6, "6 season mismatch out", tg)
  x7 <- keep_years(x6, min_years);         F[[length(F) + 1]] <- count(x7, sprintf("7 series >= %d years", min_years), tg)
}
FN <- do.call(rbind, F)
utils::write.csv(FN, file.path(odir, "funnel.csv"), row.names = FALSE)
for (i in seq_len(nrow(FN))) if (FN$group[i] %in% c("all", "NEON", "other", "other_managed"))
  if (FN$group[i] == "all") log_msg(sprintf("%-8s %-30s all %5d py %4d sites | NEON %5d py %3d sites | other %4d py %3d sites | managed %4d py %3d sites",
    FN$target[i], FN$step[i], FN$plant_years[i], FN$sites[i],
    FN$plant_years[i + 1], FN$sites[i + 1], FN$plant_years[i + 2], FN$sites[i + 2], FN$plant_years[i + 3], FN$sites[i + 3]))

# 2. min-readings sensitivity (after crops out, PRISM sites and years; season mismatches out)
base <- x2[x2$year %in% yrs & (is.null(prism_sites) | x2$site_id %in% prism_sites), ]
R <- list()
for (tg in names(targets)) for (k in c(1, 4, 6, 8, 10, 12, 16, 20)) {
  col <- targets[[tg]]
  x <- drop_mismatch(base[base$n_obs >= k & !is.na(base[[col]]), ], col)
  for (g in c("all", "NEON", "other")) {
    y <- if (g == "all") x else if (g == "other") x[x$group %in% c("other", "other_managed"), ] else x[x$group == g, ]
    R[[length(R) + 1]] <- data.frame(target = tg, min_readings = k, group = g, plant_years = nrow(y), sites = length(unique(y$site_id)),
      series_1yr = length(unique(paste(y$site_id, y$species_id))),
      series_2yr = length(unique(paste(keep_years(y, 2)$site_id, keep_years(y, 2)$species_id))),
      series_3yr = length(unique(paste(keep_years(y, 3)$site_id, keep_years(y, 3)$species_id))),
      sites_3yr = length(unique(keep_years(y, 3)$site_id)))
  }
}
RD <- do.call(rbind, R)
utils::write.csv(RD, file.path(odir, "readings.csv"), row.names = FALSE)
for (tg in names(targets)) {
  z <- RD[RD$target == tg & RD$group == "other", ]
  log_msg(sprintf("%-8s non-NEON by min readings: %s", tg, paste(sprintf("%d: %d sites (%d py, %d series)", z$min_readings, z$sites, z$plant_years, z$series_1yr), collapse = " | ")))
}

# 3. by region and state: current analysis set versus relaxed Q2 / Q4 set
B <- list()
for (tg in names(targets)) {
  col <- targets[[tg]]
  cur <- keep_years(drop_mismatch(base[base$n_obs >= min_n & !is.na(base[[col]]), ], col), min_years)
  rel <- drop_mismatch(base[base$n_obs >= relaxed_n & !is.na(base[[col]]), ], col)
  for (by in c("region", "state")) {
    lv <- sort(unique(c(cur[[by]], rel[[by]])))
    for (l in lv) {
      a <- cur[cur[[by]] == l, ]; b <- rel[rel[[by]] == l, ]
      B[[length(B) + 1]] <- data.frame(target = tg, by = by, level = as.character(l),
        current_sites = length(unique(a$site_id)), current_neon_sites = length(unique(a$site_id[a$group == "NEON"])), current_py = nrow(a),
        relaxed_sites = length(unique(b$site_id)), relaxed_neon_sites = length(unique(b$site_id[b$group == "NEON"])), relaxed_py = nrow(b),
        relaxed_series = length(unique(paste(b$site_id, b$species_id))))
    }
  }
  bt <- do.call(rbind, B[vapply(B, function(d) d$target == tg && d$by == "region", TRUE)])
  log_msg(sprintf("%-8s sites current -> relaxed (>= %d readings, >= 1 yr): %d -> %d (non-NEON %d -> %d) | regions with >= 3 sites: %d -> %d",
    tg, relaxed_n, sum(bt$current_sites), sum(bt$relaxed_sites), sum(bt$current_sites - bt$current_neon_sites), sum(bt$relaxed_sites - bt$relaxed_neon_sites),
    sum(bt$current_sites >= 3 & bt$level != "none"), sum(bt$relaxed_sites >= 3 & bt$level != "none")))
}
utils::write.csv(do.call(rbind, B), file.path(odir, "by_region.csv"), row.names = FALSE)

# 4. reading density, NEON versus other (usable, crops out, study years)
gapc <- intersect(c("g50_down_gap", "leaves_off_gap"), names(base))
DN <- do.call(rbind, lapply(c("NEON", "other", "other_managed"), function(g) {
  y <- base[base$group == g, ]; if (!nrow(y)) return(NULL)
  q <- stats::quantile(y$n_obs, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
  o <- data.frame(group = g, plant_years = nrow(y), n_obs_p10 = q[1], n_obs_p25 = q[2], n_obs_median = q[3], n_obs_p75 = q[4], n_obs_p90 = q[5],
                  share_below_min = mean(y$n_obs < min_n))
  for (gc in gapc) o[[paste0(gc, "_median")]] <- stats::median(y[[gc]], na.rm = TRUE)
  o
}))
rownames(DN) <- NULL
utils::write.csv(DN, file.path(odir, "density.csv"), row.names = FALSE)
for (i in seq_len(nrow(DN))) log_msg(sprintf("readings %-13s %5d py | median %.0f (IQR %.0f-%.0f) | below %d: %.0f %%",
  DN$group[i], DN$plant_years[i], DN$n_obs_median[i], DN$n_obs_p25[i], DN$n_obs_p75[i], min_n, 100 * DN$share_below_min[i]))

log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
