#!/usr/bin/env Rscript
# Step 51 (NPN grass curing study): responses and climate metrics.
#   series_years  one row per site x species x year (median over that site-year's plants; plant-
#                 years with >= npn_study: min_readings readings): green-up (50 % up), curing
#                 (50 % down), cured date; C3/C4 etc. from data/npn_species_pathway.csv; working
#                 region (A2 k13 map number) where known
#   site_years    per site and year: thermal (GDD, seasonal means, chill days, freezes), moisture
#                 (seasonal rain, VPD, events, dry spells, rain onset), P - PET memory on fixed dates
#   normals       per site: medians of the site-year metrics over npn_study: normals years, aridity,
#                 warm-season share, day-length gate dates
#   anchor        curing drivers in windows that START at the series' typical green-up (never end
#                 on the observed curing date): rain in 30 / 60 / 90 days, dry spell, events
# Inputs: <run_dir>/npn_study/npn_site_daily.rds (50), npn_curing_metrics.rds (GSI_OBS_DIR or
# config study_dir), data/npn_species_pathway.csv; region from <run_dir>/eval/target_check/
# site_years.csv when present (36). Runs anywhere (small).
# Output <run_dir>/npn_study/npn_metrics.rds = list(series_years, site_years, normals, anchor).
#
# Usage: GSI_OBS_DIR=<folder> Rscript scripts/51_npn_metrics.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "features.R"))
source(here::here("R", "block2.R"))
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
source(here::here("R", "npn.R"))

cfg <- load_cfg(profile)
log_file <- log_init("51_npn_metrics", profile)
ns <- cfg$npn_study %||% list()
out_dir <- file.path(cfg$run_dir, "npn_study")
t0 <- Sys.time()
log_msg("=== 51_npn_metrics | profile: ", profile)

D <- readRDS(file.path(out_dir, "npn_site_daily.rds"))
ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
cm <- readRDS(file.path(ddir, "npn_curing_metrics.rds"))
sp <- utils::read.csv(here::here("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)

# responses
SY <- npn_series_years(cm, sp, ns$min_readings %||% 12)
SY <- SY[SY$site_id %in% D$sites$site_id & SY$year %in% as.integer(format(D$dates, "%Y")), ]
tc <- file.path(cfg$run_dir, "eval", "target_check", "site_years.csv")
if (file.exists(tc)) {
  r <- utils::read.csv(tc, stringsAsFactors = FALSE)
  r$site_id <- as.integer(sub("^npn(\\d+)\\|.*$", "\\1", r$site)); r <- r[!is.na(r$site_id) & !is.na(r$region), ]
  SY$region <- r$region[match(SY$site_id, r$site_id)]
} else SY$region <- NA_integer_
log_msg(sprintf("series-years: %d (%d series, %d sites, %d species) | green-up %d, curing %d, cured %d | region known for %.0f %%",
                nrow(SY), length(unique(SY$series)), length(unique(SY$site_id)), length(unique(SY$species_id)),
                sum(!is.na(SY$greenup)), sum(!is.na(SY$curing)), sum(!is.na(SY$cured)), 100 * mean(!is.na(SY$region))))

# site-year metrics and normals
YM <- npn_site_year_metrics(D, cfg$features$block2, unlist(ns$memory_days %||% list(15, 30, 60)))
nyr <- ns$normals %||% list(start = 2001, end = 2025)
N <- npn_site_normals(YM, D, seq(nyr$start, nyr$end))
log_msg(sprintf("site-year metrics: %d sites x %d years, %d metrics; normals over %d-%d",
                length(unique(YM$site_id)), length(unique(YM$year)), ncol(YM) - 2, nyr$start, nyr$end))

# green-up-anchored curing drivers (typical green-up of the series, else of the site)
gu_series <- tapply(SY$greenup, SY$series, function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE))
gu_site <- tapply(SY$greenup, SY$site_id, function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE))
SY$gu_typical <- as.numeric(gu_series[SY$series])
miss <- is.na(SY$gu_typical); SY$gu_typical[miss] <- as.numeric(gu_site[as.character(SY$site_id[miss])])
A <- npn_anchor_metrics(D, match(SY$site_id, D$sites$site_id), SY$year, SY$gu_typical)
A <- cbind(SY[, c("site_id", "species_id", "year")], A)
log_msg(sprintf("anchor metrics: typical green-up known for %.0f %% of series-years", 100 * mean(!is.na(SY$gu_typical))))

saveRDS(list(series_years = SY, site_years = YM, normals = N, anchor = A), file.path(out_dir, "npn_metrics.rds"))
utils::write.csv(SY, file.path(out_dir, "npn_series_years.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(out_dir, "npn_metrics.rds"), log_file))
