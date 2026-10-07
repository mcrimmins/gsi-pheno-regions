#!/usr/bin/env Rscript
# Step 06 (Cut B): seasonal curves and phenology metrics from the 4 km MOD13A2 + MYD13A2
# stacks (04_cutb_download.R). Config `cutb_curves:` (defaults agreed Oct 7 2026).
# Pass 1 (per year): Terra - Aqua offset per cell (Terra interpolated to Aqua composite
#   midpoints; median over valid pairs), NDVI and NDII7 -> cutb/offsets/offset_<Y>.tif; then
#   median across years -> cutb/offsets/offset_median.tif (IQR shows drift).
# Pass 2 (per year Y): composites from Oct (Y-1) to Mar (Y+1), duplicates at year boundaries
#   removed, composites with < min_valid valid pixels or > max_snow_frac snow dropped, offset
#   split between the sensors (Terra - off/2, Aqua + off/2), observation dates from the
#   composite day-of-year layer (midpoint when that is unusable), weighted Whittaker smoothing
#   on a 4-day grid with one upper-envelope pass. Outputs:
#     cutb/curves/ndvi_<Y>.tif, ndii_<Y>.tif   smoothed curves at DOY 1, 9, ..., 361
#     cutb/metrics/metrics_<Y>.tif             per-year metrics (DOY relative to Jan 1 of Y;
#                                              < 1 or > 365 = previous / next year):
#       ndvi_max, ndvi_base, ndvi_amp, ndvi_peak_doy; sos20, sos50 (green-up: last crossing of
#       base + 20 / 50 % of amplitude before the peak); eos50, eos20 (after the peak); gsl20;
#       n_green (green spells >= 16 days above 50 %); ndii_max, ndii_base, ndii_amp,
#       ndii_peak_doy; cure50 (NDII7 falls below 50 % of its range after its peak);
#       cure_days (75 % -> 25 % decline); n_obs (observations in year Y); max_gap (days).
# Summaries: summaries/cutb_metrics_{median,iqr}.tif, cutb_ndvi_curve_{median,iqr}.tif,
#   cutb_ndii_curve_{median,iqr}.tif. Resume-safe per year and pass.
#
# Usage: Rscript scripts/06_cutb_curves.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/06_cutb_curves.R")

suppressPackageStartupMessages({ library(future); library(furrr) })

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cutb.R"))
source(here::here("R", "cutb_curves.R"))

cfg <- load_cfg(profile)
log_file <- log_init("06_cutb_curves", profile)
notify_init(cfg)
cc <- cfg$cutb_curves
chunk_rows <- cfg$features$chunk_rows
grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
prods <- unlist(cfg$cutb$products)
have <- function(p, y) all(file.exists(vapply(cc_vars, function(v) cutb_out_file(cfg, p, v, y), "")))
years <- cfg$years[vapply(cfg$years, function(y) have(prods[1], y), TRUE)]
if (!length(years)) stop("no Cut B years found; run 04_cutb_download.R first")
log_msg("=== 06_cutb_curves | profile: ", profile, " | years ", min(years), "-", max(years),
        " | lambda ", cc$lambda, " | workers ", cfg$workers)
t0 <- Sys.time()
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
memfrac <- cfg$memfrac
run_jobs <- function(jobs, fun) {
  if (!length(jobs)) return(list())
  plan(multisession, workers = min(cfg$workers, length(jobs)))
  on.exit(plan(sequential))
  future_map(jobs, fun, .options = furrr_options(seed = NULL, packages = "terra",
             globals = c("cfg", "cc", "grid_file", "chunk_rows", "memfrac", "cc_vars",
                         ls(envir = globalenv(), pattern = "^(cc_|cutb_|whit_|read_rows|first_true|last_true|write_feature_matrix)"))))
}

# ---- Pass 1: offsets -----------------------------------------------------------------------
off_years <- years[vapply(years, function(y) have(prods[2], y), TRUE)]
jobs <- lapply(off_years, function(Y) list(Y = Y, out = cc_year_file(cfg, "offsets", "offset", Y)))
jobs <- Filter(function(j) !file.exists(j$out), jobs)
log_msg("pass 1 (offsets): ", length(jobs), " of ", length(off_years), " year(s) to compute")
for (r in run_jobs(jobs, function(j) cutb_offset_year(j$Y, cfg, grid_file, j$out, cc, chunk_rows, memfrac))) {
  if (isTRUE(r$ok)) log_msg("  offset ", r$year, ": ", r$n_cells, " cells") else log_err("offset ", r$year, " failed: ", r$msg)
}
off_files <- vapply(off_years, function(Y) cc_year_file(cfg, "offsets", "offset", Y), "")
off_med <- file.path(cc_dir(cfg, "offsets"), "offset_median.tif")
if (length(off_files) && all(file.exists(off_files))) {
  summarise_years(off_files, off_med, file.path(cc_dir(cfg, "offsets"), "offset_iqr.tif"), chunk_rows)
  g <- terra::global(terra::rast(off_med)[[1:2]], function(v) stats::quantile(v, c(.05, .5, .95), na.rm = TRUE))
  log_msg(sprintf("  Terra - Aqua median offset, 5/50/95 %%: NDVI %.4f / %.4f / %.4f; NDII7 %.4f / %.4f / %.4f",
                  g[1, 1], g[1, 2], g[1, 3], g[2, 1], g[2, 2], g[2, 3]))
} else log_warn("no complete offsets: curves are merged without offset correction")

# ---- Pass 2: curves and metrics ------------------------------------------------------------
outs <- function(Y) list(ndvi = cc_year_file(cfg, "curves", "ndvi", Y), ndii = cc_year_file(cfg, "curves", "ndii", Y),
                         met = cc_year_file(cfg, "metrics", "metrics", Y))
jobs <- lapply(years, function(Y) list(Y = Y, out = outs(Y)))
jobs <- Filter(function(j) !all(file.exists(unlist(j$out))), jobs)
log_msg("pass 2 (curves, metrics): ", length(jobs), " of ", length(years), " year(s) to compute")
of <- if (file.exists(off_med)) off_med else NULL
jobs <- lapply(jobs, function(j) c(j, list(off = of)))
for (r in run_jobs(jobs, function(j) cutb_curves_year(j$Y, cfg, grid_file, j$off, j$out, cc, chunk_rows, memfrac))) {
  if (isTRUE(r$ok)) log_msg("  curves ", r$year, ": ", r$n_cells, " cells") else log_err("curves ", r$year, " failed: ", r$msg)
}

# ---- Summaries ---------------------------------------------------------------------------
done <- vapply(years, function(Y) all(file.exists(unlist(outs(Y)))), TRUE)
if (all(done)) {
  for (nm in c("met", "ndvi", "ndii")) {
    lab <- c(met = "cutb_metrics", ndvi = "cutb_ndvi_curve", ndii = "cutb_ndii_curve")[[nm]]
    summarise_years(vapply(years, function(Y) outs(Y)[[nm]], ""), summary_file(cfg, lab, "median"),
                    summary_file(cfg, lab, "iqr"), chunk_rows)
  }
  med <- terra::rast(summary_file(cfg, "cutb_metrics", "median"))
  g <- terra::global(med, c("min", "mean", "max"), na.rm = TRUE)
  for (k in seq_len(nrow(g))) log_msg(sprintf("  median %-14s min %9.3f  mean %9.3f  max %9.3f",
                                              rownames(g)[k], g[k, 1], g[k, 2], g[k, 3]))
} else log_warn("summaries skipped: ", sum(!done), " year(s) missing")

msg <- sprintf("%d years, %.1f min", sum(done), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("06_cutb_curves finished", msg, priority = 3, tags = "white_check_mark")
