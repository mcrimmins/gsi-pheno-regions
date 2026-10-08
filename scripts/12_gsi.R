#!/usr/bin/env Rscript
# Step 12 (Cut A, Block 3 / GSI): daily Growing Season Index per year on the PRISM grid, for
# every variant in config `features: gsi` (decided Oct 8 2026): parameter sets nfdrs (NFDRS2016
# literature: VPD 900-4100 Pa, daylength 10-11 h, precipitation 0-10 mm / 28 d, 21-day
# smoothing, green-up 0.5) and fems (FEMS operational: VPD 1956-3882 Pa, 11-12 h,
# 10.16-20.32 mm / 28 d, 28-day smoothing, green-up 0.3); moisture none / precip (28-day sum)
# / kbdi (KBDI ramp 200-600, decreasing); plus no-VPD sensitivity runs. Tmin ramp -2 to 5 C for
# both. Model code copied from gsi-scout (R/gsi.R Part 1) and vectorized (Part 2);
# scripts/12b_gsi_check.R checks the grid code against the original at sample points.
#
# Per year Y the daily series runs from Jan 1 of Y-1 (spin-up for smoothing, the precipitation
# sum and KBDI, which starts saturated) to Mar 31 of Y+1 (or Dec 31 of Y for the last year).
# Per-year bands per variant (DOY relative to Jan 1 of Y; can be < 1 or > 365):
#   gsi_mean, gsi_max      mean and maximum of the smoothed GSI in year Y
#   peak_doy               day of the maximum in Y
#   sos20, eos50           start (20 % of the seasonal amplitude, before the peak) and end
#                          (50 %, after the peak) of the GSI season, searched within 183 days
#                          of the peak, as for the satellite curves (NA if amplitude < 0.05)
#   gu_doy, dorm_doy       first switch to green-up in Y (2-phase classifier, threshold = the
#                          set's green-up value, 3-day persistence) and the next switch back to
#                          dormant (NA if none)
#   green_days, n_pulses   days in green-up in Y; number of switches to green-up in Y
#   lim_tmin/vpd/photo/moist  share of the days in Y on which that sub-index is the smallest
#                          and below 1 (NA when the variant has no such sub-index)
# KBDI uses each cell's mean annual precipitation over the PRISM years available
# (static/kbdi_map.tif, from the Block 1 per-year p_ann).
# Outputs features/gsi_<variant>/gsi_<variant>_<Y>.tif; summaries/gsi_<variant>_{median,iqr}.tif.
# Resume-safe per year (a year is skipped when all variant files exist).
#
# Usage: Rscript scripts/12_gsi.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/12_gsi.R")

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
source(here::here("R", "gsi.R"))

cfg <- load_cfg(profile)
sw <- step_workers(cfg, "gsi")
log_file <- log_init("12_gsi", profile)
notify_init(cfg)
g <- cfg$features$gsi
variants <- gsi_variants(g)
grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
if (!file.exists(grid_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
years <- feature_years(cfg)
log_msg("=== 12_gsi | profile: ", profile, " | years ", min(years), "-", max(years), " | variants: ",
        paste(names(variants), collapse = ", "), " | workers ", sw$n)
t0 <- Sys.time()

# Mean annual precipitation for KBDI (mm), from the Block 1 per-year water-year totals.
map_file <- file.path(static_dir(cfg), "kbdi_map.tif")
if (!file.exists(map_file)) {
  b1 <- vapply(years, function(Y) feature_year_file(cfg, "block1", Y), "")
  b1 <- b1[file.exists(b1)]
  if (!length(b1)) stop("no Block 1 per-year files (run 10_block1_seasonal.R) for the KBDI mean annual precipitation")
  pa <- terra::rast(lapply(b1, function(f) terra::rast(f)[["p_ann"]]))
  m <- terra::app(pa, mean, na.rm = TRUE); names(m) <- "map_mm"
  terra::writeRaster(m, map_file, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  log_msg("KBDI mean annual precipitation from ", length(b1), " water years -> ", map_file)
}

prism_file <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
vars <- c("ppt", "tmin", "tmax", "vpdmax")
jobs <- list()
for (Y in years) {
  outs <- lapply(names(variants), function(nm) feature_year_file(cfg, paste0("gsi_", nm), Y))
  names(outs) <- names(variants)
  if (all(file.exists(unlist(outs)))) { log_msg(Y, ": exists, skip"); next }
  ys <- c(Y - 1, Y, if (file.exists(prism_file("ppt", Y + 1)) && Y + 1 <= max(cfg$years)) Y + 1)
  fl <- lapply(setNames(vars, vars), function(v) setNames(lapply(ys, function(y) prism_file(v, y)), ys))
  miss <- unlist(fl)[!file.exists(unlist(fl))]
  if (length(miss)) { log_warn(Y, ": missing PRISM files, skip: ", paste(basename(miss), collapse = ", ")); next }
  for (o in outs) dir.create(dirname(o), showWarnings = FALSE, recursive = TRUE)
  jobs[[length(jobs) + 1]] <- list(Y = Y, files = fl, outs = outs)
}
log_msg(length(jobs), " year(s) to compute")

if (length(jobs)) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  plan(multisession, workers = min(sw$n, length(jobs)))
  memfrac <- sw$memfrac; chunk_rows <- g$chunk_rows %||% 15
  res <- future_map(jobs, function(j) gsi_year(j$Y, j$files, grid_file, map_file, j$outs, g, variants,
                                               chunk_rows, memfrac),
                    .options = furrr_options(seed = NULL, packages = "terra"))
  plan(sequential)
  for (r in res) {
    if (isTRUE(r$ok)) {
      log_msg(r$year, ": ", r$n_cells, " cells, series to ", r$to)
      notify(sprintf("12_gsi %d done", r$year), sprintf("%d cells", r$n_cells), priority = 2)
    } else log_err(r$year, " failed: ", r$msg)
  }
}

for (nm in names(variants)) {
  blk <- paste0("gsi_", nm)
  f <- vapply(years, function(Y) feature_year_file(cfg, blk, Y), "")
  if (!all(file.exists(f))) { log_warn(nm, ": summaries skipped, ", sum(!file.exists(f)), " year(s) missing"); next }
  summarise_years(f, summary_file(cfg, blk, "median"), summary_file(cfg, blk, "iqr"), cfg$features$chunk_rows)
  med <- terra::rast(summary_file(cfg, blk, "median"))
  gl <- terra::global(med[[c("gsi_max", "peak_doy", "sos20", "eos50", "gu_doy", "green_days")]],
                      "mean", na.rm = TRUE)
  log_msg(sprintf("  %-18s median: %s", nm, paste(sprintf("%s %.2f", rownames(gl), gl[, 1]), collapse = ", ")))
}

msg <- sprintf("%d years x %d variants, %.1f min", length(years), length(variants),
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("12_gsi finished", msg, priority = 3, tags = "white_check_mark")
