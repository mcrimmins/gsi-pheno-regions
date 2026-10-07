#!/usr/bin/env Rscript
# Step 11 (Cut A, Block 2): daily timing and variability features per year from daily PRISM,
# then across-year median / IQR and an across-year variability raster. Decided Oct 7 2026
# (config features: block2): rain event >= 5 mm/day (2.5 and 10 mm as sensitivity); warm-season
# onset = first wet day from May 1 starting a 10-day total >= 20 mm, with no dry spell
# >= 20 days in the next 30 days. Dry day = ppt < 1 mm. Hard freeze tmin <= -2.2 C.
# Year label Y: water year Oct(Y-1)-Sep(Y) for events, dry spells and limiting days;
# calendar year Y for onset and the spring hard-freeze GDD.
# Per-year bands (FLT4S):
#   onset_doy          warm-season rain onset DOY (end_doy + 1 = no onset that year)
#   ev_warm_<t>        days with ppt >= t mm in Apr-Sep (t = 5, 2p5, 10)
#   ev_cool_<t>        days with ppt >= t mm in Oct-Mar
#   dsl_max_warm       longest Apr-Sep dry spell (days with ppt < 1 mm)
#   n_dry20_warm       number of Apr-Sep dry spells >= 20 days
#   gdd5_lhf           GDD (base 5, from Jan 1) accumulated by the last spring tmin <= -2.2 C
#                      (0 = no hard freeze): how much warmth a late freeze can undo
#   cold_frac          share of water-year days with daily mean < 5 C (cold-limited)
#   dry_frac           share of water-year days that are not cold but have trailing 30-day
#                      P < 0.5 x 30-day Hargreaves PET (moisture-limited)
# Output: features/block2/block2_<Y>.tif; summaries/block2_median.tif, block2_iqr.tif, and
# summaries/block2_var.tif (cv_p_ann, cv_p_warm, sd_t_amj from Block 1 per-year files;
# onset_frac = share of years with an onset). Resume-safe per year.
#
# Usage: Rscript scripts/11_block2_timing.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/11_block2_timing.R")

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
source(here::here("R", "block2.R"))

cfg <- load_cfg(profile)
log_file <- log_init("11_block2_timing", profile)
notify_init(cfg)
b2 <- cfg$features$block2
chunk_rows <- cfg$features$chunk_rows
grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
if (!file.exists(grid_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
years <- feature_years(cfg)
log_msg("=== 11_block2_timing | profile: ", profile, " | years ", min(years), "-", max(years),
        " | ", length(block2_names(b2)), " features | workers ", cfg$workers)
t0 <- Sys.time()

prism_file <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
jobs <- list()
for (Y in years) {
  out <- feature_year_file(cfg, "block2", Y)
  if (file.exists(out)) { log_msg(Y, ": exists, skip"); next }
  fl <- lapply(setNames(c("ppt", "tmin", "tmax"), c("ppt", "tmin", "tmax")), function(v)
    setNames(list(prism_file(v, Y - 1), prism_file(v, Y)), c(Y - 1, Y)))
  miss <- unlist(fl)[!file.exists(unlist(fl))]
  if (length(miss)) { log_warn(Y, ": missing PRISM files, skip: ", paste(basename(miss), collapse = ", ")); next }
  jobs[[length(jobs) + 1]] <- list(Y = Y, files = fl, out = out)
}
log_msg(length(jobs), " year(s) to compute")

if (length(jobs)) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  plan(multisession, workers = min(cfg$workers, length(jobs)))
  memfrac <- cfg$memfrac
  res <- future_map(jobs, function(j) block2_year(j$Y, j$files, grid_file, j$out, b2,
                                                  chunk_rows, memfrac),
                    .options = furrr_options(seed = NULL, packages = "terra"))
  plan(sequential)
  for (r in res) {
    if (isTRUE(r$ok)) log_msg(r$year, ": ", r$n_cells, " cells") else log_err(r$year, " failed: ", r$msg)
  }
}

done <- vapply(years, function(Y) file.exists(feature_year_file(cfg, "block2", Y)), TRUE)
if (all(done)) {
  b2_files <- vapply(years, function(Y) feature_year_file(cfg, "block2", Y), "")
  log_msg("across-year median and IQR over ", length(years), " years")
  summarise_years(b2_files, summary_file(cfg, "block2", "median"), summary_file(cfg, "block2", "iqr"),
                  chunk_rows)
  med <- terra::rast(summary_file(cfg, "block2", "median"))
  g <- terra::global(med, c("min", "mean", "max"), na.rm = TRUE)
  for (k in seq_len(nrow(g))) log_msg(sprintf("  median %-13s min %9.2f  mean %9.2f  max %9.2f",
                                              rownames(g)[k], g[k, 1], g[k, 2], g[k, 3]))
  b1_files <- vapply(years, function(Y) feature_year_file(cfg, "block1", Y), "")
  if (all(file.exists(b1_files))) {
    if (length(years) < 3) log_warn("only ", length(years), " years: across-year variability is not meaningful")
    block2_variability(b1_files, b2_files, summary_file(cfg, "block2", "var"), b2, chunk_rows)
    v <- terra::rast(summary_file(cfg, "block2", "var"))
    g <- terra::global(v, c("min", "mean", "max"), na.rm = TRUE)
    for (k in seq_len(nrow(g))) log_msg(sprintf("  var    %-13s min %9.3f  mean %9.3f  max %9.3f",
                                                rownames(g)[k], g[k, 1], g[k, 2], g[k, 3]))
  } else log_warn("block2_var skipped: Block 1 per-year files missing (run 10_block1_seasonal.R)")
} else log_warn("summaries skipped: ", sum(!done), " year(s) missing")

msg <- sprintf("%d years, %.1f min", sum(done), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("11_block2_timing finished", msg, priority = 3, tags = "white_check_mark")
