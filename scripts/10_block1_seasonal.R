#!/usr/bin/env Rscript
# Step 10 (Cut A, Block 1): seasonal climate features per year from daily PRISM, then
# across-year median and IQR. Decided Oct 6 2026 (config features: block1):
#   year Y = water year Oct(Y-1)-Sep(Y) for temperature, precipitation and aridity;
#   calendar year Y for freeze dates and GDD. Seasons OND / JFM / AMJ / JAS.
# Per-year bands (FLT4S):
#   t_ond t_jfm t_amj t_jas  mean daily (tmin+tmax)/2 by season (deg C)
#   t_ann, t_range           water-year mean; warmest minus coolest season
#   p_ann                    water-year precipitation (mm)
#   pf_ond ... pf_jas        fraction of p_ann by season
#   ai_ann, ai_warm          P / PET (Hargreaves), water year and AMJ+JAS
#   lsf_<t> fff_<t> ffs_<t>  last spring freeze DOY (0 = none), first fall freeze DOY
#                            (days-in-year + 1 = none), freeze-free days; t = 0 and m2p2
#                            (tmin <= 0 and <= -2.2 deg C)
#   gdd5, gdd10              calendar-year growing degree days
# Output: features/block1/block1_<Y>.tif; summaries/block1_median.tif, block1_iqr.tif
# Features are computed for every land cell (wall-to-wall). Resume-safe per year.
#
# Usage: Rscript scripts/10_block1_seasonal.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/10_block1_seasonal.R")

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
source(here::here("R", "block1.R"))

cfg <- load_cfg(profile)
sw <- step_workers(cfg, "block1")
log_file <- log_init("10_block1_seasonal", profile)
notify_init(cfg)
b1 <- cfg$features$block1
chunk_rows <- cfg$features$chunk_rows
grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
if (!file.exists(grid_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
years <- feature_years(cfg)
log_msg("=== 10_block1_seasonal | profile: ", profile, " | years ", min(years), "-", max(years),
        " | ", length(block1_names(b1)), " features | workers ", sw$n)
t0 <- Sys.time()

prism_file <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
jobs <- list()
for (Y in years) {
  out <- feature_year_file(cfg, "block1", Y)
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
  plan(multisession, workers = min(sw$n, length(jobs)))
  memfrac <- sw$memfrac
  res <- future_map(jobs, function(j) block1_year(j$Y, j$files, grid_file, j$out, b1,
                                                  chunk_rows, memfrac),
                    .options = furrr_options(seed = NULL, packages = "terra"))
  plan(sequential)
  for (r in res) {
    if (isTRUE(r$ok)) log_msg(r$year, ": ", r$n_cells, " cells") else log_err(r$year, " failed: ", r$msg)
  }
}

done <- vapply(years, function(Y) file.exists(feature_year_file(cfg, "block1", Y)), TRUE)
if (all(done)) {
  log_msg("across-year median and IQR over ", length(years), " years")
  summarise_years(vapply(years, function(Y) feature_year_file(cfg, "block1", Y), ""),
                  summary_file(cfg, "block1", "median"), summary_file(cfg, "block1", "iqr"),
                  chunk_rows)
  med <- terra::rast(summary_file(cfg, "block1", "median"))
  g <- terra::global(med, c("min", "mean", "max"), na.rm = TRUE)
  for (k in seq_len(nrow(g))) log_msg(sprintf("  median %-8s min %9.2f  mean %9.2f  max %9.2f",
                                              rownames(g)[k], g[k, 1], g[k, 2], g[k, 3]))
} else log_warn("summaries skipped: ", sum(!done), " year(s) missing")

msg <- sprintf("%d years, %.1f min", sum(done), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("10_block1_seasonal finished", msg, priority = 3, tags = "white_check_mark")
