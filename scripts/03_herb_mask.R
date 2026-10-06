#!/usr/bin/env Rscript
# Step 03: herbaceous cover layers from Annual NLCD (decided Oct 4: grass + shrub/scrub,
# years 2001/2012/2024 averaged; strict mask herb >= 50 % of land and crops < 25 %; all in
# config.yml). Oct 5: these describe/weight cells; the analysis domain is all land cells.
#
# 1. Download NLCD land cover for each configured year to <data_root>/_shared/raw/nlcd/
#    (sequential; ~1.35 GB zip per year; only the .tif is kept).
# 2. For each year, the share of every class in every 4 km cell, from one read of the 30 m
#    data (parallel, one year per worker): static/nlcd/share_<class>_<year>.tif   [cached]
# 3. Average shares across years and apply the rule (rebuilt every run, so changing the
#    rule in config only needs this step, not the downloads or shares):
#      static/nlcd_shares.tif  one band per class, mean share of valid NLCD pixels
#      static/herb_share.tif   herb_share (core: grass + shrub) and crop_share of land
#                              area, herb_range (max - min across years), open_herb_share
#                              (broad: + pasture/hay)
#      static/herb_mask.tif    INT1U 1 = strict herbaceous cell (Cut B / sensitivity),
#                              NA elsewhere. NOT the analysis domain: regions are
#                              wall-to-wall over grid_mask.tif (decided Oct 5).
#
# Usage: Rscript scripts/03_herb_mask.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/03_herb_mask.R")
# Needs static/grid_mask.tif from 02_static.R.

suppressPackageStartupMessages({
  library(future)
  library(furrr)
})

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) {
  args[1]
} else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else {
  Sys.getenv("R_CONFIG_ACTIVE", "dev")
}

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "prism.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "nlcd.R"))

cfg <- load_cfg(profile)
log_file <- log_init("03_herb_mask", profile)
notify_init(cfg)
hm <- cfg$herb_mask
classes <- unlist(hm$classes)
years <- as.integer(unlist(hm$nlcd_years))
broad <- unlist(hm$broad_classes %||% hm$herb_classes)
stopifnot(all(unlist(hm$herb_classes) %in% names(classes)), all(broad %in% names(classes)),
          "water" %in% names(classes), "crops" %in% names(classes))

sdir <- static_dir(cfg)
mask_file <- file.path(sdir, "grid_mask.tif")
if (!file.exists(mask_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
share_dir <- file.path(sdir, "nlcd")
dir.create(share_dir, showWarnings = FALSE, recursive = TRUE)
share_file <- function(cl, yr) file.path(share_dir, sprintf("share_%s_%d.tif", cl, yr))

log_msg("=== 03_herb_mask | profile: ", profile, " | host: ", Sys.info()[["nodename"]])
log_msg("NLCD years: ", paste(years, collapse = ","), " | classes: ",
        paste(names(classes), classes, sep = "=", collapse = " "))
log_msg("rule: herb (", paste(unlist(hm$herb_classes), collapse = "+"), ") >= ",
        hm$min_herb_share, " of land, crops < ", hm$max_crop_share)
t0 <- Sys.time()

# ---- 1. Downloads (sequential) ------------------------------------------------------
nlcd_tifs <- setNames(character(length(years)), years)
for (yr in years) {
  need <- !all(file.exists(vapply(names(classes), share_file, "", yr = yr)))
  nlcd_tifs[as.character(yr)] <- if (need) nlcd_fetch_year(cfg, yr) else NA_character_
}

# ---- 2. Class shares per year (parallel: one year per worker, all classes in one pass)
jobs <- list()
for (yr in years) {
  outs <- vapply(names(classes), share_file, "", yr = yr)
  if (!all(file.exists(outs))) {
    jobs[[length(jobs) + 1]] <- list(year = yr, tif = nlcd_tifs[[as.character(yr)]], outs = outs)
  }
}
log_msg(length(jobs), " NLCD year(s) to aggregate")
if (length(jobs)) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  plan(multisession, workers = min(cfg$workers, length(jobs)))
  fact <- hm$agg_factor; memfrac <- cfg$memfrac
  res <- future_map(jobs, function(j) nlcd_year_shares(j$tif, classes, j$outs, mask_file,
                                                       fact, memfrac),
                    .options = furrr_options(seed = NULL, packages = "terra"))
  plan(sequential)
  for (k in seq_along(jobs)) {
    j <- jobs[[k]]; r <- res[[k]]
    if (isTRUE(r$ok)) {
      log_msg(sprintf("NLCD %d shares (mean over land): %s", j$year,
                      paste(sprintf("%s %.3f", names(r$means), r$means), collapse = ", ")))
    } else log_err("NLCD ", j$year, " aggregation failed: ", r$msg)
  }
  if (!all(vapply(res, function(r) isTRUE(r$ok), TRUE))) stop("some NLCD years failed")
}

# ---- 3. Combine across years and apply the rule -------------------------------------
terra::terraOptions(memfrac = 0.5, progress = 0)
mask <- terra::rast(mask_file)
mean_share <- function(cl) terra::app(terra::rast(vapply(years, share_file, "", cl = cl)),
                                      "mean", na.rm = TRUE)
shares <- terra::rast(lapply(names(classes), mean_share))
names(shares) <- names(classes)

land <- 1 - shares[["water"]]
land <- terra::ifel(land <= 0.01, NA, land)            # (nearly) all-water cells drop out
herb_of_land <- function(s) sum(s[[unlist(hm$herb_classes)]]) / land
herb <- herb_of_land(shares)
crop <- shares[["crops"]] / land
# Stability: herb share of land in each year separately, then the range.
herb_yr <- terra::rast(lapply(years, function(yr) {
  s <- terra::rast(vapply(names(classes), share_file, "", yr = yr)); names(s) <- names(classes)
  sum(s[[unlist(hm$herb_classes)]]) / terra::ifel((1 - s[["water"]]) <= 0.01, NA, 1 - s[["water"]])
}))
herb_range <- max(herb_yr) - min(herb_yr)

open_herb <- sum(shares[[broad]]) / land
hs <- c(herb, crop, herb_range, open_herb)
names(hs) <- c("herb_share", "crop_share", "herb_range", "open_herb_share")
herb_mask <- terra::ifel(herb >= hm$min_herb_share & crop < hm$max_crop_share, 1L, NA)
herb_mask <- terra::mask(herb_mask, mask); names(herb_mask) <- "herb_mask"

wr <- function(x, name, dt = "FLT4S") {
  f <- file.path(sdir, name); tmp <- paste0(f, ".tmp.tif")
  terra::writeRaster(x, tmp, overwrite = TRUE, datatype = dt,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  unlink(f); invisible(file.rename(tmp, f))
}
wr(shares, "nlcd_shares.tif"); wr(hs, "herb_share.tif"); wr(herb_mask, "herb_mask.tif", "INT1U")

# ---- Summary ------------------------------------------------------------------------
n <- function(x) as.integer(terra::global(x, "notNA")[1, 1])
n_land <- n(mask); n_mask <- n(herb_mask)
log_msg(sprintf("land cells %d | herb mask %d cells (%.1f %% of land)", n_land, n_mask,
                100 * n_mask / n_land))
for (t in c(0.25, 0.5, 0.75)) {
  log_msg(sprintf("  herb share >= %.2f: %d cells | with crops < %.2f: %d cells", t,
                  n(terra::ifel(herb >= t, 1, NA)), hm$max_crop_share,
                  n(terra::ifel(herb >= t & crop < hm$max_crop_share, 1, NA))))
}
for (k in seq_along(years)) {
  log_msg(sprintf("  %d alone: %d cells pass the rule", years[k],
                  n(terra::ifel(herb_yr[[k]] >= hm$min_herb_share & crop < hm$max_crop_share, 1, NA))))
}
for (t in c(0.25, 0.5)) {
  log_msg(sprintf("  open herb share (%s) >= %.2f: %d cells", paste(broad, collapse = "+"), t,
                  n(terra::ifel(open_herb >= t, 1, NA))))
}
mean_cls <- terra::global(terra::mask(shares, mask), "mean", na.rm = TRUE)[, 1]
log_msg("mean class shares over land cells: ",
        paste(sprintf("%s %.3f", names(classes), mean_cls), collapse = ", "))
log_msg(sprintf("cells whose herb share changes by > 0.2 across years: %d",
                n(terra::ifel(herb_range > 0.2, 1, NA))))
msg <- sprintf("herb mask: %d cells (%.1f %% of land) in %.1f min", n_mask,
               100 * n_mask / n_land, as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("03_herb_mask finished", msg, priority = 3, tags = "white_check_mark")
