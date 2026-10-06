#!/usr/bin/env Rscript
# One-off repair: Cut B outputs written before Oct 6 2026 had the MODIS scale factor
# applied twice (values 10,000x too small; NDVI ~4e-5 instead of ~0.4). Multiplies the
# affected variables back by 10,000. Idempotent: a file whose median |value| already
# looks physical is skipped, so it is safe to run more than once.
#
# Usage: Rscript scripts/04b_cutb_fix_scale.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/04b_cutb_fix_scale.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "cutb.R"))
cfg <- load_cfg(profile)
log_init("04b_cutb_fix_scale", profile)
terra::terraOptions(progress = 0)

vars <- intersect(c("ndvi", "evi", "red", "nir", "mir"), unlist(cfg$cutb$value_vars))
n_fixed <- 0L
for (p in unlist(cfg$cutb$products)) for (v in vars) for (y in cfg$years) {
  f <- cutb_out_file(cfg, p, v, y)
  if (!file.exists(f)) next
  r <- terra::rast(f)
  m <- stats::median(abs(terra::values(r[[1]])), na.rm = TRUE)
  if (!is.finite(m) || m > 0.001) { log_msg(basename(f), ": looks physical (median |v| ", signif(m, 3), "), skip"); next }
  tmp <- paste0(f, ".tmp.tif")
  terra::writeRaster(r * 10000, tmp, overwrite = TRUE, datatype = "FLT4S", names = names(r),
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "INTERLEAVE=BAND"))
  rm(r); gc(); unlink(f); file.rename(tmp, f); n_fixed <- n_fixed + 1L
  log_msg(basename(f), ": rescaled x10000 (median |v| was ", signif(m, 3), ")")
}
log_msg("=== done | ", n_fixed, " file(s) rescaled")
