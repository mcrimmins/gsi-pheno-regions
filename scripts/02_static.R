#!/usr/bin/env Rscript
# Step 02: static layers on the analysis grid (the PRISM 4 km grid written by step 01).
#
# Outputs in <data_root>/<run_name>/static/:
#   grid_mask.tif  INT1U  1 = PRISM land cell, NA elsewhere (the analysis domain)
#   elev.tif       FLT4S  elevation (m), PRISM 4 km DEM
#   lat.tif        FLT4S  cell-centre latitude (deg N, NAD83)
#   lon.tif        FLT4S  cell-centre longitude (deg E, NAD83)
# Daylength is computed from latitude where it's needed (Block 3), not stored.
# Soils: deferred (plan doc).
#
# DEM source: PRISM supporting datasets, PRISM_us_dem_4km_bil.zip (downloaded once to
# <data_root>/_shared/raw/). It predates the Oct 2025 PRISM format change, so its grid is
# checked against the daily grid: cropped if aligned, resampled (bilinear) if not, and the
# log says which.
#
# Usage: Rscript scripts/02_static.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/02_static.R")
# Resume-safe: existing outputs are kept; delete a file to rebuild it.

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

cfg <- load_cfg(profile)
log_file <- log_init("02_static", profile)
notify_init(cfg)
terra::terraOptions(memfrac = 0.5, progress = 0)

out_dir <- static_dir(cfg)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
gdal_opts <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
out <- function(name) file.path(out_dir, name)

write_layer <- function(x, name, datatype = "FLT4S", opts = gdal_opts) {
  tmp <- paste0(out(name), ".tmp.tif")
  terra::writeRaster(x, tmp, overwrite = TRUE, datatype = datatype, gdal = opts)
  unlink(out(name))
  if (!file.rename(tmp, out(name))) stop("could not move ", tmp, " into place")
  r <- terra::rast(out(name))
  rng <- terra::global(r, "range", na.rm = TRUE)[1, ]
  log_msg(sprintf("wrote %s | %d x %d | %s | range %.3f..%.3f | %d non-NA cells",
                  name, terra::nrow(r), terra::ncol(r), datatype, rng[[1]], rng[[2]],
                  as.integer(terra::global(r, "notNA")[1, 1])))
}

log_msg("=== 02_static | profile: ", profile, " | host: ", Sys.info()[["nodename"]])
log_msg("out_dir: ", out_dir)

# ---- 1. Grid template and land mask from the PRISM output of this run -------------
ref <- grid_reference_file(cfg)
log_msg("grid reference: ", ref)
tmpl <- terra::rast(ref, lyrs = 1)
if (!file.exists(out("grid_mask.tif"))) {
  m <- terra::ifel(is.na(tmpl), NA, 1L)
  names(m) <- "grid_mask"
  write_layer(m, "grid_mask.tif", "INT1U", c("COMPRESS=DEFLATE", "TILED=YES"))
} else log_msg("grid_mask.tif exists, skip")
mask <- terra::rast(out("grid_mask.tif"))
log_msg(sprintf("grid: %d x %d, res %.6f, extent %s, crs %s, %d land cells",
                terra::nrow(mask), terra::ncol(mask), terra::res(mask)[1],
                paste(round(as.vector(terra::ext(mask)), 5), collapse = ","),
                terra::crs(mask, describe = TRUE)$name,
                as.integer(terra::global(mask, "notNA")[1, 1])))

# ---- 2. Latitude / longitude --------------------------------------------------------
if (!file.exists(out("lat.tif"))) {
  lat <- terra::mask(terra::init(mask, "y"), mask); names(lat) <- "lat"
  write_layer(lat, "lat.tif")
} else log_msg("lat.tif exists, skip")
if (!file.exists(out("lon.tif"))) {
  lon <- terra::mask(terra::init(mask, "x"), mask); names(lon) <- "lon"
  write_layer(lon, "lon.tif")
} else log_msg("lon.tif exists, skip")

# ---- 3. Elevation from the PRISM 4 km DEM -------------------------------------------
if (!file.exists(out("elev.tif"))) {
  raw_dir <- shared_raw_dir(cfg)
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  zip <- file.path(raw_dir, "PRISM_us_dem_4km_bil.zip")
  url <- "https://prism.oregonstate.edu/downloads/data/PRISM_us_dem_4km_bil.zip"
  if (!file.exists(zip)) {
    log_msg("downloading ", url)
    resp <- httr2::request(url) |>
      httr2::req_user_agent(cfg$prism$user_agent) |>
      httr2::req_timeout(300) |>
      httr2::req_retry(max_tries = 3) |>
      httr2::req_perform(path = paste0(zip, ".part"))
    file.rename(paste0(zip, ".part"), zip)
  }
  dem_dir <- file.path(raw_dir, "PRISM_us_dem_4km")
  if (!dir.exists(dem_dir)) utils::unzip(zip, exdir = dem_dir)
  bil <- list.files(dem_dir, pattern = "\\.bil$", full.names = TRUE, recursive = TRUE)
  if (length(bil) != 1) stop("expected one .bil in ", dem_dir, ", found ", length(bil))
  dem <- terra::rast(bil)
  log_msg(sprintf("DEM: %s | %d x %d | res %.6f | extent %s | crs %s", basename(bil),
                  terra::nrow(dem), terra::ncol(dem), terra::res(dem)[1],
                  paste(round(as.vector(terra::ext(dem)), 5), collapse = ","),
                  terra::crs(dem, describe = TRUE)$name))

  elev <- to_grid(dem, mask, method = "bilinear", label = "DEM")
  how <- attr(elev, "how")
  if (startsWith(how, "resampled")) log_warn("DEM ", how) else log_msg("DEM ", how)
  elev <- terra::mask(elev, mask); names(elev) <- "elev"

  # Land cells with no elevation (edge-of-domain mismatches) are reported, not filled.
  n_gap <- as.integer(terra::global(terra::ifel(is.na(elev) & !is.na(mask), 1, NA),
                                    "notNA")[1, 1])
  if (n_gap > 0) log_warn(n_gap, " land cells have no DEM value (left NA)")
  write_layer(elev, "elev.tif")
} else log_msg("elev.tif exists, skip")

log_msg("=== done | log: ", log_file)
notify("02_static finished", paste("static layers in", out_dir), priority = 3,
       tags = "white_check_mark")
