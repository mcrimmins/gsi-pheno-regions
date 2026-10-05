# Annual NLCD land cover -> class shares on the 4 km analysis grid.

nlcd_dir <- function(cfg) file.path(shared_raw_dir(cfg), "nlcd")

nlcd_tif_path <- function(cfg, year) {
  file.path(nlcd_dir(cfg), sub("\\.zip$", ".tif", basename(
    gsub("{year}", year, cfg$herb_mask$nlcd_url, fixed = TRUE))))
}

# Download + unzip one NLCD year (sequential, resume-safe). Keeps only the land cover .tif.
nlcd_fetch_year <- function(cfg, year) {
  tif <- nlcd_tif_path(cfg, year)
  if (file.exists(tif)) return(tif)
  dir.create(nlcd_dir(cfg), showWarnings = FALSE, recursive = TRUE)
  url <- gsub("{year}", year, cfg$herb_mask$nlcd_url, fixed = TRUE)
  zip <- file.path(nlcd_dir(cfg), basename(url))
  if (!file.exists(zip)) {
    log_msg("NLCD ", year, ": downloading ", url)
    t0 <- Sys.time()
    httr2::request(url) |>
      httr2::req_user_agent(cfg$prism$user_agent) |>
      httr2::req_timeout(3600) |>
      httr2::req_retry(max_tries = 3) |>
      httr2::req_perform(path = paste0(zip, ".part"))
    file.rename(paste0(zip, ".part"), zip)
    log_msg(sprintf("NLCD %d: %.0f MB in %.1f min", year, file.size(zip) / 1e6,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  inside <- utils::unzip(zip, list = TRUE)$Name
  want <- grep("\\.tif$", inside, value = TRUE)
  if (length(want) != 1) stop("expected one .tif in ", basename(zip), ": ", paste(inside, collapse = ", "))
  utils::unzip(zip, files = want, exdir = nlcd_dir(cfg), junkpaths = TRUE)
  if (basename(want) != basename(tif)) file.rename(file.path(nlcd_dir(cfg), basename(want)), tif)
  if (!file.exists(tif)) stop("unzip did not produce ", tif)
  unlink(zip)
  tif
}

# Shares of all configured NLCD classes per analysis-grid cell for one year, in a single
# read of the 30 m data (no full-resolution temporary files; CONUS is ~17e9 pixels).
# Runs in a worker: paths and plain values in, paths out.
#   1. read window = grid extent projected to Albers + margin (no copy)
#   2. aggregate by `fact` with a function returning every class share at once
#      (share of valid pixels; NLCD 0/250 treated as nodata)
#   3. area-average onto the 4 km lon/lat grid, mask to PRISM land cells, one file per class
nlcd_year_shares <- function(nlcd_tif, codes, out_files, mask_file, fact, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    lc <- terra::rast(nlcd_tif)
    mask <- terra::rast(mask_file)
    e <- terra::ext(terra::project(terra::as.polygons(terra::ext(mask), crs = terra::crs(mask)),
                                   terra::crs(lc)))
    margin <- 5000   # m
    e <- terra::ext(e[1] - margin, e[2] + margin, e[3] - margin, e[4] + margin)
    terra::window(lc) <- terra::intersect(e, terra::ext(lc))
    k <- length(codes)
    share_fun <- function(v, ...) {
      v[v == 0 | v == 250] <- NA
      ok <- !is.na(v); n <- sum(ok)
      if (n == 0) return(rep(NA_real_, k))
      v <- v[ok]
      vapply(codes, function(cc) sum(v == cc) / n, numeric(1))
    }
    a <- terra::aggregate(lc, fact = fact, fun = share_fun)
    if (terra::nlyr(a) != k) stop("aggregate returned ", terra::nlyr(a), " layers, expected ", k)
    s <- terra::mask(terra::project(a, mask, method = "average"), mask)
    names(s) <- names(codes)
    means <- numeric(k)
    for (i in seq_len(k)) {
      tmp <- paste0(out_files[i], ".tmp.tif")
      terra::writeRaster(s[[i]], tmp, overwrite = TRUE, datatype = "FLT4S",
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
      unlink(out_files[i]); file.rename(tmp, out_files[i])
      means[i] <- terra::global(terra::rast(out_files[i]), "mean", na.rm = TRUE)[1, 1]
    }
    list(ok = TRUE, means = setNames(means, names(codes)))
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}
