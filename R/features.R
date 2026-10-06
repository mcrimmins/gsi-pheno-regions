# Shared helpers for Cut A feature blocks (per-year rasters on the 4 km grid).

features_dir <- function(cfg, block) file.path(cfg$run_dir, "features", block)
feature_year_file <- function(cfg, block, year) {
  file.path(features_dir(cfg, block), sprintf("%s_%d.tif", block, year))
}
summary_file <- function(cfg, block, stat) {
  file.path(cfg$run_dir, "summaries", sprintf("%s_%s.tif", block, stat))
}
# Feature years: water years need Oct-Dec of the previous PRISM year.
feature_years <- function(cfg) seq(min(cfg$years) + 1L, max(cfg$years))

# Extraterrestrial radiation Ra (MJ m-2 d-1), FAO-56 eq. 21. lat in degrees (vector),
# doy vector; returns length(lat) x length(doy) matrix.
ra_matrix <- function(lat, doy) {
  phi <- lat * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
  dec <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmin(pmax(-tan(outer(phi, rep(1, length(doy)))) *
                        tan(outer(rep(1, length(phi)), dec)), -1), 1))
  sp <- outer(sin(phi), sin(dec)); cp <- outer(cos(phi), cos(dec))
  (24 * 60 / pi) * 0.0820 * matrix(dr, length(phi), length(doy), byrow = TRUE) *
    (ws * sp + cp * sin(ws))
}

# Hargreaves PET (mm/day): 0.0023 * 0.408 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin)
hargreaves <- function(tmin, tmax, ra) {
  pmax(0.0023 * 0.408 * ra * ((tmin + tmax) / 2 + 17.8) * sqrt(pmax(tmax - tmin, 0)), 0)
}

# Read rows [row, row + nrows) of selected bands as a cells x bands matrix (scale applied).
read_rows <- function(r, row, nrows, bands) {
  x <- r[[bands]]
  terra::readStart(x); on.exit(terra::readStop(x))
  terra::readValues(x, row = row, nrows = nrows, mat = TRUE)
}

# Last column index where m is TRUE per row (0 if none); first index (ncol + 1 if none).
last_true <- function(m) {
  if (!ncol(m)) return(rep(0, nrow(m)))
  Reduce(pmax, lapply(seq_len(ncol(m)), function(j) ifelse(m[, j] %in% TRUE, j, 0)))
}
first_true <- function(m) {
  if (!ncol(m)) return(rep(1, nrow(m)))
  Reduce(pmin, lapply(seq_len(ncol(m)), function(j) ifelse(m[, j] %in% TRUE, j, ncol(m) + 1)))
}

# Write a cells x layers matrix (full grid order) as a multi-band GeoTIFF on the grid.
write_feature_matrix <- function(vals, grid, names, out) {
  r <- terra::rast(grid, nlyrs = ncol(vals))
  terra::values(r) <- vals
  names(r) <- names
  r <- terra::mask(r, grid)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(out, ".tmp.tif")
  terra::writeRaster(r, tmp, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "INTERLEAVE=BAND"))
  unlink(out); file.rename(tmp, out)
  invisible(out)
}

# Across-year summaries per feature: median and IQR (q75 - q25) over the yearly files.
summarise_years <- function(files, out_median, out_iqr, chunk_rows = 60) {
  rs <- lapply(files, terra::rast)
  nm <- names(rs[[1]]); nl <- length(nm)
  grid <- rs[[1]][[1]]
  med <- matrix(NA_real_, terra::ncell(grid), nl); iqr <- med
  for (row in seq(1, terra::nrow(grid), by = chunk_rows)) {
    nr <- min(chunk_rows, terra::nrow(grid) - row + 1)
    cells <- terra::cellFromRowCol(grid, row, 1):terra::cellFromRowCol(grid, row + nr - 1, terra::ncol(grid))
    arr <- vapply(rs, function(r) read_rows(r, row, nr, seq_len(nl)),
                  matrix(0, length(cells), nl))           # cells x layers x years
    for (k in seq_len(nl)) {
      q <- apply(arr[, k, , drop = FALSE], 1, stats::quantile, probs = c(0.25, 0.5, 0.75),
                 na.rm = TRUE, names = FALSE)
      med[cells, k] <- q[2, ]; iqr[cells, k] <- q[3, ] - q[1, ]
    }
  }
  mask <- terra::rast(files[1])[[1]]
  write_feature_matrix(med, mask, nm, out_median)
  write_feature_matrix(iqr, mask, nm, out_iqr)
}
