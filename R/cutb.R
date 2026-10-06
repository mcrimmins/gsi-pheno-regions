# Cut B: paths and stacking of AppEEARS MODIS VI composites.

cutb_short <- function(product) sub("\\..*$", "", product)          # MOD13C1.061 -> MOD13C1
cutb_dir <- function(cfg) file.path(cfg$run_dir, "cutb")
cutb_raw_dir <- function(cfg, product, year) {
  file.path(cfg$scratch_dir, "raw", "cutb", cutb_short(product), year)
}
cutb_out_file <- function(cfg, product, var, year) {
  file.path(cutb_dir(cfg), cutb_short(product), var,
            sprintf("%s_%s_%d.tif", cutb_short(product), var, year))
}
cutb_bbox <- function(cfg) if (is.null(cfg$bbox)) unlist(cfg$cutb$conus_bbox) else cfg$bbox

# Save the AppEEARS layer metadata (scale factors, fill values) next to the stacks.
# Stacks keep the raw integer values; later steps apply ScaleFactor from this file.
cutb_write_meta <- function(cfg, product, resolved) {
  d <- file.path(cutb_dir(cfg), cutb_short(product))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  m <- lapply(names(resolved$layers), function(nm) {
    x <- resolved$meta[[resolved$layers[[nm]]]]
    list(var = nm, layer = resolved$layers[[nm]], scale_factor = x$ScaleFactor %||% 1,
         fill_value = x$FillValue %||% NA, data_type = x$DataType %||% NA,
         units = x$Units %||% NA, description = x$Description %||% NA)
  })
  jsonlite::write_json(m, file.path(d, "layers.json"), auto_unbox = TRUE, pretty = TRUE,
                       na = "null")
  invisible(m)
}

cutb_read_meta <- function(cfg, product) {
  m <- jsonlite::read_json(file.path(cutb_dir(cfg), cutb_short(product), "layers.json"))
  setNames(m, vapply(m, `[[`, "", "var"))
}

# Output variables written per product-year: the configured value_vars plus QA summaries.
cutb_out_vars <- function(cfg) c(unlist(cfg$cutb$value_vars), "n_valid", "snow_frac")

# Aggregate one product-year of 1 km composites to the 4 km analysis grid.
# For each composite: 1 km pixels whose pixel reliability is in `valid_rel` are scaled to
# physical units (ScaleFactor) and area-averaged onto the grid; n_valid = number of valid
# 1 km pixels per cell; snow_frac = share of non-fill pixels flagged snow/ice (reliability 2).
# Output: one multi-band file per variable (bands = composites, named by composite start
# date), FLT4S (n_valid INT2U), masked to grid_mask. Runs in the main process.
cutb_aggregate_year <- function(raw_dir, layers, meta, value_vars, valid_rel, grid_file,
                                out_files) {
  grid <- terra::rast(grid_file)
  f_all <- list.files(raw_dir, pattern = "\\.tif$", full.names = TRUE)
  files_for <- function(nm) {
    f <- f_all[grepl(paste0(layers[[nm]], "_doy\\d{7}"), basename(f_all))]
    setNames(f, sub(".*_doy(\\d{7}).*", "\\1", basename(f)))
  }
  rel_f <- files_for("reliability")
  if (!length(rel_f)) stop("no pixel reliability files in ", raw_dir)
  vf <- lapply(setNames(value_vars, value_vars), files_for)
  doys <- sort(Reduce(intersect, c(list(names(rel_f)), lapply(vf, names))))
  if (!length(doys)) stop("no composites with all layers in ", raw_dir)
  dates <- as.Date(as.integer(substr(doys, 5, 7)) - 1,
                   origin = as.Date(sprintf("%s-01-01", substr(doys, 1, 4))))
  scale_of <- function(nm) {
    sf <- meta[[layers[[nm]]]]$ScaleFactor
    if (is.null(sf) || is.na(sf) || nm == "doy") 1 else as.numeric(sf)
  }
  res <- setNames(vector("list", length(value_vars) + 2), c(value_vars, "n_valid", "snow_frac"))
  for (k in seq_along(res)) res[[k]] <- vector("list", length(doys))
  # AppEEARS "geographic" output is WGS84 lon/lat; the grid is NAD83 lon/lat. At 4 km the
  # ~1 m datum difference is irrelevant, so label inputs with the grid CRS and resample
  # (no datum transformation, which can make PROJ fetch shift grids over the network).
  rd <- function(f) {
    x <- terra::rast(f)
    if (terra::is.lonlat(x) && terra::is.lonlat(grid)) terra::crs(x) <- terra::crs(grid)
    x
  }
  to_grid4 <- function(x, method) {
    if (terra::same.crs(x, grid)) terra::resample(x, grid, method = method)
    else terra::project(x, grid, method = method)
  }
  for (d in seq_along(doys)) {
    rel <- rd(rel_f[[doys[d]]])
    # Comparisons, not %in%: terra's %in% method is only used when terra is attached, and
    # these scripts call terra:: without attaching it (base %in% fails on a SpatRaster).
    ok <- Reduce(`|`, lapply(valid_rel, function(v) rel == v))
    good <- terra::ifel(ok, 1, 0)                              # NA (fill) -> NA
    res$n_valid[[d]] <- to_grid4(terra::ifel(is.na(good), 0, good), "sum")
    res$snow_frac[[d]] <- to_grid4(terra::ifel(rel == 2, 1, 0), "average")
    for (nm in value_vars) {
      x <- rd(vf[[nm]][[doys[d]]]) * scale_of(nm)
      x <- terra::mask(x, good, maskvalues = c(0, NA))
      res[[nm]][[d]] <- to_grid4(x, "average")
    }
  }
  info <- list()
  for (nm in names(res)) {
    r <- terra::mask(terra::rast(res[[nm]]), grid)
    names(r) <- format(dates, "%Y-%m-%d")
    out <- out_files[[nm]]
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    tmp <- paste0(out, ".tmp.tif")
    dt <- if (nm == "n_valid") "INT2U" else "FLT4S"
    terra::writeRaster(r, tmp, overwrite = TRUE, datatype = dt,
                       gdal = c("COMPRESS=DEFLATE", if (dt == "FLT4S") "PREDICTOR=3" else "PREDICTOR=2",
                                "TILED=YES", "INTERLEAVE=BAND"))
    if (terra::nlyr(terra::rast(tmp)) != length(doys)) stop(nm, ": band count mismatch")
    unlink(out); file.rename(tmp, out)
    info[[nm]] <- list(n = length(doys), first = min(dates), last = max(dates))
  }
  # share of grid cells with at least one valid pixel, per composite (for the log)
  nv <- terra::rast(out_files[["n_valid"]])
  info$valid_share <- as.numeric(terra::global(nv > 0, "mean", na.rm = TRUE)[, 1])
  rm(grid, rel, good, res, nv); gc()
  try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE)
  info
}
