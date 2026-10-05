# The analysis grid is the PRISM 4 km grid as written by step 01 for this run (cropped to
# the profile's bbox). Every other layer is matched to it.

static_dir <- function(cfg) file.path(cfg$run_dir, "static")
shared_raw_dir <- function(cfg) file.path(cfg$data_root, "_shared", "raw")

# Path of a PRISM variable-year file to use as the grid reference (first one found).
grid_reference_file <- function(cfg) {
  for (v in cfg$prism$vars) for (yr in cfg$years) {
    f <- prism_out_file(cfg, v, yr)
    if (file.exists(f)) return(f)
  }
  stop("No PRISM output found for profile '", cfg$profile, "'; run 01_prism_daily.R first")
}

# Land mask on the analysis grid: 1 where PRISM has data, NA elsewhere (ocean, outside CONUS).
grid_mask <- function(cfg) {
  f <- file.path(static_dir(cfg), "grid_mask.tif")
  if (!file.exists(f)) stop("grid_mask.tif missing; run 02_static.R first")
  terra::rast(f)
}

# Does `x` share the template's CRS, resolution and cell alignment (origin)?
grid_aligned <- function(x, tmpl, tol = 1e-6) {
  same_crs <- terra::same.crs(x, tmpl)
  same_res <- all(abs(terra::res(x) - terra::res(tmpl)) < tol)
  # origin difference modulo the cell size (0 when cell edges line up)
  d <- (terra::origin(x) - terra::origin(tmpl)) %% terra::res(tmpl)
  d <- pmin(d, terra::res(tmpl) - d)
  list(ok = same_crs && same_res && all(d < tol), same_crs = same_crs,
       same_res = same_res, origin_offset = d)
}

# Bring a layer onto the analysis grid: crop if aligned, else resample (and say so).
to_grid <- function(x, tmpl, method = "bilinear", label = "layer") {
  a <- grid_aligned(x, tmpl)
  if (a$ok) {
    y <- terra::crop(terra::extend(x, tmpl), tmpl, snap = "near")
    if (!terra::compareGeom(y, tmpl, stopOnError = FALSE)) {
      y <- terra::resample(y, tmpl, method = "near")
    }
    attr(y, "how") <- "cropped (grids aligned)"
  } else {
    if (!a$same_crs) x <- terra::project(x, terra::crs(tmpl))
    y <- terra::resample(x, tmpl, method = method)
    attr(y, "how") <- sprintf("resampled (%s): crs %s, res %s, origin offset %s",
                              method, if (a$same_crs) "same" else "differs",
                              if (a$same_res) "same" else "differs",
                              paste(signif(a$origin_offset, 3), collapse = ","))
  }
  y
}
