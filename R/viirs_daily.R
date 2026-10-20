# Daily VIIRS surface reflectance (VNP09GA) for the fire study's monitoring test (Oct 20 2026).
# Used by scripts/05_viirs_daily.R (download + aggregation) and 65c_fire_viirs_daily.R (weekly states).
# Config: fire_study: viirs_daily. AppEEARS client: R/appeears.R.

vd_conf <- function(cfg) {
  v <- cfg$fire_study$viirs_daily %||% list()
  v$product <- v$product %||% "VNP09GA.002"
  v$regions <- as.integer(unlist(v$regions %||% list(3, 10)))
  v$years <- v$years %||% list(start = 2012, end = 2024)
  v$clear <- v$clear %||% list()
  v
}
vd_short <- function(product) sub("\\..*$", "", product)
vd_dir <- function(cfg, ...) file.path(cfg$run_dir, if (nzchar(Sys.getenv("GSI_VD_TEST"))) "viirs_daily_test" else "viirs_daily", ...)
vd_out_file <- function(cfg, product, region, year) vd_dir(cfg, vd_short(product), sprintf("r%02d", region), sprintf("daily_%d.rds", year))
vd_raw_dir <- function(cfg, product, region, year) file.path(cfg$scratch_dir, "raw", "viirs_daily", vd_short(product), sprintf("r%02d", region), year)

# Herbaceous phenology cells of a pilot region, as 65 / 65b: rules `herb` cells, or in pheno_source
# regions (10, 13) cells with fuel share >= min_fuel_share. Returns cell numbers on the full grid.
vd_cells <- function(cfg, region) {
  fs <- cfg$fire_study; ps <- fs$pheno_source %||% list()
  TD <- terra::rast(file.path(cfg$run_dir, "rules", "typical_dates.tif"))
  tdv <- terra::values(TD[[c("herb", "region")]]); reg <- as.integer(tdv[, "region"])
  if (region %in% as.integer(unlist(ps$regions))) {
    fuel <- fire_fuel_share(cfg)
    which(reg == region & !is.na(fuel) & fuel >= as.numeric(ps$min_fuel_share %||% 0.35))
  } else which(reg == region & tdv[, "herb"] %in% 1)
}

# Bounding box of a set of grid cells (cell edges, plus a small pad in degrees).
vd_bbox <- function(grid, cells, pad = 0.05) {
  xy <- terra::xyFromCell(grid, cells); h <- terra::res(grid) / 2
  c(xmin = min(xy[, 1]) - h[1] - pad, ymin = min(xy[, 2]) - h[2] - pad,
    xmax = max(xy[, 1]) + h[1] + pad, ymax = max(xy[, 2]) + h[2] + pad)
}

# Resolve configured short names (regex) to exact AppEEARS layer names; stop on 0 or > 1 matches.
vd_resolve_layers <- function(cfg, product, pats) {
  meta <- ae_product_layers(cfg, product)
  out <- vapply(names(pats), function(nm) {
    hit <- grep(pats[[nm]], names(meta), value = TRUE)
    if (length(hit) != 1) stop(sprintf("%s: pattern '%s' (%s) matched %d layers. Available: %s", product,
                                       pats[[nm]], nm, length(hit), paste(names(meta), collapse = ", ")))
    hit
  }, "")
  list(layers = out, meta = meta[out])
}

# One AppEEARS area task for a date range (Date objects).
vd_submit <- function(cfg, token, product, layers, start, end, bbox, task_name) {
  body <- list(task_type = "area", task_name = task_name,
    params = list(dates = list(list(startDate = format(start, "%m-%d-%Y"), endDate = format(end, "%m-%d-%Y"))),
                  layers = lapply(unname(layers), function(l) list(product = product, layer = l)),
                  output = list(format = list(type = "geotiff"), projection = "geographic"),
                  geo = ae_bbox_geojson(bbox)))
  resp <- ae_req(cfg, "task", token) |> httr2::req_body_json(body, auto_unbox = TRUE) |> httr2::req_perform()
  httr2::resp_body_json(resp)$task_id
}

# Bits k .. k + n - 1 of an integer raster (no bitwAnd on SpatRasters: integer arithmetic).
vd_bits <- function(x, k, n = 1) (x %/% (2^k)) %% (2^n)

# Clear-land mask (500 m, on the I1 grid) from the VNP09GA quality layers. Rules (config clear:):
# day (QF1 bit 4 = 0); cloud confidence QF1 bits 2-3 <= max_cloud_conf (0 confident clear,
# 1 probably clear); land background (QF2 bits 0-2 in 0-1); no cloud shadow (QF2 bit 3), snow/ice
# (QF2 bit 5; QF7 bit 0), thin cirrus (QF2 bits 6-7, if no_cirrus); not adjacent to cloud (QF7 bit 1,
# if no_adjacent); I1 and I2 good quality (QF6 bits 3-4); sensor zenith <= max_vza degrees.
vd_clear <- function(q, i1, rules) {
  near <- function(x) if (is.null(x)) NULL else terra::resample(x, i1, method = "near")
  q1 <- near(q$qf1); q2 <- near(q$qf2); q6 <- near(q$qf6); q7 <- near(q$qf7); vz <- near(q$vza)
  ok <- vd_bits(q1, 4) == 0 & vd_bits(q1, 2, 2) <= as.integer(rules$max_cloud_conf %||% 0) &
    vd_bits(q2, 0, 3) <= 1 & vd_bits(q2, 3) == 0 & vd_bits(q2, 5) == 0
  if (isTRUE(rules$no_cirrus %||% TRUE)) ok <- ok & vd_bits(q2, 6, 2) == 0
  if (!is.null(q7)) {
    ok <- ok & vd_bits(q7, 0) == 0
    if (isTRUE(rules$no_adjacent %||% TRUE)) ok <- ok & vd_bits(q7, 1) == 0
  }
  if (!is.null(q6)) ok <- ok & vd_bits(q6, 3, 2) == 0
  if (!is.null(vz)) ok <- ok & vz <= as.numeric(rules$max_vza %||% 55)
  land <- vd_bits(q2, 0, 3) <= 1
  list(clear = ok, land = land)
}

# Read a layer: label lon/lat with the grid CRS (as cutb_aggregate_year: no datum shift); set fill
# to NA; apply the scale factor when the file holds raw integers.
vd_read <- function(f, grid, fill = NULL, scale = 1) {
  x <- terra::rast(f)
  if (terra::is.lonlat(x) && terra::is.lonlat(grid)) terra::crs(x) <- terra::crs(grid)
  if (!is.null(fill) && !is.na(fill)) x <- terra::classify(x, cbind(as.numeric(fill), NA))
  already <- grepl("FLT", terra::datatype(x)[1]) || terra::scoff(x)[1, 1] != 1
  if (!already && scale != 1) x <- x * scale
  x
}

# One day -> per-cell mean NDVI of clear 500 m pixels, count of clear pixels, count of observed land
# pixels, for the region's cells (cell numbers on the full grid). files: named list of paths for the
# day (i1, i2, qf1, qf2 required; qf6, qf7, vza optional). meta: AppEEARS layer metadata by short name.
vd_day <- function(files, meta, gridc, cidx, rules) {
  m <- function(nm, key) { x <- meta[[nm]][[key]]; if (is.null(x)) NA else as.numeric(x) }
  i1 <- vd_read(files$i1, gridc, m("i1", "FillValue"), m("i1", "ScaleFactor"))
  i2 <- vd_read(files$i2, gridc, m("i2", "FillValue"), m("i2", "ScaleFactor"))
  rq <- function(nm) if (is.null(files[[nm]])) NULL else vd_read(files[[nm]], gridc, m(nm, "FillValue"),
                                                                   if (nm == "vza") m(nm, "ScaleFactor") else 1)
  q <- list(qf1 = rq("qf1"), qf2 = rq("qf2"), qf6 = rq("qf6"), qf7 = rq("qf7"), vza = rq("vza"))
  if (!terra::compareGeom(i1, i2, stopOnError = FALSE)) i2 <- terra::resample(i2, i1, method = "near")
  cl <- vd_clear(q, i1, rules)
  okr <- i1 >= -0.01 & i1 <= 1.6 & i2 >= -0.01 & i2 <= 1.6
  clear <- cl$clear & okr
  nd <- terra::values((i2 - i1) / (i2 + i1), mat = FALSE)
  cv <- terra::values(clear, mat = FALSE) %in% c(1, TRUE)
  lv <- terra::values(cl$land & !is.na(i1), mat = FALSE) %in% c(1, TRUE)
  # 500 m pixel -> 4 km cell by pixel centre (exact; terra's "average" warp mishandles NA blocks)
  idx <- terra::cellFromXY(gridc, terra::xyFromCell(i1, seq_len(terra::ncell(i1))))
  nb <- terra::ncell(gridc)
  sumf <- function(val, keep) {
    k <- keep & !is.na(idx); out <- numeric(nb)
    if (any(k)) { r <- rowsum(val[k], idx[k]); out[as.integer(rownames(r))] <- r[, 1] }
    out
  }
  cv <- cv & is.finite(nd)
  n_c <- sumf(rep(1, length(cv)), cv); s_nd <- sumf(nd, cv); n_l <- sumf(rep(1, length(lv)), lv)
  ndvi <- ifelse(n_c > 0, s_nd / n_c, NA_real_)
  list(ndvi = ndvi[cidx], n_clear = n_c[cidx], n_land = n_l[cidx])
}

# Files of one AppEEARS bundle by layer short name and day (yyyyddd), for days with i1, i2, qf1, qf2.
# layers: named vector short name -> AppEEARS layer name.
vd_files_by_day <- function(raw_dir, layers) {
  f_all <- list.files(raw_dir, pattern = "\\.tif$", full.names = TRUE)
  by <- lapply(layers, function(l) {
    f <- f_all[grepl(paste0(l, "_doy\\d{7}"), basename(f_all))]
    stats::setNames(f, sub(".*_doy(\\d{7}).*", "\\1", basename(f)))
  })
  days <- sort(Reduce(intersect, lapply(by[c("i1", "i2", "qf1", "qf2")], names)))
  files <- lapply(stats::setNames(days, days), function(d) {
    l <- lapply(by, function(x) if (d %in% names(x)) unname(x[[d]]) else NULL)
    l[!vapply(l, is.null, TRUE)]
  })
  list(days = days, files = files)
}

vd_doy_date <- function(yyyyddd) as.Date(as.integer(substr(yyyyddd, 5, 7)) - 1, origin = as.Date(sprintf("%s-01-01", substr(yyyyddd, 1, 4))))
