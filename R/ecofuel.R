# Fuel phenology by ecoregion x fuel type (scripts 70-79; config `ecofuel:`; scope in the project doc
# claude/gsi-ecofuel-climatology-handoff.md). New work, separate from the NPN study: it may READ the
# NPN study's site tables (npn_study/data/) and the fire study's tables, and never writes to either.
# Outputs live under <run_dir>/ecofuel/ (never in git); the page report/ecofuel.qmd reads small
# copies in report/tables/ecofuel_*.csv (made on the laptop by 79).
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

ecofuel_dir <- function(cfg, ...) file.path(cfg$run_dir, "ecofuel", ...)

# ---- Scott & Burgan (2005) 40 fuel models as mapped by LANDFIRE FBFM40 ---------------------------
# One row per code: model label (band name in fbfm_shares.tif) and group (band name in
# fuel_groups.tif). Nonburnable models keep their own groups so the agriculture question
# (cropland is NB3 in FBFM40) stays open.
fbfm_table <- function() {
  m <- data.frame(
    code = c(91, 92, 93, 98, 99, 101:109, 121:124, 141:149, 161:165, 181:189, 201:204),
    model = c("NB1", "NB2", "NB3", "NB8", "NB9", paste0("GR", 1:9), paste0("GS", 1:4),
              paste0("SH", 1:9), paste0("TU", 1:5), paste0("TL", 1:9), paste0("SB", 1:4)),
    stringsAsFactors = FALSE)
  nb <- c(NB1 = "NB_urban", NB2 = "NB_snow", NB3 = "NB_ag", NB8 = "NB_water", NB9 = "NB_barren")
  m$group <- ifelse(m$model %in% names(nb), nb[m$model], substr(m$model, 1, 2))
  m
}
# Group order everywhere (bands of fuel_groups.tif, rows of the tables)
ecofuel_groups <- function() c("GR", "GS", "SH", "TU", "TL", "SB",
                               "NB_ag", "NB_urban", "NB_barren", "NB_snow", "NB_water")
ecofuel_burnable <- function() c("GR", "GS", "SH", "TU", "TL", "SB")
ecofuel_group_labels <- function() c(
  GR = "Grass (GR)", GS = "Grass-shrub (GS)", SH = "Shrub (SH)", TU = "Timber-understory (TU)",
  TL = "Timber litter (TL)", SB = "Slash-blowdown (SB)", NB_ag = "Agriculture (NB3)",
  NB_urban = "Urban (NB1)", NB_barren = "Barren (NB9)", NB_snow = "Snow / ice (NB2)",
  NB_water = "Water (NB8)")

# ---- inputs ----------------------------------------------------------------------------------
# LANDFIRE FBFM40 GeoTIFF: config ecofuel: landfire_fbfm40 (file or folder, per OS), default
# <data_root>/landfire; in a folder, the first .tif whose name or folder says FBFM40 / F40 (as 59b).
ecofuel_landfire_file <- function(cfg) {
  os <- if (.Platform$OS.type == "windows") "windows" else "linux"
  lf <- cfg$ecofuel$landfire_fbfm40; if (is.list(lf)) lf <- lf[[os]]
  if (is.null(lf) || !nzchar(lf)) lf <- file.path(cfg$data_root, "landfire")
  lf <- path.expand(lf)
  if (dir.exists(lf)) {
    f <- list.files(lf, "\\.tif$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    f <- f[grepl("f40|fbfm40", basename(f), ignore.case = TRUE) | grepl("fbfm40", dirname(f), ignore.case = TRUE)]
    want <- cfg$ecofuel$landfire_version %||% ""
    if (nzchar(want) && any(grepl(want, basename(f), fixed = TRUE))) f <- f[grepl(want, basename(f), fixed = TRUE)]
    if (!length(f)) stop("no FBFM40 .tif under ", lf, " (copy it from the laptop or run npn_study/scripts/59a)")
    lf <- f[1]
  }
  if (!file.exists(lf)) stop("LANDFIRE FBFM40 raster not found: ", lf)
  lf
}

# Ecoregion polygons (config ecofuel: ecoregions). Looks in the shared cache
# <data_root>/_shared/raw/ecoregions/l<level>/, then (read only) in the NPN study's copy
# <run_dir>/npn_study/ecoregions/ for Level II, else downloads the zip into the cache.
ecofuel_ecoregion_shp <- function(cfg) {
  eo <- cfg$ecofuel$ecoregions
  lev <- as.integer(eo$level %||% 2)
  find_shp <- function(d) {
    if (!dir.exists(d)) return(NA_character_)
    s <- list.files(d, "\\.shp$", recursive = TRUE, full.names = TRUE)
    pat <- eo$shp_pattern %||% sprintf("l%d|level%d", lev, lev)
    s <- s[grepl(pat, basename(s), ignore.case = TRUE)]
    if (length(s)) s[1] else NA_character_
  }
  cache <- file.path(cfg$data_root, "_shared", "raw", "ecoregions", paste0("l", lev))
  s <- find_shp(cache)
  if (!is.na(s)) return(s)
  if (lev == 2) {
    s <- find_shp(file.path(cfg$run_dir, "npn_study", "ecoregions"))
    if (is.na(s)) s <- find_shp(file.path(cfg$data_root, "full", "npn_study", "ecoregions"))
    if (!is.na(s)) return(s)
  }
  url <- eo$url %||% stop("config ecofuel: ecoregions: url is missing")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  zf <- file.path(cache, basename(url))
  if (!file.exists(zf)) {
    log_msg("downloading ecoregions: ", url)
    tmp <- paste0(zf, ".part")
    st <- utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
    if (st != 0 || !file.exists(tmp)) stop("ecoregion download failed: ", url)
    file.rename(tmp, zf)
  }
  utils::unzip(zf, exdir = cache)
  s <- find_shp(cache)
  if (is.na(s)) stop("no Level ", lev, " shapefile after unzipping ", basename(zf))
  s
}

# ---- 70 worker: one band of LANDFIRE rows -> ~1 km class shares ---------------------------------
# Paths and plain values in, a status list out (SpatRasters don't serialize). `win` is the band's
# extent (xmin, xmax, ymin, ymax) in the LANDFIRE CRS, aligned to its cells and to `fact`.
# Output bands: one share per model (of mapped pixels in the block) + `mapped` (share of the
# block's pixels with a valid FBFM40 code: 0 outside the map, e.g. ocean or Canada).
ecofuel_band_shares <- function(lf_file, win, codes, labels, out_file, fact, memfrac, tmpdir) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0, tempdir = tmpdir)
    lc <- terra::rast(lf_file)[[1]]
    if (any(terra::is.factor(lc))) levels(lc) <- NULL       # LANDFIRE ships a category table
    terra::window(lc) <- terra::ext(win)
    k <- length(codes)
    lut <- integer(max(codes)); lut[codes] <- seq_len(k)
    mx <- max(codes)
    share_fun <- function(v, ...) {
      n_all <- length(v)
      v <- v[!is.na(v)]
      v <- v[v >= 1 & v <= mx]
      i <- lut[v]; i <- i[i > 0L]
      n <- length(i)
      if (n == 0L) return(c(rep(NA_real_, k), 0))
      c(tabulate(i, k) / n, n / n_all)
    }
    tmp <- paste0(out_file, ".tmp.tif")
    a <- terra::aggregate(lc, fact = fact, fun = share_fun, filename = tmp, overwrite = TRUE,
                          wopt = list(datatype = "FLT4S", names = c(labels, "mapped"),
                                      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")))
    if (terra::nlyr(a) != k + 1) stop("aggregate returned ", terra::nlyr(a), " layers, expected ", k + 1)
    rm(a); gc()
    unlink(out_file); file.rename(tmp, out_file)
    list(ok = TRUE)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
}

# Bands of LANDFIRE rows covering the analysis grid: window = grid extent projected to the
# LANDFIRE CRS + margin, snapped to LANDFIRE cells; split into bands of rows_per_band * fact rows
# (multiples of fact, so the ~1 km blocks of all bands line up).
# Extent of a lon/lat raster's bounding box in another CRS, from n points along every edge. Projecting
# only the 4 corners (as 03 / R/nlcd.R do) misses the bulge of the southern parallel in Albers: on the
# CONUS grid that cut off land south of ~28 N (south Florida, south Texas; found Oct 26 2026).
ext_in_crs <- function(r, crs_to, n = 500) {
  e <- terra::ext(r)
  xs <- seq(e[1], e[2], length.out = n); ys <- seq(e[3], e[4], length.out = n)
  pts <- rbind(cbind(xs, e[3]), cbind(xs, e[4]), cbind(e[1], ys), cbind(e[2], ys))
  p <- terra::project(pts, terra::crs(r), crs_to)
  terra::ext(min(p[, 1]), max(p[, 1]), min(p[, 2]), max(p[, 2]))
}

ecofuel_bands <- function(lf_file, mask_file, fact, rows_per_band, margin_m = 5000) {
  lc <- terra::rast(lf_file)[[1]]
  mask <- terra::rast(mask_file)
  e <- ext_in_crs(mask, terra::crs(lc))
  e <- terra::ext(e[1] - margin_m, e[2] + margin_m, e[3] - margin_m, e[4] + margin_m)
  e <- terra::intersect(e, terra::ext(lc))
  if (is.null(e)) stop("the analysis grid does not overlap the LANDFIRE map")
  e <- terra::align(e, lc, snap = "out")
  r <- terra::res(lc)[2]
  h <- fact * rows_per_band * r                                 # band height in map units
  n <- ceiling((e[4] - e[3]) / h)
  lapply(seq_len(n), function(b) {
    top <- e[4] - (b - 1) * h
    bot <- max(e[3], top - h)
    c(xmin = e[1], xmax = e[2], ymin = bot, ymax = top)
  })
}

# ---- units (72+) ------------------------------------------------------------------------------------
# Cells of every climatology unit (ecoregion x fuel group) at a purity threshold: a land cell belongs to
# group g's unit in its ecoregion when g holds >= purity of the cell's land (share / (1 - water), as 71),
# and LANDFIRE maps at least min_mapped of the cell. Groups from config ecofuel: climatology (groups
# everywhere + extra_units by ecoregion). Shares >= 0.5 sum to <= 1, so a cell is in at most one unit.
# Returns a data.frame: row (index into the land cells), cell (grid cell number), eco, code, group,
# share, unit ("<code> <group>").
ecofuel_unit_cells <- function(cfg, purity) {
  ec <- cfg$ecofuel; cl <- ec$climatology
  lev <- as.integer(ec$ecoregions$level %||% 2)
  mask <- grid_mask(cfg)
  land <- which(!is.na(terra::values(mask, mat = FALSE)))
  G <- terra::values(terra::rast(ecofuel_dir(cfg, "fuel_groups.tif")))[land, , drop = FALSE]
  E <- terra::values(terra::rast(ecofuel_dir(cfg, sprintf("eco_l%d.tif", lev)))[["eco"]], mat = FALSE)[land]
  L <- utils::read.csv(ecofuel_dir(cfg, sprintf("eco_l%d_lookup.csv", lev)), stringsAsFactors = FALSE, colClasses = c(code = "character"))
  w <- G[, "NB_water"]; w[is.na(w)] <- 0
  mp <- G[, "mapped"]; mp[is.na(mp)] <- 0
  lw <- ifelse(1 - w > 0.01, 1 - w, NA)
  ok_map <- mp >= as.numeric(ec$min_mapped %||% 0.5)
  code <- L$code[match(E, L$eco)]
  extra <- cl$extra_units %||% list()
  out <- list()
  for (g in unique(c(unlist(cl$groups), names(extra)))) {
    sl <- G[, g] / lw
    inside <- !is.na(sl) & sl >= purity & ok_map & !is.na(code)
    if (!(g %in% unlist(cl$groups))) inside <- inside & code %in% unlist(extra[[g]])
    i <- which(inside)
    if (length(i)) out[[g]] <- data.frame(row = i, cell = land[i], eco = E[i], code = code[i], group = g,
                                          share = sl[i], stringsAsFactors = FALSE)
  }
  U <- do.call(rbind, out)
  U$unit <- paste(U$code, U$group)
  U
}

# Circular day difference a - b in (-182.5, 182.5]
cdiff <- function(a, b) ((a - b + 182.5) %% 365) - 182.5
