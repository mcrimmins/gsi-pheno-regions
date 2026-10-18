# Fire activity study helpers (scripts 60-69; config `fire_study:`; scope in the project doc
# claude/gsi-fire-study-scope.md). Data live under <run_dir>/fire_study/ (never in git); raw
# downloads in <data_root>/_shared/raw/fire/ so every profile shares one copy.
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

fire_dir <- function(cfg, ...) file.path(cfg$run_dir, "fire_study", ...)
fire_cache_dir <- function(cfg, ...) file.path(cfg$data_root, "_shared", "raw", "fire", ...)

# Download url to dest once (temp file, renamed on success). Never re-fetches a file on disk.
fire_download <- function(url, dest, timeout_sec = 3600, ua = "gsi-pheno-regions (R httr2)") {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  httr2::request(url) |>
    httr2::req_user_agent(ua) |>
    httr2::req_timeout(timeout_sec) |>
    httr2::req_retry(max_tries = 2, backoff = function(i) 60) |>
    httr2::req_perform(path = tmp)
  if (!file.exists(tmp) || file.size(tmp) == 0) stop("empty download: ", url)
  file.rename(tmp, dest)
  invisible(dest)
}

# Tables: Parquet when arrow is installed (FIRMS needs it anyway), else .rds. `stem` has no
# extension; the reader takes whichever exists (Parquet first).
fire_write_table <- function(df, stem) {
  dir.create(dirname(stem), showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("arrow", quietly = TRUE)) {
    f <- paste0(stem, ".parquet"); tmp <- paste0(stem, ".tmp.parquet")
    arrow::write_parquet(df, tmp)
  } else {
    f <- paste0(stem, ".rds"); tmp <- paste0(stem, ".tmp.rds")
    saveRDS(df, tmp)
  }
  if (file.exists(f)) unlink(f)
  file.rename(tmp, f)
  f
}
fire_table_file <- function(stem) {
  f <- paste0(stem, c(".parquet", ".rds")); f <- f[file.exists(f)]
  if (length(f)) f[1] else NA_character_
}
fire_read_table <- function(stem) {
  f <- fire_table_file(stem)
  if (is.na(f)) stop("no table ", stem, ".parquet / .rds")
  if (grepl("\\.parquet$", f)) as.data.frame(arrow::read_parquet(f)) else readRDS(f)
}

# Dates in FPA-FOD can arrive as Date, POSIXct or text depending on the format and driver.
# Parse what we can; the caller cross-checks against FIRE_YEAR + DISCOVERY_DOY.
fire_parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.character(x) && all(grepl("^-?[0-9]+(\\.[0-9]+)?$", x[!is.na(x) & nzchar(x)]))) x <- as.numeric(x)
  if (is.numeric(x)) {          # epoch days, seconds or milliseconds (MTBS Ig_Date differs between copies)
    d <- ifelse(abs(x) > 1e11, x / 86400000, ifelse(abs(x) > 1e7, x / 86400, x))
    return(as.Date(floor(d), origin = "1970-01-01"))
  }
  x <- trimws(as.character(x))
  out <- as.Date(substr(x, 1, 10), format = "%Y-%m-%d")
  for (fm in c("%Y/%m/%d", "%m/%d/%Y")) {
    i <- is.na(out) & !is.na(x) & nzchar(x)
    if (!any(i)) break
    out[i] <- as.Date(sub("[ T].*$", "", x[i]), format = fm)
  }
  out
}

# Region (A2 k13 map number) for every grid cell. The rules package (40) stores it with the
# typical dates; without it, rebuild from the cluster rasters exactly as 40 does.
fire_region_values <- function(cfg, grid, partition = "A2_groups", k = 13) {
  td <- file.path(cfg$run_dir, "rules", "typical_dates.tif")
  if (file.exists(td)) {
    r <- terra::rast(td)
    if ("region" %in% names(r)) return(list(v = terra::values(r[["region"]], mat = FALSE), src = "rules/typical_dates.tif"))
  }
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  D <- region_display(cfg, grid, list(list(partition = partition, k = k)))[[1]]
  v <- rep(NA_real_, terra::ncell(grid)); v[land] <- D$num
  list(v = v, src = sprintf("region_display(%s, k%d)", partition, k))
}

# Herbaceous fuel share of land per grid cell (decided Oct 17 2026): sum of the NLCD class shares
# in config fire_study: fuel_classes (grass, shrub, pasture, wetland_herb) from
# static/nlcd_shares.tif, divided by the land share (1 - water). Returns a vector over all cells.
fire_fuel_share <- function(cfg) {
  f <- file.path(cfg$run_dir, "static", "nlcd_shares.tif")
  if (!file.exists(f)) stop("no ", f, " (run 03_herb_mask.R)")
  r <- terra::rast(f)
  cls <- unlist(cfg$fire_study$fuel_classes %||% list("grass", "shrub", "pasture", "wetland_herb"))
  miss <- setdiff(c(cls, "water"), names(r))
  if (length(miss)) stop("nlcd_shares.tif lacks band(s) ", paste(miss, collapse = ", "))
  v <- terra::values(r[[c(cls, "water")]])
  land <- 1 - v[, "water"]; land[land <= 0.01] <- NA     # as 03_herb_mask.R
  out <- rowSums(v[, cls, drop = FALSE]) / land
  out[!is.finite(out)] <- NA
  pmin(pmax(out, 0), 1)
}

# Discovery / analysis week (decided Oct 18 2026; config fire_study: week): fixed 7-day blocks
# from Jan 1, days 365-366 folded into week 52, so every year has the same 52 weeks.
fire_week <- function(date) {
  doy <- as.integer(format(as.Date(date), "%j"))
  pmin(52L, (doy - 1L) %/% 7L + 1L)
}
