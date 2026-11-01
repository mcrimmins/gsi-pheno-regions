# NPN x fire study: small shared helpers (tables, argument parsing).

# Write a table: Parquet with arrow, else gzipped CSV. Writes to a temp name, then renames.
# Returns the path written.
nf_write_table <- function(df, path_no_ext) {
  dir.create(dirname(path_no_ext), showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("arrow", quietly = TRUE)) {
    out <- paste0(path_no_ext, ".parquet"); tmp <- paste0(out, ".tmp")
    arrow::write_parquet(df, tmp, compression = "zstd")
  } else {
    out <- paste0(path_no_ext, ".csv.gz"); tmp <- paste0(out, ".tmp")
    con <- gzfile(tmp, "w"); utils::write.csv(df, con, row.names = FALSE); close(con)
  }
  if (file.exists(out)) file.remove(out)
  file.rename(tmp, out)
  out
}
# Path of an existing table written by nf_write_table (parquet or csv.gz), or NA.
nf_table_path <- function(path_no_ext) {
  p <- paste0(path_no_ext, c(".parquet", ".csv.gz"))
  p <- p[file.exists(p)]
  if (length(p)) p[1] else NA_character_
}
nf_read_table <- function(path) {
  if (grepl("\\.parquet$", path)) as.data.frame(arrow::read_parquet(path))
  else utils::read.csv(path, stringsAsFactors = FALSE)
}

# Tables written by the parent pipeline (fire_write_table: Parquet, else .rds). `stem` has no extension.
nf_read_parent <- function(stem) {
  f <- paste0(stem, c(".parquet", ".rds")); f <- f[file.exists(f)]
  if (!length(f)) stop("no parent table ", stem, ".parquet / .rds")
  if (grepl("\\.parquet$", f[1])) as.data.frame(arrow::read_parquet(f[1])) else readRDS(f[1])
}

# Week of year: fixed 7-day blocks from Jan 1, days 365-366 folded into week 52 (same as the fire
# study's fire_week(), so the two studies' weeks line up).
nf_week <- function(date) {
  doy <- as.integer(format(as.Date(date), "%j"))
  pmin(52L, (doy - 1L) %/% 7L + 1L)
}

# Unique 375 m pixel-day key: detections snapped to a 375 m sinusoidal grid (sphere R = 6371007.181 m),
# once per pixel and solar date (same definition as the fire study's fw_pixel_key()).
nf_pixel_key <- function(lon, lat, date, size = 375, R = 6371007.181) {
  y <- R * lat * pi / 180; x <- R * lon * pi / 180 * cos(lat * pi / 180)
  paste(floor(x / size), floor(y / size), format(date))
}

# Read a CSV written by rnpn: everything as text, then the known numeric columns typed, so every
# state-year file has the same column types (a column that is all NA in one file would otherwise
# come back logical there and text elsewhere).
nf_npn_numeric <- c("observation_id", "site_id", "latitude", "longitude", "elevation_in_meters",
                    "species_id", "individual_id", "phenophase_id", "day_of_year", "phenophase_status",
                    "intensity_category_id", "abundance_value", "observedby_person_id", "dataset_id",
                    "observed_status_conflict_flag_id")
nf_read_npn_csv <- function(path) {
  df <- utils::read.csv(path, colClasses = "character", na.strings = c("", "NA", "-9999"), check.names = FALSE)
  names(df) <- tolower(names(df))
  for (k in intersect(names(df), nf_npn_numeric)) df[[k]] <- suppressWarnings(as.numeric(df[[k]]))
  df
}

# NPN intensity values are text classes. Type: "percent" ("25-49%", "Less than 5%", "95% or more"),
# "count" ("3 to 10", "Less than 3", "More than 10,000"), else "qualitative" ("Little", "Some", "Lots").
nf_intensity_type <- function(v) {
  v <- as.character(v)
  ifelse(grepl("%", v), "percent",
  ifelse(grepl("^\\s*(less than|more than)?\\s*[0-9,]+(\\s*to\\s*[0-9,]+)?\\s*$", v, ignore.case = TRUE), "count",
  "qualitative"))
}
# Midpoint (%) of a percent class: "a-b%" -> (a + b) / 2; "Less than a%" -> a / 2;
# "a% or more" -> (a + 100) / 2; a bare "a%" -> a. NA for anything else.
nf_pct_mid <- function(v) {
  v <- tolower(gsub("\\s+", " ", trimws(as.character(v))))
  num <- function(re) suppressWarnings(as.numeric(sub(re, "\\1", v)))
  out <- rep(NA_real_, length(v))
  r <- grepl("^[0-9.]+ ?- ?[0-9.]+ ?%$", v)
  if (any(r)) { a <- as.numeric(sub("^([0-9.]+).*", "\\1", v[r])); b <- as.numeric(sub("^[0-9.]+ ?- ?([0-9.]+).*", "\\1", v[r])); out[r] <- (a + b) / 2 }
  r <- grepl("^less than [0-9.]+ ?%$", v); if (any(r)) out[r] <- num("^less than ([0-9.]+).*")[r] / 2
  r <- grepl("^[0-9.]+ ?% or more$", v); if (any(r)) out[r] <- (num("^([0-9.]+).*")[r] + 100) / 2
  r <- grepl("^[0-9.]+ ?%$", v); if (any(r)) out[r] <- num("^([0-9.]+).*")[r]
  out
}

# Parse "2015", "2015:2018", "2015,2017" -> integer years; "AZ,NM" -> states. NULL if empty.
nf_parse_years <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  if (!grepl("^[0-9]{4}(:[0-9]{4})?(,[0-9]{4}(:[0-9]{4})?)*$", x)) stop("years must look like 2015, 2015:2018 or 2015,2017: ", x)
  as.integer(unlist(lapply(strsplit(x, ",")[[1]], function(p) {
    r <- as.integer(strsplit(p, ":")[[1]]); if (length(r) == 2) r[1]:r[2] else r
  })))
}
nf_parse_states <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  toupper(trimws(strsplit(x, ",")[[1]]))
}
