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
