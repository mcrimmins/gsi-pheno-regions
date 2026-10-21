#!/usr/bin/env Rscript
# Step 59a (NPN grass curing study; laptop or P720 with internet): download the LANDFIRE
# 40 Scott & Burgan fire behavior fuel models (FBFM40) for CONUS and unzip the GeoTIFF where
# 59b looks for it (config npn_study: landfire_fbfm40, default <data_root>/landfire).
# The zip is large (LF2023: 3.27 GB). Re-running skips the download when the zip is already
# complete (size checked against the server), and skips unzipping when the GeoTIFF is there.
# Only the raster files (.tif, .tfw, .ovr, .aux.xml, .vat.dbf) are extracted; add "keep_zip" as a
# second argument to keep the zip afterwards (default: deleted once the GeoTIFF is extracted).
#
# Usage: Rscript npn_study/scripts/59a_npn_landfire_download.R [profile] [keep_zip]
#   then: Rscript npn_study/scripts/59b_npn_fuelmodels.R full

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
keep_zip <- "keep_zip" %in% args
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("59a_npn_landfire_download", profile)
ns <- cfg$npn_study %||% list()
url <- ns$landfire_fbfm40_url %||% "https://landfire.gov/data-downloads/CONUS_LF2023/LF2023_FBFM40_CONUS.zip"
os <- if (.Platform$OS.type == "windows") "windows" else "linux"
dest <- ns$landfire_fbfm40; if (is.list(dest)) dest <- dest[[os]]
if (is.null(dest) || !nzchar(dest)) dest <- file.path(cfg$data_root, "landfire")
dest <- path.expand(dest)
if (grepl("\\.tif$", dest, ignore.case = TRUE)) dest <- dirname(dest)   # config names the file: use its folder
dir.create(dest, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 59a_npn_landfire_download | ", url, " -> ", dest)

have_tif <- function() { f <- list.files(dest, "\\.tif$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f[grepl("f40|fbfm40", basename(f), ignore.case = TRUE) | grepl("fbfm40", dirname(f), ignore.case = TRUE)] }
if (length(have_tif())) {
  log_msg("GeoTIFF already there: ", have_tif()[1], " (delete it to download again)")
  log_msg("=== done | log: ", log_file); quit(save = "no")
}

zf <- file.path(dest, basename(url))
remote_size <- tryCatch({
  h <- curlGetHeaders(url); cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  if (length(cl)) as.numeric(trimws(sub("^[^:]*:", "", cl[length(cl)]))) else NA_real_
}, error = function(e) NA_real_)
log_msg(sprintf("server size: %s", if (is.na(remote_size)) "unknown" else sprintf("%.2f GB", remote_size / 1e9)))

complete <- function() file.exists(zf) && (is.na(remote_size) || file.size(zf) == remote_size)
if (!complete()) {
  if (file.exists(zf)) { log_warn(sprintf("partial zip (%.2f GB): downloading again", file.size(zf) / 1e9)); unlink(zf) }
  options(timeout = max(getOption("timeout"), 4 * 3600))               # 4 h for a slow connection
  t0 <- Sys.time()
  st <- tryCatch(utils::download.file(url, zf, mode = "wb", method = "libcurl", quiet = FALSE), error = function(e) { log_err(conditionMessage(e)); 1 })
  if (st != 0 || !complete()) stop("download failed or incomplete (", if (file.exists(zf)) sprintf("%.2f GB", file.size(zf) / 1e9) else "no file",
                                   "); run again, or download ", url, " by hand into ", dest)
  log_msg(sprintf("downloaded %.2f GB in %.0f min", file.size(zf) / 1e9, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
} else log_msg("zip already complete: ", zf)

# extract only the raster files
L <- utils::unzip(zf, list = TRUE)
keep <- L$Name[grepl("\\.(tif|tfw|ovr|tif\\.aux\\.xml|vat\\.dbf|tif\\.vat\\.dbf|tif\\.vat\\.cpg)$", L$Name, ignore.case = TRUE)]
if (!length(keep)) stop("no GeoTIFF in ", zf, "; contents: ", paste(utils::head(L$Name, 20), collapse = ", "))
log_msg(sprintf("extracting %d of %d files (%.2f GB)", length(keep), nrow(L), sum(L$Length[L$Name %in% keep]) / 1e9))
utils::unzip(zf, files = keep, exdir = dest)
tif <- have_tif()
if (!length(tif)) stop("extracted, but no FBFM40 .tif found under ", dest)
log_msg("GeoTIFF: ", tif[1])
if (!keep_zip) { unlink(zf); log_msg("zip deleted (add keep_zip to keep it)") }
log_msg("=== done | next: Rscript npn_study/scripts/59b_npn_fuelmodels.R ", profile, " | log: ", log_file)
