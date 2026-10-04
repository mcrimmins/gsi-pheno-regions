# PRISM daily download and conversion.
#
# Source: PRISM Group web service (services.nacse.org), documented in
#   https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf (26 Mar 2025)
#   URL: <base>/<region>/<res>/<element>/<YYYYMMDD>  -> zip containing a COG .tif
# Limits: a given file may be downloaded at most twice per 24 h (then refused for that
# period); excessive activity can get the IP blocked. So: sequential requests, a pause
# after each, never re-download a day already on disk, at most one retry per day.
# Data older than ~6 months are stable until a new time-series version is released.
# Citation: PRISM Group, Oregon State University, https://prism.oregonstate.edu,
#   accessed <date>.

prism_url <- function(cfg, var, date) {
  p <- cfg$prism
  sprintf("%s/%s/%s/%s/%s", p$base_url, p$region, p$res, var, format(date, "%Y%m%d"))
}

# Dates to fetch for a year: complete days only (never today or the future).
prism_days <- function(year, max_days = NULL) {
  d <- seq(as.Date(sprintf("%d-01-01", year)), as.Date(sprintf("%d-12-31", year)),
           by = "day")
  d <- d[d <= Sys.Date() - 1]
  if (!is.null(max_days)) d <- utils::head(d, max_days)
  d
}

prism_day_file <- function(raw_dir, var, date) {
  file.path(raw_dir, sprintf("%s_%s.tif", var, format(date, "%Y%m%d")))
}

# Download one day; keep only the .tif. Returns a list with status:
#   "downloaded" | "error" | "refused" (server answered but not with a zip: don't retry)
prism_download_day <- function(cfg, var, date, raw_dir) {
  out_tif <- prism_day_file(raw_dir, var, date)
  tmp_dir <- file.path(raw_dir, ".tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  zip <- tempfile(pattern = paste0(var, "_"), tmpdir = tmp_dir, fileext = ".zip")
  on.exit(unlink(zip), add = TRUE)

  req <- httr2::request(prism_url(cfg, var, date)) |>
    httr2::req_user_agent(cfg$prism$user_agent) |>
    httr2::req_timeout(cfg$prism$timeout_sec) |>
    httr2::req_error(is_error = function(resp) FALSE)

  resp <- tryCatch(httr2::req_perform(req, path = zip), error = function(e) e)
  if (inherits(resp, "error")) {
    return(list(status = "error", msg = conditionMessage(resp)))
  }
  code <- httr2::resp_status(resp)
  if (code != 200) {
    return(list(status = "error", msg = sprintf("HTTP %d", code)))
  }

  files <- tryCatch(utils::unzip(zip, list = TRUE)$Name, error = function(e) NULL)
  if (is.null(files)) {
    txt <- tryCatch(paste(readLines(zip, n = 5, warn = FALSE), collapse = " "),
                    error = function(e) "<unreadable>")
    return(list(status = "refused", msg = paste("response is not a zip:", substr(txt, 1, 300))))
  }
  tif <- grep("\\.tif$", files, value = TRUE)
  if (length(tif) != 1) {
    return(list(status = "refused",
                msg = paste("expected one .tif in zip, found:", paste(files, collapse = ", "))))
  }

  utils::unzip(zip, files = tif, exdir = tmp_dir, junkpaths = TRUE)
  extracted <- file.path(tmp_dir, basename(tif))
  if (!file.rename(extracted, out_tif)) {
    unlink(extracted)
    return(list(status = "error", msg = "could not move extracted tif into place"))
  }
  ok <- tryCatch(terra::nlyr(terra::rast(out_tif)) == 1, error = function(e) FALSE)
  if (!ok) {
    unlink(out_tif)
    return(list(status = "error", msg = "extracted tif failed to open"))
  }
  list(status = "downloaded", src = basename(tif), bytes = file.size(zip))
}

# Download all missing days for one variable-year, sequentially and politely.
# Returns the dates still missing afterwards.
prism_download_var_year <- function(cfg, var, year, dates) {
  raw_dir <- prism_raw_dir(cfg, var, year)
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  have <- file.exists(prism_day_file(raw_dir, var, dates))
  todo <- dates[!have]
  log_msg(sprintf("%s %d: %d days, %d on disk, %d to download",
                  var, year, length(dates), sum(have), length(todo)))
  if (!length(todo)) return(as.Date(character()))

  t0 <- Sys.time()
  n_ok <- 0L; bytes <- 0
  for (i in seq_along(todo)) {
    d <- todo[i]
    res <- prism_download_day(cfg, var, d, raw_dir)
    Sys.sleep(cfg$prism$sleep_sec)
    if (res$status == "error") {
      log_warn(sprintf("%s %s: %s; retrying once in %d s", var, d, res$msg,
                       cfg$prism$retry_wait_sec))
      Sys.sleep(cfg$prism$retry_wait_sec)
      res <- prism_download_day(cfg, var, d, raw_dir)
      Sys.sleep(cfg$prism$sleep_sec)
    }
    if (res$status == "downloaded") {
      n_ok <- n_ok + 1L; bytes <- bytes + res$bytes
      if (n_ok == 1L) log_msg(sprintf("%s %d: first file from server: %s", var, year, res$src))
    } else {
      log_err(sprintf("%s %s: %s: %s", var, d, res$status, res$msg))
      if (res$status == "refused" && grepl("limit|block|exceed", res$msg, ignore.case = TRUE)) {
        log_err("Server reports a download limit; stopping this variable-year. ",
                "Wait 24 h before rerunning these days.")
        notify(sprintf("PRISM download limit: %s %d", var, year),
               paste0("Server refused ", var, " ", d, ": ", res$msg,
                      "\nStopped this variable-year. Wait 24 h before rerunning."),
               priority = 5, tags = "rotating_light")
        break
      }
    }
    if (i %% 50 == 0 || i == length(todo)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      log_msg(sprintf("%s %d: %d/%d requested, %d ok, %.1f MB, %.1f s/file",
                      var, year, i, length(todo), n_ok, bytes / 1e6, el / i))
    }
  }
  unlink(file.path(raw_dir, ".tmp"), recursive = TRUE)
  missing <- dates[!file.exists(prism_day_file(raw_dir, var, dates))]
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (length(missing)) {
    fail_file <- file.path(prism_paths(cfg)$raw_root, sprintf("_missing_%s_%d.txt", var, year))
    writeLines(format(missing), fail_file)
    log_warn(sprintf("%s %d: %d days still missing (listed in %s)", var, year,
                     length(missing), fail_file))
    notify_error(sprintf("%s %d: %d days missing", var, year, length(missing)),
                 sprintf("%d of %d days downloaded in %s. Year not converted; rerun later.\nMissing list: %s",
                         length(dates) - length(missing), length(dates), fmt_dur(el), fail_file))
  } else if (notify_downloads_on()) {
    notify(sprintf("%s %d downloaded", var, year),
           sprintf("%d new files, %.0f MB, %s (%.1f s/file)", n_ok, bytes / 1e6, fmt_dur(el),
                   el / length(todo)),
           priority = 2, tags = "arrow_down")
  }
  missing
}

# Stack daily tifs into one cropped, compressed multi-band GeoTIFF (one band per day,
# band names = ISO dates). Self-contained so it can run in a future worker: it gets
# file paths and plain values only, never terra objects.
prism_convert_year <- function(tifs, dates, out_file, bbox, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    r <- terra::rast(tifs)
    if (!is.null(bbox)) {
      e <- terra::ext(bbox[["xmin"]], bbox[["xmax"]], bbox[["ymin"]], bbox[["ymax"]])
      if (!terra::is.lonlat(r)) {
        e <- terra::ext(terra::project(terra::as.polygons(e, crs = "EPSG:4269"),
                                       terra::crs(r)))
      }
      r <- terra::crop(r, e, snap = "out")
    }
    names(r) <- format(dates, "%Y-%m-%d")
    terra::time(r) <- dates

    tmp_dir <- file.path(dirname(out_file), ".tmp")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    tmp <- file.path(tmp_dir, basename(out_file))
    terra::writeRaster(r, tmp, overwrite = TRUE, datatype = "FLT4S",
                       gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=6",
                                "TILED=YES", "INTERLEAVE=BAND", "BIGTIFF=IF_SAFER"))
    chk <- terra::rast(tmp)
    stopifnot(terra::nlyr(chk) == length(dates))
    info <- list(nrow = terra::nrow(chk), ncol = terra::ncol(chk),
                 ext = as.vector(terra::ext(chk)), crs = terra::crs(chk, describe = TRUE)$name,
                 res = terra::res(chk),
                 range = as.vector(terra::global(chk[[1]], "range", na.rm = TRUE)[1, ]))
    rm(chk, r); gc()

    for (suffix in c("", ".aux.xml")) {
      if (file.exists(paste0(tmp, suffix))) {
        if (!file.rename(paste0(tmp, suffix), paste0(out_file, suffix))) {
          stop("could not move ", basename(tmp), suffix, " into place")
        }
      }
    }
    unlink(tmp_dir, recursive = TRUE)
    c(list(ok = TRUE, out_file = out_file, nlyr = length(dates),
           mb = file.size(out_file) / 1e6), info)
  }, error = function(e) list(ok = FALSE, out_file = out_file, msg = conditionMessage(e)))
}
