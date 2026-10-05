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

  # Ctrl-C during a transfer reaches us as a curl error ("aborted by an application
  # callback"), not an interrupt; re-raise it so the run stops instead of retrying.
  resp <- tryCatch(httr2::req_perform(req, path = zip),
                   interrupt = function(e) stop("interrupted by user", call. = FALSE),
                   error = function(e) e)
  if (inherits(resp, "error")) {
    msg <- conditionMessage(resp)
    if (grepl("aborted by an application callback|interrupt", msg, ignore.case = TRUE)) {
      stop("interrupted by user", call. = FALSE)
    }
    return(list(status = "error", msg = msg))
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

# Write a daily stack to its final GeoTIFF: temp file, full QA against the source, then
# rename into place. `scale` NULL -> FLT4S; a number (e.g. 0.01) -> INT2S with that GDAL
# scale factor (terra applies it on read, so values come back in physical units).
# QA: every band, max |written - source| <= scale/2 (0 for float) and identical NA cells.
prism_write_stack <- function(r, out_file, dates, scale = NULL) {
  names(r) <- format(dates, "%Y-%m-%d")
  terra::time(r) <- dates
  src <- r
  if (is.null(scale)) {
    dt <- "FLT4S"; pred <- "PREDICTOR=3"; wargs <- list()
  } else {
    # terra truncates toward zero when scaling to integers: round to the storage step,
    # then nudge away from zero by a tenth of a step so truncation hits the right integer.
    # One app() pass: chained raster arithmetic would write a full temp copy per operation
    # (~1.3 GB each for a CONUS year).
    r <- terra::app(r, function(v) {
      x <- round(v / scale) * scale
      x + sign(x) * (scale / 10)
    })
    names(r) <- format(dates, "%Y-%m-%d")
    dt <- "INT2S"; pred <- "PREDICTOR=2"; wargs <- list(scale = scale, offset = 0)
  }
  tmp_dir <- file.path(dirname(out_file), ".tmp")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  tmp <- file.path(tmp_dir, basename(out_file))
  do.call(terra::writeRaster, c(list(
    x = r, filename = tmp, overwrite = TRUE, datatype = dt,
    gdal = c("COMPRESS=DEFLATE", pred, "ZLEVEL=6", "TILED=YES", "INTERLEAVE=BAND",
             "BIGTIFF=IF_SAFER")), wargs))

  chk <- terra::rast(tmp)
  stopifnot(terra::nlyr(chk) == length(dates))
  # QA in one pass over written + source bands: per cell, max |difference| and the number of
  # days where only one of them is NA. Output is 2 layers, so no large temp files.
  k <- length(dates)
  qa <- terra::app(c(chk, src), function(v) {
    a <- v[seq_len(k)]; b <- v[k + seq_len(k)]
    d <- abs(a - b)
    c(if (all(is.na(d))) 0 else max(d, na.rm = TRUE), sum(is.na(a) != is.na(b)))
  })
  max_err <- terra::global(qa[[1]], "max", na.rm = TRUE)[1, 1]
  na_mismatch <- terra::global(qa[[2]], "sum", na.rm = TRUE)[1, 1]
  # Rounding error is at most scale/2, but source values are float32: near 100 (e.g. VPD in
  # hPa) one float32 step is ~8e-6, so an exact scale/2 bound fails on legitimate values.
  # A 0.2 % margin covers that and still catches real faults (truncation errors reach a
  # full step, overflow gives NA or huge errors).
  tol <- if (is.null(scale)) 1e-4 else scale * 0.501
  if (!is.finite(max_err) || max_err > tol || na_mismatch != 0) {
    stop(sprintf("QA failed: max error %.7f (tol %.7f), %.0f cell-days NA in only one of written/source",
                 max_err, tol, na_mismatch))
  }
  info <- list(nrow = terra::nrow(chk), ncol = terra::ncol(chk),
               ext = as.vector(terra::ext(chk)),
               crs = terra::crs(chk, describe = TRUE)$name, res = terra::res(chk),
               datatype = terra::datatype(chk)[1], max_err = max_err,
               range = as.vector(terra::global(chk[[1]], "range", na.rm = TRUE)[1, ]))
  rm(chk, r, src, qa); gc()

  # Replace any existing output (and its stale .aux.xml) with the checked temp file.
  for (suffix in c("", ".aux.xml")) unlink(paste0(out_file, suffix))
  for (suffix in c("", ".aux.xml")) {
    if (file.exists(paste0(tmp, suffix)) &&
        !file.rename(paste0(tmp, suffix), paste0(out_file, suffix))) {
      stop("could not move ", basename(tmp), suffix, " into place")
    }
  }
  unlink(tmp_dir, recursive = TRUE)
  c(list(ok = TRUE, out_file = out_file, nlyr = length(dates),
         mb = file.size(out_file) / 1e6), info)
}

# Stack daily tifs into one cropped, compressed multi-band GeoTIFF (one band per day,
# band names = ISO dates). Self-contained so it can run in a future worker: it gets
# file paths and plain values only, never terra objects.
prism_convert_year <- function(tifs, dates, out_file, bbox, memfrac, scale = NULL) {
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
    prism_write_stack(r, out_file, dates, scale)
  }, error = function(e) list(ok = FALSE, out_file = out_file, msg = conditionMessage(e)),
  finally = prism_clean_tmp())
}

# Remove this R session's terra temp files. Workers live for the whole run, so without
# this their temp files pile up year after year. current = TRUE only: other workers'
# files (same parent temp dir) are left alone.
prism_clean_tmp <- function() {
  gc()
  try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE),
      silent = TRUE)
  invisible(NULL)
}

# Re-encode an existing variable-year file to the configured storage (no download).
prism_recode_file <- function(f, scale, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    r <- terra::rast(f)
    dates <- as.Date(names(r))
    if (anyNA(dates)) stop("band names are not ISO dates")
    # Read into the temp-file path from a copy so the source can be replaced afterwards.
    src_copy <- file.path(tempdir(), paste0("src_", basename(f)))
    file.copy(f, src_copy, overwrite = TRUE)
    if (file.exists(paste0(f, ".aux.xml"))) file.copy(paste0(f, ".aux.xml"),
                                                    paste0(src_copy, ".aux.xml"), overwrite = TRUE)
    rm(r); gc()
    r <- terra::rast(src_copy)
    res <- prism_write_stack(r, f, dates, scale)
    rm(r); gc()
    unlink(c(src_copy, paste0(src_copy, ".aux.xml")))
    res
  }, error = function(e) list(ok = FALSE, out_file = f, msg = conditionMessage(e)),
  finally = prism_clean_tmp())
}

# Storage scale for a variable from config (NULL = float32).
prism_scale <- function(cfg, var) {
  s <- cfg$prism$scale[[var]]
  if (is.null(s)) NULL else as.numeric(s)
}
