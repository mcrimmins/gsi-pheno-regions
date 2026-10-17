#!/usr/bin/env Rscript
# Step 62 (fire activity study; P720): VIIRS 375 m active fire detections from NASA FIRMS for the
# lower 48 -> one Parquet table per sensor and year, on the 4 km grid, with region and cover.
#   1. Which dates each source holds: FIRMS data_availability (falls back to config dates).
#      Science-quality "SP" data are used where they exist; the near-real-time "NRT" feed only for
#      later dates (FIRMS replaces NRT with SP about 5 months after acquisition).
#   2. Download in pieces of window_days (the FIRMS area API now takes at most 5 days per request) over
#      the CONUS box (config fire_study: firms: bbox), one gzipped CSV per piece in
#      <data_root>/_shared/raw/fire/firms/<SOURCE>/<year>/. Pieces on disk are never fetched again,
#      so the step resumes where it stopped. Sequential, with a pause, and a check of the key's
#      transaction count (5,000 per 10 minutes) before each request.
#      The key comes ONLY from the FIRMS_MAP_KEY environment variable (~/.Renviron); it is never
#      logged (URLs in logs and errors are masked).
#   3. Per source-year: combine pieces; drop non-fire detection types where the data say so
#      (SP `type`: 0 vegetation fire kept; 1 volcano, 2 static land source such as gas flares and
#      industry, 3 offshore dropped); keep every confidence class with a `conf` column
#      (low / nominal / high); drop points off the grid's land cells. Attach cell, region, fuel
#      share, crop share, strict herb mask; acq date in UTC and the local solar date (UTC + lon/15 h),
#      since evening and night passes can fall on the next UTC day.
#   4. Gaps: days per month with no detections at all over CONUS (an orbit or ground-system outage,
#      not a fire-free day) are listed in firms_gaps.csv.
# Primary series VIIRS S-NPP (from 2012-01-20); NOAA-20 (from 2018-04-01) as the check series
# (config fire_study: firms: primary / check). Linking to fires and the cropland screen are 63.
# Outputs <run_dir>/fire_study/firms/:
#   <SENSOR>/firms_<SENSOR>_<year>.parquet   detections (lat, lon, date_utc, time_utc, date_solar,
#       daynight, conf, frp, bright_ti4, bright_ti5, scan, track, type, source (SP / NRT), cell,
#       region, fuel_share, crop_share, herb_mask, pilot)
#   firms_summary.csv  sensor x year: rows downloaded, by type, by confidence, kept, pilot regions
#   firms_gaps.csv     sensor x year x month: days with no CONUS detections
#   firms_requests.csv one row per downloaded piece: source, start, days, rows, bytes, seconds
#
# Usage: Rscript scripts/62_fire_firms.R [profile] [sensors] [years] [force]
#   sensors: VIIRS_SNPP, VIIRS_NOAA20 or both (default: primary and check from config)
#   years:   2012-2025 (default config fire_study: firms: years), or one year, e.g. 2020
#   force:   rebuild the per-year tables (downloads are still never repeated)
# Needs FIRMS_MAP_KEY in ~/.Renviron (free key: https://firms.modaps.eosdis.nasa.gov/api/map_key/).

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
force <- "force" %in% args
rest <- setdiff(args[-1], "force")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("62_fire_firms", profile)
notify_init(cfg)
fs <- cfg$fire_study; fc <- fs$firms
sensors <- c(fc$primary, fc$check)
if (length(rest) && grepl("VIIRS", rest[1])) { sensors <- strsplit(rest[1], ",")[[1]]; rest <- rest[-1] }
yrs <- seq(as.integer(fc$years$start), as.integer(fc$years$end))
if (length(rest)) { r <- as.integer(strsplit(rest[1], "-")[[1]]); yrs <- seq(r[1], r[length(r)]) }
out_dir <- fire_dir(cfg, "firms"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
key <- Sys.getenv("FIRMS_MAP_KEY", "")
mask_key <- function(x) if (nzchar(key)) gsub(key, "<MAP_KEY>", x, fixed = TRUE) else x
log_msg("=== 62_fire_firms | profile: ", profile, " | sensors ", paste(sensors, collapse = ", "), " | years ",
        min(yrs), "-", max(yrs), " | key ", if (nzchar(key)) "set" else "NOT SET")
bbox <- unlist(fc$bbox); area <- paste(bbox[c("xmin", "ymin", "xmax", "ymax")], collapse = ",")
api <- sub("/+$", "", fc$api)

# ---- helpers ---------------------------------------------------------------------------------------
http_text <- function(url, timeout = 300, tries = 1) {
  req <- httr2::request(url) |> httr2::req_user_agent(cfg$prism$user_agent %||% "gsi-pheno-regions (R httr2)") |>
    httr2::req_timeout(timeout) |> httr2::req_error(is_error = function(resp) FALSE)
  if (tries > 1) req <- httr2::req_retry(req, max_tries = tries, backoff = function(i) 30 * i)
  resp <- tryCatch(httr2::req_perform(req), error = function(e) stop(mask_key(conditionMessage(e)), call. = FALSE))
  list(status = httr2::resp_status(resp), text = httr2::resp_body_string(resp))
}
# transactions used in the current 10-minute window (fail-soft: NA if the status call fails)
key_status <- function() {
  r <- tryCatch(http_text(paste0(fc$status_url, "?MAP_KEY=", key), timeout = 30), error = function(e) NULL)
  if (is.null(r) || r$status != 200) return(c(used = NA, limit = NA))
  j <- tryCatch(jsonlite::fromJSON(r$text), error = function(e) NULL)
  c(used = as.numeric(j$current_transactions %||% NA), limit = as.numeric(j$transaction_limit %||% NA))
}
wait_for_quota <- function() {
  for (i in 1:15) {
    s <- key_status()
    if (is.na(s["used"]) || s["used"] < s["limit"] * 0.8) return(invisible(TRUE))
    log_msg(sprintf("  key at %s of %s transactions; waiting 60 s", s["used"], s["limit"]))
    Sys.sleep(60)
  }
}
# FIRMS date coverage per source id (e.g. VIIRS_SNPP_SP), from the API, else from config
availability <- function() {
  av <- NULL
  if (nzchar(key)) {
    r <- tryCatch(http_text(sprintf("%s/%s/ALL", sub("/area/csv$", "/data_availability/csv", api), key), timeout = 60),
                  error = function(e) { log_warn("data_availability: ", conditionMessage(e)); NULL })
    if (!is.null(r) && r$status == 200 && grepl("^data_id", r$text))
      av <- utils::read.csv(text = r$text, stringsAsFactors = FALSE)
  }
  if (is.null(av)) {
    log_warn("FIRMS data_availability not reachable: using config dates (fire_study: firms: fallback_dates)")
    fb <- fc$fallback_dates
    av <- data.frame(data_id = names(fb), min_date = vapply(fb, function(x) x$min, ""), max_date = vapply(fb, function(x) x$max, ""))
  }
  av$min_date <- as.Date(av$min_date); av$max_date <- as.Date(av$max_date)
  av
}
# request windows of config window_days (FIRMS max: 5 days as of Oct 2026) covering [from, to]
wdays <- as.integer(fc$window_days %||% 5)
windows <- function(from, to) {
  if (from > to) return(data.frame(start = as.Date(character()), days = integer()))
  s <- seq(from, to, by = wdays)
  data.frame(start = s, days = as.integer(pmin(wdays, as.numeric(to - s) + 1)))
}
piece_file <- function(src, start, days) file.path(fire_cache_dir(cfg, "firms", src, format(start, "%Y")),
                                                  sprintf("%s_%s_%02d.csv.gz", src, format(start, "%Y%m%d"), days))
REQ <- list()
fetch_piece <- function(src, start, days) {
  f <- piece_file(src, start, days)
  if (file.exists(f)) return(f)
  if (!nzchar(key)) stop("FIRMS_MAP_KEY is not set (add it to ~/.Renviron) and piece ", basename(f), " is not on disk")
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  wait_for_quota()
  url <- sprintf("%s/%s/%s/%s/%d/%s", api, key, src, area, days, format(start, "%Y-%m-%d"))
  t1 <- Sys.time()
  r <- http_text(url, timeout = as.numeric(fc$timeout_sec %||% 300), tries = 3)
  if (r$status != 200 || !grepl("^latitude", r$text)) {
    msg <- substr(gsub("\\s+", " ", r$text), 1, 200)
    if (grepl("exceed|limit", msg, ignore.case = TRUE)) {
      log_warn("  transaction limit reached; waiting 10 min"); Sys.sleep(600)
      return(fetch_piece(src, start, days))
    }
    stop(sprintf("FIRMS %s %s: HTTP %d: %s%s", src, format(start), r$status, mask_key(msg),
                 if (r$status == 400) " (request rejected: stopping the run; check fire_study: firms: window_days)" else ""))
  }
  tmp <- paste0(f, ".part")
  con <- gzfile(tmp, "w"); writeLines(r$text, con, sep = ""); close(con)
  file.rename(tmp, f)
  n <- max(0L, length(strsplit(r$text, "\n", fixed = TRUE)[[1]]) - 1L)
  REQ[[length(REQ) + 1]] <<- data.frame(source = src, start = start, days = days, rows = n, bytes = nchar(r$text),
                                        seconds = round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1))
  Sys.sleep(as.numeric(fc$sleep_sec %||% 1))
  f
}
read_piece <- function(f) {
  d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, colClasses = c(acq_time = "character")), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d
}
conf_class <- function(x) {
  x <- tolower(as.character(x))
  ifelse(x %in% c("l", "low"), "low", ifelse(x %in% c("n", "nominal"), "nominal", ifelse(x %in% c("h", "high"), "high", NA)))
}

# ---- grid layers ---------------------------------------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
gv <- terra::values(grid, mat = FALSE)
reg <- fire_region_values(cfg, grid, fs$partition %||% "A2_groups", fs$k %||% 13)
fuel <- fire_fuel_share(cfg)
hs <- terra::rast(file.path(static_dir(cfg), "herb_share.tif"))
crop <- terra::values(hs[["crop_share"]], mat = FALSE)
hmv <- terra::values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), mat = FALSE)
herb <- as.integer(!is.na(hmv) & hmv > 0)
pilot_regions <- as.integer(unlist(fs$pilot_regions))
log_msg("regions from ", reg$src, " | request box ", area)

av <- availability()
for (i in seq_len(nrow(av))) if (grepl("VIIRS_(SNPP|NOAA20)", av$data_id[i]))
  log_msg(sprintf("  %-18s %s to %s", av$data_id[i], format(av$min_date[i]), format(av$max_date[i])))

# request plan for one sensor-year (NULL if nothing to do)
plan_year <- function(sen, y) {
  sp <- av[av$data_id == paste0(sen, "_SP"), ]; nrt <- av[av$data_id == paste0(sen, "_NRT"), ]
  if (!nrow(sp)) return(NULL)
  y0 <- max(as.Date(sprintf("%d-01-01", y)), sp$min_date); y1 <- min(as.Date(sprintf("%d-12-31", y)), Sys.Date() - 1)
  if (y0 > y1) return(NULL)
  rbind(cbind(windows(y0, min(y1, sp$max_date)), src = paste0(sen, "_SP")),
        if (nrow(nrt) && y1 > sp$max_date) cbind(windows(max(y0, sp$max_date + 1), y1), src = paste0(sen, "_NRT")))
}
# whole-run plan: pieces still to download, so progress lines can show an overall count and ETA
todo_total <- 0L; pieces_total <- 0L
for (sen in sensors) for (y in yrs) {
  if (!force && !is.na(fire_table_file(file.path(out_dir, sen, sprintf("firms_%s_%d", sen, y))))) next
  W0 <- plan_year(sen, y); if (is.null(W0)) next
  pieces_total <- pieces_total + nrow(W0)
  todo_total <- todo_total + sum(!file.exists(mapply(piece_file, W0$src, W0$start, W0$days)))
}
log_msg(sprintf("plan: %s pieces, %s already on disk, %s to download (%d-day windows)", format(pieces_total, big.mark = ","),
                format(pieces_total - todo_total, big.mark = ","), format(todo_total, big.mark = ","), wdays))
done_total <- 0L; sec_total <- 0; every <- as.integer(fc$progress_every %||% 1)
eta <- function(left, per) if (!is.finite(per) || left <= 0) "-" else fmt_dur(left * per)

SUM <- list(); GAP <- list()
for (sen in sensors) {
  sp <- av[av$data_id == paste0(sen, "_SP"), ]; nrt <- av[av$data_id == paste0(sen, "_NRT"), ]
  if (!nrow(sp)) { log_err("no availability for ", sen, "_SP"); next }
  for (y in yrs) {
    of <- file.path(out_dir, sen, sprintf("firms_%s_%d", sen, y))
    if (!force && !is.na(fire_table_file(of))) { log_msg(sprintf("%s %d: table exists", sen, y)); next }
    y0 <- max(as.Date(sprintf("%d-01-01", y)), sp$min_date); y1 <- min(as.Date(sprintf("%d-12-31", y)), Sys.Date() - 1)
    if (y0 > y1) { log_msg(sprintf("%s %d: before the sensor record", sen, y)); next }
    W <- plan_year(sen, y)
    t1 <- Sys.time(); n_new <- 0
    files <- character(nrow(W))
    have <- file.exists(mapply(piece_file, W$src, W$start, W$days))
    log_msg(sprintf("%s %d: %d pieces, %d on disk, %d to download", sen, y, nrow(W), sum(have), sum(!have)))
    for (j in seq_len(nrow(W))) {
      had <- have[j]; tp <- Sys.time()
      files[j] <- tryCatch(fetch_piece(W$src[j], W$start[j], W$days[j]), error = function(e) {
        m <- conditionMessage(e)
        if (grepl("MAP_KEY|HTTP 400", m, ignore.case = TRUE)) stop(m, call. = FALSE)   # bad key or rejected request: stop the run
        log_err(m); NA })
      if (!had && !is.na(files[j])) {
        n_new <- n_new + 1; done_total <- done_total + 1L
        dt <- as.numeric(difftime(Sys.time(), tp, units = "secs")); sec_total <- sec_total + dt
        if (n_new %% every == 0 || n_new == sum(!have)) {
          rq <- REQ[[length(REQ)]]
          log_msg(sprintf("  %s %d [%d/%d] %s +%dd: %s rows, %.1f MB, %.0f s | year left %s | run %s/%s, left %s",
                          sub("VIIRS_", "", sen), y, n_new, sum(!have), format(W$start[j], "%m-%d"), W$days[j],
                          format(rq$rows, big.mark = ","), rq$bytes / 2^20, dt, eta(sum(!have) - n_new, sec_total / done_total),
                          format(done_total, big.mark = ","), format(todo_total, big.mark = ","), eta(todo_total - done_total, sec_total / done_total)))
        }
      }
    }
    if (any(is.na(files))) { log_err(sprintf("%s %d: %d pieces failed; year left incomplete (rerun to resume)", sen, y, sum(is.na(files)))); next }
    log_msg(sprintf("%s %d: %d pieces (%d downloaded now, %s)", sen, y, nrow(W), n_new, fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))))
    D <- lapply(seq_along(files), function(j) { d <- read_piece(files[j]); if (!is.null(d)) d$source <- sub(".*_", "", W$src[j]); d })
    D <- D[!vapply(D, is.null, TRUE)]
    if (!length(D)) { log_warn(sprintf("%s %d: no detections at all", sen, y)); next }
    cols <- unique(unlist(lapply(D, names)))
    D <- do.call(rbind, lapply(D, function(d) { for (cn in setdiff(cols, names(d))) d[[cn]] <- NA; d[, cols] }))
    n_raw <- nrow(D)
    D$date_utc <- as.Date(D$acq_date)
    D <- D[!is.na(D$date_utc) & format(D$date_utc, "%Y") == as.character(y), ]
    typ <- if ("type" %in% names(D)) suppressWarnings(as.integer(D$type)) else rep(NA_integer_, nrow(D))
    by_type <- table(factor(ifelse(is.na(typ), "na", typ), levels = c("0", "1", "2", "3", "na")))
    D$conf <- conf_class(D$confidence)
    keep_t <- is.na(typ) | typ %in% as.integer(unlist(fc$keep_types %||% list(0)))
    cell <- terra::cellFromXY(grid, cbind(D$longitude, D$latitude))
    on <- !is.na(cell) & !is.na(gv[cell])
    D <- D[keep_t & on, ]; cell <- cell[keep_t & on]
    hhmm <- sprintf("%04d", as.integer(D$acq_time))
    utc <- as.POSIXct(paste(D$date_utc, hhmm), format = "%Y-%m-%d %H%M", tz = "UTC")
    X <- data.frame(
      lat = D$latitude, lon = D$longitude, date_utc = D$date_utc, time_utc = hhmm,
      date_solar = as.Date(utc + D$longitude / 15 * 3600), daynight = D$daynight, conf = D$conf,
      frp = as.numeric(D$frp), bright_ti4 = as.numeric(D$bright_ti4), bright_ti5 = as.numeric(D$bright_ti5),
      scan = as.numeric(D$scan), track = as.numeric(D$track), type = if ("type" %in% names(D)) as.integer(D$type) else NA_integer_,
      source = D$source, cell = as.integer(cell), region = as.integer(reg$v[cell]), fuel_share = fuel[cell],
      crop_share = crop[cell], herb_mask = herb[cell], stringsAsFactors = FALSE)
    X$pilot <- X$region %in% pilot_regions
    fo <- fire_write_table(X, of)
    days_seen <- unique(as.Date(substr(paste(D$acq_date), 1, 10)))
    months <- seq(as.Date(sprintf("%d-01-01", y)), y1, by = "month")
    for (m in seq_along(months)) {
      dm <- seq(max(months[m], y0), min(seq(months[m], by = "month", length.out = 2)[2] - 1, y1), by = 1)
      miss <- sum(!dm %in% days_seen)
      if (miss) GAP[[length(GAP) + 1]] <- data.frame(sensor = sen, year = y, month = as.integer(format(months[m], "%m")), days_missing = miss,
                                                   first_missing = format(min(dm[!dm %in% days_seen])))
    }
    cf <- table(factor(X$conf, levels = c("low", "nominal", "high")))
    SUM[[length(SUM) + 1]] <- data.frame(sensor = sen, year = y, downloaded = n_raw, type0 = by_type[["0"]], type1 = by_type[["1"]],
      type2 = by_type[["2"]], type3 = by_type[["3"]], type_na = by_type[["na"]], kept = nrow(X),
      low = cf[["low"]], nominal = cf[["nominal"]], high = cf[["high"]], pilot = sum(X$pilot),
      nrt_rows = sum(X$source == "NRT"), days_with_data = length(days_seen))
    s <- SUM[[length(SUM)]]
    log_msg(sprintf("  %s %d: %s downloaded | static/volcano/offshore %s | kept on grid %s (low %.0f %%) | pilot %s | days with data %d%s",
                    sen, y, format(n_raw, big.mark = ","), format(s$type1 + s$type2 + s$type3, big.mark = ","), format(s$kept, big.mark = ","),
                    100 * s$low / max(1, s$kept), format(s$pilot, big.mark = ","), s$days_with_data,
                    if (s$nrt_rows) sprintf(" | NRT rows %s", format(s$nrt_rows, big.mark = ",")) else ""))
    notify(sprintf("62 FIRMS %s %d", sen, y), sprintf("%s detections kept (%s pilot)", format(s$kept, big.mark = ","), format(s$pilot, big.mark = ",")), priority = 2)
  }
}

# ---- logs ------------------------------------------------------------------------------------------
app <- function(df, f) { if (is.null(df) || !nrow(df)) return(invisible())
  old <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) { key_cols <- intersect(c("sensor", "year", "month", "source", "start"), names(df))
    old <- old[!do.call(paste, old[key_cols]) %in% do.call(paste, lapply(df[key_cols], as.character)), , drop = FALSE]
    df <- rbind(old[, names(df)], df) }
  utils::write.csv(df, f, row.names = FALSE) }
app(if (length(SUM)) do.call(rbind, SUM), file.path(out_dir, "firms_summary.csv"))
app(if (length(GAP)) do.call(rbind, GAP), file.path(out_dir, "firms_gaps.csv"))
app(if (length(REQ)) do.call(rbind, REQ), file.path(out_dir, "firms_requests.csv"))
if (length(GAP)) {
  G <- do.call(rbind, GAP); G <- G[G$days_missing >= 3, ]
  if (nrow(G)) log_msg("months with >= 3 days without any CONUS detection: ",
                       paste(sprintf("%s %d-%02d (%d)", sub("VIIRS_", "", G$sensor), G$year, G$month, G$days_missing), collapse = ", "))
}
msg <- sprintf("%d sensor-years written; %s", length(SUM), fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("62 FIRMS done", msg, priority = 3, tags = "white_check_mark")
