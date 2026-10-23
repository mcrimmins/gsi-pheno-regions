#!/usr/bin/env Rscript
# Step 05 (fire study monitoring test, Oct 20 2026): daily VIIRS surface reflectance (VNP09GA, S-NPP,
# version 2) for the herbaceous phenology cells of the test regions, via AppEEARS area requests, cut to
# the 4 km grid as it arrives. The delay diagnostic (65b GSI_RT_LAG) showed the fire-growth signal comes
# back only with about a week or less of delay; this pulls the data a near-daily monitor would use.
#
# One AppEEARS task = one region x year x chunk of the year (config fire_study: viirs_daily:
# chunks_per_year, default 4), bbox = the region's phenology cells (as 65 / 65b) plus 0.05 deg; layers
# I1 (red), I2 (NIR) at 500 m and the quality layers QF1, QF2, QF6, QF7 and sensor zenith at 1 km. Up
# to max_pending_tasks are queued; finished bundles are downloaded (size + SHA-256 checked), each day
# is reduced to the region's cells (forked over days, workers_by_step viirs_daily), and the raw files
# are deleted. Clear-land rules: R/viirs_daily.R vd_clear() (config clear:).
#
# Output: <run_dir>/viirs_daily/VNP09GA/r<region>/daily_<Y>_c<k>.rds: list(region, year, chunk, cells
#   (full-grid cell numbers), dates, ndvi (cells x days: mean NDVI of clear 500 m pixels), n_clear,
#   n_land (observed land pixels), layers, clear rules). State: viirs_daily/tasks.csv. Layer metadata:
#   viirs_daily/VNP09GA/layers.json.
# Needs rules/typical_dates.tif (40), static/grid_mask.tif, static/nlcd_shares.tif.
# Credentials: Earthdata login in a netrc file (R/appeears.R). Never in code.
#
# Smoke test: GSI_VD_TEST=1 -> first region only, 10 days (config test_dates), outputs under
#   viirs_daily_test/, raw files kept for inspection, a per-day summary in the log.
# Usage: Rscript scripts/05_viirs_daily.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "appeears.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "viirs_daily.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
test <- nzchar(Sys.getenv("GSI_VD_TEST"))
log_file <- log_init(if (test) "05_viirs_daily_test" else "05_viirs_daily", profile)
notify_init(cfg)
terra::terraOptions(progress = 0)
vd <- vd_conf(cfg); product <- vd$product
sw <- step_workers(cfg, "viirs_daily")
pending_max <- as.integer(vd$max_pending_tasks %||% cfg$cutb$max_pending_tasks %||% 4)
poll_min <- as.numeric(vd$poll_minutes %||% cfg$cutb$poll_minutes %||% 5)
n_chunk <- as.integer(vd$chunks_per_year %||% 4)
dir.create(vd_dir(cfg), showWarnings = FALSE, recursive = TRUE)
state_file <- vd_dir(cfg, "tasks.csv")
grid <- terra::rast(file.path(cfg$run_dir, "static", "grid_mask.tif"))

regs <- vd$regions; yrs <- seq(as.integer(vd$years$start), as.integer(vd$years$end))
if (test) {
  td <- as.Date(unlist(vd$test_dates %||% list("2023-06-01", "2023-06-10")))
  regs <- regs[1]; yrs <- as.integer(format(td[1], "%Y"))
}
log_msg("=== 05_viirs_daily | profile: ", profile, if (test) " | TEST" else "", " | ", product, " | regions ", paste(regs, collapse = ","),
        " | years ", min(yrs), "-", max(yrs), " | chunks/yr ", if (test) 1 else n_chunk, " | workers ", sw$n)

# cells and bbox per region
CELLS <- lapply(stats::setNames(regs, regs), function(r) vd_cells(cfg, r))
BBOX <- lapply(CELLS, function(cc) vd_bbox(grid, cc))
for (r in regs) log_msg(sprintf("region %d: %s cells, bbox %s", r, format(length(CELLS[[as.character(r)]]), big.mark = ","),
                                paste(sprintf("%.2f", BBOX[[as.character(r)]]), collapse = " ")))

# layers (resolved against AppEEARS, metadata saved)
pats <- unlist(vd$layers %||% list(i1 = "SurfReflect_I1", i2 = "SurfReflect_I2", qf1 = "SurfReflect_QF1",
                                   qf2 = "SurfReflect_QF2", qf6 = "SurfReflect_QF6", qf7 = "SurfReflect_QF7",
                                   vza = "SensorZenith"))
res <- vd_resolve_layers(cfg, product, pats)
LAY <- res$layers; META <- stats::setNames(res$meta, names(LAY))
dir.create(vd_dir(cfg, vd_short(product)), showWarnings = FALSE, recursive = TRUE)
jsonlite::write_json(lapply(META, function(x) x[intersect(names(x), c("Description", "DataType", "FillValue", "ScaleFactor", "Units", "ValidMin", "ValidMax"))]),
                     vd_dir(cfg, vd_short(product), "layers.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
log_msg("layers: ", paste(names(LAY), LAY, sep = "=", collapse = " | "))
for (nm in c("i1", "i2")) log_msg(sprintf("  %s scale %s fill %s type %s", nm, META[[nm]]$ScaleFactor %||% "-", META[[nm]]$FillValue %||% "-", META[[nm]]$DataType %||% "-"))

# jobs: region x year x chunk (date ranges)
chunk_dates <- function(Y) {
  if (test) return(list(td))
  b <- seq(as.Date(sprintf("%d-01-01", Y)), as.Date(sprintf("%d-01-01", Y + 1)), length.out = n_chunk + 1)
  lapply(seq_len(n_chunk), function(k) c(b[k], if (k < n_chunk) b[k + 1] - 1 else as.Date(sprintf("%d-12-31", Y))))
}
jobs <- do.call(rbind, lapply(regs, function(r) do.call(rbind, lapply(yrs, function(Y) {
  cd <- chunk_dates(Y)
  data.frame(region = r, year = Y, chunk = seq_along(cd), start = as.Date(vapply(cd, function(x) format(x[1]), "")),
             end = as.Date(vapply(cd, function(x) format(x[2]), "")))
}))))
jobs$key <- sprintf("r%02d_%d_c%d", jobs$region, jobs$year, jobs$chunk)
out_of <- function(i) sub("\\.rds$", sprintf("_c%d.rds", jobs$chunk[i]), vd_out_file(cfg, product, jobs$region[i], jobs$year[i]))
raw_of <- function(i) file.path(vd_raw_dir(cfg, product, jobs$region[i], jobs$year[i]), sprintf("c%d", jobs$chunk[i]))
prev <- if (file.exists(state_file)) utils::read.csv(state_file, stringsAsFactors = FALSE, colClasses = c(task_id = "character")) else
  data.frame(key = character(), task_id = character(), status = character(), attempts = integer())
st <- merge(jobs, prev[, c("key", "task_id", "status", "attempts")], by = "key", all.x = TRUE, sort = FALSE)
st <- st[order(st$year, st$region, st$chunk), ]
st$status[is.na(st$status)] <- "new"; st$attempts[is.na(st$attempts)] <- 0L
jobs <- st[, c("region", "year", "chunk", "start", "end", "key")]
for (i in seq_len(nrow(st))) if (file.exists(out_of(i))) st$status[i] <- "complete"
st$status[st$status == "stack_failed"] <- "done"
save_state <- function() { st$updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  utils::write.csv(st[, c("key", "region", "year", "chunk", "task_id", "status", "attempts", "updated")], state_file, row.names = FALSE) }
save_state()
log_msg(sum(st$status == "complete"), " of ", nrow(st), " region-year-chunks already complete")

rules <- vd$clear
# Reduce one downloaded bundle to the region's cells.
reduce_job <- function(i) {
  r <- st$region[i]; cells <- CELLS[[as.character(r)]]
  gridc <- terra::crop(grid, terra::ext(BBOX[[as.character(r)]][c("xmin", "xmax", "ymin", "ymax")]), snap = "out")
  cidx <- terra::cellFromXY(gridc, terra::xyFromCell(grid, cells))
  fb <- vd_files_by_day(raw_of(i), LAY)
  if (!length(fb$days)) stop("no complete days (i1, i2, qf1, qf2) in ", raw_of(i))
  gf <- terra::sources(grid); ec <- as.vector(terra::ext(gridc))
  one <- function(d) tryCatch({
    terra::terraOptions(memfrac = sw$memfrac, progress = 0)
    vd_day(fb$files[[d]], META, terra::crop(terra::rast(gf), terra::ext(ec)), cidx, rules)
  }, error = function(e) paste("ERROR", d, conditionMessage(e)))
  out <- if (.Platform$OS.type == "unix" && sw$n > 1) parallel::mclapply(fb$days, one, mc.cores = sw$n, mc.preschedule = FALSE) else lapply(fb$days, one)
  bad <- vapply(out, is.character, TRUE)
  if (any(bad)) stop(paste(unlist(out[bad])[1:min(3, sum(bad))], collapse = "; "))
  mk <- function(v, type) { m <- do.call(cbind, lapply(out, `[[`, v)); storage.mode(m) <- type; m }
  O <- list(region = r, year = st$year[i], chunk = st$chunk[i], cells = cells, dates = vd_doy_date(fb$days),
            ndvi = mk("ndvi", "double"), n_clear = mk("n_clear", "integer"), n_land = mk("n_land", "integer"),
            product = product, layers = LAY, clear_rules = rules)
  O$ndvi[O$n_clear == 0] <- NA
  f <- out_of(i); dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
  saveRDS(O, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
  if (test) {
    for (k in seq_along(O$dates)) log_msg(sprintf("  %s: cells with clear pixels %.0f %%, clear share of observed land %.0f %%, median NDVI %.3f",
      format(O$dates[k]), 100 * mean(O$n_clear[, k] > 0), 100 * sum(O$n_clear[, k]) / max(1, sum(O$n_land[, k])), stats::median(O$ndvi[, k], na.rm = TRUE)))
  }
  list(days = length(O$dates), clear_cells = mean(O$n_clear > 0), clear_px = sum(O$n_clear) / max(1, sum(O$n_land)),
       ndvi_med = stats::median(O$ndvi, na.rm = TRUE))
}

active <- function() st$status %in% c("submitted", "queued", "pending", "processing", "done")
todo <- function() st$status %in% c("new", "retry")
if (!any(active() | todo())) {
  log_msg("nothing to do")
} else {
  token <- ae_login(cfg); t_login <- Sys.time(); log_msg("AppEEARS login ok")
  notify("05_viirs_daily started", sprintf("%d region-year-chunks to fetch", sum(active() | todo())), priority = 2, tags = "arrow_forward")
  t0 <- Sys.time(); n_fail <- 0L; failed_this_run <- character()
  dl_fail <- list(); dl_max_rounds <- as.integer(vd$download_rounds %||% 4)   # failed download rounds per job before resubmitting
  repeat {
    if (difftime(Sys.time(), t_login, units = "hours") > 24) { token <- ae_login(cfg); t_login <- Sys.time() }
    while (sum(active()) < pending_max && any(todo())) {
      i <- which(todo())[1]
      tid <- vd_submit(cfg, token, product, LAY, st$start[i], st$end[i], BBOX[[as.character(st$region[i])]],
                       sprintf("gsi_%s_vd_%s", profile, st$key[i]))
      st$task_id[i] <- tid; st$status[i] <- "submitted"; st$attempts[i] <- st$attempts[i] + 1L
      save_state(); log_msg("submitted ", st$key[i], " (", st$start[i], " .. ", st$end[i], "; task ", tid, ")")
    }
    for (i in which(active())) {
      s <- tryCatch(ae_task_status(cfg, token, st$task_id[i]), error = function(e) { log_warn("status ", st$key[i], ": ", conditionMessage(e)); NA })
      if (is.na(s)) next
      if (s != st$status[i]) { st$status[i] <- s; save_state(); log_msg(st$key[i], ": ", s) }
      if (s == "done" && st$key[i] %in% failed_this_run) next
      if (s == "done") {
        raw <- raw_of(i); dir.create(raw, showWarnings = FALSE, recursive = TRUE)
        files <- ae_bundle_files(cfg, token, st$task_id[i])
        files <- files[grepl("\\.tif$", files$file_name), , drop = FALSE]
        if (!nrow(files)) { log_warn(st$key[i], ": no GeoTIFFs in the bundle; marked empty"); st$status[i] <- "empty"; save_state(); next }
        t1 <- Sys.time()
        # Oct 23: a 403 on one file (AppEEARS under load before the S-NPP cutoff) stopped the whole run. Now a
        # failed download keeps the files already fetched, logs in again, and retries the job at the next poll;
        # after dl_max_rounds failed rounds the job is resubmitted (fresh bundle), then marked failed.
        dl <- tryCatch({
          for (k in seq_len(nrow(files))) {
            dest <- file.path(raw, files$file_name[k])
            if (!file.exists(dest)) ae_download_file(cfg, token, st$task_id[i], files[k, ], dest)
          }
          TRUE
        }, error = function(e) conditionMessage(e))
        if (!isTRUE(dl)) {
          k0 <- st$key[i]; dl_fail[[k0]] <- (dl_fail[[k0]] %||% 0L) + 1L
          log_warn(sprintf("%s: download round %d failed (%s); %d of %d files on disk", k0, dl_fail[[k0]], dl,
                           sum(file.exists(file.path(raw, files$file_name))), nrow(files)))
          token <- tryCatch({ tk <- ae_login(cfg); t_login <- Sys.time(); tk }, error = function(e) token)
          if (dl_fail[[k0]] >= dl_max_rounds) {
            if (st$attempts[i] < 2) {
              log_warn(k0, ": resubmitting the task after ", dl_fail[[k0]], " failed download rounds")
              unlink(raw, recursive = TRUE); st$status[i] <- "retry"; dl_fail[[k0]] <- 0L
            } else {
              log_err(k0, ": downloads failed on two tasks (task ", st$task_id[i], "); marked failed")
              st$status[i] <- "failed"; n_fail <- n_fail + 1L
            }
            save_state()
          }
          next
        }
        mb <- sum(files$file_size, na.rm = TRUE) / 1e6
        log_msg(sprintf("%s: downloaded %d files, %.0f MB, %.1f min", st$key[i], nrow(files), mb, as.numeric(difftime(Sys.time(), t1, units = "mins"))))
        t1 <- Sys.time()
        rr <- tryCatch(reduce_job(i), error = function(e) e)
        if (inherits(rr, "error")) {
          log_err(st$key[i], ": reduction failed (raw kept): ", conditionMessage(rr))
          st$status[i] <- "stack_failed"; n_fail <- n_fail + 1L; failed_this_run <- c(failed_this_run, st$key[i])
        } else {
          if (!test) unlink(raw, recursive = TRUE)
          st$status[i] <- "complete"
          log_msg(sprintf("%s: %d days reduced in %.1f min; cell-days with clear pixels %.0f %%, clear share of land pixels %.0f %%, median NDVI %.3f",
                          st$key[i], rr$days, as.numeric(difftime(Sys.time(), t1, units = "mins")), 100 * rr$clear_cells, 100 * rr$clear_px, rr$ndvi_med))
          notify(sprintf("VIIRS daily %s done", st$key[i]), sprintf("%.0f MB; %d left", mb, sum(st$status != "complete")), priority = 2, tags = "arrow_down")
        }
        save_state()
      } else if (s == "error") {
        if (st$attempts[i] < 2) { log_warn(st$key[i], ": AppEEARS task error; resubmitting"); st$status[i] <- "retry" }
        else { log_err(st$key[i], ": AppEEARS task failed twice (task ", st$task_id[i], ")"); st$status[i] <- "failed"; n_fail <- n_fail + 1L }
        save_state()
      }
    }
    if (!any((active() | todo()) & !(st$key %in% failed_this_run))) break
    Sys.sleep(60 * poll_min)
  }
  ae_logout(cfg, token)
  msg <- sprintf("%d complete, %d empty, %d failed of %d | %.1f h", sum(st$status == "complete"), sum(st$status == "empty"),
                 sum(st$status %in% c("failed", "stack_failed")), nrow(st), as.numeric(difftime(Sys.time(), t0, units = "hours")))
  log_msg("=== done | ", msg, " | state: ", state_file)
  notify(if (n_fail) "05_viirs_daily finished WITH PROBLEMS" else "05_viirs_daily finished", msg, priority = if (n_fail) 4 else 3,
         tags = if (n_fail) "warning" else "tada")
  if (n_fail > 0 && !interactive()) quit(status = 1)
}
log_msg("log: ", log_file)
