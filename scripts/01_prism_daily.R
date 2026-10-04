#!/usr/bin/env Rscript
# Step 01: download daily PRISM (4 km) and convert to one compressed GeoTIFF per
# variable-year, cropped to the configured bbox. Raw daily files are dropped afterwards.
#
# Usage (from the project root):
#   Rscript scripts/01_prism_daily.R [profile]     # profile: test | dev | full
# From an RStudio console (project open), as a background job so the console stays free:
#   gsi_profile <- "test"
#   rstudioapi::jobRunScript("scripts/01_prism_daily.R", workingDir = getwd(), importEnv = TRUE)
# or in the console itself: gsi_profile <- "test"; source("scripts/01_prism_daily.R")
#
# Resume-safe: finished variable-years are skipped; daily files already on disk are not
# re-downloaded. Downloads are sequential with a pause; only conversion runs in parallel
# (future multisession via furrr, one variable per worker, file paths in, list out).
#
# Output: <data_root>/<run_name>/prism_daily/<var>/prism_<var>_<year>.tif
#   one band per day, band names = ISO dates, FLT4S, DEFLATE + float predictor, tiled.

suppressPackageStartupMessages({
  library(future)
  library(furrr)
})

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) {
  args[1]
} else if (exists("gsi_profile", envir = globalenv())) {
  get("gsi_profile", envir = globalenv())
} else {
  Sys.getenv("R_CONFIG_ACTIVE", "dev")
}

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "prism.R"))
source(here::here("R", "notify.R"))

cfg <- load_cfg(profile)
log_file <- log_init("01_prism_daily", profile)
paths <- prism_paths(cfg)
dir.create(paths$out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$raw_root, showWarnings = FALSE, recursive = TRUE)

bbox_txt <- if (is.null(cfg$bbox)) "none (full CONUS grid)" else
  paste(names(cfg$bbox), cfg$bbox, sep = "=", collapse = " ")
log_msg("=== 01_prism_daily | profile: ", profile, " | host: ", Sys.info()[["nodename"]],
        " | ", R.version.string, " | terra ", as.character(utils::packageVersion("terra")))
log_msg("run_dir: ", cfg$run_dir, " | scratch: ", cfg$scratch_dir)
log_msg("vars: ", paste(cfg$prism$vars, collapse = ","), " | years: ",
        min(cfg$years), "-", max(cfg$years), " | bbox: ", bbox_txt,
        " | workers: ", cfg$workers, " | memfrac: ", cfg$memfrac)
log_msg("source: ", cfg$prism$base_url, "/", cfg$prism$region, "/", cfg$prism$res,
        " | sleep ", cfg$prism$sleep_sec, " s between requests")

notify_init(cfg)
notify("01_prism_daily started",
       sprintf("vars %s | years %d-%d | bbox %s | %d workers",
               paste(cfg$prism$vars, collapse = ","), min(cfg$years), max(cfg$years),
               if (is.null(cfg$bbox)) "CONUS" else "set", cfg$workers),
       priority = 2, tags = "arrow_forward")

# Fatal errors: notify, then let R stop as usual. Restored at the end of a normal run.
old_error_opt <- getOption("error")
options(error = function() {
  notify("01_prism_daily FAILED", paste("Stopped with error:", geterrmessage()),
         priority = 5, tags = "rotating_light")
  options(error = old_error_opt)
  if (!interactive()) quit(status = 1)
})

# Single-threaded BLAS/OpenMP in workers (inherited by multisession workers).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
plan(multisession, workers = cfg$workers)

t_start <- Sys.time()
n_skip <- 0L; n_done <- 0L; n_fail <- 0L

year_secs <- numeric()   # elapsed time of years that did work, for the ETA

for (yr in cfg$years) {
  t_year <- Sys.time()
  dates <- prism_days(yr, cfg$prism$max_days_per_year)
  if (!length(dates)) { log_warn(yr, ": no complete days, skipping"); next }
  if (max(dates) > Sys.Date() - 183) {
    log_warn(yr, ": includes days < 6 months old (may be provisional and change later)")
  }

  todo <- character()
  for (v in cfg$prism$vars) {
    if (file.exists(prism_out_file(cfg, v, yr))) {
      log_msg(v, " ", yr, ": output exists, skip"); n_skip <- n_skip + 1L
    } else todo <- c(todo, v)
  }
  if (!length(todo)) next

  # 1. Sequential, polite downloads for every variable this year.
  ready <- character()
  for (v in todo) {
    missing <- prism_download_var_year(cfg, v, yr, dates)
    if (length(missing)) n_fail <- n_fail + 1L else ready <- c(ready, v)
  }
  if (!length(ready)) next

  # 2. Parallel conversion: one variable per worker; pass paths and plain values only.
  jobs <- lapply(ready, function(v) list(
    var = v,
    tifs = prism_day_file(prism_raw_dir(cfg, v, yr), v, dates),
    out_file = prism_out_file(cfg, v, yr)
  ))
  for (j in jobs) dir.create(dirname(j$out_file), showWarnings = FALSE, recursive = TRUE)
  log_msg(yr, ": converting ", paste(ready, collapse = ","), " on ", nbrOfWorkers(), " worker(s)")
  bbox <- cfg$bbox; memfrac <- cfg$memfrac
  results <- future_map(
    jobs,
    function(j) prism_convert_year(j$tifs, dates, j$out_file, bbox, memfrac),
    .options = furrr_options(seed = NULL, packages = "terra")
  )

  ok_txt <- character(); bad_txt <- character()
  for (k in seq_along(jobs)) {
    r <- results[[k]]; v <- jobs[[k]]$var
    if (isTRUE(r$ok)) {
      ok_txt <- c(ok_txt, sprintf("%s %.0f MB", v, r$mb))
      n_done <- n_done + 1L
      log_msg(sprintf("%s %d: wrote %s | %d bands | %d x %d | res %.5f | %.1f MB | crs %s | day1 range %.2f..%.2f",
                      v, yr, basename(r$out_file), r$nlyr, r$nrow, r$ncol, r$res[1], r$mb,
                      r$crs, r$range[1], r$range[2]))
      if (isTRUE(cfg$prism$drop_raw)) {
        unlink(prism_raw_dir(cfg, v, yr), recursive = TRUE)
      }
    } else {
      n_fail <- n_fail + 1L
      bad_txt <- c(bad_txt, v)
      log_err(v, " ", yr, ": conversion failed: ", r$msg, " (raw kept)")
    }
  }

  # Year summary + ETA from the average time of years that did work.
  year_secs <- c(year_secs, as.numeric(difftime(Sys.time(), t_year, units = "secs")))
  left <- sum(cfg$years > yr)
  eta <- if (left > 0) sprintf(" | %d years left, ETA ~%s", left,
                               fmt_dur(left * mean(year_secs))) else ""
  notify(sprintf("%d processed", yr),
         paste0(if (length(ok_txt)) paste("written:", paste(ok_txt, collapse = ", ")) else "",
                if (length(bad_txt)) paste0("\nFAILED: ", paste(bad_txt, collapse = ", ")) else "",
                "\nyear took ", fmt_dur(tail(year_secs, 1)), eta),
         priority = if (length(bad_txt)) 4 else 3,
         tags = if (length(bad_txt)) "warning" else "white_check_mark")
}

plan(sequential)

# Summary of everything on disk for this run.
outs <- list.files(paths$out_root, pattern = "^prism_.*\\.tif$", recursive = TRUE,
                   full.names = TRUE)
log_msg(sprintf("=== done in %.1f min | written %d | skipped %d | failed/incomplete %d",
                as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                n_done, n_skip, n_fail))
log_msg(sprintf("%d variable-year files in %s (%.1f MB total)", length(outs),
                paths$out_root, sum(file.size(outs)) / 1e6))
for (f in outs) log_msg("  ", sub(paste0(paths$out_root, "/"), "", f, fixed = TRUE),
                        sprintf("  %.1f MB", file.size(f) / 1e6))
log_msg("log: ", log_file)

notify(if (n_fail > 0) "01_prism_daily finished WITH PROBLEMS" else "01_prism_daily finished",
       sprintf("%.1f h | written %d, skipped %d, failed/incomplete %d\n%d files, %.1f GB in %s",
               as.numeric(difftime(Sys.time(), t_start, units = "hours")), n_done, n_skip,
               n_fail, length(outs), sum(file.size(outs)) / 1e9, paths$out_root),
       priority = if (n_fail > 0) 4 else 3,
       tags = if (n_fail > 0) "warning" else "tada")
options(error = old_error_opt)
if (n_fail > 0) {
  if (interactive()) warning(n_fail, " variable-year(s) failed or incomplete; see log") else quit(status = 1)
}
