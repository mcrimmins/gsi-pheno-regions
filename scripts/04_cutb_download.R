#!/usr/bin/env Rscript
# Step 04 (Cut B): download MODIS 16-day VI composites (MOD13A2 + MYD13A2, 1 km) for the
# profile bbox (CONUS when null) and years via AppEEARS area requests, then aggregate them
# to the 4 km analysis grid (mean of valid-quality 1 km pixels).
#
# One AppEEARS task = one product-year (all configured layers). Up to max_pending_tasks
# are queued at once; the script polls, downloads finished bundles (size + SHA-256
# checked), aggregates them, deletes the raw 1 km files, and submits the next.
#
# Output: <run_dir>/cutb/<MOD13A2|MYD13A2>/<var>/<product>_<var>_<year>.tif on the 4 km
#   grid; bands = composites (named by composite start date); physical units (scale factor
#   applied); vars = value_vars + n_valid (valid 1 km pixels) + snow_frac.
#   <run_dir>/cutb/<product>/layers.json keeps the AppEEARS layer metadata.
# Needs static/grid_mask.tif (02_static.R).
# State: <run_dir>/cutb/tasks.csv (task ids), so a restart never resubmits finished work.
#
# Credentials: NASA Earthdata login in a netrc file (see R/appeears.R). Never in code.
# Usage: Rscript scripts/04_cutb_download.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/04_cutb_download.R")

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
source(here::here("R", "notify.R"))
source(here::here("R", "appeears.R"))
source(here::here("R", "cutb.R"))
source(here::here("R", "grid.R"))

cfg <- load_cfg(profile)
log_file <- log_init("04_cutb_download", profile)
notify_init(cfg)
terra::terraOptions(progress = 0)
cb <- cfg$cutb
products <- unlist(cb$products)
bbox <- cutb_bbox(cfg)
dir.create(cutb_dir(cfg), showWarnings = FALSE, recursive = TRUE)
state_file <- file.path(cutb_dir(cfg), "tasks.csv")
grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
if (!file.exists(grid_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
out_vars <- cutb_out_vars(cfg)

log_msg("=== 04_cutb_download | profile: ", profile, " | host: ", Sys.info()[["nodename"]])
log_msg("products: ", paste(products, collapse = ", "), " | years: ", min(cfg$years), "-",
        max(cfg$years), " | bbox: ", paste(names(bbox), bbox, sep = "=", collapse = " "))

old_error_opt <- getOption("error")
options(error = function() {
  notify("04_cutb_download FAILED", paste("Stopped with error:", geterrmessage()),
         priority = 5, tags = "rotating_light")
  options(error = old_error_opt)
  if (!interactive()) quit(status = 1)
})

# ---- layers per product (checked against AppEEARS, saved with scale/fill) ----------
resolved <- list()
for (p in products) {
  resolved[[p]] <- ae_resolve_layers(cfg, p)
  cutb_write_meta(cfg, p, resolved[[p]])
  log_msg(p, " layers: ", paste(names(resolved[[p]]$layers), resolved[[p]]$layers,
                                sep = "=", collapse = " | "))
}

# ---- jobs and state -------------------------------------------------------------------
jobs <- expand.grid(product = products, year = cfg$years, stringsAsFactors = FALSE)
jobs$key <- paste(jobs$product, jobs$year, sep = "_")
outs_done <- function(p, y) all(file.exists(vapply(out_vars,
                                                   function(v) cutb_out_file(cfg, p, v, y), "")))
prev <- if (file.exists(state_file)) {
  utils::read.csv(state_file, stringsAsFactors = FALSE, colClasses = c(task_id = "character"))
} else data.frame(key = character(), task_id = character(), status = character(),
                  attempts = integer(), stringsAsFactors = FALSE)
st <- merge(jobs, prev[, c("key", "task_id", "status", "attempts")], by = "key", all.x = TRUE,
            sort = FALSE)
st <- st[order(st$year, st$product), ]
st$status[is.na(st$status)] <- "new"; st$attempts[is.na(st$attempts)] <- 0L
for (i in seq_len(nrow(st))) if (outs_done(st$product[i], st$year[i])) st$status[i] <- "complete"
# A task that finished at NASA but wasn't downloaded/stacked yet is picked up again; failed
# stacking is retried once per run (raw files were kept).
st$status[st$status == "stack_failed"] <- "done"
save_state <- function() {
  st$updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  utils::write.csv(st[, c("key", "product", "year", "task_id", "status", "attempts", "updated")],
                   state_file, row.names = FALSE)
}
save_state()
log_msg(sum(st$status == "complete"), " of ", nrow(st), " product-years already complete")

active <- function() st$status %in% c("submitted", "queued", "pending", "processing", "done")
todo <- function() st$status %in% c("new", "retry")
if (!any(active() | todo())) {
  log_msg("nothing to do"); options(error = old_error_opt)
} else {
  token <- ae_login(cfg); t_login <- Sys.time()
  log_msg("AppEEARS login ok")
  notify("04_cutb_download started", sprintf("%d product-years to fetch",
                                             sum(active() | todo())), priority = 2,
         tags = "arrow_forward")
  t0 <- Sys.time(); n_done <- 0L; n_fail <- 0L; failed_this_run <- character()

  repeat {
    if (difftime(Sys.time(), t_login, units = "hours") > 24) {
      token <- ae_login(cfg); t_login <- Sys.time()
    }
    # submit while slots are free
    while (sum(active()) < cb$max_pending_tasks && any(todo())) {
      i <- which(todo())[1]
      tid <- ae_submit_year(cfg, token, st$product[i], resolved[[st$product[i]]]$layers,
                            st$year[i], bbox, sprintf("gsi_%s_%s", profile, st$key[i]))
      st$task_id[i] <- tid; st$status[i] <- "submitted"; st$attempts[i] <- st$attempts[i] + 1L
      save_state(); log_msg("submitted ", st$key[i], " (task ", tid, ")")
    }
    # check active tasks
    for (i in which(active())) {
      s <- tryCatch(ae_task_status(cfg, token, st$task_id[i]),
                    error = function(e) { log_warn("status ", st$key[i], ": ", conditionMessage(e)); NA })
      if (is.na(s)) next
      if (s != st$status[i]) { st$status[i] <- s; save_state(); log_msg(st$key[i], ": ", s) }
      if (s == "done" && st$key[i] %in% failed_this_run) next
      if (s == "done") {
        p <- st$product[i]; y <- st$year[i]
        raw <- cutb_raw_dir(cfg, p, y); dir.create(raw, showWarnings = FALSE, recursive = TRUE)
        files <- ae_bundle_files(cfg, token, st$task_id[i])
        files <- files[grepl("\\.tif$", files$file_name), , drop = FALSE]
        for (k in seq_len(nrow(files))) {
          dest <- file.path(raw, files$file_name[k])
          if (!file.exists(dest)) ae_download_file(cfg, token, st$task_id[i], files[k, ], dest)
        }
        res <- tryCatch({
          outs <- setNames(lapply(out_vars, function(v) cutb_out_file(cfg, p, v, y)), out_vars)
          cutb_aggregate_year(raw, resolved[[p]]$layers, resolved[[p]]$meta,
                              unlist(cb$value_vars), as.integer(unlist(cb$valid_reliability)),
                              grid_file, outs)
        }, error = function(e) e)
        if (inherits(res, "error")) {
          log_err(st$key[i], ": aggregation failed (raw kept): ", conditionMessage(res))
          st$status[i] <- "stack_failed"; n_fail <- n_fail + 1L
          failed_this_run <- c(failed_this_run, st$key[i])
        } else {
          unlink(raw, recursive = TRUE)
          st$status[i] <- "complete"; n_done <- n_done + 1L
          n <- res[[1]]$n
          log_msg(sprintf("%s: aggregated %d vars x %d composites (%s .. %s), %.0f MB downloaded; cells with valid data per composite: median %.0f %%, min %.0f %%",
                          st$key[i], length(out_vars), n, res[[1]]$first, res[[1]]$last,
                          sum(files$file_size, na.rm = TRUE) / 1e6,
                          100 * stats::median(res$valid_share), 100 * min(res$valid_share)))
          left <- sum(st$status != "complete")
          notify(sprintf("Cut B %s done", st$key[i]),
                 sprintf("%d composites; %d product-years left", n, left - 0L),
                 priority = 2, tags = "arrow_down")
        }
        save_state()
      } else if (s == "error") {
        if (st$attempts[i] < 2) {
          log_warn(st$key[i], ": AppEEARS task error; resubmitting")
          st$status[i] <- "retry"
        } else {
          log_err(st$key[i], ": AppEEARS task failed twice (task ", st$task_id[i], ")")
          st$status[i] <- "failed"; n_fail <- n_fail + 1L
        }
        save_state()
      }
    }
    if (!any((active() | todo()) & !(st$key %in% failed_this_run))) break
    Sys.sleep(60 * cb$poll_minutes)
  }

  ae_logout(cfg, token)
  msg <- sprintf("%d complete, %d failed, %d total | %.1f h", sum(st$status == "complete"),
                 sum(st$status %in% c("failed", "stack_failed")), nrow(st),
                 as.numeric(difftime(Sys.time(), t0, units = "hours")))
  log_msg("=== done | ", msg, " | state: ", state_file)
  notify(if (n_fail) "04_cutb_download finished WITH PROBLEMS" else "04_cutb_download finished",
         msg, priority = if (n_fail) 4 else 3, tags = if (n_fail) "warning" else "tada")
  options(error = old_error_opt)
  if (n_fail > 0 && !interactive()) quit(status = 1)
}
log_msg("log: ", log_file)
