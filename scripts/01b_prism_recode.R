#!/usr/bin/env Rscript
# Step 01b: re-encode existing PRISM variable-year files to the storage set in config
# (prism$scale), without downloading anything. Used once to convert the float32 dev files
# to scaled INT2S; also handles any later change of storage settings.
#
# Usage: Rscript scripts/01b_prism_recode.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/01b_prism_recode.R")
#
# Resume-safe: files already in the target storage are skipped. Each file is written to a
# temp path, checked band-by-band against the original, then swapped in.

suppressPackageStartupMessages({
  library(future)
  library(furrr)
})

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) {
  args[1]
} else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else {
  Sys.getenv("R_CONFIG_ACTIVE", "dev")
}

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "prism.R"))
source(here::here("R", "notify.R"))

cfg <- load_cfg(profile)
log_file <- log_init("01b_prism_recode", profile)
notify_init(cfg)
log_msg("=== 01b_prism_recode | profile: ", profile, " | run_dir: ", cfg$run_dir)

jobs <- list()
for (v in cfg$prism$vars) {
  scale <- prism_scale(cfg, v)
  target <- if (is.null(scale)) "FLT4S" else "INT2S"
  for (yr in cfg$years) {
    f <- prism_out_file(cfg, v, yr)
    if (!file.exists(f)) next
    current <- terra::datatype(terra::rast(f))[1]
    if (identical(current, target)) {
      log_msg(v, " ", yr, ": already ", target, ", skip")
    } else {
      jobs[[length(jobs) + 1]] <- list(var = v, year = yr, f = f, scale = scale,
                                       from = current, to = target,
                                       mb_before = file.size(f) / 1e6)
    }
  }
}
log_msg(length(jobs), " file(s) to re-encode")

if (length(jobs)) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  plan(multisession, workers = cfg$workers)
  memfrac <- cfg$memfrac
  results <- future_map(jobs, function(j) prism_recode_file(j$f, j$scale, memfrac),
                        .options = furrr_options(seed = NULL, packages = "terra"))
  plan(sequential)

  n_fail <- 0L; mb_before <- 0; mb_after <- 0
  for (k in seq_along(jobs)) {
    j <- jobs[[k]]; r <- results[[k]]
    if (isTRUE(r$ok)) {
      mb_before <- mb_before + j$mb_before; mb_after <- mb_after + r$mb
      log_msg(sprintf("%s %d: %s -> %s | %.1f -> %.1f MB | QA max err %.4f",
                      j$var, j$year, j$from, r$datatype, j$mb_before, r$mb, r$max_err))
    } else {
      n_fail <- n_fail + 1L
      log_err(j$var, " ", j$year, ": re-encode failed (original kept): ", r$msg)
    }
  }
  msg <- sprintf("%d re-encoded, %d failed | %.0f -> %.0f MB", length(jobs) - n_fail,
                 n_fail, mb_before, mb_after)
  log_msg("=== done | ", msg)
  notify("01b_prism_recode finished", msg, priority = if (n_fail) 4 else 3,
         tags = if (n_fail) "warning" else "white_check_mark")
  if (n_fail > 0 && !interactive()) quit(status = 1)
}
log_msg("log: ", log_file)
