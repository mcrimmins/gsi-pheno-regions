# NPN grass curing study: paths, config, logging and the small climate helpers the scripts use.
# Self-contained (Oct 18 2026): nothing here reads the parent gsi-pheno-regions repo. Works with
# npn_study/ as a subfolder of that repo or, after export, as the root of its own repo.
# Each script starts with:
#   npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
#   source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

if (!exists("npn_root")) npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
npn_path <- function(...) file.path(npn_root, ...)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- config ----------------------------------------------------------------------------------
# profile from the first argument; npn_study/config.yml. Adds run_dir (outputs go to
# <run_dir>/npn_study/), study_dir, prism_daily_dir; resolved per OS.
load_npn_cfg <- function(profile = NULL) {
  if (is.null(profile) || !nzchar(profile)) profile <- Sys.getenv("R_CONFIG_ACTIVE", "full")
  cfg <- config::get(config = profile, file = npn_path("config.yml"), use_parent = FALSE)
  os <- if (.Platform$OS.type == "windows") "windows" else "linux"
  pick <- function(x) if (is.list(x)) x[[os]] else x
  root <- Sys.getenv("NPN_DATA_ROOT", "")
  if (!nzchar(root)) root <- pick(cfg$data_root)
  cfg$data_root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
  cfg$run_dir <- file.path(cfg$data_root, cfg$run_name)
  cfg$study_dir <- path.expand(pick(cfg$study_dir) %||% "")
  pd <- pick(cfg$prism_daily_dir)
  cfg$prism_daily_dir <- if (is.null(pd) || !nzchar(pd)) file.path(cfg$run_dir, "prism_daily") else path.expand(pd)
  cfg$profile <- profile
  cfg
}
# Folder with npn_curing_metrics.rds and npn_green_obs.rds (GSI_OBS_DIR overrides)
npn_obs_dir <- function(cfg) Sys.getenv("GSI_OBS_DIR", file.path(cfg$study_dir, "data", "derived"))

# ---- logging: timestamped lines to the console and npn_study/logs/<script>_<profile>_<time>.log
log_init <- function(script, profile) {
  dir.create(npn_path("logs"), showWarnings = FALSE)
  f <- npn_path("logs", sprintf("%s_%s_%s.log", script, profile, format(Sys.time(), "%Y%m%d-%H%M%S")))
  options(npn.log_file = f)
  f
}
log_msg <- function(..., level = "INFO") {
  line <- sprintf("%s [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(...))
  message(line)
  f <- getOption("npn.log_file")
  if (!is.null(f)) cat(line, "\n", file = f, append = TRUE, sep = "")
  invisible(line)
}
log_warn <- function(...) log_msg(..., level = "WARN")
log_err <- function(...) log_msg(..., level = "ERROR")

# ---- climate helpers (same definitions as gsi-pheno-regions R/features.R and R/block2.R) --------
# Extraterrestrial radiation Ra (MJ m-2 d-1), FAO-56 eq. 21. lat in degrees (vector),
# doy vector; returns length(lat) x length(doy) matrix.
ra_matrix <- function(lat, doy) {
  phi <- lat * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
  dec <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmin(pmax(-tan(outer(phi, rep(1, length(doy)))) *
                        tan(outer(rep(1, length(phi)), dec)), -1), 1))
  sp <- outer(sin(phi), sin(dec)); cp <- outer(cos(phi), cos(dec))
  (24 * 60 / pi) * 0.0820 * matrix(dr, length(phi), length(doy), byrow = TRUE) *
    (ws * sp + cp * sin(ws))
}
# Hargreaves PET (mm/day): 0.0023 * 0.408 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin)
hargreaves <- function(tmin, tmax, ra) {
  pmax(0.0023 * 0.408 * ra * ((tmin + tmax) / 2 + 17.8) * sqrt(pmax(tmax - tmin, 0)), 0)
}
# Last column index where m is TRUE per row (0 if none); first index (ncol + 1 if none).
last_true <- function(m) {
  if (!ncol(m)) return(rep(0, nrow(m)))
  Reduce(pmax, lapply(seq_len(ncol(m)), function(j) ifelse(m[, j] %in% TRUE, j, 0)))
}
first_true <- function(m) {
  if (!ncol(m)) return(rep(1, nrow(m)))
  Reduce(pmin, lapply(seq_len(ncol(m)), function(j) ifelse(m[, j] %in% TRUE, j, ncol(m) + 1)))
}
# Dry-run lengths ending at each day (0 on wet days), for a logical sites x days matrix.
run_lengths <- function(D) {
  R <- matrix(0L, nrow(D), ncol(D)); cur <- integer(nrow(D))
  for (j in seq_len(ncol(D))) { cur <- ifelse(D[, j], cur + 1L, 0L); R[, j] <- cur }
  R
}
row_cumsum <- function(X) { for (j in seq_len(ncol(X))[-1]) X[, j] <- X[, j - 1] + X[, j]; X }
# Forward window sums: S[, d] = sum(X[, d .. d + w - 1]) for d in idx (window must fit).
fwd_sum <- function(X, w, idx) {
  C <- cbind(0, row_cumsum(X))
  C[, idx + w, drop = FALSE] - C[, idx, drop = FALSE]
}
# Warm-season rain onset (calendar year): first day d in [start_doy, end_doy] that is wet,
# starts a window_days total >= total_mm, and is not followed by a dry spell >= max_dry_spell
# within check_days. Returns DOY, or end_doy + 1 when no onset (censored). b2 = cfg$climate.
onset_doy <- function(P, cal_doy, b2) {
  o <- b2$onset
  idx <- which(cal_doy >= o$start_doy & cal_doy <= o$end_doy)
  need <- max(idx) + max(o$window_days - 1, o$check_days)
  if (need > ncol(P)) stop("onset window runs past the end of the year")
  wet <- P[, idx, drop = FALSE] >= b2$dry_day_mm
  S <- fwd_sum(P, o$window_days, idx)
  RL <- run_lengths(P[, seq_len(need), drop = FALSE] < b2$dry_day_mm)
  M <- matrix(0L, nrow(P), length(idx))
  for (k in seq_len(o$check_days)) M <- pmax(M, RL[, idx + k, drop = FALSE])
  cand <- wet & S >= o$total_mm & M < o$max_dry_spell
  f <- first_true(cand)
  ifelse(f > length(idx), o$end_doy + 1, cal_doy[idx][pmin(f, length(idx))])
}
