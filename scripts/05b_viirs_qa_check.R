#!/usr/bin/env Rscript
# Step 05b (Oct 21; P720, seconds): which clear-sky rule removes the most pixels in a 05 smoke test.
# Reads the raw VNP09GA files the GSI_VD_TEST=1 run kept (scratch raw/viirs_daily/VNP09GA/r<RR>/<Y>/c1)
# and, over observed land pixels (I1 present, QF2 land background), reports per day the share failing
# each rule on its own and the clear share under relaxed rule sets -> viirs_daily_test/qa_check.csv.
# Usage: GSI_VD_TEST=1 Rscript scripts/05b_viirs_qa_check.R [profile]
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "dev")
Sys.setenv(GSI_VD_TEST = "1")
source(here::here("R", "config.R")); source(here::here("R", "log.R")); source(here::here("R", "viirs_daily.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
cfg <- load_cfg(profile); log_file <- log_init("05b_viirs_qa_check", profile)
vd <- vd_conf(cfg); td <- as.Date(unlist(vd$test_dates %||% list("2023-06-01")))
raw <- file.path(vd_raw_dir(cfg, vd$product, vd$regions[1], as.integer(format(td[1], "%Y"))), "c1")
lay <- c(unlist(vd$layers), qf4 = "SurfReflect_QF4"); lay <- vapply(lay, function(p) {
  f <- list.files(raw, pattern = "\\.tif$"); h <- unique(sub("_doy\\d{7}.*$", "", f[grepl(p, f)])); h <- sub("^[^_]+\\.\\d{3}_", "", h); h[1] }, "")
log_msg("raw: ", raw, " | layers: ", paste(names(lay), lay, sep = "=", collapse = " "))
fb <- vd_files_by_day(raw, lay)
b <- function(x, k, n = 1) (x %/% 2^k) %% 2^n
OUT <- list()
for (d in fb$days) {
  f <- fb$files[[d]]
  i1 <- terra::rast(f$i1); v1 <- terra::values(i1, mat = FALSE); v2 <- terra::values(terra::rast(f$i2), mat = FALSE)
  if (stats::quantile(abs(v1[is.finite(v1) & v1 != -28672]), 0.99) > 2) { v1 <- v1 * 1e-4; v2 <- v2 * 1e-4 }
  nd <- (v2 - v1) / (v2 + v1)
  rq <- function(nm) terra::values(terra::resample(terra::rast(f[[nm]]), i1, method = "near"), mat = FALSE)
  q1 <- rq("qf1"); q2 <- rq("qf2"); q6 <- rq("qf6"); q7 <- rq("qf7"); vz <- rq("vza"); if (stats::median(vz, na.rm = TRUE) > 180) vz <- vz / 100
  q4 <- if (!is.null(f$qf4)) rq("qf4") else rep(0, length(v1))
  obs <- !is.na(v1) & v1 > -1 & v1 < 2 & !is.na(q2) & b(q2, 0, 3) <= 1 & !is.na(q1)
  cc <- b(q1, 2, 2)
  fail <- list(night = b(q1, 4) == 1, cloud_conf_gt0 = cc > 0, cloud_conf_gt1 = cc > 1, shadow = b(q2, 3) == 1,
               snow_qf2 = b(q2, 5) == 1, cirrus_refl = b(q2, 6) == 1, cirrus_emis = b(q2, 7) == 1, heavy_aerosol = b(q2, 4) == 1,
               snow_qf7 = b(q7, 0) == 1, adjacent = b(q7, 1) == 1, cirrus_qf7 = b(q7, 4) == 1,
               i1_bad = b(q6, 3) == 1, i2_bad = b(q6, 4) == 1, sdr_bad_qf4 = b(q4, 1, 2) > 0, vza_gt55 = vz > 55, vza_gt65 = vz > 65)
  base_ok <- function(conf, cirrus, adj, vmax, qual = TRUE, sdr = FALSE) {
    ok <- obs & b(q1, 4) == 0 & cc <= conf & b(q2, 3) == 0 & b(q2, 5) == 0 & b(q7, 0) == 0
    if (cirrus) ok <- ok & b(q2, 6, 2) == 0
    if (adj) ok <- ok & b(q7, 1) == 0
    if (qual) ok <- ok & b(q6, 3, 2) == 0
    if (sdr) ok <- ok & b(q4, 1, 2) == 0
    if (is.finite(vmax)) ok <- ok & vz <= vmax
    ok
  }
  sets <- list(current = base_ok(0, TRUE, TRUE, 55),
               noq = base_ok(0, TRUE, TRUE, 55, FALSE),
               noq_conf1 = base_ok(1, TRUE, TRUE, 55, FALSE),
               noq_conf1_noadj = base_ok(1, TRUE, FALSE, 55, FALSE),
               noq_conf1_noadj_vza65 = base_ok(1, TRUE, FALSE, 65, FALSE),
               noq_conf1_noadj_nocirr_vza65 = base_ok(1, FALSE, FALSE, 65, FALSE),
               noq_conf1_noadj_vza65_sdr = base_ok(1, TRUE, FALSE, 65, FALSE, TRUE),
               all_obs_vza65 = obs & b(q1, 4) == 0 & vz <= 65)
  n <- sum(obs, na.rm = TRUE)
  row <- c(date = d, obs_px = n, vapply(fail, function(x) round(mean(x[obs], na.rm = TRUE), 3), 0),
           stats::setNames(vapply(sets, function(x) round(sum(x, na.rm = TRUE) / n, 3), 0), paste0("clear_", names(sets))),
           stats::setNames(vapply(sets, function(x) round(stats::median(nd[x %in% TRUE], na.rm = TRUE), 3), 0), paste0("ndvi_", names(sets))))
  OUT[[d]] <- row
  log_msg(d, ": obs ", n, " | fail ", paste(names(fail), round(100 * vapply(fail, function(x) mean(x[obs], na.rm = TRUE), 0)), sep = " ", collapse = ", "),
          " % | clear ", paste(names(sets), round(100 * vapply(sets, function(x) sum(x, na.rm = TRUE) / n, 0)), sep = " ", collapse = ", "), " % | median NDVI ",
          paste(names(sets), round(vapply(sets, function(x) stats::median(nd[x %in% TRUE], na.rm = TRUE), 0), 3), sep = " ", collapse = ", "))
}
O <- as.data.frame(do.call(rbind, OUT)); dir.create(vd_dir(cfg), showWarnings = FALSE, recursive = TRUE); out <- vd_dir(cfg, "qa_check.csv"); utils::write.csv(O, out, row.names = FALSE)
log_msg("=== done | ", out, " | log: ", log_file)
