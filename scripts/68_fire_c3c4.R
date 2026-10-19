#!/usr/bin/env Rscript
# Step 68 (fire activity study; P720, minutes): Cheryl's C3/C4 question (scope decision 23e). In mixed
# prairie, does fire activity dip in the gap between the cool-season (C3) and warm-season (C4) growth
# pulses?
# Two-pulse cells: Plains herbaceous cells (config fire_study: c3c4: regions; cells as 65: strict herb
# cells, or fuel share >= min_fuel_share in pheno_source regions) whose typical (across-year median)
# satellite NDVI curve has two peaks between Mar and Oct: the two highest local maxima at least
# min_sep_days apart, a trough between them at least dip_frac of the curve's amplitude below the lower
# peak, and the lower peak at least peak_frac of the amplitude above the curve's minimum. Others in the
# same regions are one-pulse cells.
# Test: fire rate (FPA-FOD wildfires >= 10 ac, no debris, 2001-2024) in each two-pulse cell's own
# windows (+/- half_window_days around its first peak, trough and second peak), divided by the rate
# one-pulse cells have on the same calendar days (their daily rate, smoothed over 15 days). index =
# two-pulse rate / one-pulse rate per window; dip = index(gap) / mean(index(peak 1), index(peak 2)),
# < 1 = relatively fewer fires in the gap. Bootstrap over years (2000) for 90 % intervals.
# Outputs <run_dir>/fire_study/c3c4/: two_pulse.tif (class 1 one-pulse / 2 two-pulse, peak1, trough,
# peak2 DOY), cells_summary.csv, windows.csv (rates, index, dip with intervals), weekly.csv (fires per
# 1000 cells per week for both classes, with the two-pulse median peak and trough weeks).
#
# Usage: Rscript scripts/68_fire_c3c4.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("68_fire_c3c4", profile)
notify_init(cfg)
fs <- cfg$fire_study; cc <- fs$c3c4 %||% list()
t0 <- Sys.time(); terra::terraOptions(progress = 0); set.seed(68)
out_dir <- fire_dir(cfg, "c3c4"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
regs <- as.integer(unlist(cc$regions %||% list(10, 2)))
yrs <- seq(as.integer(fs$fod$years$start), as.integer(fs$fod$years$end))
hw <- as.integer(cc$half_window_days %||% 14)
log_msg("=== 68_fire_c3c4 | profile: ", profile, " | regions ", paste(regs, collapse = ", "))

grid <- terra::rast(file.path(cfg$run_dir, "static", "grid_mask.tif"))
fuel <- fire_fuel_share(cfg)
TD <- terra::rast(file.path(cfg$run_dir, "rules", "typical_dates.tif"))
tdv <- terra::values(TD[[c("herb", "region")]]); region <- as.integer(tdv[, "region"])
ps <- fs$pheno_source %||% list()
cells <- unlist(lapply(regs, function(r)
  if (r %in% as.integer(unlist(ps$regions))) which(region == r & !is.na(fuel) & fuel >= as.numeric(ps$min_fuel_share %||% 0.35))
  else which(region == r & tdv[, "herb"] %in% 1)))
log_msg(format(length(cells), big.mark = ","), " herbaceous cells")

# ---- typical NDVI curve (46 x 8-day values) -------------------------------------------------------------
sm <- file.path(cfg$run_dir, "summaries", "cutb_ndvi_curve_median.tif")
if (file.exists(sm)) {
  V <- as.matrix(terra::extract(terra::rast(sm), cells))
  log_msg("typical curves from ", basename(sm))
} else {
  cf <- file.path(cfg$run_dir, "cutb", "curves", sprintf("ndvi_%d.tif", yrs))
  cf <- cf[file.exists(cf)]
  if (!length(cf)) stop("no summaries/cutb_ndvi_curve_median.tif and no cutb/curves/ndvi_<Y>.tif")
  A <- simplify2array(lapply(cf, function(f) as.matrix(terra::extract(terra::rast(f), cells))))
  V <- apply(A, c(1, 2), stats::median, na.rm = TRUE); rm(A)
  log_msg("typical curves = median of ", length(cf), " yearly curve files")
}
doy <- seq(1, by = 8, length.out = ncol(V))
two_pulse <- function(v) {
  if (any(is.na(v))) return(c(NA, NA, NA, NA))
  s <- (c(v[length(v)], v[-length(v)]) + v + c(v[-1], v[1])) / 3            # light circular smoothing
  amp <- max(s) - min(s); if (amp <= 0) return(c(1, NA, NA, NA))
  k <- which(doy >= 60 & doy <= 300)
  pk <- k[s[k] >= s[k - 1] & s[k] >= s[k + 1]]
  if (length(pk) < 2) return(c(1, NA, NA, NA))
  best <- NULL
  for (a in pk) for (b in pk) if (b > a && (doy[b] - doy[a]) >= (cc$min_sep_days %||% 48)) {
    tr <- a - 1 + which.min(s[a:b]); low <- min(s[a], s[b])
    if (low - s[tr] >= (cc$dip_frac %||% 0.15) * amp && low - min(s) >= (cc$peak_frac %||% 0.4) * amp) {
      sc <- s[a] + s[b]
      if (is.null(best) || sc > best$sc) best <- list(sc = sc, a = a, b = b, t = tr)
    }
  }
  if (is.null(best)) c(1, NA, NA, NA) else c(2, doy[best$a] + 3.5, doy[best$t] + 3.5, doy[best$b] + 3.5)
}
CL <- t(apply(V, 1, two_pulse)); colnames(CL) <- c("class", "peak1", "trough", "peak2")
out <- matrix(NA_real_, terra::ncell(grid), 4); out[cells, ] <- CL
write_feature_matrix(out, grid, colnames(CL), file.path(out_dir, "two_pulse.tif"))
T2 <- which(CL[, "class"] %in% 2); O1 <- which(CL[, "class"] %in% 1)
SM <- do.call(rbind, lapply(regs, function(r) { i <- region[cells] == r
  data.frame(region = r, cells = sum(i), two_pulse = sum(CL[i, "class"] %in% 2), one_pulse = sum(CL[i, "class"] %in% 1),
             peak1_doy = stats::median(CL[i & CL[, "class"] %in% 2, "peak1"]), trough_doy = stats::median(CL[i & CL[, "class"] %in% 2, "trough"]),
             peak2_doy = stats::median(CL[i & CL[, "class"] %in% 2, "peak2"])) }))
utils::write.csv(SM, file.path(out_dir, "cells_summary.csv"), row.names = FALSE)
for (k in seq_len(nrow(SM))) with(SM[k, ], log_msg(sprintf("region %2d: %d cells, %d two-pulse (%.0f %%); median peaks DOY %.0f and %.0f, trough %.0f",
  region, cells, two_pulse, 100 * two_pulse / max(cells, 1), peak1_doy, peak2_doy, trough_doy)))
if (length(T2) < 25) { log_warn("fewer than 25 two-pulse cells: no test"); if (!interactive()) quit(save = "no") else stop("too few two-pulse cells") }

# ---- fires -------------------------------------------------------------------------------------------------
FOD <- fire_read_table(fire_dir(cfg, "fod", "fod_fires"))
FOD <- FOD[FOD$c_plus & !FOD$debris & FOD$year %in% yrs & FOD$cell %in% cells, ]
FOD$doy <- as.integer(format(as.Date(FOD$date), "%j")); FOD$doy[FOD$doy > 365] <- 365
FOD$k <- match(FOD$cell, cells)
cls <- CL[FOD$k, "class"]
log_msg(sprintf("fires: %d in two-pulse cells, %d in one-pulse cells", sum(cls %in% 2), sum(cls %in% 1)))

win_stats <- function(F) {     # F: fires (with k, doy, year); returns rate ratios per window
  ny <- length(yrs)
  f1 <- F[CL[F$k, "class"] %in% 2, ]; f0 <- F[CL[F$k, "class"] %in% 1, ]
  d0 <- tabulate(f0$doy, 365) / (length(O1) * ny)                 # one-pulse daily rate per cell-year
  d0 <- as.numeric(stats::filter(d0, rep(1 / 15, 15), circular = TRUE))
  res <- c()
  for (w in c("peak1", "trough", "peak2")) {
    c0 <- round(CL[T2, w])
    lo <- c0 - hw; hi <- c0 + hw
    inwin <- f1$doy >= lo[match(f1$k, T2)] & f1$doy <= hi[match(f1$k, T2)]
    r1 <- sum(inwin) / (length(T2) * ny * (2 * hw + 1))
    r0 <- mean(vapply(seq_along(T2), function(j) mean(d0[((lo[j]:hi[j]) - 1) %% 365 + 1]), 0))
    res[w] <- r1 / r0
  }
  c(res, dip = unname(res["trough"] / mean(res[c("peak1", "peak2")])))
}
obs <- win_stats(FOD)
B <- replicate(as.integer(cc$n_boot %||% 2000), {
  ys <- sample(yrs, replace = TRUE)
  win_stats(do.call(rbind, lapply(ys, function(y) FOD[FOD$year == y, ])))
})
WN <- data.frame(window = names(obs), index = unname(obs), lo = apply(B, 1, stats::quantile, 0.05, na.rm = TRUE),
                 hi = apply(B, 1, stats::quantile, 0.95, na.rm = TRUE))
utils::write.csv(WN, file.path(out_dir, "windows.csv"), row.names = FALSE)
log_msg(sprintf("two-pulse / one-pulse fire rate: first peak %.2f [%.2f, %.2f], gap %.2f [%.2f, %.2f], second peak %.2f [%.2f, %.2f]; dip %.2f [%.2f, %.2f]",
                WN$index[1], WN$lo[1], WN$hi[1], WN$index[2], WN$lo[2], WN$hi[2], WN$index[3], WN$lo[3], WN$hi[3], WN$index[4], WN$lo[4], WN$hi[4]))

wk <- fire_week(as.Date("2001-01-01") + FOD$doy - 1)
WK <- data.frame(week = 1:52,
  two_pulse_per_1000 = tabulate(wk[cls %in% 2], 52) / length(T2) * 1000 / length(yrs),
  one_pulse_per_1000 = tabulate(wk[cls %in% 1], 52) / length(O1) * 1000 / length(yrs))
WK$peak1_week <- fire_week(as.Date("2001-01-01") + stats::median(CL[T2, "peak1"]) - 1)
WK$trough_week <- fire_week(as.Date("2001-01-01") + stats::median(CL[T2, "trough"]) - 1)
WK$peak2_week <- fire_week(as.Date("2001-01-01") + stats::median(CL[T2, "peak2"]) - 1)
utils::write.csv(WK, file.path(out_dir, "weekly.csv"), row.names = FALSE)
msg <- sprintf("%d two-pulse cells, dip %.2f, %.1f min", length(T2), obs[["dip"]], as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("68_fire_c3c4 finished", msg, priority = 3, tags = "white_check_mark")
