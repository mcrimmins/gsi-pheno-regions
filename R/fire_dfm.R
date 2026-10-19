# Dead-fuel dryness from daily PRISM for the fire study (script 64; config fire_study: dfm).
# All functions work on cells x days matrices (one row per grid cell, one column per day).
#
# 1000-hour fuel moisture, NFDRS 1978 (Fosberg et al. 1981; Bradshaw et al. 1984 eqs. 38-42;
# Cohen and Deeming 1985). Checked Oct 18 2026 against Bradshaw et al. (1984, INT-GTR-169),
# Fosberg et al. (1981) and Ottmar and Sandberg (1985, PNW-RP-336):
#   EMC (%), T in F, RH in % (Simard 1968):
#     RH <= 10:      0.03229 + 0.281073 RH - 0.000578 RH T
#     10 < RH <= 50: 2.22749 + 0.160107 RH - 0.014784 T
#     RH > 50:       21.0606 + 0.005565 RH^2 - 0.00035 RH T - 0.483199 RH
#   daily EMC = [n EMC(Tmax, RHmin) + (24 - n) EMC(Tmin, RHmax)] / 24, n = daylength (h)
#   daily boundary D = [(24 - pd) EMC + (2.7 pd + 76) pd] / 24, pd = rain duration (h)
#   7-day mean boundary Dbar(t) = mean of D over days t-6 .. t
#   MC(t) = MC(t-7) + (Dbar(t) - MC(t-7)) (1 - 0.82 exp(-168 / 1000))
#   (INT-GTR-169 eq. 42 prints the bracket as (MC0 - D), a typo: the 100-h analogue, Fosberg
#   et al. and PNW-RP-336 all have (D - MC0).)
#   Rain duration when not observed (Bradshaw eq. 39): pd = round((P_in + 0.02) / rate), at
#   most 8 h; rate 0.25 in/h for NFDRS climate classes 1-2 (arid, semiarid, subhumid with dry
#   summers), 0.05 in/h for classes 3-4.
# PRISM has no humidity: vapour pressure ea = es(Tmax) - VPDmax (PRISM VPDmax is es(Tmax) minus
# the vapour pressure of the mean dewpoint), RHmin = ea / es(Tmax), RHmax = ea / es(Tmin)
# (constant vapour pressure through the day), capped at 100 %.

fdm_es_kpa <- function(t_c) 0.6108 * exp(17.27 * t_c / (t_c + 237.3))

# EMC (%) from temperature (F) and RH (%); keeps matrix dimensions.
fdm_emc <- function(Tf, H) {
  E <- 21.0606 + 0.005565 * H^2 - 0.00035 * H * Tf - 0.483199 * H
  i <- H <= 50; E[i] <- (2.22749 + 0.160107 * H - 0.014784 * Tf)[i]
  i <- H <= 10; E[i] <- (0.03229 + 0.281073 * H - 0.000578 * H * Tf)[i]
  E
}

# Relative humidity extremes (%) from PRISM tmin, tmax (C) and vpdmax (hPa).
fdm_rh <- function(TN, TX, VPD_hpa, ea_min_kpa = 0.01) {
  ea <- fdm_es_kpa(TX) - VPD_hpa / 10
  ea[ea < ea_min_kpa] <- ea_min_kpa
  rh_min <- 100 * ea / fdm_es_kpa(TX); rh_max <- 100 * ea / fdm_es_kpa(TN)
  rh_min[rh_min > 100] <- 100; rh_max[rh_max > 100] <- 100
  list(rh_min = rh_min, rh_max = rh_max)
}

# Rain duration (h) from daily precipitation (mm); rate_in_h per cell (vector) or scalar.
fdm_rain_hours <- function(P, rate_in_h, max_h = 8) {
  pin <- P / 25.4
  h <- floor((pin + 0.02) / rate_in_h + 0.5)          # round half up (R's round() is half-even)
  h[!is.na(P) & P <= 0] <- 0
  h[h > max_h] <- max_h
  h
}

# Daily boundary condition D (%), cells x days. photo_h = daylength in hours (cells x days).
fdm_boundary <- function(P, TN, TX, VPD_hpa, photo_h, rate_in_h, max_h = 8) {
  rh <- fdm_rh(TN, TX, VPD_hpa)
  e_day <- fdm_emc(TX * 9 / 5 + 32, rh$rh_min)
  e_night <- fdm_emc(TN * 9 / 5 + 32, rh$rh_max)
  emc <- (photo_h * e_day + (24 - photo_h) * e_night) / 24
  pd <- fdm_rain_hours(P, rate_in_h, max_h)
  ((24 - pd) * emc + (2.7 * pd + 76) * pd) / 24
}

# 1000-h moisture from the daily boundary D (cells x days). The first 7 days start at their
# own 7-day mean boundary (equilibrium start); 64 runs a spin-up year first.
fdm_1000h <- function(D, fr = 1 - 0.82 * exp(-168 / 1000)) {
  n <- ncol(D)
  C <- cbind(0, D); for (j in 2:ncol(C)) C[, j] <- C[, j - 1] + C[, j]
  idx <- seq_len(n); lo <- pmax(0L, idx - 7L)
  Db <- sweep(C[, idx + 1L, drop = FALSE] - C[, lo + 1L, drop = FALSE], 2, idx - lo, "/")
  M <- Db
  if (n > 7) for (t in 8:n) M[, t] <- M[, t - 7] + (Db[, t] - M[, t - 7]) * fr
  M
}

# Degree-day snowpack (mm SWE): precipitation falls as snow when mean temperature <= t_snow;
# melt = ddf x (Tmean - t_melt)+ per day. Starts snow-free.
fdm_snow <- function(P, TN, TX, t_snow = 1, t_melt = 0, ddf = 3) {
  TM <- (TN + TX) / 2
  S <- matrix(0, nrow(P), ncol(P)); s <- numeric(nrow(P))
  for (t in seq_len(ncol(P))) {
    p <- P[, t]; p[is.na(p)] <- 0; tm <- TM[, t]; tm[is.na(tm)] <- 0
    s <- s + ifelse(tm <= t_snow, p, 0)
    s <- s - pmin(s, ddf * pmax(tm - t_melt, 0))
    S[, t] <- s
  }
  S
}

# Snow-off day of year per cell for one year: the day after the last day from Jan 1 to
# `by` with SWE > swe_min; NA when the cell had no such day. cols = that year's columns.
fdm_snowoff <- function(S, dates, year, swe_min = 5, by = "07-31") {
  j <- which(dates >= as.Date(sprintf("%d-01-01", year)) & dates <= as.Date(sprintf("%d-%s", year, by)))
  W <- S[, j, drop = FALSE] > swe_min
  last <- rep(0L, nrow(W))
  for (k in seq_along(j)) last[W[, k]] <- k
  out <- ifelse(last > 0, as.integer(format(dates[j][1], "%j")) + last, NA_integer_)  # day after last snowy day
  as.integer(out)
}

# NFDRS climate-class rain rate per cell from aridity: 0.25 in/h (classes 1-2: arid, semiarid,
# dry summers) where Apr-Sep P/PET < ai_warm_lt (or annual P/PET < ai_ann_lt, off by default),
# else 0.05 in/h (classes 3-4). Decided Oct 18: warm season only (an annual cutoff of 0.8 split
# the upper Midwest cell by cell, its median being 0.81).
fdm_rain_rate <- function(P, PET, dates, ai_ann_lt = -Inf, ai_warm_lt = 0.5, dry = 0.25, wet = 0.05) {
  m <- as.integer(format(dates, "%m")); warm <- m >= 4 & m <= 9
  ai <- rowSums(P, na.rm = TRUE) / pmax(rowSums(PET, na.rm = TRUE), 1e-6)
  aw <- rowSums(P[, warm, drop = FALSE], na.rm = TRUE) / pmax(rowSums(PET[, warm, drop = FALSE], na.rm = TRUE), 1e-6)
  list(rate = ifelse(ai < ai_ann_lt | aw < ai_warm_lt, dry, wet), ai_ann = ai, ai_warm = aw)
}
