# Candidate curing cues per cell and year (37_curing_cues.R; worked examples). Each cue is
# computed from the daily series for year Y and expressed as a day of year relative to Jan 1 of
# Y, a temperature (C), a precipitation total (mm) or a count of days. Needs R/features.R
# (first_true, last_true) and R/block2.R (onset_doy, fwd_sum, run_lengths).
#   onset          warm-season rain onset, Block 2 rule (May 1 - Sep 30; none -> Oct 1)
#   rain_end       last day Jul 1 - Oct 31 that starts a 10-day total >= 20 mm (summer rains end)
#   spring_rain    last day Jan 1 - Jun 30 with >= 10 mm (last soaking spring rain)
#   fall_frz0      first day from Jul 1 with tmin <= 0 C (through Mar 31 of Y+1)
#   fall_frz2      the same for tmin <= -2.2 C (hard freeze)
#   tmax_spring    mean daily maximum temperature, Mar - May
#   tmax_summer    mean daily maximum temperature, Jun - Aug
#   p_cool         precipitation Oct (Y-1) - Mar (Y)
#   p_spring       precipitation Apr - Jun
#   p_summer       precipitation Jul - Sep
#   p_fall         precipitation Sep - Nov
#   dry_spell      longest run of days under 1 mm, Jun - Sep
# Daylength is not a cue: it is the same every year, so it cannot explain early and late years.

cue_names <- c("onset", "rain_end", "spring_rain", "fall_frz0", "fall_frz2", "tmax_spring", "tmax_summer",
               "p_cool", "p_spring", "p_summer", "p_fall", "dry_spell")

# P, TN, TX: cells x days; dates: Date per column. b2: config features: block2 (onset rule).
cue_year <- function(P, TN, TX, dates, Y, b2) {
  m <- nrow(P); out <- matrix(NA_real_, m, length(cue_names), dimnames = list(NULL, cue_names))
  d0 <- as.Date(sprintf("%d-01-01", Y))
  rel <- as.numeric(dates - d0) + 1                         # day relative to Jan 1 of Y
  iy <- which(format(dates, "%Y") == as.character(Y))
  if (length(iy) < 365) return(out)
  between <- function(a, b) which(dates >= as.Date(sprintf(a, Y)) & dates <= as.Date(sprintf(b, Y)))
  Py <- P[, iy, drop = FALSE]
  out[, "onset"] <- onset_doy(Py, seq_along(iy), b2)
  # end of summer rains: last start of a 10-day total >= 20 mm, Jul 1 - Oct 31
  w <- between("%d-07-01", "%d-10-31"); S <- fwd_sum(P[, c(w, max(w) + 1:9), drop = FALSE], 10, seq_along(w))
  k <- last_true(S >= 20 & P[, w, drop = FALSE] >= 1); out[, "rain_end"] <- ifelse(k > 0, rel[w][pmax(k, 1)], NA)
  w <- between("%d-01-01", "%d-06-30")
  k <- last_true(P[, w, drop = FALSE] >= 10); out[, "spring_rain"] <- ifelse(k > 0, rel[w][pmax(k, 1)], NA)
  w <- which(dates >= as.Date(sprintf("%d-07-01", Y)) & dates <= as.Date(sprintf("%d-03-31", Y + 1)))
  for (th in list(c("fall_frz0", 0), c("fall_frz2", -2.2))) {
    f <- first_true(TN[, w, drop = FALSE] <= as.numeric(th[2]))
    out[, th[1]] <- ifelse(f <= length(w), rel[w][pmin(f, length(w))], NA)
  }
  mm <- function(a, b) rowMeans(TX[, between(a, b), drop = FALSE], na.rm = TRUE)
  out[, "tmax_spring"] <- mm("%d-03-01", "%d-05-31"); out[, "tmax_summer"] <- mm("%d-06-01", "%d-08-31")
  ps <- function(idx) rowSums(P[, idx, drop = FALSE], na.rm = TRUE)
  wc <- which(dates >= as.Date(sprintf("%d-10-01", Y - 1)) & dates <= as.Date(sprintf("%d-03-31", Y)))
  out[, "p_cool"] <- if (length(wc) >= 180) ps(wc) else NA
  out[, "p_spring"] <- ps(between("%d-04-01", "%d-06-30")); out[, "p_summer"] <- ps(between("%d-07-01", "%d-09-30"))
  out[, "p_fall"] <- ps(between("%d-09-01", "%d-11-30"))
  w <- between("%d-06-01", "%d-09-30")
  out[, "dry_spell"] <- apply(run_lengths(P[, w, drop = FALSE] < 1), 1, max)
  out
}

# Anomalies: each cell's values minus its own mean over the years it has (cells x years).
cue_anom <- function(X) X - rowMeans(X, na.rm = TRUE)

# Screening scores for one target (cells x years anomalies `y`) and a list of cue anomaly
# matrices `C` (same shape), pooled over cells and years:
#   r     correlation of target and cue anomalies
#   cv_r2 share of the target's anomaly variance a small linear model of the cues predicts in
#         held-out years (fold A fits, fold B scores, and the reverse); cues picked forward on
#         the fitting fold (up to max_cues, each must raise the fitting-fold R^2 by min_gain)
cue_screen <- function(y, C, folds, max_cues = 3, min_gain = 0.01) {
  r <- vapply(C, function(x) { u <- !is.na(y) & !is.na(x); if (sum(u) < 30 || stats::sd(x[u]) == 0) NA else stats::cor(y[u], x[u]) }, 0)
  n <- vapply(C, function(x) sum(!is.na(y) & !is.na(x)), 0)
  sse <- sst <- 0; picked <- list()
  for (f in names(folds)) {
    tr <- folds[[f]]; te <- setdiff(seq_len(ncol(y)), tr)
    yt <- as.vector(y[, tr]); Xt <- vapply(C, function(x) as.vector(x[, tr]), numeric(length(yt)))
    sel <- character(); best <- 0
    for (step in seq_len(max_cues)) {
      cand <- setdiff(colnames(Xt), sel); gain <- c()
      for (cn in cand) {
        v <- c(sel, cn); u <- stats::complete.cases(yt, Xt[, v])
        if (sum(u) < 30) next
        fit <- stats::lm.fit(cbind(Xt[u, v, drop = FALSE]), yt[u])
        gain[cn] <- 1 - sum(fit$residuals^2) / sum(yt[u]^2)
      }
      gain <- gain[is.finite(gain)]
      if (!length(gain) || max(gain) - best < min_gain) break
      sel <- c(sel, names(which.max(gain))); best <- max(gain)
    }
    picked[[f]] <- sel
    ye <- as.vector(y[, te])
    if (!length(sel)) { u <- !is.na(ye); sse <- sse + sum(ye[u]^2); sst <- sst + sum(ye[u]^2); next }
    u <- stats::complete.cases(yt, Xt[, sel])
    b <- stats::lm.fit(Xt[u, sel, drop = FALSE], yt[u])$coefficients
    Xe <- vapply(C[sel], function(x) as.vector(x[, te]), numeric(length(ye)))
    if (is.null(dim(Xe))) Xe <- matrix(Xe, ncol = length(sel))
    u <- stats::complete.cases(ye, Xe)
    sse <- sse + sum((ye[u] - Xe[u, , drop = FALSE] %*% b)^2); sst <- sst + sum(ye[u]^2)
  }
  list(r = r, n = n, cv_r2 = if (sst > 0) 1 - sse / sst else NA,
       cues = vapply(picked, function(s) paste(s, collapse = "+"), ""))
}

# ---- Typical date plus a weather adjustment (38_curing_adjust.R) ------------------------------
# Lead-time cues: each cell's cues use only weather up to a cutoff day `cut` (day relative to
# Jan 1 of Y; the cell's typical curing date minus the lead time), so the adjustment could be
# made before the date is needed. All windows end at the cutoff:
#   p30, p90, p180   precipitation in the 30, 90, 180 days before (mm)
#   tx30, tx60       mean daily maximum temperature in the 30, 60 days before (C)
#   dry60            days under 1 mm in the 60 days before
#   since_rain       days since the last day with >= 10 mm (at most 180)
#   greenup          this year's satellite green-up date if it is before the cutoff, else NA
adj_cue_names <- c("p30", "p90", "p180", "tx30", "tx60", "dry60", "since_rain", "greenup")

adj_cues <- function(P, TX, dates, Y, cut, sos_y = NULL) {
  m <- nrow(P); out <- matrix(NA_real_, m, length(adj_cue_names), dimnames = list(NULL, adj_cue_names))
  ti <- match(as.Date(sprintf("%d-01-01", Y)) + round(cut) - 1, dates)
  ok <- which(!is.na(ti) & ti > 180)
  if (!length(ok)) return(out)
  lo <- min(ti[ok]) - 180; hi <- max(ti[ok]); cols <- lo:hi; tt <- ti[ok] - lo + 1   # position in block
  rr <- seq_along(tt)
  blk <- function(X) X[ok, cols, drop = FALSE]
  cs <- function(X) { X[is.na(X)] <- 0; if (length(cols) == 1) cbind(0, X) else cbind(0, t(apply(X, 1, cumsum))) }  # col j + 1 = sum 1..j
  wsum <- function(C, w) C[cbind(rr, tt + 1)] - C[cbind(rr, tt + 1 - w)]
  Pb <- blk(P); CP <- cs(Pb); CT <- cs(blk(TX)); CD <- cs((Pb < 1) * 1)
  out[ok, "p30"] <- wsum(CP, 30); out[ok, "p90"] <- wsum(CP, 90); out[ok, "p180"] <- wsum(CP, 180)
  out[ok, "tx30"] <- wsum(CT, 30) / 30; out[ok, "tx60"] <- wsum(CT, 60) / 60
  out[ok, "dry60"] <- wsum(CD, 60)
  since <- rep(180, length(tt)); found <- rep(FALSE, length(tt))
  for (k in 0:179) {
    h <- !found & Pb[cbind(rr, tt - k)] >= 10
    h[is.na(h)] <- FALSE; since[h] <- k; found <- found | h
  }
  out[ok, "since_rain"] <- since
  if (!is.null(sos_y)) out[ok, "greenup"] <- ifelse(!is.na(sos_y[ok]) & sos_y[ok] <= cut[ok], sos_y[ok], NA)
  out
}

# Fit and score "typical date + adjustment" for one region and one fold. y_tr / y_te: departures
# from the typical date (cells x years, NA where unknown; season jumps already NA in y_tr).
# X_tr / X_te: lists of cue anomaly matrices (same shapes; NA cue -> 0, i.e. "typical").
# Forward selection on the fitting years (no intercept, so no cue = the typical date), up to
# max_cues, each raising the fitting R^2 by min_gain. Returns the cues, coefficients and the
# predicted adjustment for the test years, clamped to +/- max_adj days.
adj_fit <- function(y_tr, X_tr, X_te, max_cues = 3, min_gain = 0.01, max_adj = 60) {
  yv <- as.vector(y_tr); u <- !is.na(yv)
  Xv <- vapply(X_tr, function(x) { v <- as.vector(x); v[is.na(v)] <- 0; v }, numeric(length(yv)))[u, , drop = FALSE]
  yv <- yv[u]; sel <- character(); best <- 0; sst <- sum(yv^2)
  if (length(yv) >= 30 && sst > 0) for (step in seq_len(max_cues)) {
    gain <- c()
    for (cn in setdiff(colnames(Xv), sel)) {
      v <- c(sel, cn); if (all(Xv[, cn] == 0)) next
      gain[cn] <- 1 - sum(stats::lm.fit(Xv[, v, drop = FALSE], yv)$residuals^2) / sst
    }
    gain <- gain[is.finite(gain)]
    if (!length(gain) || max(gain) - best < min_gain) break
    sel <- c(sel, names(which.max(gain))); best <- max(gain)
  }
  te_dim <- dim(X_te[[1]])
  if (!length(sel)) return(list(cues = character(), b = numeric(), r2_fit = 0, adj = matrix(0, te_dim[1], te_dim[2])))
  b <- stats::lm.fit(Xv[, sel, drop = FALSE], yv)$coefficients; b[is.na(b)] <- 0
  adj <- Reduce(`+`, lapply(sel, function(cn) { x <- X_te[[cn]]; x[is.na(x)] <- 0; b[[cn]] * x }))
  list(cues = sel, b = b, r2_fit = best, adj = pmax(pmin(adj, max_adj), -max_adj))
}
