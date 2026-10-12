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
