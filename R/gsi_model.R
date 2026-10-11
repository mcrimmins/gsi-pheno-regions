# GSI-per-region model test (scripts/34_gsi_model.R): fit GSI ramp thresholds to satellite
# green-up and curing dates, per region, CONUS-wide, or not at all (defaults), and score each on
# held-out years. Works on a sample of cells with their full daily series (cells x days).
# Moisture options (gm_moist) each have their own pair of fitted thresholds; crossing forms
# (gm_cross) set how green-up and curing dates are read from the smoothed index.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Moisture options: driver (cells x days, in drv), whether the ramp decreases, search bounds.
#   precip / precip60 / precip90  trailing 28 / 60 / 90-day precipitation sum (mm)
#   kbdi                          Keetch-Byram drought index
#   kbdi30                        30-day trailing mean of KBDI (longer memory)
#   none                          no moisture ramp
gm_moist <- list(
  precip   = list(drv = "prec28", decreasing = FALSE, m_lo = c(0, 40),  m_w = c(1, 60)),
  precip60 = list(drv = "prec60", decreasing = FALSE, m_lo = c(0, 90),  m_w = c(1, 130)),
  precip90 = list(drv = "prec90", decreasing = FALSE, m_lo = c(0, 130), m_w = c(1, 200)),
  kbdi     = list(drv = "kbdi",   decreasing = TRUE,  m_lo = c(0, 700), m_w = c(50, 800)),
  kbdi30   = list(drv = "kbdi30", decreasing = TRUE,  m_lo = c(0, 700), m_w = c(50, 800)),
  # soil-water bucket (Oct 11): a different model form, not a GSI moisture ramp. Daily bucket of
  # capacity wmax (mm, fitted) filled by rain and drained by evapotranspiration at PET x the
  # bucket's relative fill; index = ramp(relative fill) x ramp(Tmin), no VPD or daylength.
  swb      = list(drv = "swb",    decreasing = FALSE, m_lo = c(0, 0.8), m_w = c(0.02, 1)),
  none     = NULL)
# How season dates are read from the smoothed GSI (test of curing forms, 35_gsi_curing.R):
#   rel         relative crossings around the yearly peak: green-up at base + 0.2 x amplitude,
#               curing at base + 0.5 x amplitude (as the satellite curves; v1-v2 default)
#   rel_split   relative, with the green-up and curing fractions fitted separately
#   abs_nfdrs4  absolute GSI levels as NFDRS4: green-up when GSI rises past GU (fitted),
#               50 % cured when it falls to (1 + GU) / 2 (NFDRS4 Cure(): cured share =
#               (1 - GSI) / (1 - GU)); the annual-herb ratchet only acts below that level
#   abs_split   absolute, green-up and curing levels fitted separately
#   rel_anchor  as rel_split, but the season's peak is searched within +/- anchor_halfwidth_days
#               of the cell's usual mid-season (midway between its typical green-up and curing
#               dates in the fitting years), not as the year's highest value: for two-season
#               climates (monsoon), the index's season nearest the cell's usual one is read
gm_cross <- list(rel = list(), rel_split = list(up = c(0.05, 0.8), down = c(0.05, 0.95)),
                 abs_nfdrs4 = list(gu = c(0.05, 0.9)), abs_split = list(up = c(0.02, 0.9), down = c(0.02, 0.95)),
                 rel_anchor = list(up = c(0.05, 0.8), down = c(0.05, 0.95)))
gm_bounds_base <- list(tmin_lo = c(-8, 5), tmin_w = c(1, 15), vpd_lo = c(0, 4000), vpd_w = c(200, 6000),
                       photo_lo_h = c(8, 13))
gm_bounds_swb <- list(tmin_lo = c(-8, 5), tmin_w = c(1, 15), wmax = c(10, 400))
gm_bounds <- function(moisture, cross = "rel") {
  mo <- gm_moist[[moisture]]
  if (!moisture %in% names(gm_moist)) stop("unknown moisture option ", moisture)
  if (!cross %in% names(gm_cross)) stop("unknown crossing form ", cross)
  cr <- gm_cross[[cross]]; if (length(cr)) names(cr) <- paste0("x_", names(cr))
  base <- if (moisture == "swb") gm_bounds_swb else gm_bounds_base
  c(base, if (!is.null(mo)) mo[c("m_lo", "m_w")], cr)
}

gm_decode <- function(th, moisture, cross = "rel") {
  b <- gm_bounds(moisture, cross)
  v <- vapply(seq_along(b), function(i) b[[i]][1] + (b[[i]][2] - b[[i]][1]) / (1 + exp(-th[i])), 0)
  names(v) <- names(b)
  par <- list(moisture = moisture, cross = cross, tmin = c(v[["tmin_lo"]], v[["tmin_lo"]] + v[["tmin_w"]]))
  if (moisture == "swb") par$wmax <- v[["wmax"]] else {
    par$vpd <- c(v[["vpd_lo"]], v[["vpd_lo"]] + v[["vpd_w"]])
    par$photo_h <- c(v[["photo_lo_h"]], v[["photo_lo_h"]] + 1)
  }
  if (moisture != "none") par$moist <- c(v[["m_lo"]], v[["m_lo"]] + v[["m_w"]])
  if (cross != "rel") par$x <- v[grep("^x_", names(v))] |> stats::setNames(sub("^x_", "", grep("^x_", names(v), value = TRUE)))
  par
}
gm_encode <- function(par) {
  cross <- par$cross %||% "rel"
  b <- gm_bounds(par$moisture, cross)
  v <- if (par$moisture == "swb") c(tmin_lo = par$tmin[1], tmin_w = diff(par$tmin), wmax = par$wmax) else
    c(tmin_lo = par$tmin[1], tmin_w = diff(par$tmin), vpd_lo = par$vpd[1], vpd_w = diff(par$vpd),
      photo_lo_h = par$photo_h[1])
  if (par$moisture != "none") v <- c(v, m_lo = par$moist[1], m_w = diff(par$moist))
  if (cross != "rel") v <- c(v, stats::setNames(unlist(par$x), paste0("x_", names(par$x))))
  vapply(seq_along(b), function(i) {
    p <- (v[[names(b)[i]]] - b[[i]][1]) / (b[[i]][2] - b[[i]][1]); p <- min(max(p, 1e-4), 1 - 1e-4)
    log(p / (1 - p))
  }, 0)
}
gm_par_vec <- function(par) {
  x <- par$x; gx <- function(nm) if (nm %in% names(x)) unname(x[[nm]]) else NA
  g1 <- function(v, i) if (is.null(v)) NA else v[i]
  c(tmin_lo = par$tmin[1], tmin_hi = par$tmin[2], vpd_lo = g1(par$vpd, 1),
    vpd_hi = g1(par$vpd, 2), photo_lo_h = g1(par$photo_h, 1), photo_hi_h = g1(par$photo_h, 2),
    wmax = g1(par$wmax, 1),
    moist_lo = if (is.null(par$moist)) NA else par$moist[1],
    moist_hi = if (is.null(par$moist)) NA else par$moist[2],
    x_up = gx("up"), x_down = gx("down"), x_gu = gx("gu"))
}

# Smoothed GSI (cells x days). drv: list(tmin, vpd_pa, photo_s, plus the moisture drivers
# named in gm_moist: prec28, prec60, prec90, kbdi, kbdi30), each cells x days (full series).
# par$vpd = NULL: no VPD ramp.
# Relative fill of a daily soil-water bucket (cells x days): W += P - PET x W / wmax, kept
# within [0, wmax] (rain beyond capacity runs off); starts half full. P, PET in mm/day.
# A C++ version (gm_bucket_init(), ~100x faster) is used when it compiled; same result.
gm_bucket_cpp <- NULL
gm_bucket_init <- function() {
  if (!is.null(gm_bucket_cpp) || !requireNamespace("Rcpp", quietly = TRUE)) return(invisible(!is.null(gm_bucket_cpp)))
  f <- tryCatch(Rcpp::cppFunction('
    NumericMatrix bucket_cpp(NumericMatrix P, NumericMatrix PET, double wmax) {
      int m = P.nrow(), n = P.ncol(); NumericMatrix R(m, n);
      for (int i = 0; i < m; i++) {
        double W = wmax / 2;
        for (int t = 0; t < n; t++) {
          double p = P(i, t), e = PET(i, t);
          if (ISNAN(p)) p = 0; if (ISNAN(e)) e = 0;
          W = W + p - e * W / wmax;
          if (W < 0) W = 0; if (W > wmax) W = wmax;
          R(i, t) = W / wmax;
        }
      }
      return R;
    }'), error = function(e) NULL)
  if (!is.null(f)) assign("gm_bucket_cpp", f, envir = environment(gm_bucket_init))
  invisible(!is.null(f))
}
m_bucket <- function(P, PET, wmax) {
  if (!is.null(gm_bucket_cpp)) return(gm_bucket_cpp(P, PET, wmax))
  n <- ncol(P); W <- rep(wmax / 2, nrow(P)); R <- matrix(0, nrow(P), n)
  P[is.na(P)] <- 0; PET[is.na(PET)] <- 0
  for (t in seq_len(n)) {
    W <- W + P[, t] - PET[, t] * W / wmax
    W[W < 0] <- 0; W[W > wmax] <- wmax
    R[, t] <- W
  }
  R / wmax
}

gm_gsi <- function(drv, par, smooth) {
  ig <- m_ramp(drv$tmin, par$tmin[1], par$tmin[2])
  if (identical(par$moisture, "swb")) {                       # bucket model: no VPD, no daylength
    ig <- ig * m_ramp(m_bucket(drv$ppt, drv$pet, par$wmax), par$moist[1], par$moist[2])
    return(m_roll_mean(ig, smooth))
  }
  if (!is.null(par$photo_h)) ig <- ig * m_ramp(drv$photo_s, par$photo_h[1] * 3600, par$photo_h[2] * 3600)
  if (!is.null(par$vpd)) ig <- ig * m_ramp(drv$vpd_pa, par$vpd[1], par$vpd[2], decreasing = TRUE)
  mo <- gm_moist[[par$moisture]]
  if (!is.null(mo)) {
    x <- drv[[mo$drv]]; if (is.null(x)) stop("driver ", mo$drv, " missing for moisture option ", par$moisture)
    ig <- ig * m_ramp(x, par$moist[1], par$moist[2], decreasing = mo$decreasing)
  }
  m_roll_mean(ig, smooth)
}

# Season dates from the smoothed GSI around the year-Y peak (generalizes gsi_season_from_G:
# cross = "rel" with up 0.2 / down 0.5 gives the same sos20 / eos50). Relative: thresholds
# base + frac x amplitude, season if amplitude >= min_amp. Absolute: thresholds are GSI levels
# (GSImax = 1); no season if the peak stays below the green-up level, no curing date if it
# stays below the curing level.
gm_season_dates <- function(G, dd, iy, g, cross = "rel", x = NULL, anchor = NULL) {
  m <- nrow(G); n <- ncol(G); rows <- seq_len(m)
  if (identical(cross, "rel_anchor") && !is.null(anchor)) {
    # peak within +/- aw days of the cell's usual mid-season; cells without one use year Y
    aw <- g$anchor_halfwidth_days %||% 75
    M <- abs(outer(anchor, dd, "-")) <= aw
    nob <- is.na(anchor) | rowSums(M, na.rm = TRUE) == 0
    M[is.na(M)] <- FALSE
    if (any(nob)) M[nob, ] <- matrix(seq_len(n) %in% iy, sum(nob), n, byrow = TRUE)
    Gm <- G; Gm[!M] <- -Inf
    pmx <- apply(Gm, 1, max)
    pcol <- first_true(Gm >= pmx - 1e-6)
  } else {
    Gy <- G[, iy, drop = FALSE]
    pmx <- apply(Gy, 1, max)
    pcol <- iy[first_true(Gy >= pmx - 1e-6)]
  }
  hw <- g$season_halfwidth_days
  if (startsWith(cross, "rel")) {
    base <- pmx
    for (k in seq_len(hw)) for (cc in list(pcol - k, pcol + k)) {
      ok <- cc >= 1 & cc <= n
      base[ok] <- pmin(base[ok], G[cbind(rows[ok], cc[ok])])
    }
    amp <- pmx - base
    fu <- if (cross == "rel") 0.2 else x[["up"]]; fd <- if (cross == "rel") 0.5 else x[["down"]]
    t_up <- base + fu * amp; t_dn <- base + fd * amp
    seas_up <- seas_dn <- amp >= g$min_amp
  } else {
    lu <- if (cross == "abs_nfdrs4") x[["gu"]] else x[["up"]]
    ld <- if (cross == "abs_nfdrs4") (1 + x[["gu"]]) / 2 else x[["down"]]
    t_up <- rep(lu, m); t_dn <- rep(ld, m)
    seas_up <- pmx >= lu; seas_dn <- seas_up & pmx > ld
  }
  sos <- eos <- rep(NA_real_, m)
  for (k in seq_len(hw)) {
    cb <- pcol - k; ok <- is.na(sos) & cb >= 1
    hit <- ok; hit[ok] <- G[cbind(rows[ok], cb[ok])] < t_up[ok]
    sos[hit] <- dd[cb[hit]] + 1
    ca <- pcol + k; ok <- is.na(eos) & ca <= n
    hit <- ok; hit[ok] <- G[cbind(rows[ok], ca[ok])] < t_dn[ok]
    eos[hit] <- dd[ca[hit]]
  }
  cbind(sos = ifelse(seas_up, sos, NA), eos = ifelse(seas_dn, eos, NA))
}

# Season dates per year (cells x years), each year on its own window Jan 1 (Y-1) .. Mar 31
# (Y+1) as in 12_gsi.R. Default (rel): the GSI's sos20 / eos50.
# anchor (optional, rel_anchor): cells-long vector of the usual mid-season day (relative to
# Jan 1 of each year), used for every year.
gm_seasons <- function(G, dates, years, g, cross = "rel", x = NULL, anchor = NULL) {
  sos <- eos <- matrix(NA_real_, nrow(G), length(years))
  for (j in seq_along(years)) {
    Y <- years[j]
    w <- which(dates >= as.Date(sprintf("%d-01-01", Y - 1)) & dates <= as.Date(sprintf("%d-03-31", Y + 1)))
    dd <- as.numeric(dates[w] - as.Date(sprintf("%d-01-01", Y))) + 1
    iy <- which(format(dates[w], "%Y") == as.character(Y))
    ss <- gm_season_dates(G[, w, drop = FALSE], dd, iy, g, cross, x, anchor)
    sos[, j] <- ss[, "sos"]; eos[, j] <- ss[, "eos"]
  }
  list(sos = sos, eos = eos)
}
# Season dates for a parameter set (GSI + the set's own crossing form).
gm_predict <- function(drv, par, smooth, dates, years, g, anchor = NULL)
  gm_seasons(gm_gsi(drv, par, smooth), dates, years, g, par$cross %||% "rel", par$x, anchor)
# Usual mid-season per cell from observed dates (cells x years): midway between the median
# green-up and the median curing date.
gm_anchor <- function(sos, cure) {
  ms <- apply(sos, 1, stats::median, na.rm = TRUE); mc <- apply(cure, 1, stats::median, na.rm = TRUE)
  ifelse(is.na(mc), ms, ifelse(is.na(ms), mc, (ms + mc) / 2))
}

# Loss: mean absolute error (days) over green-up and curing, each error capped at `cap` days
# (a satellite season that jumps to another time of year, e.g. spring vs monsoon, would
# otherwise dominate); a missing GSI date where the satellite has one counts `penalty` days.
# obs$fill_sos / obs$fill_cure (optional, cells x years): dates used where the GSI gives none
# (e.g. the cell's typical date), so that skipping a date never scores better than the
# fallback; the penalty then only applies where the fallback is missing too.
gm_fill <- function(p, f) { if (is.null(f)) return(p); i <- is.na(p); p[i] <- f[i]; p }
gm_loss <- function(pred, obs, penalty, cap = 90) {
  pred$sos <- gm_fill(pred$sos, obs$fill_sos); pred$eos <- gm_fill(pred$eos, obs$fill_cure)
  e <- c(abs(pred$sos - obs$sos)[!is.na(obs$sos)], abs(pred$eos - obs$cure)[!is.na(obs$cure)])
  e[is.na(e)] <- penalty
  mean(pmin(e, cap))
}

# Fit thresholds for one moisture option on the given cells and years: Nelder-Mead from the
# start values, then `restarts` more searches from the best point so far (a restart re-forms
# the simplex, which often gets Nelder-Mead out of a stall).
gm_fit <- function(drv, obs, dates, years, g, gmc, start) {
  mo <- start$moisture; cr <- start$cross %||% "rel"
  anc <- if (cr == "rel_anchor") gm_anchor(obs$sos, obs$cure) else NULL     # fitting years only
  fn <- function(th) {
    pr <- gm_predict(drv, gm_decode(th, mo, cr), gmc$smooth, dates, years, g, anc)
    gm_loss(pr, obs, gmc$penalty_days, gmc$cap_days)
  }
  th <- gm_encode(start); best <- Inf; evals <- 0
  for (r in 0:(gmc$restarts %||% 0)) {
    o <- stats::optim(th, fn, method = "Nelder-Mead", control = list(maxit = gmc$maxit))
    evals <- evals + o$counts[[1]]
    if (o$value < best - 1e-6) { best <- o$value; th <- o$par } else break
  }
  list(par = gm_decode(th, mo, cr), loss = best, evals = evals)
}

# Scores for predictions vs observations on test cell-years: MAE (days; errors capped at `cap`,
# missing GSI dates counted as the penalty), median absolute error, share within 30 days, share
# of observed dates with a prediction, and the correlation of yearly anomalies (each cell's
# dates minus that cell's mean over the years; anom = FALSE for predictions with no year-to-year
# variation, e.g. climatology).
# fill (optional): list(sos, eos) of fallback dates (cells x years), used where the model gives
# no date (coverage is counted before the fill). With fill, also the paired comparison on the
# cell-years where the model gives a date: model MAE and the fallback's MAE there.
gm_score <- function(pred, obs, penalty, cap = 90, anom = TRUE, fill = NULL) {
  one <- function(p, o, f) {
    has <- !is.na(o)
    if (!any(has)) return(c(mae = NA, mdae = NA, within30 = NA, coverage = NA, anom_r = NA, paired_mae = NA, paired_fill_mae = NA))
    cov <- mean(!is.na(p[has]))
    pm <- pf <- NA
    if (!is.null(f)) {
      u <- has & !is.na(p) & !is.na(f)
      if (any(u)) { pm <- mean(pmin(abs(p - o)[u], cap)); pf <- mean(pmin(abs(f - o)[u], cap)) }
    }
    p0 <- p; p <- gm_fill(p, f)
    e <- abs(p - o)[has]; e[is.na(e)] <- penalty
    mdae <- stats::median(e); w30 <- mean(e <= 30); e <- pmin(e, cap)
    ok <- !is.na(p0) & !is.na(o)                 # anomaly skill on the model's own dates
    pa <- p0; oa <- o; pa[!ok] <- NA; oa[!ok] <- NA
    pa <- pa - rowMeans(pa, na.rm = TRUE); oa <- oa - rowMeans(oa, na.rm = TRUE)
    u <- !is.na(pa) & !is.na(oa)
    r <- if (anom && sum(u) > 10 && stats::sd(pa[u]) > 0) stats::cor(pa[u], oa[u]) else NA
    c(mae = mean(e), mdae = mdae, within30 = w30, coverage = cov, anom_r = r, paired_mae = pm, paired_fill_mae = pf)
  }
  s <- one(pred$sos, obs$sos, fill$sos); c_ <- one(pred$eos, obs$cure, fill$eos)
  out <- c(s, c_); names(out) <- c(paste0("sos_", names(s)), paste0("cure_", names(c_)))
  out
}
