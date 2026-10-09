# GSI-per-region model test (scripts/34_gsi_model.R): fit GSI ramp thresholds to satellite
# green-up and curing dates, per region, CONUS-wide, or not at all (defaults), and score each on
# held-out years. Works on a sample of cells with their full daily series (cells x days).
# Moisture options: "precip" (28-day precipitation sum ramp), "kbdi" (KBDI ramp, decreasing),
# "none". Each has its own pair of fitted moisture thresholds (none: no moisture ramp).

gm_bounds_base <- list(tmin_lo = c(-8, 5), tmin_w = c(1, 15), vpd_lo = c(0, 4000), vpd_w = c(200, 6000),
                       photo_lo_h = c(8, 13))
gm_bounds_moist <- list(precip = list(m_lo = c(0, 40), m_w = c(1, 60)),
                        kbdi = list(m_lo = c(0, 700), m_w = c(50, 800)),
                        none = list())
gm_bounds <- function(moisture) c(gm_bounds_base, gm_bounds_moist[[moisture]])

gm_decode <- function(th, moisture) {
  b <- gm_bounds(moisture)
  v <- vapply(seq_along(b), function(i) b[[i]][1] + (b[[i]][2] - b[[i]][1]) / (1 + exp(-th[i])), 0)
  names(v) <- names(b)
  par <- list(moisture = moisture,
              tmin = c(v[["tmin_lo"]], v[["tmin_lo"]] + v[["tmin_w"]]),
              vpd = c(v[["vpd_lo"]], v[["vpd_lo"]] + v[["vpd_w"]]),
              photo_h = c(v[["photo_lo_h"]], v[["photo_lo_h"]] + 1))
  if (moisture != "none") par$moist <- c(v[["m_lo"]], v[["m_lo"]] + v[["m_w"]])
  par
}
gm_encode <- function(par) {
  b <- gm_bounds(par$moisture)
  v <- c(tmin_lo = par$tmin[1], tmin_w = diff(par$tmin), vpd_lo = par$vpd[1], vpd_w = diff(par$vpd),
         photo_lo_h = par$photo_h[1])
  if (par$moisture != "none") v <- c(v, m_lo = par$moist[1], m_w = diff(par$moist))
  vapply(seq_along(b), function(i) {
    p <- (v[[names(b)[i]]] - b[[i]][1]) / (b[[i]][2] - b[[i]][1]); p <- min(max(p, 1e-4), 1 - 1e-4)
    log(p / (1 - p))
  }, 0)
}
gm_par_vec <- function(par) c(tmin_lo = par$tmin[1], tmin_hi = par$tmin[2], vpd_lo = par$vpd[1],
                              vpd_hi = par$vpd[2], photo_lo_h = par$photo_h[1], photo_hi_h = par$photo_h[2],
                              moist_lo = if (is.null(par$moist)) NA else par$moist[1],
                              moist_hi = if (is.null(par$moist)) NA else par$moist[2])

# Smoothed GSI (cells x days). drv: list(tmin, vpd_pa, photo_s, prec28 = 28-day precipitation
# sum, kbdi), each cells x days (full series). par$vpd = NULL: no VPD ramp.
gm_gsi <- function(drv, par, smooth) {
  ig <- m_ramp(drv$tmin, par$tmin[1], par$tmin[2]) *
    m_ramp(drv$photo_s, par$photo_h[1] * 3600, par$photo_h[2] * 3600)
  if (!is.null(par$vpd)) ig <- ig * m_ramp(drv$vpd_pa, par$vpd[1], par$vpd[2], decreasing = TRUE)
  if (identical(par$moisture, "precip")) ig <- ig * m_ramp(drv$prec28, par$moist[1], par$moist[2])
  if (identical(par$moisture, "kbdi")) ig <- ig * m_ramp(drv$kbdi, par$moist[1], par$moist[2], decreasing = TRUE)
  m_roll_mean(ig, smooth)
}

# GSI sos20 / eos50 per year (cells x years), each year on its own window Jan 1 (Y-1) .. Mar 31
# (Y+1) as in 12_gsi.R.
gm_seasons <- function(G, dates, years, g) {
  sos <- eos <- matrix(NA_real_, nrow(G), length(years))
  for (j in seq_along(years)) {
    Y <- years[j]
    w <- which(dates >= as.Date(sprintf("%d-01-01", Y - 1)) & dates <= as.Date(sprintf("%d-03-31", Y + 1)))
    dd <- as.numeric(dates[w] - as.Date(sprintf("%d-01-01", Y))) + 1
    iy <- which(format(dates[w], "%Y") == as.character(Y))
    ss <- gsi_season_from_G(G[, w, drop = FALSE], dd, iy, g)
    sos[, j] <- ss[, "sos20"]; eos[, j] <- ss[, "eos50"]
  }
  list(sos = sos, eos = eos)
}

# Loss: mean absolute error (days) over green-up and curing, each error capped at `cap` days
# (a satellite season that jumps to another time of year, e.g. spring vs monsoon, would
# otherwise dominate); a missing GSI date where the satellite has one counts `penalty` days.
gm_loss <- function(pred, obs, penalty, cap = 90) {
  e <- c(abs(pred$sos - obs$sos)[!is.na(obs$sos)], abs(pred$eos - obs$cure)[!is.na(obs$cure)])
  e[is.na(e)] <- penalty
  mean(pmin(e, cap))
}

# Fit thresholds for one moisture option on the given cells and years: Nelder-Mead from the
# start values, then `restarts` more searches from the best point so far (a restart re-forms
# the simplex, which often gets Nelder-Mead out of a stall).
gm_fit <- function(drv, obs, dates, years, g, gmc, start) {
  mo <- start$moisture
  fn <- function(th) {
    pr <- gm_seasons(gm_gsi(drv, gm_decode(th, mo), gmc$smooth), dates, years, g)
    gm_loss(pr, obs, gmc$penalty_days, gmc$cap_days)
  }
  th <- gm_encode(start); best <- Inf; evals <- 0
  for (r in 0:(gmc$restarts %||% 0)) {
    o <- stats::optim(th, fn, method = "Nelder-Mead", control = list(maxit = gmc$maxit))
    evals <- evals + o$counts[[1]]
    if (o$value < best - 1e-6) { best <- o$value; th <- o$par } else break
  }
  list(par = gm_decode(th, mo), loss = best, evals = evals)
}

# Scores for predictions vs observations on test cell-years: MAE (days; errors capped at `cap`,
# missing GSI dates counted as the penalty), median absolute error, share within 30 days, share
# of observed dates with a prediction, and the correlation of yearly anomalies (each cell's
# dates minus that cell's mean over the years; anom = FALSE for predictions with no year-to-year
# variation, e.g. climatology).
gm_score <- function(pred, obs, penalty, cap = 90, anom = TRUE) {
  one <- function(p, o) {
    has <- !is.na(o); if (!any(has)) return(c(mae = NA, mdae = NA, within30 = NA, coverage = NA, anom_r = NA))
    e <- abs(p - o)[has]; cov <- mean(!is.na(e)); e[is.na(e)] <- penalty
    mdae <- stats::median(e); w30 <- mean(e <= 30); e <- pmin(e, cap)
    ok <- !is.na(p) & !is.na(o)
    pa <- p; oa <- o; pa[!ok] <- NA; oa[!ok] <- NA
    pa <- pa - rowMeans(pa, na.rm = TRUE); oa <- oa - rowMeans(oa, na.rm = TRUE)
    u <- !is.na(pa) & !is.na(oa)
    r <- if (anom && sum(u) > 10 && stats::sd(pa[u]) > 0) stats::cor(pa[u], oa[u]) else NA
    c(mae = mean(e), mdae = mdae, within30 = w30, coverage = cov, anom_r = r)
  }
  s <- one(pred$sos, obs$sos); c_ <- one(pred$eos, obs$cure)
  out <- c(s, c_); names(out) <- c(paste0("sos_", names(s)), paste0("cure_", names(c_)))
  out
}
