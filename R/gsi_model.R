# GSI-per-region model test (scripts/34_gsi_model.R): fit GSI ramp thresholds to satellite
# green-up and curing dates, per region, CONUS-wide, or not at all (defaults), and score each on
# held-out years. Works on a sample of cells with their full daily series (cells x days).

# Fitted parameters and their bounds (logistic transform keeps the search inside them).
gm_bounds <- list(tmin_lo = c(-8, 5), tmin_w = c(1, 15), vpd_lo = c(0, 4000), vpd_w = c(200, 6000),
                  photo_lo_h = c(8, 13), prec_lo = c(0, 40), prec_w = c(1, 60))
gm_decode <- function(th) {
  b <- gm_bounds
  v <- vapply(seq_along(b), function(i) b[[i]][1] + (b[[i]][2] - b[[i]][1]) / (1 + exp(-th[i])), 0)
  names(v) <- names(b)
  list(tmin = c(v[["tmin_lo"]], v[["tmin_lo"]] + v[["tmin_w"]]),
       vpd = c(v[["vpd_lo"]], v[["vpd_lo"]] + v[["vpd_w"]]),
       photo_h = c(v[["photo_lo_h"]], v[["photo_lo_h"]] + 1),
       prec = c(v[["prec_lo"]], v[["prec_lo"]] + v[["prec_w"]]))
}
gm_encode <- function(par) {
  v <- c(tmin_lo = par$tmin[1], tmin_w = diff(par$tmin), vpd_lo = par$vpd[1], vpd_w = diff(par$vpd),
         photo_lo_h = par$photo_h[1], prec_lo = par$prec[1], prec_w = diff(par$prec))
  b <- gm_bounds
  vapply(seq_along(b), function(i) {
    p <- (v[[i]] - b[[i]][1]) / (b[[i]][2] - b[[i]][1]); p <- min(max(p, 1e-4), 1 - 1e-4)
    log(p / (1 - p))
  }, 0)
}
gm_par_vec <- function(par) c(tmin_lo = par$tmin[1], tmin_hi = par$tmin[2], vpd_lo = par$vpd[1],
                              vpd_hi = par$vpd[2], photo_lo_h = par$photo_h[1], photo_hi_h = par$photo_h[2],
                              prec_lo = par$prec[1], prec_hi = par$prec[2])

# Smoothed GSI (cells x days) for one parameter set. drv: list(tmin, vpd_pa, ppt, photo_s),
# each cells x days (full series). prec = NULL: no moisture ramp.
gm_gsi <- function(drv, par, smooth, prec_window = 28) {
  ig <- m_ramp(drv$tmin, par$tmin[1], par$tmin[2]) *
    m_ramp(drv$photo_s, par$photo_h[1] * 3600, par$photo_h[2] * 3600)
  if (!is.null(par$vpd)) ig <- ig * m_ramp(drv$vpd_pa, par$vpd[1], par$vpd[2], decreasing = TRUE)
  if (!is.null(par$prec)) ig <- ig * m_ramp(m_roll_sum(drv$ppt, prec_window), par$prec[1], par$prec[2])
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

# Fit thresholds on the given cells (rows) and years (columns of the obs matrices).
gm_fit <- function(drv, obs, dates, years, g, gmc, start) {
  fn <- function(th) {
    par <- gm_decode(th)
    pr <- gm_seasons(gm_gsi(drv, par, gmc$smooth), dates, years, g)
    gm_loss(pr, obs, gmc$penalty_days, gmc$cap_days)
  }
  o <- stats::optim(gm_encode(start), fn, method = "Nelder-Mead", control = list(maxit = gmc$maxit))
  list(par = gm_decode(o$par), loss = o$value, evals = o$counts[[1]])
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
