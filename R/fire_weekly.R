# Helpers for the region x week x year table (script 65; fire study). Weeks: fire_week() in R/fire.R.

# Mid-week day of year for weeks 1..52 (week 52 runs to Dec 31).
fw_mid_doy <- function(n_weeks = 52) (seq_len(n_weeks) - 1) * 7 + 4

# Is each cell green on each day d (DOY relative to Jan 1 of Y; matrix cells x days)? Seasons
# from year Y and its neighbours, all as DOY relative to their own Jan 1: green when d falls in
# [sos, eos) of Y, of Y-1 (shifted back by the length of Y-1) or of Y+1 (shifted forward by the
# length of Y). A missing season counts as no green.
fw_green <- function(d, S, E, len_prev = 365, len_this = 365) {
  # S, E: cells x 3 matrices (columns Y-1, Y, Y+1)
  shift <- c(-len_prev, 0, len_this)
  G <- matrix(FALSE, nrow(S), length(d))
  for (k in 1:3) {
    s <- S[, k] + shift[k]; e <- E[, k] + shift[k]
    ok <- !is.na(s) & !is.na(e) & e > s
    if (!any(ok)) next
    G[ok, ] <- G[ok, ] | (outer(s[ok], d, "<=") & outer(e[ok], d, ">"))
  }
  G
}

# Weighted mean and weighted median (NA dropped).
fw_wmean <- function(x, w) { i <- !is.na(x) & !is.na(w) & w > 0; if (!any(i)) NA_real_ else sum(x[i] * w[i]) / sum(w[i]) }
fw_wmedian <- function(x, w) {
  i <- !is.na(x) & !is.na(w) & w > 0
  if (!any(i)) return(NA_real_)
  o <- order(x[i]); xs <- x[i][o]; cw <- cumsum(w[i][o]) / sum(w[i])
  xs[which(cw >= 0.5)[1]]
}

# Unique 375 m pixel-days: detections snapped to a 375 m grid on the sinusoidal projection
# (the MODIS / VIIRS land grid, sphere R = 6371007.181 m; computed directly, no PROJ datum grids),
# counted once per pixel and solar date. Returns a key per detection.
fw_pixel_key <- function(lon, lat, date, size = 375, R = 6371007.181) {
  y <- R * lat * pi / 180; x <- R * lon * pi / 180 * cos(lat * pi / 180)
  paste(floor(x / size), floor(y / size), format(date))
}
