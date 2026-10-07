# Cut B curve step: merge Terra + Aqua 16-day composites into one series per cell, correct
# the sensor offset, drop poor composites, smooth (weighted Whittaker), and derive per-year
# greenness and curing metrics. Config `cutb_curves:`. Matrices are cells x composites/nodes.

cc_dir <- function(cfg, sub) file.path(cutb_dir(cfg), sub)
cc_year_file <- function(cfg, sub, name, year) file.path(cc_dir(cfg, sub), sprintf("%s_%d.tif", name, year))
cc_vars <- c("ndvi", "nir", "mir", "doy", "n_valid", "snow_frac")

# Composites available to a window: one row per (product, composite start date), taken from
# the file of the start date's own year when possible (AppEEARS repeats the previous
# December's composite at the start of each year's file).
cc_composites <- function(cfg, years) {
  out <- list()
  for (p in unlist(cfg$cutb$products)) for (y in years) {
    f <- cutb_out_file(cfg, p, "ndvi", y)
    if (!file.exists(f)) next
    st <- as.Date(names(terra::rast(f)))
    out[[length(out) + 1]] <- data.frame(product = p, year = y, band = seq_along(st), start = st,
                                         own = as.integer(format(st, "%Y")) == y)
  }
  if (!length(out)) return(NULL)
  d <- do.call(rbind, out)
  d <- d[order(d$product, d$start, !d$own), ]
  d[!duplicated(d[, c("product", "start")]), ]
}

# Observation date per cell: composite start + (doy - start doy) when that falls inside the
# 16-day window; otherwise (doy averaged across the year boundary, or missing) the midpoint.
cc_obs_offset <- function(doy, start) {
  sd <- as.integer(format(start, "%j"))
  ny <- ifelse(as.integer(format(start, "%Y")) %% 4 == 0, 366, 365)
  off <- (doy - sd) %% ny
  ifelse(!is.na(off) & off <= 16, off, 8)
}

# Read all composites of `comp` (rows of cc_composites) for one product, for grid rows
# [row, row + nr) and land cells; returns list of cells x composites matrices.
cc_read <- function(cfg, comp, row, nr, land) {
  res <- lapply(setNames(cc_vars, cc_vars), function(v) NULL)
  for (y in unique(comp$year)) {
    k <- comp$year == y
    for (v in cc_vars) {
      r <- terra::rast(cutb_out_file(cfg, comp$product[k][1], v, y))
      res[[v]] <- cbind(res[[v]], read_rows(r, row, nr, comp$band[k])[land, , drop = FALSE])
    }
  }
  res
}

# Quality mask and indices for one product's composites.
cc_prepare <- function(x, cc) {
  ok <- !is.na(x$ndvi) & !is.na(x$n_valid) & x$n_valid >= cc$min_valid &
    (is.na(x$snow_frac) | x$snow_frac <= cc$max_snow_frac)
  ndii <- (x$nir - x$mir) / (x$nir + x$mir)
  ok <- ok & !is.na(ndii)
  w <- pmin(x$n_valid / cc$full_valid, 1)
  list(ok = ok, ndvi = x$ndvi, ndii = ndii, w = ifelse(ok, w, 0))
}

# Vectorized Whittaker smoother (second differences) on a regular node grid: solves
# (W + lambda D'D) z = W y for every row at once with a banded Cholesky factorization.
# Exact against a dense solve (max difference ~1e-13 in tests).
whit_fit <- function(Y, W, lambda) {
  n <- ncol(Y); m <- nrow(Y)
  if (n < 4) stop("need at least 4 nodes")
  d0 <- c(1, 5, rep(6, n - 4), 5, 1); d1 <- c(0, -2, rep(-4, n - 3), -2); d2 <- c(0, 0, rep(1, n - 2))
  L0 <- L1 <- L2 <- U <- matrix(0, m, n)
  B <- W * ifelse(W > 0, Y, 0)
  for (i in seq_len(n)) {
    a0 <- W[, i] + lambda * d0[i] + 1e-8
    l2 <- if (i >= 3) lambda * d2[i] / L0[, i - 2] else 0
    l1 <- if (i >= 2) (lambda * d1[i] - (if (i >= 3) l2 * L1[, i - 1] else 0)) / L0[, i - 1] else 0
    L2[, i] <- l2; L1[, i] <- l1
    L0[, i] <- sqrt(a0 - l1^2 - l2^2)
    U[, i] <- (B[, i] - (if (i >= 2) l1 * U[, i - 1] else 0) - (if (i >= 3) l2 * U[, i - 2] else 0)) / L0[, i]
  }
  Z <- matrix(0, m, n)
  for (i in n:1) {
    Z[, i] <- (U[, i] - (if (i + 1 <= n) L1[, i + 1] * Z[, i + 1] else 0) -
                 (if (i + 2 <= n) L2[, i + 2] * Z[, i + 2] else 0)) / L0[, i]
  }
  Z
}

# Smooth with an upper-envelope pass: observations below the first fit get their weight
# multiplied by below_weight (residual cloud / aerosol bias is negative), then refit.
whit_envelope <- function(Y, W, lambda, iter = 1, below_weight = 0.5) {
  Z <- whit_fit(Y, W, lambda)
  for (k in seq_len(iter)) {
    W2 <- W; low <- W > 0 & Y < Z
    W2[low] <- W2[low] * below_weight
    Z <- whit_fit(Y, W2, lambda)
  }
  Z
}

# ---- Pass 1: Terra - Aqua offset for one year ------------------------------------------
# Terra is linearly interpolated (nominal composite midpoints) to each Aqua composite midpoint;
# the per-cell median difference over the year's valid pairs is written per index.
cutb_offset_year <- function(Y, cfg, grid_file, out_file, cc, chunk_rows, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    grid <- terra::rast(grid_file)
    comp <- cc_composites(cfg, Y)
    comp <- comp[comp$own, ]
    prods <- unlist(cfg$cutb$products)
    ct <- comp[comp$product == prods[1], ]; ca <- comp[comp$product == prods[2], ]
    if (!nrow(ct) || !nrow(ca)) stop("need both products for year ", Y)
    tm <- as.numeric(ct$start) + 8; am <- as.numeric(ca$start) + 8
    j <- findInterval(am, tm)                                 # Terra bracket for each Aqua composite
    use <- j >= 1 & j < length(tm)
    j <- j[use]; f <- (am[use] - tm[j]) / (tm[j + 1] - tm[j])
    nms <- c("ndvi_off", "ndii_off", "n_pairs")
    out <- matrix(NA_real_, terra::ncell(grid), 3, dimnames = list(NULL, nms))
    nc <- terra::ncol(grid)
    for (row in seq(1, terra::nrow(grid), by = chunk_rows)) {
      nr <- min(chunk_rows, terra::nrow(grid) - row + 1)
      cells <- (row - 1) * nc + seq_len(nr * nc)
      land <- !is.na(read_rows(grid, row, nr, 1)[, 1])
      if (!any(land)) next
      t <- cc_prepare(cc_read(cfg, ct, row, nr, land), cc)
      a <- cc_prepare(cc_read(cfg, ca, row, nr, land), cc)
      a_ok <- a$ok[, use, drop = FALSE]
      pair <- t$ok[, j, drop = FALSE] & t$ok[, j + 1, drop = FALSE] & a_ok
      o <- matrix(NA_real_, sum(land), 3)
      for (k in 1:2) {
        ix <- c("ndvi", "ndii")[k]
        ti <- t[[ix]][, j, drop = FALSE] * rep(1 - f, each = sum(land)) + t[[ix]][, j + 1, drop = FALSE] * rep(f, each = sum(land))
        dd <- ti - a[[ix]][, use, drop = FALSE]; dd[!pair] <- NA
        o[, k] <- apply(dd, 1, stats::median, na.rm = TRUE)
      }
      o[, 3] <- rowSums(pair)
      o[o[, 3] < cc$offset_min_pairs, 1:2] <- NA
      out[cells[land], ] <- o
    }
    write_feature_matrix(out, grid, nms, out_file)
    list(ok = TRUE, year = Y, n_cells = sum(!is.na(out[, 2])))
  }, error = function(e) list(ok = FALSE, year = Y, msg = conditionMessage(e)),
  finally = try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE))
}

# ---- Pass 2: smoothed curves and metrics for one year ----------------------------------
cc_metric_names <- c("ndvi_max", "ndvi_base", "ndvi_amp", "ndvi_peak_doy",
                     "sos20", "sos50", "eos50", "eos20", "gsl20", "n_green",
                     "ndii_max", "ndii_base", "ndii_amp", "ndii_peak_doy", "cure50", "cure_days",
                     "n_obs", "max_gap")

# Threshold crossings around a per-row peak on a node series Z (rows = cells).
# before = TRUE: last node before the peak below thr (within win), interpolated upward;
# before = FALSE: first node after the peak below thr. Returns day (node_day units) or NA.
cc_crossing <- function(Z, peak, thr, node_day, win, before = TRUE) {
  n <- ncol(Z); m <- nrow(Z)
  colm <- matrix(seq_len(n), m, n, byrow = TRUE)
  rel <- colm - peak
  cand <- Z < thr & (if (before) rel < 0 else rel > 0) & abs(rel) <= win
  cand[is.na(cand)] <- FALSE
  k <- if (before) last_true(cand) else first_true(cand)
  found <- if (before) k > 0 else k <= n
  out <- rep(NA_real_, m)
  i <- which(found); if (!length(i)) return(out)
  kk <- k[i]; nb <- if (before) kk + 1 else kk - 1       # neighbour on the peak side (above thr)
  z1 <- Z[cbind(i, kk)]; z2 <- Z[cbind(i, nb)]
  fr <- (thr[i] - z1) / (z2 - z1); fr[!is.finite(fr)] <- 0
  out[i] <- node_day[kk] + fr * (node_day[nb] - node_day[kk])
  out
}

cc_row_max <- function(Z, keep) { Z[!keep] <- -Inf; i <- max.col(Z, ties.method = "first"); list(i = i, v = Z[cbind(seq_len(nrow(Z)), i)]) }
cc_row_min_win <- function(Z, centre, win) {
  colm <- matrix(seq_len(ncol(Z)), nrow(Z), ncol(Z), byrow = TRUE)
  Z[abs(colm - centre) > win] <- Inf
  do.call(pmin, c(lapply(seq_len(ncol(Z)), function(j) Z[, j]), list(na.rm = TRUE)))
}

cutb_curves_year <- function(Y, cfg, grid_file, offset_file, out_files, cc, chunk_rows, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    grid <- terra::rast(grid_file)
    w0 <- as.Date(sprintf("%d-%02d-01", Y - 1, cc$window$start_month))
    w1 <- as.Date(sprintf("%d-%02d-01", Y + 1, cc$window$end_month + 1)) - 1
    comp <- cc_composites(cfg, (Y - 1):(Y + 1))
    comp <- comp[comp$start >= w0 - 16 & comp$start <= w1, ]
    prods <- unlist(cfg$cutb$products)
    nodes <- seq(w0, w1, by = cc$node_days); nn <- length(nodes)
    node_day <- as.numeric(nodes - as.Date(sprintf("%d-01-01", Y))) + 1   # DOY of Y (can be < 1 or > 365)
    in_y <- node_day >= 1 & node_day <= as.numeric(as.Date(sprintf("%d-12-31", Y)) - as.Date(sprintf("%d-01-01", Y))) + 1
    out_doy <- seq(1, 361, by = cc$curve_step)                             # output curve DOYs
    jx <- findInterval(out_doy, node_day); fx <- (out_doy - node_day[jx]) / cc$node_days
    win <- ceiling(cc$season_halfwidth_days / cc$node_days)
    offr <- if (!is.null(offset_file) && file.exists(offset_file)) terra::rast(offset_file) else NULL
    ncv <- length(out_doy)
    OUT <- list(ndvi = matrix(NA_real_, terra::ncell(grid), ncv), ndii = matrix(NA_real_, terra::ncell(grid), ncv),
                met = matrix(NA_real_, terra::ncell(grid), length(cc_metric_names), dimnames = list(NULL, cc_metric_names)))
    nc <- terra::ncol(grid)
    for (row in seq(1, terra::nrow(grid), by = chunk_rows)) {
      nr <- min(chunk_rows, terra::nrow(grid) - row + 1)
      cells <- (row - 1) * nc + seq_len(nr * nc)
      land <- !is.na(read_rows(grid, row, nr, 1)[, 1])
      if (!any(land)) next
      m <- sum(land)
      off <- if (!is.null(offr)) read_rows(offr, row, nr, 1:2)[land, , drop = FALSE] else matrix(0, m, 2)
      off[is.na(off)] <- 0
      YS <- list(ndvi = matrix(0, m, nn), ndii = matrix(0, m, nn)); WS <- matrix(0, m, nn)
      n_obs <- numeric(m); obs_days <- vector("list", 0)
      lastd <- rep(NA_real_, m); maxgap <- rep(0, m); firstd <- rep(NA_real_, m)
      for (p in prods) {
        cp <- comp[comp$product == p, ]
        if (!nrow(cp)) next
        x <- cc_read(cfg, cp, row, nr, land)
        q <- cc_prepare(x, cc)
        sgn <- if (p == prods[1]) -0.5 else 0.5                            # split the offset between sensors
        for (k in seq_len(nrow(cp))) {
          ok <- q$ok[, k]
          if (!any(ok)) next
          od <- as.numeric(cp$start[k] - w0) + cc_obs_offset(x$doy[, k], cp$start[k])
          ok <- ok & od >= 0 & od <= as.numeric(w1 - w0)
          if (!any(ok)) next
          node <- pmin(round(od / cc$node_days) + 1, nn)
          i <- which(ok); ij <- cbind(i, node[i]); wk <- q$w[i, k]
          YS$ndvi[ij] <- YS$ndvi[ij] + wk * (q$ndvi[i, k] + sgn * off[i, 1])
          YS$ndii[ij] <- YS$ndii[ij] + wk * (q$ndii[i, k] + sgn * off[i, 2])
          WS[ij] <- WS[ij] + wk
          dY <- node_day[node[i]]
          iny <- dY >= 1 & dY <= 366
          n_obs[i[iny]] <- n_obs[i[iny]] + 1
        }
      }
      for (ix in c("ndvi", "ndii")) { h <- WS > 0; YS[[ix]][h] <- YS[[ix]][h] / WS[h] }
      # largest gap (days) between nodes with observations, within year Y (edges included)
      has <- WS > 0
      for (j in which(in_y)) {
        gap_now <- ifelse(is.na(lastd), node_day[j] - 1, node_day[j] - lastd)
        upd <- has[, j]
        maxgap[upd] <- pmax(maxgap[upd], gap_now[upd])
        lastd[upd] <- node_day[j]
      }
      endY <- max(node_day[in_y])
      maxgap <- pmax(maxgap, ifelse(is.na(lastd), endY, endY - lastd))
      enough <- rowSums(has) >= cc$min_obs_window
      Z <- lapply(YS, function(Yv) whit_envelope(Yv, WS, cc$lambda, cc$envelope$iter, cc$envelope$below_weight))
      met <- matrix(NA_real_, m, length(cc_metric_names), dimnames = list(NULL, cc_metric_names))
      # NDVI season
      pk <- cc_row_max(Z$ndvi, matrix(in_y, m, nn, byrow = TRUE))
      base <- cc_row_min_win(Z$ndvi, pk$i, win)
      amp <- pk$v - base
      met[, "ndvi_max"] <- pk$v; met[, "ndvi_base"] <- base; met[, "ndvi_amp"] <- amp
      met[, "ndvi_peak_doy"] <- node_day[pk$i]
      seas <- amp >= cc$min_amp
      for (fr in c(20, 50)) {
        thr <- base + fr / 100 * amp
        s <- cc_crossing(Z$ndvi, pk$i, thr, node_day, win, TRUE)
        e <- cc_crossing(Z$ndvi, pk$i, thr, node_day, win, FALSE)
        met[, paste0("sos", fr)] <- ifelse(seas, s, NA); met[, paste0("eos", fr)] <- ifelse(seas, e, NA)
      }
      met[, "gsl20"] <- met[, "eos20"] - met[, "sos20"]
      # green spells: runs of >= min_spell_days above base + green_frac x amp within year Y
      above <- (Z$ndvi > base + cc$green_frac * amp)[, in_y, drop = FALSE]
      need <- ceiling(cc$min_spell_days / cc$node_days)
      run <- integer(m); cnt <- integer(m)
      for (j in seq_len(ncol(above))) { run <- ifelse(above[, j], run + 1L, 0L); cnt <- cnt + (run == need) }
      met[, "n_green"] <- ifelse(seas, cnt, 0)
      # NDII7 curing: peak within year Y, then decline after the peak
      pk2 <- cc_row_max(Z$ndii, matrix(in_y, m, nn, byrow = TRUE))
      base2 <- cc_row_min_win(Z$ndii, pk2$i, win)
      amp2 <- pk2$v - base2
      met[, "ndii_max"] <- pk2$v; met[, "ndii_base"] <- base2; met[, "ndii_amp"] <- amp2
      met[, "ndii_peak_doy"] <- node_day[pk2$i]
      seas2 <- amp2 >= cc$min_amp
      c50 <- cc_crossing(Z$ndii, pk2$i, base2 + 0.5 * amp2, node_day, win, FALSE)
      c75 <- cc_crossing(Z$ndii, pk2$i, base2 + 0.75 * amp2, node_day, win, FALSE)
      c25 <- cc_crossing(Z$ndii, pk2$i, base2 + 0.25 * amp2, node_day, win, FALSE)
      met[, "cure50"] <- ifelse(seas2, c50, NA); met[, "cure_days"] <- ifelse(seas2, c25 - c75, NA)
      met[, "n_obs"] <- n_obs; met[, "max_gap"] <- maxgap
      met[!enough, ] <- NA
      cv <- function(Zm) { o <- Zm[, jx, drop = FALSE] * rep(1 - fx, each = m) + Zm[, jx + 1, drop = FALSE] * rep(fx, each = m); o[!enough, ] <- NA; o }
      OUT$ndvi[cells[land], ] <- cv(Z$ndvi); OUT$ndii[cells[land], ] <- cv(Z$ndii)
      OUT$met[cells[land], ] <- met
    }
    cn <- sprintf("doy%03d", out_doy)
    write_feature_matrix(OUT$ndvi, grid, cn, out_files$ndvi)
    write_feature_matrix(OUT$ndii, grid, cn, out_files$ndii)
    write_feature_matrix(OUT$met, grid, cc_metric_names, out_files$met)
    list(ok = TRUE, year = Y, n_cells = sum(!is.na(OUT$met[, 1])))
  }, error = function(e) list(ok = FALSE, year = Y, msg = conditionMessage(e)),
  finally = try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE))
}
