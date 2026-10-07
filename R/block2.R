# Block 2: daily timing and variability features for one year label Y (config features: block2).
# Runs in a worker: gets file paths and plain values, returns a small status list.
# Matrices are cells x days for one row chunk; helpers work row-wise without apply().

b2_tag <- function(t) sub("-", "m", sub(".", "p", as.character(t), fixed = TRUE))   # 2.5 -> "2p5" (per element)

block2_names <- function(b2) {
  ev <- unlist(b2$event_mm)
  c("onset_doy",
    paste0("ev_warm_", b2_tag(ev)), paste0("ev_cool_", b2_tag(ev)),
    "dsl_max_warm", paste0("n_dry", b2$dry_spell_days, "_warm"),
    paste0("gdd", b2$gdd_base, "_lhf"),
    "cold_frac", "dry_frac")
}

# Dry-run lengths ending at each day (0 on wet days), for a logical cells x days matrix.
run_lengths <- function(D) {
  R <- matrix(0L, nrow(D), ncol(D)); cur <- integer(nrow(D))
  for (j in seq_len(ncol(D))) { cur <- ifelse(D[, j], cur + 1L, 0L); R[, j] <- cur }
  R
}
# Row-wise cumulative sum along days.
row_cumsum <- function(X) { for (j in seq_len(ncol(X))[-1]) X[, j] <- X[, j - 1] + X[, j]; X }
# Forward window sums: S[, d] = sum(X[, d .. d + w - 1]) for d in idx (window must fit).
fwd_sum <- function(X, w, idx) {
  C <- cbind(0, row_cumsum(X))
  C[, idx + w, drop = FALSE] - C[, idx, drop = FALSE]
}
# Trailing window sums: S[, t] = sum(X[, t - w + 1 .. t]) for t in idx (t >= w).
trail_sum <- function(X, w, idx) {
  C <- cbind(0, row_cumsum(X))
  C[, idx + 1, drop = FALSE] - C[, idx - w + 1, drop = FALSE]
}

# Warm-season rain onset (calendar year): first day d in [start_doy, end_doy] that is wet,
# starts a window_days total >= total_mm, and is not followed by a dry spell >= max_dry_spell
# within check_days. Returns DOY, or end_doy + 1 when no onset (censored).
onset_doy <- function(P, cal_doy, b2) {
  o <- b2$onset
  idx <- which(cal_doy >= o$start_doy & cal_doy <= o$end_doy)
  need <- max(idx) + max(o$window_days - 1, o$check_days)
  if (need > ncol(P)) stop("onset window runs past the end of the year")
  wet <- P[, idx, drop = FALSE] >= b2$dry_day_mm
  S <- fwd_sum(P, o$window_days, idx)
  RL <- run_lengths(P[, seq_len(need), drop = FALSE] < b2$dry_day_mm)
  M <- matrix(0L, nrow(P), length(idx))
  for (k in seq_len(o$check_days)) M <- pmax(M, RL[, idx + k, drop = FALSE])
  cand <- wet & S >= o$total_mm & M < o$max_dry_spell
  f <- first_true(cand)
  ifelse(f > length(idx), o$end_doy + 1, cal_doy[idx][pmin(f, length(idx))])
}

block2_year <- function(Y, files, grid_file, out_file, b2, chunk_rows, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    grid <- terra::rast(grid_file)
    R <- lapply(files, function(v) lapply(v, terra::rast))   # files[[var]][[as.character(year)]]
    d_prev <- as.Date(names(R$ppt[[as.character(Y - 1)]]))
    d_cur <- as.Date(names(R$ppt[[as.character(Y)]]))
    lw <- b2$limit$window_days
    wy_start <- as.Date(sprintf("%d-10-01", Y - 1))
    # water year plus (lw - 1) lead days from September of Y-1 for the trailing windows
    i_prev <- which(d_prev >= wy_start - (lw - 1))
    i_cur_wy <- which(as.integer(format(d_cur, "%m")) <= 9)
    ext_dates <- c(d_prev[i_prev], d_cur[i_cur_wy])
    wy_pos <- which(ext_dates >= wy_start)                    # water-year days within ext
    if (length(wy_pos) < 360) stop("water year ", Y, " has only ", length(wy_pos), " days")
    ext_month <- as.integer(format(ext_dates, "%m"))
    ext_doy <- as.integer(format(ext_dates, "%j"))
    warm <- wy_pos[ext_month[wy_pos] %in% unlist(b2$warm_months)]
    cool <- wy_pos[ext_month[wy_pos] %in% unlist(b2$cool_months)]
    if (any(diff(warm) != 1)) stop("warm_months must be contiguous within the water year")
    cal_doy <- as.integer(format(d_cur, "%j"))
    spring <- which(cal_doy <= 181)
    ev <- unlist(b2$event_mm); nms <- block2_names(b2)

    out <- matrix(NA_real_, terra::ncell(grid), length(nms), dimnames = list(NULL, nms))
    nr_all <- terra::nrow(grid); nc <- terra::ncol(grid)
    for (row in seq(1, nr_all, by = chunk_rows)) {
      nr <- min(chunk_rows, nr_all - row + 1)
      cells <- (row - 1) * nc + seq_len(nr * nc)
      land <- !is.na(read_rows(grid, row, nr, 1)[, 1])
      if (!any(land)) next
      ext <- function(v) cbind(read_rows(R[[v]][[as.character(Y - 1)]], row, nr, i_prev),
                               read_rows(R[[v]][[as.character(Y)]], row, nr, i_cur_wy))[land, , drop = FALSE]
      cal <- function(v) read_rows(R[[v]][[as.character(Y)]], row, nr, seq_along(d_cur))[land, , drop = FALSE]
      P <- ext("ppt"); TN <- ext("tmin"); TX <- ext("tmax"); TM <- (TN + TX) / 2
      CP <- cal("ppt"); CN <- cal("tmin"); CX <- cal("tmax")
      o <- matrix(NA_real_, sum(land), length(nms), dimnames = list(NULL, nms))

      # onset (calendar year)
      o[, "onset_doy"] <- onset_doy(CP, cal_doy, b2)
      # rain events (water year, by season)
      for (t in ev) {
        o[, paste0("ev_warm_", b2_tag(t))] <- rowSums(P[, warm, drop = FALSE] >= t)
        o[, paste0("ev_cool_", b2_tag(t))] <- rowSums(P[, cool, drop = FALSE] >= t)
      }
      # warm-season dry spells (runs cut at the season edges)
      RL <- run_lengths(P[, warm, drop = FALSE] < b2$dry_day_mm)
      o[, "dsl_max_warm"] <- do.call(pmax, c(lapply(seq_len(ncol(RL)), function(j) RL[, j]), list(0L)))
      o[, paste0("n_dry", b2$dry_spell_days, "_warm")] <- rowSums(RL == b2$dry_spell_days)
      # GDD accumulated (from Jan 1) by the last spring hard freeze; 0 if none
      G <- row_cumsum(pmax((CN[, spring, drop = FALSE] + CX[, spring, drop = FALSE]) / 2 - b2$gdd_base, 0))
      lhf <- last_true(CN[, spring, drop = FALSE] <= b2$hard_freeze)
      g <- numeric(sum(land)); has <- lhf > 0
      g[has] <- G[cbind(which(has), lhf[has])]
      o[, paste0("gdd", b2$gdd_base, "_lhf")] <- g
      # limiting days (water year): cold, or warm with trailing P < dry_ratio x trailing PET
      lat <- terra::yFromRow(grid, row + (seq_len(nr) - 1))
      ra <- ra_matrix(lat, ext_doy)[rep(seq_len(nr), each = nc)[land], , drop = FALSE]
      PET <- hargreaves(TN, TX, ra)
      P30 <- trail_sum(P, lw, wy_pos); E30 <- trail_sum(PET, lw, wy_pos)
      cold <- TM[, wy_pos, drop = FALSE] < b2$limit$cold_tmean
      dry <- !cold & P30 < b2$limit$dry_ratio * E30
      o[, "cold_frac"] <- rowMeans(cold)
      o[, "dry_frac"] <- rowMeans(dry)

      bad <- rowSums(is.na(P)) > 0 | rowSums(is.na(TN)) > 0 | rowSums(is.na(TX)) > 0 |
        rowSums(is.na(CP)) > 0 | rowSums(is.na(CN)) > 0 | rowSums(is.na(CX)) > 0
      o[bad, ] <- NA
      out[cells[land], ] <- o
    }
    write_feature_matrix(out, grid, nms, out_file)
    list(ok = TRUE, year = Y, n_cells = sum(!is.na(out[, 1])))
  }, error = function(e) list(ok = FALSE, year = Y, msg = conditionMessage(e)),
  finally = try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE))
}

# Across-year variability (one raster): CV of water-year and warm-season precipitation and SD
# of spring (AMJ) temperature from the Block 1 per-year files; share of years with a
# warm-season onset from the Block 2 per-year files.
block2_variability <- function(b1_files, b2_files, out, b2, chunk_rows = 60) {
  r1 <- lapply(b1_files, terra::rast); r2 <- lapply(b2_files, terra::rast)
  grid <- r1[[1]][[1]]
  nms <- c("cv_p_ann", "cv_p_warm", "sd_t_amj", "onset_frac")
  res <- matrix(NA_real_, terra::ncell(grid), length(nms))
  ws <- paste0("pf_", c("amj", "jas"))
  b1n <- names(r1[[1]]); need <- c("p_ann", ws, "t_amj")
  if (!all(need %in% b1n)) stop("Block 1 files lack: ", paste(setdiff(need, b1n), collapse = ", "))
  ib1 <- match(need, b1n); ion <- match("onset_doy", names(r2[[1]]))
  cv <- function(m) { s <- apply(m, 1, stats::sd); mu <- rowMeans(m); ifelse(mu > 0, s / mu, NA) }
  for (row in seq(1, terra::nrow(grid), by = chunk_rows)) {
    nr <- min(chunk_rows, terra::nrow(grid) - row + 1)
    cells <- terra::cellFromRowCol(grid, row, 1):terra::cellFromRowCol(grid, row + nr - 1, terra::ncol(grid))
    a <- vapply(r1, function(r) read_rows(r, row, nr, ib1), matrix(0, length(cells), length(need)))
    on <- vapply(r2, function(r) read_rows(r, row, nr, ion)[, 1], numeric(length(cells)))
    on <- matrix(on, length(cells))
    pann <- matrix(a[, 1, ], length(cells)); pwarm <- pann * matrix(a[, 2, ] + a[, 3, ], length(cells))
    res[cells, 1] <- cv(pann); res[cells, 2] <- cv(pwarm)
    res[cells, 3] <- apply(matrix(a[, 4, ], length(cells)), 1, stats::sd)
    res[cells, 4] <- rowMeans(on <= b2$onset$end_doy)
  }
  write_feature_matrix(res, terra::rast(b1_files[1])[[1]], nms, out)
}
