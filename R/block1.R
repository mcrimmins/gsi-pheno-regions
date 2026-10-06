# Block 1: seasonal climate features for one year label Y (see config features: block1).
# Runs in a worker: gets file paths and plain values, returns a small status list.

block1_names <- function(b1) {
  s <- names(b1$seasons)
  thr <- unlist(b1$freeze_thresholds)
  tag <- function(t) sub("-", "m", sub("\\.", "p", format(t)))      # 0 -> "0", -2.2 -> "m2p2"
  c(paste0("t_", s), "t_ann", "t_range",
    "p_ann", paste0("pf_", s),
    "ai_ann", "ai_warm",
    unlist(lapply(thr, function(t) paste0(c("lsf_", "fff_", "ffs_"), tag(t)))),
    paste0("gdd", unlist(b1$gdd_bases)))
}

block1_year <- function(Y, files, grid_file, out_file, b1, chunk_rows, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    grid <- terra::rast(grid_file)
    R <- lapply(files, function(v) lapply(v, terra::rast))   # files[[var]][[as.character(year)]]
    d_prev <- as.Date(names(R$ppt[[as.character(Y - 1)]]))
    d_cur <- as.Date(names(R$ppt[[as.character(Y)]]))
    # water year: Oct-Dec of Y-1 + Jan-Sep of Y; calendar year: all of Y
    i_prev <- which(as.integer(format(d_prev, "%m")) >= 10)
    i_cur_wy <- which(as.integer(format(d_cur, "%m")) <= 9)
    wy_dates <- c(d_prev[i_prev], d_cur[i_cur_wy])
    if (length(wy_dates) < 360) stop("water year ", Y, " has only ", length(wy_dates), " days")
    wy_month <- as.integer(format(wy_dates, "%m"))
    wy_doy <- as.integer(format(wy_dates, "%j"))
    cal_doy <- as.integer(format(d_cur, "%j"))
    seas <- lapply(b1$seasons, function(m) which(wy_month %in% unlist(m)))
    warm <- unlist(seas[unlist(b1$warm_seasons)])
    spring <- which(cal_doy <= b1$spring_end_doy); fall <- which(cal_doy > b1$spring_end_doy)
    thr <- unlist(b1$freeze_thresholds); bases <- unlist(b1$gdd_bases)
    nms <- block1_names(b1)

    out <- matrix(NA_real_, terra::ncell(grid), length(nms), dimnames = list(NULL, nms))
    nr_all <- terra::nrow(grid); nc <- terra::ncol(grid)
    for (row in seq(1, nr_all, by = chunk_rows)) {
      nr <- min(chunk_rows, nr_all - row + 1)
      cells <- (row - 1) * nc + seq_len(nr * nc)
      land <- !is.na(read_rows(grid, row, nr, 1)[, 1])
      if (!any(land)) next
      wy <- function(v) cbind(read_rows(R[[v]][[as.character(Y - 1)]], row, nr, i_prev),
                              read_rows(R[[v]][[as.character(Y)]], row, nr, i_cur_wy))[land, , drop = FALSE]
      P <- wy("ppt"); TN <- wy("tmin"); TX <- wy("tmax"); TM <- (TN + TX) / 2
      lat <- terra::yFromRow(grid, row + (seq_len(nr) - 1))
      ra <- ra_matrix(lat, wy_doy)[rep(seq_len(nr), each = nc)[land], , drop = FALSE]
      PET <- hargreaves(TN, TX, ra)
      o <- matrix(NA_real_, sum(land), length(nms), dimnames = list(NULL, nms))
      st <- sapply(seas, function(ix) rowMeans(TM[, ix, drop = FALSE]))
      st <- matrix(st, nrow = sum(land))
      o[, paste0("t_", names(seas))] <- st
      o[, "t_ann"] <- rowMeans(TM)
      o[, "t_range"] <- apply(st, 1, max) - apply(st, 1, min)
      p_ann <- rowSums(P)
      o[, "p_ann"] <- p_ann
      for (s in names(seas)) {
        o[, paste0("pf_", s)] <- ifelse(p_ann > 0, rowSums(P[, seas[[s]], drop = FALSE]) / p_ann, NA)
      }
      o[, "ai_ann"] <- p_ann / pmax(rowSums(PET), 1e-6)
      o[, "ai_warm"] <- rowSums(P[, warm, drop = FALSE]) / pmax(rowSums(PET[, warm, drop = FALSE]), 1e-6)
      # calendar year Y: freeze dates and GDD
      CN <- read_rows(R$tmin[[as.character(Y)]], row, nr, seq_along(d_cur))[land, , drop = FALSE]
      CX <- read_rows(R$tmax[[as.character(Y)]], row, nr, seq_along(d_cur))[land, , drop = FALSE]
      CM <- (CN + CX) / 2
      tag <- function(t) sub("-", "m", sub("\\.", "p", format(t)))
      for (t in thr) {
        lsf <- cal_doy[spring][pmax(last_true(CN[, spring, drop = FALSE] <= t), 1)]
        lsf[last_true(CN[, spring, drop = FALSE] <= t) == 0] <- 0          # no spring freeze
        ff_i <- first_true(CN[, fall, drop = FALSE] <= t)
        fff <- ifelse(ff_i > length(fall), max(cal_doy) + 1, cal_doy[fall][pmin(ff_i, length(fall))])
        o[, paste0("lsf_", tag(t))] <- lsf
        o[, paste0("fff_", tag(t))] <- fff
        o[, paste0("ffs_", tag(t))] <- fff - lsf - 1
      }
      for (b in bases) o[, paste0("gdd", b)] <- rowSums(pmax(CM - b, 0))
      # rows with any missing input stay NA
      bad <- rowSums(is.na(P)) > 0 | rowSums(is.na(TN)) > 0 | rowSums(is.na(CN)) > 0
      o[bad, ] <- NA
      out[cells[land], ] <- o
    }
    write_feature_matrix(out, grid, nms, out_file)
    list(ok = TRUE, year = Y, n_cells = sum(!is.na(out[, 1])))
  }, error = function(e) list(ok = FALSE, year = Y, msg = conditionMessage(e)),
  finally = try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE))
}
