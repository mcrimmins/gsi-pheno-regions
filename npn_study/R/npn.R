# NPN grass curing study helpers (scripts 51, 52). Needs R/features.R (ra_matrix, hargreaves,
# first_true, last_true) and R/block2.R (run_lengths, onset_doy).

# Region of each site: EPA / CEC Level II ecoregion code (e.g. "10.1"), from
# data/npn_site_ecoregions.csv (scripts/59_npn_ecoregions.R); "none" when unknown. The NPN study
# uses this published regionalization, not the fuel phenology regions (decided 18 Oct 2026).
npn_site_regions <- function(site_id) {
  f <- npn_path("data", "npn_site_ecoregions.csv")
  if (!file.exists(f)) stop("needs data/npn_site_ecoregions.csv: run scripts/59_npn_ecoregions.R first")
  e <- utils::read.csv(f, stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character", eco_l1_code = "character"))
  r <- e$eco_l2_code[match(site_id, e$site_id)]
  r[is.na(r) | !nzchar(r)] <- "none"
  r
}
# Level II code -> name ("10.1" -> "Cold Deserts"), title case
npn_region_names <- function() {
  f <- npn_path("data", "npn_site_ecoregions.csv")
  if (!file.exists(f)) return(character())
  e <- unique(utils::read.csv(f, stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))[, c("eco_l2_code", "eco_l2_name")])
  nm <- tools::toTitleCase(tolower(e$eco_l2_name))
  nm <- gsub("\\bUsa\\b", "USA", nm); nm <- gsub("/", "/ ", nm, fixed = TRUE)   # wrap-friendly
  stats::setNames(nm, e$eco_l2_code)
}

# Sensitivity variant for 52-54 from the second command-line argument (config npn_study:
# sensitivity: <name>). Returns name, years to exclude, and the output sub-folder under npn_study/
# ("" for the main run, "sens/<name>" otherwise).
npn_sensitivity <- function(ns, args) {
  v <- if (length(args) >= 2) args[2] else ""
  if (!nzchar(v)) return(list(name = "", exclude_years = integer(), sub = ""))
  s <- ns$sensitivity[[v]]
  if (is.null(s)) stop("unknown sensitivity '", v, "'; config npn_study: sensitivity has: ", paste(names(ns$sensitivity), collapse = ", "))
  list(name = v, exclude_years = as.integer(unlist(s$exclude_years)), sub = file.path("sens", v))
}

# Series-years: one row per site x species x year, the median over that site-year's plants
# (NPN rows are plants; several plants of a species per site are common). Dates are day of year.
#   cm   npn_curing_metrics.rds;  sp  data/npn_species_pathway.csv;  min_readings per plant-year
npn_series_years <- function(cm, sp, min_readings = 12) {
  u <- cm[cm$usable %in% TRUE & cm$n_obs >= min_readings, ]
  u$key <- paste(u$site_id, u$species_id, u$year)
  med <- function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)
  g <- split(u, u$key)
  out <- do.call(rbind, lapply(g, function(d) data.frame(
    site_id = d$site_id[1], species_id = d$species_id[1], year = d$year[1], lat = d$lat[1], lon = d$lon[1],
    state = d$state[1], genus = d$genus[1], species = d$species[1], common_name = d$common_name[1],
    n_plants = length(unique(d$individual_id)), readings = stats::median(d$n_obs),
    greenup = med(d$g50_up_doy), curing = med(d$g50_down_doy), cured = med(d$cured_doy))))
  rownames(out) <- NULL
  out$series <- paste(out$site_id, out$species_id)
  m <- match(out$species_id, sp$species_id)
  out$pathway <- sp$pathway[m]; out$life_cycle <- sp$life_cycle[m]; out$origin <- sp$origin[m]; out$crop <- sp$crop[m]
  out
}

# Exponentially weighted moving average along days (rows = sites): m_t = m_{t-1} + (x_t - m_{t-1}) / tau
ewma_rows <- function(X, tau) {
  X[is.na(X)] <- 0
  t(apply(X, 1, function(x) as.numeric(stats::filter(x / tau, 1 - 1 / tau, method = "recursive"))))
}

# Per site and calendar year: thermal, freeze, moisture and memory metrics from daily PRISM.
# D: npn_site_daily.rds (sites, dates, ppt, tmin, tmax, vpdmax); b2: config features block2;
# taus: memory time constants (days). Returns a long data.frame site_id x year x metrics.
npn_site_year_metrics <- function(D, b2, taus = c(15, 30, 60)) {
  dates <- D$dates; yr <- as.integer(format(dates, "%Y")); doy <- as.integer(format(dates, "%j"))
  TM <- (D$tmin + D$tmax) / 2
  PET <- hargreaves(D$tmin, D$tmax, ra_matrix(D$sites$lat, doy))
  W <- D$ppt - PET
  EW <- lapply(stats::setNames(taus, paste0("ewma", taus)), function(tau) ewma_rows(W, tau))
  md <- function(Y, m, d) which(dates == as.Date(sprintf("%d-%02d-%02d", Y, m, d)))
  win <- function(a, b) which(dates >= a & dates <= b)
  dt <- function(Y, m, d) as.Date(sprintf("%d-%02d-%02d", Y, m, d))
  rows <- list()
  for (Y in sort(unique(yr))[-1]) {                      # first year only as the "year before"
    iy <- which(yr == Y); if (length(iy) < 365) next
    n <- nrow(D$ppt); o <- data.frame(site_id = D$sites$site_id, year = Y)
    cs <- function(X, base) t(apply(pmax(X[, iy, drop = FALSE] - base, 0), 1, cumsum))
    G5 <- cs(TM, 5); G10 <- cs(TM, 10)
    at <- function(C, m, d) C[, which(dates[iy] == dt(Y, m, d))]
    o$gdd5_mar31 <- at(G5, 3, 31); o$gdd5_apr30 <- at(G5, 4, 30); o$gdd5_may31 <- at(G5, 5, 31)
    o$gdd10_apr30 <- at(G10, 4, 30); o$gdd10_may31 <- at(G10, 5, 31); o$gdd10_jun30 <- at(G10, 6, 30)
    o$t_marapr <- rowMeans(TM[, win(dt(Y, 3, 1), dt(Y, 4, 30)), drop = FALSE])
    o$t_mayjun <- rowMeans(TM[, win(dt(Y, 5, 1), dt(Y, 6, 30)), drop = FALSE])
    o$t_julaug <- rowMeans(TM[, win(dt(Y, 7, 1), dt(Y, 8, 31)), drop = FALSE])
    o$t_sepoct <- rowMeans(TM[, win(dt(Y, 9, 1), dt(Y, 10, 31)), drop = FALSE])
    o$chill_days <- rowSums(TM[, win(dt(Y - 1, 11, 1), dt(Y, 2, 28)), drop = FALSE] < 5)
    sp_w <- win(dt(Y, 1, 1), dt(Y, 7, 15)); fa_w <- win(dt(Y, 7, 16), dt(Y, 12, 31))
    for (th in list(c("0", 0), c("2", -2.2))) {
      lt <- last_true(D$tmin[, sp_w, drop = FALSE] <= as.numeric(th[2]))
      o[[paste0("last_frz", th[1])]] <- ifelse(lt > 0, doy[sp_w][pmax(lt, 1)], NA)
      ft <- first_true(D$tmin[, fa_w, drop = FALSE] <= as.numeric(th[2]))
      o[[paste0("first_frz", th[1])]] <- ifelse(ft <= length(fa_w), doy[fa_w][pmin(ft, length(fa_w))], 367)
    }
    o$gdd5_at_last_frz2 <- vapply(seq_len(n), function(i) if (is.na(o$last_frz2[i])) 0 else G5[i, o$last_frz2[i]], 0)
    ps <- function(a, b) rowSums(D$ppt[, win(a, b), drop = FALSE])
    o$p_son_prev <- ps(dt(Y - 1, 9, 1), dt(Y - 1, 11, 30)); o$p_djf <- ps(dt(Y - 1, 12, 1), dt(Y, 2, 28))
    o$p_mam <- ps(dt(Y, 3, 1), dt(Y, 5, 31)); o$p_jja <- ps(dt(Y, 6, 1), dt(Y, 8, 31)); o$p_son <- ps(dt(Y, 9, 1), dt(Y, 11, 30))
    o$p_amj <- ps(dt(Y, 4, 1), dt(Y, 6, 30)); o$p_jas <- ps(dt(Y, 7, 1), dt(Y, 9, 30)); o$p_ann <- ps(dt(Y, 1, 1), dt(Y, 12, 31))
    o$pet_ann <- rowSums(PET[, iy, drop = FALSE])
    vm <- function(a, b) rowMeans(D$vpdmax[, win(a, b), drop = FALSE])
    o$vpd_mam <- vm(dt(Y, 3, 1), dt(Y, 5, 31)); o$vpd_jja <- vm(dt(Y, 6, 1), dt(Y, 8, 31))
    w49 <- win(dt(Y, 4, 1), dt(Y, 9, 30))
    o$events_amjjas <- rowSums(D$ppt[, w49, drop = FALSE] >= (b2$event_mm[[1]] %||% 5))
    o$dryspell_amjjas <- apply(run_lengths(D$ppt[, w49, drop = FALSE] < (b2$dry_day_mm %||% 1)), 1, max)
    o$dryspell_jjas <- apply(run_lengths(D$ppt[, win(dt(Y, 6, 1), dt(Y, 9, 30)), drop = FALSE] < (b2$dry_day_mm %||% 1)), 1, max)
    iy2 <- c(iy, max(iy) + seq_len(40))
    iy2 <- iy2[iy2 <= length(dates)]
    o$onset <- tryCatch(onset_doy(D$ppt[, iy2, drop = FALSE], doy[iy2], b2), error = function(e) rep(NA_real_, n))
    for (nm in names(EW)) for (m in c(5, 6, 7, 8, 9, 10))
      o[[sprintf("%s_%s1", nm, tolower(month.abb[m]))]] <- EW[[nm]][, md(Y, m, 1)]
    rows[[length(rows) + 1]] <- o
  }
  do.call(rbind, rows)
}

# Site normals: medians of the site-year metrics over years, plus a few derived ones.
npn_site_normals <- function(SY, D, years) {
  x <- SY[SY$year %in% years, ]
  num <- setdiff(names(x), c("site_id", "year"))
  N <- stats::aggregate(x[, num], list(site_id = x$site_id), stats::median, na.rm = TRUE)
  N$aridity <- N$p_ann / N$pet_ann
  N$p_jas_share <- N$p_jas / N$p_ann
  N <- merge(D$sites[, c("site_id", "lat", "lon")], N, by = "site_id")
  # day length gates (geometric): first day after the summer solstice day length < 12 h, < 11 h
  dl <- function(lat, d) { dec <- 0.409 * sin(2 * pi * d / 365 - 1.39); ws <- acos(pmin(pmax(-tan(lat * pi / 180) * tan(dec), -1), 1)); 24 / pi * ws }
  for (h in c(12, 11)) N[[sprintf("daylen_lt%dh", h)]] <- vapply(N$lat, function(la) { d <- 172:365; d[which(dl(la, d) < h)[1]] }, 0)
  N
}

# Green-up-anchored metrics per series-year (curing drivers): windows start at the series' typical
# green-up (median over its years; the site's median across species where the series has none),
# never at the observed curing date. gu: typical green-up doy per row; site_row: row in D$sites.
npn_anchor_metrics <- function(D, site_row, Y, gu, taus = c(30)) {
  dates <- D$dates; doy0 <- as.Date(sprintf("%d-01-01", Y)) - 1
  PET <- NULL
  out <- data.frame(p_gu30 = NA_real_, p_gu60 = NA_real_, p_gu90 = NA_real_, dry_gu90 = NA_real_, events_gu90 = NA_real_)[rep(1, length(Y)), ]
  for (k in seq_along(Y)) {
    if (is.na(gu[k]) || is.na(site_row[k])) next
    a <- match(doy0[k] + round(gu[k]), dates); if (is.na(a)) next
    p <- D$ppt[site_row[k], ]
    sw <- function(w) { b <- min(a + w - 1, length(p)); sum(p[a:b]) }
    out$p_gu30[k] <- sw(30); out$p_gu60[k] <- sw(60); out$p_gu90[k] <- sw(90)
    seg <- p[a:min(a + 89, length(p))] < 1
    r <- rle(seg); out$dry_gu90[k] <- if (any(r$values)) max(r$lengths[r$values]) else 0
    out$events_gu90[k] <- sum(p[a:min(a + 89, length(p))] >= 5)
  }
  out
}
