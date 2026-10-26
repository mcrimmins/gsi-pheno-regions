#!/usr/bin/env Rscript
# Step 58i (NPN grass curing study; laptop, ~2-5 min): does NPN percent green track in-season fuel
# status, not only the season's dates? Two parts, on every usable plant-year (crops out; C3 / C4):
#
# 1. Reversible dips. Within the green season (first reading >= 50 % to the sustained curing date
#    of 58f, i.e. the first fall below 50 % after the peak with no reading >= 50 % within
#    npn_study: sustained_days, default 21), a dip is a fall of >= inseason: dip_depth points
#    (default 25) from the running maximum, followed by a rise of >= dip_depth points from the low
#    within inseason: dip_recover_days (default 42) of the low. Readings are class midpoints.
#    Reported: share of plant-years with a dip, dips per plant-year, depth and length, how many dips
#    cross 50 % (those set the reference curing date early when they come after the peak), by grass
#    type, ecoregion and month; and the weather during each dip (VPD and rain anomalies).
#
# 2. Weekly greenness change vs weather. Pairs of consecutive readings of one plant 5-14 days apart
#    inside the green season. Rate = change in percent green per week. Anomaly = rate minus the same
#    series' mean rate in that week of the year in the other years (needs >= 2 other years). Weather
#    for the pair: mean daily VPDmax over the 7 days before the first reading to the second, and rain
#    over the 21 days to the second reading, each as an anomaly from the site's daily
#    climatology (npn_study: inseason: clim_years, default 2001-2024; 15-day smoothed), standardized. Fitted at two scales, plant pairs and series-weeks
#    (mean over the plants of a series in one interval week), by grass type and by season half
#    (before / after the series' typical 50 % green-up + 60 days): standardized slopes, R^2, and 95 %
#    intervals from a bootstrap over sites (B = inseason: boot, default 200).
#
# Inputs: npn_curing_metrics.rds and npn_green_obs.rds (GSI_OBS_DIR or config study_dir),
# data/npn_species_pathway.csv, <run_dir>/npn_study/npn_site_daily.rds (50).
# Outputs <run_dir>/npn_study/inseason/: dips.csv (one row per dip), dips_series.csv (dips on series
# curves, the mean of >= 2 plants per date), dips_summary.csv,
# pairs_fit.csv; report/tables/npn_inseason_dips.csv, npn_inseason_fit.csv.
#
# Usage: Rscript npn_study/scripts/58i_npn_inseason.R [profile]   (after 50, 59)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("58i_npn_inseason", profile)
ns <- cfg$npn_study %||% list(); ip <- ns$inseason %||% list()
base <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(base, "inseason"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
sustained_days <- as.numeric(ns$sustained_days %||% 21)
depth <- as.numeric(ip$dip_depth %||% 25); max_rec <- as.numeric(ip$dip_recover_days %||% 42); B <- as.integer(ip$boot %||% 200)
t0 <- Sys.time(); set.seed(58)
log_msg(sprintf("=== 58i_npn_inseason | profile: %s | dip >= %g points, recovery within %g d; sustained curing %g d; B = %d", profile, depth, max_rec, sustained_days, B))
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id); cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$pathway %in% c("C3", "C4"), ]
u$py <- paste(u$individual_id, u$year); u$region <- npn_site_regions(u$site_id); u$series <- paste(u$site_id, u$species_id)
go <- go[!is.na(go$pct_green), ]; go$py <- paste(go$individual_id, go$year); go <- go[go$py %in% u$py, ]
go <- go[order(go$py, go$doy), ]
idx <- split(seq_len(nrow(go)), go$py)
log_msg(sprintf("plant-years: %d with readings (%d readings)", length(idx), nrow(go)))

# ---- green season per plant-year: first reading >= 50 % to the sustained curing date ---------------
season <- function(d, p) {
  n <- length(p); s0 <- which(p >= 50)[1]; if (is.na(s0)) return(c(NA, NA))
  st <- which.max(p); idn <- which(p[-n] >= 50 & p[-1] < 50); idn <- idn[idn >= st]
  for (i in idn) { x <- d[i] + (p[i] - 50) / (p[i] - p[i + 1]) * (d[i + 1] - d[i])
    j <- which(d > x & d <= x + sustained_days); if (!any(p[j] >= 50)) return(c(d[s0], x)) }
  c(d[s0], d[n])                                         # no sustained fall: season runs to the last reading
}
# ---- dips -------------------------------------------------------------------------------------------
dips_one <- function(d, p) {
  out <- list(); n <- length(p); i <- 1
  while (i < n) {
    M <- p[i]; iM <- i; j <- i + 1; hit <- FALSE
    while (j <= n) {
      if (p[j] > M) { M <- p[j]; iM <- j }
      if (M - p[j] >= depth) {
        k <- j; mn <- p[j]; imn <- j; rec <- NA
        while (k <= n && d[k] - d[imn] <= max_rec) {
          if (p[k] < mn) { mn <- p[k]; imn <- k }
          if (p[k] - mn >= depth) { rec <- k; break }
          k <- k + 1 }
        if (!is.na(rec)) { nl <- sum(p[iM:rec] <= (M + mn) / 2); out[[length(out) + 1]] <- c(d[iM], d[imn], M, mn, d[rec], p[rec], nl); i <- rec } else i <- n
        hit <- TRUE; break }
      j <- j + 1 }
    if (!hit) break }
  if (length(out)) do.call(rbind, out) else NULL
}
SE <- t(vapply(idx, function(i) season(go$doy[i], go$pct_green[i]), c(0, 0)))
DP <- list(); n_season <- 0; season_py <- character()
for (k in seq_along(idx)) {
  if (is.na(SE[k, 1])) next
  i <- idx[[k]]; d <- go$doy[i]; p <- go$pct_green[i]; w <- d >= SE[k, 1] & d <= SE[k, 2] + sustained_days
  if (sum(w) < 4) next
  n_season <- n_season + 1; season_py <- c(season_py, names(idx)[k])
  r <- dips_one(d[w], p[w]); if (is.null(r)) next
  DP[[length(DP) + 1]] <- data.frame(py = names(idx)[k], start_doy = r[, 1], low_doy = r[, 2], high = r[, 3], low = r[, 4], recover_doy = r[, 5], recover = r[, 6], low_readings = r[, 7])
}
# the same on series curves: mean percent green over the plants of a series on each date with >= 2
# plants read (a dip there is shared by the plants, not one plant or one reading)
go$series <- paste(go$site_id, go$species_id)
SC <- stats::aggregate(pct_green ~ series + year + doy, go, mean); SC$n <- stats::aggregate(pct_green ~ series + year + doy, go, length)$pct_green
SC <- SC[SC$n >= 2, ]; SC <- SC[order(SC$series, SC$year, SC$doy), ]
sidx <- split(seq_len(nrow(SC)), paste(SC$series, SC$year))
SD <- list(); n_sseason <- 0; sy_ok <- character()
for (k in seq_along(sidx)) {
  i <- sidx[[k]]; d <- SC$doy[i]; p <- SC$pct_green[i]; se <- season(d, p); if (is.na(se[1])) next
  w <- d >= se[1] & d <= se[2] + sustained_days; if (sum(w) < 4) next
  n_sseason <- n_sseason + 1; sy_ok <- c(sy_ok, names(sidx)[k])
  r <- dips_one(d[w], p[w]); if (is.null(r)) next
  SD[[length(SD) + 1]] <- data.frame(sy = names(sidx)[k], series = SC$series[i[1]], year = SC$year[i[1]], start_doy = r[, 1], low_doy = r[, 2], high = r[, 3], low = r[, 4], recover_doy = r[, 5], low_readings = r[, 7])
}
SD <- if (length(SD)) do.call(rbind, SD) else NULL
if (!is.null(SD)) { SD$pathway <- u$pathway[match(SD$series, u$series)]; SD$crosses50 <- SD$high >= 50 & SD$low < 50 }
D1 <- if (length(DP)) do.call(rbind, DP) else data.frame(py = character(), start_doy = numeric(), low_doy = numeric(), high = numeric(), low = numeric(), recover_doy = numeric(), recover = numeric())
mm <- match(D1$py, u$py)
D1 <- cbind(D1, u[mm, c("site_id", "species_id", "series", "year", "region", "pathway")])
D1$depth <- D1$high - D1$low; D1$days <- D1$recover_doy - D1$start_doy; D1$crosses50 <- D1$high >= 50 & D1$low < 50; D1$single_reading <- D1$low_readings <= 1
D1$after_peak <- D1$start_doy >= u$peak_doy[mm]
D1$month <- as.integer(format(as.Date(D1$low_doy - 1, origin = sprintf("%d-01-01", D1$year)), "%m"))
log_msg(sprintf("dips: %d in %d of %d plant-years with a green season (%.0f %%); %d cross 50 %% (%d of them after the peak); %.0f %% rest on a single low reading",
                nrow(D1), length(unique(D1$py)), n_season, 100 * length(unique(D1$py)) / n_season, sum(D1$crosses50), sum(D1$crosses50 & D1$after_peak), 100 * mean(D1$single_reading)))
sy_pw <- u$pathway[match(sub(" [0-9]+$", "", sy_ok), u$series)]
if (!is.null(SD)) log_msg(sprintf("series curves (>= 2 plants per date): %d dips in %d of %d series-years (%.0f %%; C3 %.0f %%, C4 %.0f %%); %d cross 50 %%; %.0f %% on a single low reading",
  nrow(SD), length(unique(SD$sy)), n_sseason, 100 * length(unique(SD$sy)) / n_sseason,
  100 * length(unique(SD$sy[SD$pathway %in% "C3"])) / max(1, sum(sy_pw %in% "C3")), 100 * length(unique(SD$sy[SD$pathway %in% "C4"])) / max(1, sum(sy_pw %in% "C4")),
  sum(SD$crosses50), 100 * mean(SD$low_readings <= 1)))

# ---- daily weather anomalies at the sites ----------------------------------------------------------
W <- readRDS(file.path(base, "npn_site_daily.rds"))
dates <- as.Date(W$dates); yr <- as.integer(format(dates, "%Y")); dy <- as.integer(format(dates, "%j"))
clim_years <- unlist(ip$clim_years %||% c(2001, 2024))
cl <- function(X) {                                       # sites x 366 smoothed daily climatology
  inb <- yr >= clim_years[1] & yr <= clim_years[2]; if (!any(inb)) inb <- rep(TRUE, length(yr))
  C <- vapply(1:366, function(k) rowMeans(X[, inb & dy == k, drop = FALSE], na.rm = TRUE), numeric(nrow(X)))
  if (is.null(dim(C))) C <- matrix(C, nrow = 1)
  C[, 366] <- ifelse(is.na(C[, 366]), C[, 365], C[, 366])
  K <- 15; h <- (K - 1) / 2
  t(apply(C, 1, function(z) { zz <- c(utils::tail(z, h), z, utils::head(z, h)); as.numeric(stats::filter(zz, rep(1 / K, K)))[(h + 1):(h + 366)] }))
}
VPa <- W$vpdmax - cl(W$vpdmax)[, dy, drop = FALSE]
P21 <- t(apply(W$ppt, 1, function(z) { z[is.na(z)] <- 0; as.numeric(stats::filter(z, rep(1, 21), sides = 1)) }))
P21a <- P21 - cl(P21)[, dy, drop = FALSE]
srow <- match(go$site_id, W$sites$site_id)
dix <- function(Y, d) match(as.Date(sprintf("%d-01-01", Y)) + round(d) - 1, dates)
win_mean <- function(M, r, a, b) vapply(seq_along(r), function(k) if (is.na(r[k]) || is.na(a[k]) || is.na(b[k])) NA_real_ else mean(M[r[k], a[k]:b[k]], na.rm = TRUE), 0)
# weather during each dip: high reading to the low
if (nrow(D1)) {
  r <- match(D1$site_id, W$sites$site_id); a <- dix(D1$year, D1$start_doy); b <- dix(D1$year, D1$low_doy)
  D1$vpd_anom <- win_mean(VPa, r, a, b); D1$rain21_anom <- vapply(seq_along(r), function(k) if (is.na(r[k]) || is.na(b[k])) NA_real_ else P21a[r[k], b[k]], 0)
}

# ---- reading pairs ------------------------------------------------------------------------------------
PR <- list()
for (k in seq_along(idx)) {
  if (is.na(SE[k, 1])) next
  i <- idx[[k]]; d <- go$doy[i]; p <- go$pct_green[i]; n <- length(i); if (n < 2) next
  g <- diff(d); ok <- g >= 5 & g <= 14 & d[-n] >= SE[k, 1] & d[-1] <= SE[k, 2]
  if (!any(ok)) next
  PR[[length(PR) + 1]] <- data.frame(py = names(idx)[k], row = i[-n][ok], d1 = d[-n][ok], d2 = d[-1][ok], rate = (diff(p)[ok]) / g[ok] * 7)
}
PR <- do.call(rbind, PR)
mm <- match(PR$py, u$py); PR$site_id <- u$site_id[mm]; PR$series <- u$series[mm]; PR$year <- u$year[mm]; PR$pathway <- u$pathway[mm]; PR$region <- u$region[mm]
PR$week <- floor(((PR$d1 + PR$d2) / 2 - 1) / 7)
# expected rate: the series' mean in that week in the other years
key <- paste(PR$series, PR$week)
s_all <- tapply(PR$rate, key, sum); n_all <- tapply(PR$rate, key, length)
ky <- paste(key, PR$year); s_y <- tapply(PR$rate, ky, sum); n_y <- tapply(PR$rate, ky, length)
yrs_k <- tapply(PR$year, key, function(z) length(unique(z)))
oth_n <- n_all[key] - n_y[ky]; PR$expect <- ifelse(yrs_k[key] >= 3 & oth_n > 0, (s_all[key] - s_y[ky]) / oth_n, NA)
PR$anom <- PR$rate - PR$expect
r <- match(PR$site_id, W$sites$site_id); a <- dix(PR$year, PR$d1 - 7); b <- dix(PR$year, PR$d2)
PR$vpd_anom <- win_mean(VPa, r, a, b)
PR$rain21_anom <- vapply(seq_along(r), function(k) if (is.na(r[k]) || is.na(b[k])) NA_real_ else P21a[r[k], b[k]], 0)
gu_typ <- tapply(u$g50_up_doy, u$series, cmed)
PR$half <- ifelse(!is.na(gu_typ[PR$series]) & PR$d1 >= gu_typ[PR$series] + 60, "late", "early")
PR <- PR[!is.na(PR$anom) & !is.na(PR$vpd_anom) & !is.na(PR$rain21_anom), ]
log_msg(sprintf("reading pairs in the green season with an expected rate and weather: %d (%d series, %d sites)", nrow(PR), length(unique(PR$series)), length(unique(PR$site_id))))

# ---- fits -------------------------------------------------------------------------------------------------
fit1 <- function(X) {
  if (nrow(X) < 50) return(NULL)
  z <- function(v) (v - mean(v)) / stats::sd(v)
  f <- stats::lm(z(anom) ~ z(vpd_anom) + z(rain21_anom), data = X); cf <- stats::coef(f)
  c(b_vpd = unname(cf[2]), b_rain = unname(cf[3]), r2 = summary(f)$r.squared, r_vpd = stats::cor(X$anom, X$vpd_anom), r_rain = stats::cor(X$anom, X$rain21_anom))
}
boot1 <- function(X) {
  s <- unique(X$site_id); est <- fit1(X); if (is.null(est)) return(NULL)
  bs <- replicate(B, { ss <- sample(s, length(s), replace = TRUE); Xb <- do.call(rbind, lapply(seq_along(ss), function(k) X[X$site_id == ss[k], ])); v <- fit1(Xb); if (is.null(v)) rep(NA, 5) else v })
  lo <- apply(bs, 1, stats::quantile, 0.025, na.rm = TRUE); hi <- apply(bs, 1, stats::quantile, 0.975, na.rm = TRUE)
  data.frame(stat = names(est), estimate = est, lo = lo, hi = hi, n = nrow(X), sites = length(s))
}
SW <- stats::aggregate(cbind(anom, vpd_anom, rain21_anom) ~ series + site_id + year + week + pathway + half, PR, mean)
FT <- list()
for (sc in c("plant_pair", "series_week")) for (pw in c("all", "C3", "C4")) for (hf in c("all", "early", "late")) {
  X <- if (sc == "plant_pair") PR else SW
  X <- X[(pw == "all" | X$pathway == pw) & (hf == "all" | X$half == hf), ]
  b <- boot1(X); if (is.null(b)) next
  FT[[length(FT) + 1]] <- cbind(scale = sc, pathway = pw, half = hf, b)
  v <- b$estimate; names(v) <- b$stat
  log_msg(sprintf("  %-11s %-3s %-5s n %6d, %3d sites | VPD slope %+.2f [%+.2f, %+.2f], rain %+.2f [%+.2f, %+.2f] | R2 %.3f [%.3f, %.3f]", sc, pw, hf, b$n[1], b$sites[1],
                  v[["b_vpd"]], b$lo[1], b$hi[1], v[["b_rain"]], b$lo[2], b$hi[2], v[["r2"]], b$lo[3], b$hi[3]))
}
FT <- do.call(rbind, FT); rownames(FT) <- NULL

# ---- summaries -------------------------------------------------------------------------------------------
pys <- data.frame(py = season_py); pys$pathway <- u$pathway[match(pys$py, u$py)]; pys$region <- u$region[match(pys$py, u$py)]
summ <- function(by, lev, P, Dd) data.frame(by = by, value = lev, plant_years = nrow(P), with_dip = length(unique(Dd$py)), share_with_dip = length(unique(Dd$py)) / max(nrow(P), 1),
  dips = nrow(Dd), dips_per_plant_year = nrow(Dd) / max(nrow(P), 1), median_depth = if (nrow(Dd)) stats::median(Dd$depth) else NA, median_days = if (nrow(Dd)) stats::median(Dd$days) else NA,
  cross50 = sum(Dd$crosses50), cross50_after_peak = sum(Dd$crosses50 & Dd$after_peak), single_reading_share = if (nrow(Dd)) mean(Dd$single_reading) else NA,
  dip_vpd_anom_median = if (nrow(Dd)) stats::median(Dd$vpd_anom, na.rm = TRUE) else NA, dip_rain21_anom_median = if (nrow(Dd)) stats::median(Dd$rain21_anom, na.rm = TRUE) else NA)
SM <- rbind(summ("all", "all", pys, D1),
            do.call(rbind, lapply(c("C3", "C4"), function(k) summ("pathway", k, pys[pys$pathway == k, ], D1[D1$pathway == k, ]))),
            do.call(rbind, lapply(sort(unique(pys$region[pys$region != "none"])), function(k) summ("region", k, pys[pys$region == k, ], D1[D1$region == k, ]))),
            do.call(rbind, lapply(1:12, function(k) { Dd <- D1[D1$month == k, ]; s <- summ("month", month.abb[k], pys, Dd); s$plant_years <- NA; s$share_with_dip <- NA; s$dips_per_plant_year <- NA; s })))
if (!is.null(SD)) {
  sy_all <- data.frame(sy = sy_ok, pathway = sy_pw)
  sx <- function(lev, P, Dd) data.frame(by = "series_curve", value = lev, plant_years = nrow(P), with_dip = length(unique(Dd$sy)), share_with_dip = length(unique(Dd$sy)) / max(nrow(P), 1),
    dips = nrow(Dd), dips_per_plant_year = nrow(Dd) / max(nrow(P), 1), median_depth = if (nrow(Dd)) stats::median(Dd$high - Dd$low) else NA,
    median_days = if (nrow(Dd)) stats::median(Dd$recover_doy - Dd$start_doy) else NA, cross50 = sum(Dd$crosses50), cross50_after_peak = NA,
    single_reading_share = if (nrow(Dd)) mean(Dd$low_readings <= 1) else NA, dip_vpd_anom_median = NA, dip_rain21_anom_median = NA)
  SM <- rbind(SM, sx("all", sy_all, SD), sx("C3", sy_all[sy_all$pathway %in% "C3", ], SD[SD$pathway %in% "C3", ]), sx("C4", sy_all[sy_all$pathway %in% "C4", ], SD[SD$pathway %in% "C4", ]))
  utils::write.csv(SD, file.path(odir, "dips_series.csv"), row.names = FALSE)
}
# weather during dips against the pairs' weather (all green-season intervals)
log_msg(sprintf("weather during dips: VPD anomaly median %+.2f hPa (all green-season intervals %+.2f); 21-day rain anomaly %+.1f mm (all %+.1f)",
                stats::median(D1$vpd_anom, na.rm = TRUE), stats::median(PR$vpd_anom), stats::median(D1$rain21_anom, na.rm = TRUE), stats::median(PR$rain21_anom)))
for (i in seq_len(min(nrow(SM), 3))) with(SM[i, ], log_msg(sprintf("  %-7s %-4s %5d plant-years: %.0f %% with a dip, %.2f dips each, depth %.0f points, %.0f days; %d cross 50 %% (%d after the peak)",
  by, value, plant_years, 100 * share_with_dip, dips_per_plant_year, median_depth, median_days, cross50, cross50_after_peak)))

utils::write.csv(D1, file.path(odir, "dips.csv"), row.names = FALSE)
utils::write.csv(SM, file.path(odir, "dips_summary.csv"), row.names = FALSE)
utils::write.csv(FT, file.path(odir, "pairs_fit.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(SM, file.path(tdir, "npn_inseason_dips.csv"), row.names = FALSE)
utils::write.csv(FT, file.path(tdir, "npn_inseason_fit.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
