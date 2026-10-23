#!/usr/bin/env Rscript
# Step 55f (NPN grass curing study; laptop): does weather explain the year-to-year departures once
# many series are averaged? The ladder (52) scored single series-years; the LTER check (55e) found
# a weather signal at the level of many plots in a year. Same series-years and site-year weather
# metrics as the ladder (51), scored at three scales:
#   series   one site x species in one year (as the ladder)
#   site     mean over the species at a site in a year
#   region   mean over the series in an EPA Level II ecoregion in a year (>= 3 series, >= 2 sites)
# For each scale, date (green-up, 50 % curing, cured) and grass type (all, C3, C4):
#   - departures: each date minus its series' circular median (|d| > 90 d dropped, as 52; series
#     need >= 3 years); drivers: site-year metric minus the site's mean over 2002-2025; both
#     averaged to the scale
#   - signal: variance of the scale's mean departures against what plant-to-plant / series noise
#     alone would give (mean within-group variance / n); ratio > 1 = a shared year signal
#   - skill: leave one year out, the single driver with the largest |Spearman r| in the training
#     years, a straight-line fit, predicted for the left-out year; skill = 1 - MAE / MAE(0)
#   - the drivers' Spearman r at each scale
#   - two driver sets: lead (windows ending >= lead_days before the target's median typical date,
#     as the ladder) and all (any window in the year, which can see past the date)
# Outputs <run_dir>/npn_study/scale/{skill,drivers,groups}.csv; report/tables/npn_scale_skill.csv.
# Sensitivity: second argument no1617 drops 2016-2017 (outputs under sens/no1617/).
#
# Usage: Rscript npn_study/scripts/55f_npn_scale.R [profile] [no1617]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
ns <- cfg$npn_study %||% list()
sens <- npn_sensitivity(ns, args)
log_file <- log_init(paste0("55f_npn_scale", if (nzchar(sens$name %||% "")) paste0("_", sens$name) else ""), profile)
base <- file.path(cfg$run_dir, "npn_study")
odir <- file.path(base, sens$sub %||% "", "scale"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
wcsv <- function(x, f) utils::write.csv(x, file.path(odir, f), row.names = FALSE)
log_msg("=== 55f_npn_scale | profile: ", profile, if (nzchar(sens$name %||% "")) paste0(" | sensitivity ", sens$name) else "", " -> ", odir)

m <- readRDS(file.path(base, "npn_metrics.rds"))
SY <- m$series_years; YM <- m$site_years
SY <- SY[SY$pathway %in% c("C3", "C4") & !(SY$crop %in% "yes"), ]
if (length(sens$exclude_years %||% integer())) SY <- SY[!SY$year %in% sens$exclude_years, ]
min_years <- ns$min_years %||% 3
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }

# drivers: site-year departures from the site mean over the normal years
nyr <- ns$normals %||% list(start = 2002, end = 2025)
drv <- setdiff(names(YM)[vapply(YM, is.numeric, TRUE)], c("site_id", "year"))
YMn <- YM[YM$year >= nyr$start & YM$year <= nyr$end, ]
mu <- stats::aggregate(YMn[, drv], list(site_id = YMn$site_id), function(z) mean(z, na.rm = TRUE))
YA <- YM[, c("site_id", "year", drv)]
YA[, drv] <- YA[, drv] - mu[match(YA$site_id, mu$site_id), drv]

# when each driver's window ends (day of year), for the lead set: drivers ending at least
# lead_days (14) before the median typical date of the target (as the ladder's lead set)
end_doy <- function(z) {
  e <- c(gdd5_mar31 = 90, gdd5_apr30 = 120, gdd5_may31 = 151, gdd10_apr30 = 120, gdd10_may31 = 151, gdd10_jun30 = 181,
         t_marapr = 120, t_mayjun = 181, t_julaug = 243, t_sepoct = 304, chill_days = 59, last_frz0 = 151, last_frz2 = 151,
         first_frz0 = 365, first_frz2 = 365, gdd5_at_last_frz2 = 151, p_son_prev = 0, p_djf = 59, p_mam = 151, p_jja = 243,
         p_son = 334, p_amj = 181, p_jas = 273, p_ann = 365, pet_ann = 365, vpd_mam = 151, vpd_jja = 243, events_amjjas = 273,
         dryspell_amjjas = 273, dryspell_jjas = 273, onset = 273)
  mo <- c(may1 = 121, jun1 = 152, jul1 = 182, aug1 = 213, sep1 = 244, oct1 = 274)
  v <- unname(e[z]); i <- is.na(v) & grepl("^ewma", z); v[i] <- mo[sub(".*_", "", z[i])]; v }
lead_days <- ns$lead_days %||% 14
targets <- c("greenup", "curing", "cured")
SK <- list(); DR <- list(); GR <- list()
loyo <- function(X, dcols) {                                          # X: y, year, drivers
  pred <- rep(NA_real_, nrow(X))
  for (yy in unique(X$year)) { tr <- X$year != yy
    if (sum(tr) < 8) next
    r_tr <- vapply(dcols, function(z) { v <- suppressWarnings(stats::cor(X$y[tr], X[[z]][tr], method = "spearman", use = "complete.obs")); if (is.na(v)) 0 else abs(v) }, 0)
    b <- names(which.max(r_tr)); ok <- tr & !is.na(X[[b]])
    if (sum(ok) < 6) next
    f <- stats::lm(X$y[ok] ~ X[[b]][ok]); xb <- X[[b]][!tr]
    pred[!tr] <- ifelse(is.na(xb), 0, stats::coef(f)[1] + stats::coef(f)[2] * xb) }
  ok <- !is.na(pred)
  c(n = sum(ok), skill = 1 - mean(abs(X$y[ok] - pred[ok])) / mean(abs(X$y[ok])))
}
for (tg in targets) for (pw in c("all", "C3", "C4")) {
  d <- SY[!is.na(SY[[tg]]) & (pw == "all" | SY$pathway == pw), c("site_id", "species_id", "year", "series", "region", "pathway", tg)]
  names(d)[7] <- "x"
  med <- c(tapply(d$x, d$series, cmed)); d$y <- cdiff(d$x, med[d$series]); d <- d[abs(d$y) <= 90, ]
  ny <- tapply(d$year, d$series, function(z) length(unique(z))); d <- d[d$series %in% names(ny)[ny >= min_years], ]
  if (nrow(d) < 30) next
  d$y <- d$y - stats::ave(d$y, d$series)                             # departures from each series' mean
  d <- merge(d, YA, by = c("site_id", "year"))
  typ <- cmed(d$x); lead_drv <- drv[!is.na(end_doy(drv)) & end_doy(drv) <= typ - lead_days]
  scales <- list(series = d$series, site = d$site_id, region = d$region)
  for (sc in names(scales)) for (set in c("lead", "all")) {
    dset <- if (set == "lead") lead_drv else drv
    if (length(dset) < 2) next
    g <- paste(scales[[sc]], d$year)
    if (sc == "series") X <- d[, c("year", "y", drv)] else {
      keep <- if (sc == "region") d$region != "none" else rep(TRUE, nrow(d))
      dd <- d[keep, ]; gg <- g[keep]
      ns_ <- tapply(dd$series, gg, function(z) length(unique(z))); nsite <- tapply(dd$site_id, gg, function(z) length(unique(z)))
      ok_g <- names(ns_)[ns_ >= (if (sc == "region") 3 else 2) & (sc != "region" | nsite[names(ns_)] >= 2)]
      dd <- dd[gg %in% ok_g, ]; gg <- gg[gg %in% ok_g]
      if (length(unique(gg)) < 10) next
      X <- stats::aggregate(dd[, c("y", drv)], list(group = gg), function(z) mean(z, na.rm = TRUE))
      X$year <- as.integer(sub(".* ", "", X$group)); X$n <- as.integer(table(gg)[X$group])
      wv <- tapply(dd$y, gg, stats::var); X$within_var <- as.numeric(wv[X$group])
      if (set == "lead") GR[[length(GR) + 1]] <- data.frame(target = tg, pathway = pw, scale = sc, X[, c("group", "year", "n", "y")])
    }
    k <- loyo(X, dset)
    rr <- vapply(dset, function(z) suppressWarnings(stats::cor(X$y, X[[z]], method = "spearman", use = "complete.obs")), 0)
    sig <- if (sc == "series") NA else stats::var(X$y) / mean(X$within_var / X$n, na.rm = TRUE)
    SK[[length(SK) + 1]] <- data.frame(target = tg, pathway = pw, scale = sc, drivers = set, typical = round(typ), n_drivers = length(dset), groups = nrow(X), years = length(unique(X$year)),
      sd_days = stats::sd(X$y), signal_ratio = sig, best_driver = names(which.max(abs(rr))), best_r = rr[which.max(abs(rr))],
      loyo_n = k[["n"]], skill = k[["skill"]])
    DR[[length(DR) + 1]] <- data.frame(target = tg, pathway = pw, scale = sc, drivers = set, driver = dset, r = rr, n = nrow(X))
  }
}
SK <- do.call(rbind, SK); DR <- do.call(rbind, DR); GR <- do.call(rbind, GR)
wcsv(SK, "skill.csv"); wcsv(DR, "drivers.csv"); wcsv(GR, "groups.csv")
if (!nzchar(sens$name %||% "")) { tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE); utils::write.csv(SK, file.path(tdir, "npn_scale_skill.csv"), row.names = FALSE) }
for (i in seq_len(nrow(SK))) with(SK[i, ], log_msg(sprintf("%-8s %-3s %-6s %-4s %4d groups %2d yrs, SD %4.1f d, signal ratio %s, best %-18s r %+.2f, LOYO skill %+.2f",
  target, pathway, scale, drivers, groups, years, sd_days, if (is.na(signal_ratio)) "  - " else sprintf("%4.1f", signal_ratio), best_driver, best_r, skill)))
log_msg("=== done | ", odir, " | log: ", log_file)
