#!/usr/bin/env Rscript
# Step 59e (NPN grass curing study; laptop, ~20-40 min with tuning): the operational Growing Season
# Index (NFDRS v4, Jolly et al. 2024) at every NPN grass site, against the ground record. EXPLORATORY:
# under consideration for the paper (Mike, 27 Oct).
# GSI code: the project's existing implementation, vendored in npn_study/R/gsi_nfdrs.R (run_gsi(),
# gsi_dates(), NFDRS v4 parameters). Weather: daily PRISM at the sites (50): Tmin, daily maximum VPD
# (vpdmax, hPa -> Pa), precipitation; daylength from latitude (geometric, as v4).
#
# GSI events per site-year (gsi_dates: onset = first upward crossing held >= 3 days; end = last day
# at or above the threshold after the annual peak), against the median NPN date of each series-year
# (site x species x year; relaxed set, crops out):
#   green-up  GSI onset vs NPN 50 % green-up
#   curing    GSI end vs NPN 50 % curing (also against the sustained curing date of 58f)
#   cured     GSI end at GU vs NPN cured date
# Operational thresholds: GU = 0.2 (herbaceous LFM leaves its 30 % minimum, Eq. 10) for green-up and
# cured; GSI 0.527 (herbaceous LFM 120 %, where NFDRS starts moving herbaceous load to dead) for curing.
# Also a sweep of thresholds 0.1-0.9 (does any single threshold work?).
# Tuned GSI: the same model with its ramps refitted to the NPN dates (Tmin, VPD and daylength ramps,
# precipitation window and limit, smoothing window, green-up and curing thresholds; Nelder-Mead on
# the mean absolute date error, a missing event counted as gsi: miss_days, default 60 d), for all
# grasses and for C3 and C4 separately. Scored on held-out years (fit on odd years, score on even, and
# the reverse); the reported parameters are fitted on all years.
# Also: (a) year to year: departures of GSI dates from their series mean vs the observed departures;
# (b) which sub-index binds around the observed dates, default vs tuned; (c) status: series percent
# green (>= 2 plants per date) vs GSI on the same day, raw and as weekly anomalies.
# Baseline in every table: the series' typical date (mean of the other years).
#
# Outputs <run_dir>/npn_study/gsi/: series_years.csv (obs and GSI dates per model), scores.csv,
# sweep.csv, params.csv, binding.csv, status.csv; report/tables/npn_gsi_scores.csv, npn_gsi_sweep.csv,
# npn_gsi_params.csv, npn_gsi_binding.csv, npn_gsi_status.csv.
# Config npn_study: gsi: {tune: true, maxit: 400, miss_days: 60}.
#
# Usage: Rscript npn_study/scripts/59e_npn_gsi.R [profile]   (after 50, 58, 59; ~20-40 min, mostly the 9 tuning fits)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R")); source(file.path(npn_root, "R", "gsi_nfdrs.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("59e_npn_gsi", profile)
ns <- cfg$npn_study %||% list(); gp <- ns$gsi %||% list()
do_tune <- !isFALSE(gp$tune); maxit <- as.integer(gp$maxit %||% 400); miss <- as.numeric(gp$miss_days %||% 60)
base <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(base, "gsi"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time(); set.seed(59)
log_msg(sprintf("=== 59e_npn_gsi | profile: %s | tuning %s (maxit %d), missing event = %g d", profile, do_tune, maxit, miss))
cmean <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }

# ---- observed series-years ---------------------------------------------------------------------
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id); cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"
W <- readRDS(file.path(base, "npn_site_daily.rds"))
dates <- as.Date(W$dates); yr <- as.integer(format(dates, "%Y")); dy <- as.integer(format(dates, "%j"))
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$pathway %in% c("C3", "C4") & cm$year %in% unique(yr) & cm$site_id %in% W$sites$site_id, ]
u$series <- paste(u$site_id, u$species_id); u$region <- npn_site_regions(u$site_id)
# sustained curing per plant-year (as 58f): first fall below 50 % after the peak with no reading >= 50 % within sustained_days
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds")); go <- go[!is.na(go$pct_green), ]
go$py <- paste(go$individual_id, go$year); u$py <- paste(u$individual_id, u$year); go <- go[go$py %in% u$py, ]; go <- go[order(go$py, go$doy), ]
sdays <- as.numeric(ns$sustained_days %||% 21)
sus <- vapply(split(seq_len(nrow(go)), go$py), function(i) { d <- go$doy[i]; p <- go$pct_green[i]; n <- length(p); if (n < 2) return(NA_real_)
  idn <- which(p[-n] >= 50 & p[-1] < 50); idn <- idn[idn >= which.max(p)]
  for (k in idn) { x <- d[k] + (p[k] - 50) / (p[k] - p[k + 1]) * (d[k + 1] - d[k]); j <- which(d > x & d <= x + sdays); if (!any(p[j] >= 50)) return(x) }
  NA_real_ }, 0)
u$sustained <- sus[u$py]
agg <- function(v) stats::aggregate(u[[v]], list(series = u$series, site_id = u$site_id, species_id = u$species_id, year = u$year, pathway = u$pathway, region = u$region), cmean)
SY <- Reduce(function(a, b) merge(a, b, all = TRUE), lapply(c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy", sustained = "sustained"), function(v) { x <- agg(v); names(x)[ncol(x)] <- v; x }))
names(SY)[names(SY) == "g50_up_doy"] <- "obs_greenup"; names(SY)[names(SY) == "g50_down_doy"] <- "obs_curing"; names(SY)[names(SY) == "cured_doy"] <- "obs_cured"; names(SY)[names(SY) == "sustained"] <- "obs_sustained"
SY$both_types <- SY$site_id %in% intersect(SY$site_id[SY$pathway == "C3"], SY$site_id[SY$pathway == "C4"])
log_msg(sprintf("series-years: %d (%d series, %d sites); with green-up %d, curing %d, cured %d", nrow(SY), length(unique(SY$series)), length(unique(SY$site_id)),
                sum(!is.na(SY$obs_greenup)), sum(!is.na(SY$obs_curing)), sum(!is.na(SY$obs_cured))))

# ---- weather matrices for the sites in the record ------------------------------------------------
srows <- match(sort(unique(SY$site_id)), W$sites$site_id); site_ids <- W$sites$site_id[srows]; lat <- W$sites$lat[srows]
fill <- function(X) t(apply(X, 1, function(z) if (all(is.na(z))) rep(0, length(z)) else stats::approx(seq_along(z), z, seq_along(z), rule = 2)$y))
TM <- fill(W$tmin[srows, , drop = FALSE]); VP <- fill(W$vpdmax[srows, , drop = FALSE]) * 100; PR <- W$ppt[srows, , drop = FALSE]; PR[is.na(PR)] <- 0
PH <- t(vapply(lat, function(la) photoperiod_seconds(dy, la, "geometric"), numeric(length(dy))))
CP <- t(apply(PR, 1, cumsum)); CP <- cbind(0, CP)
rowshift_sum <- function(C, n) { D <- ncol(C) - 1; j <- seq_len(D); C[, j + 1, drop = FALSE] - C[, pmax(0, j - n) + 1, drop = FALSE] }
trail_mean <- function(X, n) { C <- cbind(0, t(apply(X, 1, cumsum))); j <- seq_len(ncol(X)); lo <- pmax(0, j - n); (C[, j + 1, drop = FALSE] - C[, lo + 1, drop = FALSE]) / rep(j - lo, each = nrow(X)) }
# fast matrix GSI with the vendored ramps; checked against run_gsi() below
gsi_mat <- function(p, keep_idx = FALSE) {
  dm <- function(v) { dim(v) <- dim(TM); v }                          # pmax() in the ramps drops matrix dimensions
  it <- dm(idx_tmin(TM, p$tmin_lo, p$tmin_hi)); iv <- dm(idx_vpd(VP, p$vpd_lo, p$vpd_hi)); ip <- dm(idx_photo(PH, p$photo_lo, p$photo_hi))
  iq <- dm(gsi_ramp(rowshift_sum(CP, round(p$precip_window)), p$precip_lo, p$precip_hi))
  ig <- it * iv * ip * iq
  g <- dm(pmin(1, trail_mean(ig, round(p$window)) / p$gsi_max))
  if (keep_idx) list(g = g, it = it, iv = iv, ip = ip, iq = iq) else g
}
p4 <- npn_gsi_params_v4()
chk <- run_gsi(data.frame(date = dates, doy = dy, tmin_c = TM[1, ], vpd_pa = VP[1, ], precip_mm = PR[1, ]), lat[1], p4)
G4 <- gsi_mat(p4)
log_msg(sprintf("check: matrix GSI vs run_gsi() at the first site, max |difference| %.2e", max(abs(G4[1, ] - chk$gsi_rel))))

# ---- events per site-year -------------------------------------------------------------------------
SY$srow <- match(SY$site_id, site_ids)
sy_keys <- unique(SY[, c("srow", "year")])
cols_of <- split(seq_along(yr), yr)
events <- function(G, thr_up, thr_dn, keys = sy_keys) {
  o <- matrix(NA_real_, nrow(keys), 2)
  for (k in seq_len(nrow(keys))) {
    j <- cols_of[[as.character(keys$year[k])]]; g <- G[keys$srow[k], j]
    a <- g >= thr_up; r <- rle(a); e <- cumsum(r$lengths); s <- e - r$lengths + 1
    ok <- r$values & r$lengths >= 3 & s > 1; if (any(ok)) o[k, 1] <- s[ok][1]
    ipk <- which.max(g); af <- which(g >= thr_dn & seq_along(g) >= ipk); if (length(af)) o[k, 2] <- max(af)
  }
  data.frame(srow = keys$srow, year = keys$year, up = o[, 1], dn = o[, 2])
}
attach_ev <- function(E, prefix) { i <- match(paste(SY$srow, SY$year), paste(E$srow, E$year)); SY[[paste0(prefix, "_up")]] <<- E$up[i]; SY[[paste0(prefix, "_dn")]] <<- E$dn[i] }
thr120 <- 0.2 + 0.8 * (120 - 30) / (250 - 30)                          # GSI at which v4 herbaceous LFM = 120 %
attach_ev(events(G4, 0.2, thr120), "v4")                                # green-up at GU, curing at LFM 120 %
SY$v4_cured <- events(G4, 0.2, 0.2)$dn[match(paste(SY$srow, SY$year), paste(sy_keys$srow, sy_keys$year))]
log_msg(sprintf("v4 default: GSI reaches GU in %.0f %% of site-years; falls through LFM 120 %% (GSI %.3f) after the peak in %.0f %%",
                100 * mean(!is.na(SY$v4_up)), thr120, 100 * mean(!is.na(SY$v4_dn))))

# ---- scoring helpers --------------------------------------------------------------------------------
typical <- function(obs, series, year) vapply(seq_along(obs), function(i) { o <- obs[series == series[i] & year != year[i]]; cmean(o) }, 0)
for (tg in c("greenup", "curing", "cured", "sustained")) SY[[paste0("typ_", tg)]] <- typical(SY[[paste0("obs_", tg)]], SY$series, SY$year)
score <- function(pred, obs, X, label, target, set) {
  ok <- !is.na(obs); if (!any(ok)) return(NULL)
  e <- cdiff(pred[ok], obs[ok]); miss_n <- sum(is.na(e)); ae <- ifelse(is.na(e), miss, abs(e))
  typ <- X[[paste0("typ_", target)]][ok]; et <- abs(cdiff(typ, obs[ok])); both <- !is.na(e) & !is.na(et)
  # departures from each series' own mean, series with >= 3 years of both
  s <- X$series[ok]; pd <- pred[ok]; od <- obs[ok]
  dep <- do.call(rbind, lapply(split(seq_along(s), s), function(k) { k <- k[!is.na(pd[k])]; if (length(k) < 3) return(NULL)
    data.frame(p = cdiff(pd[k], cmean(pd[k])), o = cdiff(od[k], cmean(od[k]))) }))
  data.frame(model = label, target = target, set = set, n = sum(ok), never = miss_n / sum(ok), mae = mean(ae), bias = stats::median(e, na.rm = TRUE),
             mae_typical = mean(et, na.rm = TRUE), mae_both = mean(abs(e[both])), mae_typical_both = mean(et[both]),
             skill_vs_typical = 1 - mean(ae) / mean(et, na.rm = TRUE), dep_r = if (!is.null(dep) && nrow(dep) >= 10) stats::cor(dep$p, dep$o) else NA, dep_n = if (is.null(dep)) 0 else nrow(dep))
}
sets <- list(all = rep(TRUE, nrow(SY)), C3 = SY$pathway == "C3", C4 = SY$pathway == "C4", C3_at_both = SY$pathway == "C3" & SY$both_types, C4_at_both = SY$pathway == "C4" & SY$both_types)
SC <- list()
add_scores <- function(lab, up, dn, cured = NULL, use = names(sets)) for (sn in use) { k <- sets[[sn]]; X <- SY[k, ]
  SC[[length(SC) + 1]] <<- score(up[k], X$obs_greenup, X, lab, "greenup", sn)
  SC[[length(SC) + 1]] <<- score(dn[k], X$obs_curing, X, lab, "curing", sn)
  SC[[length(SC) + 1]] <<- score(dn[k], X$obs_sustained, X, lab, "sustained", sn)
  if (!is.null(cured)) SC[[length(SC) + 1]] <<- score(cured[k], X$obs_cured, X, lab, "cured", sn) }
add_scores("v4 default", SY$v4_up, SY$v4_dn, SY$v4_cured)
for (r in Filter(Negate(is.null), SC)) if (r$set %in% c("all", "C3", "C4")) with(r, log_msg(sprintf("  %-12s %-9s %-10s n %4d | MAE %5.1f d (typical %5.1f) | bias %+6.1f d | never %3.0f %% | departure r %+.2f (%d)",
  model, target, set, n, mae, mae_typical, bias, 100 * never, dep_r, dep_n)))

# ---- threshold sweep (default inputs) ----------------------------------------------------------------
thr <- seq(0.1, 0.9, by = 0.05)
SW <- do.call(rbind, lapply(thr, function(t) { E <- events(G4, t, t); i <- match(paste(SY$srow, SY$year), paste(E$srow, E$year))
  do.call(rbind, lapply(c("all", "C3", "C4"), function(sn) { k <- sets[[sn]]
    do.call(rbind, lapply(list(c("greenup", "up"), c("curing", "dn"), c("cured", "dn")), function(z) { o <- SY[[paste0("obs_", z[1])]][k]; p <- E[[z[2]]][i][k]; ok <- !is.na(o)
      e <- cdiff(p[ok], o[ok]); data.frame(thr = t, set = sn, target = z[1], n = sum(ok), never = mean(is.na(e)), mae = mean(ifelse(is.na(e), miss, abs(e))), bias = stats::median(e, na.rm = TRUE)) })) })) }))
best <- do.call(rbind, lapply(split(SW, paste(SW$set, SW$target)), function(d) d[which.min(d$mae), ]))
for (i in seq_len(nrow(best))) with(best[i, ], log_msg(sprintf("  best single threshold %-3s %-8s GSI %.2f: MAE %5.1f d, bias %+5.1f d, never %3.0f %%", set, target, thr, mae, bias, 100 * never)))
# v4 inputs with the best thresholds per target (all grasses): separates the threshold from the inputs
bu <- best$thr[best$set == "all" & best$target == "greenup"]; bd <- best$thr[best$set == "all" & best$target == "curing"]
Eb <- events(G4, bu, bd); ib <- match(paste(SY$srow, SY$year), paste(Eb$srow, Eb$year)); SY$v4b_up <- Eb$up[ib]; SY$v4b_dn <- Eb$dn[ib]
add_scores(sprintf("v4, best thresholds (%.2f / %.2f)", bu, bd), SY$v4b_up, SY$v4b_dn)

# ---- tuning -------------------------------------------------------------------------------------------
# theta: tmin_lo, log tmin width, vpd_lo (kPa), log vpd width, photo_lo (h), log photo width, log precip_hi,
#        log precip window, log smoothing window, logit thr_up, logit thr_dn
lg <- function(x) log(x / (1 - x)); il <- function(z) 1 / (1 + exp(-z))
to_p <- function(th) modifyList(p4, list(tmin_lo = th[1], tmin_hi = th[1] + exp(th[2]), vpd_lo = 1000 * th[3], vpd_hi = 1000 * (th[3] + exp(th[4])),
  photo_lo = 3600 * th[5], photo_hi = 3600 * (th[5] + exp(th[6])), precip_hi = exp(th[7]), precip_window = min(180, max(3, round(exp(th[8])))),
  window = min(90, max(3, round(exp(th[9])))), thr_up = il(th[10]), thr_dn = il(th[11])))
th0 <- c(-2, log(7), 0.9, log(3.2), 10, log(1), log(10), log(28), log(28), lg(0.2), lg(thr120))
obj <- function(th, k) {
  p <- to_p(th); if (p$precip_hi < 0.5 || p$tmin_hi - p$tmin_lo < 0.2) return(1e6)
  G <- gsi_mat(p); keys <- unique(SY[k, c("srow", "year")]); E <- events(G, p$thr_up, p$thr_dn, keys)
  i <- match(paste(SY$srow[k], SY$year[k]), paste(E$srow, E$year))
  e1 <- cdiff(E$up[i], SY$obs_greenup[k]); e2 <- cdiff(E$dn[i], SY$obs_curing[k])
  a <- c(ifelse(is.na(e1), miss, abs(e1))[!is.na(SY$obs_greenup[k])], ifelse(is.na(e2), miss, abs(e2))[!is.na(SY$obs_curing[k])])
  mean(a)
}
fit <- function(k) {                                                    # Nelder-Mead, restarted once from its own best point
  f <- stats::optim(th0, obj, k = k, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-4))
  g <- stats::optim(f$par, obj, k = k, method = "Nelder-Mead", control = list(maxit = maxit, reltol = 1e-4))
  if (g$value < f$value) { g$counts[1] <- g$counts[1] + f$counts[1]; f <- g } else f$counts[1] <- f$counts[1] + g$counts[1]
  list(th = f$par, value = f$value, evals = f$counts[1]) }
PAR <- list(); BIND <- list()
pdesc <- function(p, lab) data.frame(model = lab, tmin_lo_C = p$tmin_lo, tmin_hi_C = p$tmin_hi, vpd_lo_kPa = p$vpd_lo / 1000, vpd_hi_kPa = p$vpd_hi / 1000,
  daylength_lo_h = p$photo_lo / 3600, daylength_hi_h = p$photo_hi / 3600, precip_window_d = p$precip_window, precip_hi_mm = p$precip_hi, smoothing_d = p$window,
  thr_greenup = p$thr_up %||% 0.2, thr_curing = p$thr_dn %||% thr120)
PAR[[1]] <- pdesc(p4, "v4 default")
if (do_tune) {
  odd <- SY$year %% 2 == 1
  for (sn in c("all", "C3", "C4")) {
    k <- sets[[sn]] & (!is.na(SY$obs_greenup) | !is.na(SY$obs_curing))
    t1 <- Sys.time()
    f_odd <- fit(k & odd); f_even <- fit(k & !odd); f_all <- fit(k)
    log_msg(sprintf("  tuned %-3s: training MAE odd %.1f / even %.1f / all %.1f d (%d, %d, %d evaluations; %.1f min)", sn, f_odd$value, f_even$value, f_all$value,
                    f_odd$evals, f_even$evals, f_all$evals, as.numeric(difftime(Sys.time(), t1, units = "mins"))))
    # held-out predictions: even years from the odd-year fit and the reverse
    pr <- function(f, rows) { p <- to_p(f$th); G <- gsi_mat(p); keys <- unique(SY[rows, c("srow", "year")]); E <- events(G, p$thr_up, p$thr_dn, keys)
      i <- match(paste(SY$srow[rows], SY$year[rows]), paste(E$srow, E$year)); cbind(E$up[i], E$dn[i]) }
    up <- rep(NA_real_, nrow(SY)); dn <- up
    r1 <- which(sets[[sn]] & !odd); r2 <- which(sets[[sn]] & odd)
    if (length(r1)) { z <- pr(f_odd, r1); up[r1] <- z[, 1]; dn[r1] <- z[, 2] }
    if (length(r2)) { z <- pr(f_even, r2); up[r2] <- z[, 1]; dn[r2] <- z[, 2] }
    lab <- sprintf("tuned on %s (held-out years)", sn)
    SY[[paste0("tuned_", sn, "_up")]] <- up; SY[[paste0("tuned_", sn, "_dn")]] <- dn
    add_scores(lab, up, dn, use = if (sn == "all") names(sets) else grep(sn, names(sets), value = TRUE))
    PAR[[length(PAR) + 1]] <- cbind(pdesc(to_p(f_all$th), paste("tuned on", sn)), train_mae = f_all$value)
    PAR[[length(PAR) + 1]] <- cbind(pdesc(to_p(f_odd$th), paste("tuned on", sn, "(odd years)")), train_mae = f_odd$value)
    PAR[[length(PAR) + 1]] <- cbind(pdesc(to_p(f_even$th), paste("tuned on", sn, "(even years)")), train_mae = f_even$value)
    assign(paste0("p_", sn), to_p(f_all$th))
  }
}
if (do_tune) {                                                          # C3 rows from the C3 fit, C4 rows from the C4 fit
  up <- ifelse(SY$pathway == "C3", SY$tuned_C3_up, SY$tuned_C4_up); dn <- ifelse(SY$pathway == "C3", SY$tuned_C3_dn, SY$tuned_C4_dn)
  add_scores("tuned by grass type (held-out years)", up, dn)
}
PAR <- do.call(rbind, lapply(PAR, function(x) { if (!"train_mae" %in% names(x)) x$train_mae <- NA; x }))

# ---- which sub-index binds around the observed dates ---------------------------------------------------
bind_tab <- function(p, lab) {
  L <- gsi_mat(p, keep_idx = TRUE)
  do.call(rbind, lapply(c("greenup", "curing"), function(tg) do.call(rbind, lapply(c("C3", "C4"), function(pw) {
    X <- SY[SY$pathway == pw & !is.na(SY[[paste0("obs_", tg)]]), ]; if (!nrow(X)) return(NULL)
    cnt <- c(tmin = 0, vpd = 0, photo = 0, precip = 0); lvl <- numeric(0)
    for (i in seq_len(nrow(X))) { j0 <- match(as.Date(sprintf("%d-01-01", X$year[i])) + round(X[[paste0("obs_", tg)]][i]) - 1, dates); if (is.na(j0)) next
      j <- max(1, j0 - 14):min(length(dates), j0 + 14); r <- X$srow[i]
      M <- cbind(tmin = L$it[r, j], vpd = L$iv[r, j], photo = L$ip[r, j], precip = L$iq[r, j]); b <- colnames(M)[max.col(-M, ties.method = "first")]
      cnt <- cnt + as.numeric(table(factor(b, names(cnt)))); lvl <- c(lvl, L$g[r, j0]) }
    data.frame(model = lab, target = tg, pathway = pw, as.list(stats::setNames(as.numeric(cnt) / sum(cnt), c("tmin", "vpd", "photo", "precip"))), gsi_at_obs_median = stats::median(lvl), gsi_at_obs_iqr = stats::IQR(lvl), check.names = FALSE) }))))
}
BIND <- rbind(bind_tab(p4, "v4 default"), if (do_tune) rbind(bind_tab(p_all, "tuned on all"), bind_tab(p_C3, "tuned on C3"), bind_tab(p_C4, "tuned on C4")))
for (i in seq_len(nrow(BIND))) with(BIND[i, ], log_msg(sprintf("  binding %-14s %-8s %s: Tmin %3.0f %%, VPD %3.0f %%, daylength %3.0f %%, precip %3.0f %% | GSI on the observed date %.2f",
  model, target, pathway, 100 * tmin, 100 * vpd, 100 * photo, 100 * precip, gsi_at_obs_median)))

# ---- status: series percent green vs GSI on the same day -------------------------------------------------
go$series <- paste(go$site_id, go$species_id)
SC2 <- stats::aggregate(pct_green ~ series + site_id + year + doy, go, mean); SC2$n <- stats::aggregate(pct_green ~ series + site_id + year + doy, go, length)$pct_green
SC2 <- SC2[SC2$n >= 2 & SC2$site_id %in% site_ids, ]; SC2$pathway <- u$pathway[match(SC2$series, u$series)]
jj <- match(as.Date(sprintf("%d-01-01", SC2$year)) + SC2$doy - 1, dates); SC2 <- SC2[!is.na(jj), ]; jj <- jj[!is.na(jj)]
SC2$g_v4 <- G4[cbind(match(SC2$site_id, site_ids), jj)]
if (do_tune) { Gt <- list(C3 = gsi_mat(p_C3), C4 = gsi_mat(p_C4)); SC2$g_tuned <- ifelse(SC2$pathway == "C3", Gt$C3[cbind(match(SC2$site_id, site_ids), jj)], Gt$C4[cbind(match(SC2$site_id, site_ids), jj)]) }
SC2$week <- (SC2$doy - 1) %/% 7
an <- function(v, X) v - stats::ave(v, X$series, X$week)
ST <- do.call(rbind, lapply(c("all", "C3", "C4"), function(pw) { X <- SC2[pw == "all" | SC2$pathway == pw, ]
  do.call(rbind, lapply(intersect(c("g_v4", "g_tuned"), names(X)), function(gc) data.frame(set = pw, model = c(g_v4 = "v4 default", g_tuned = "tuned by grass type")[[gc]], n = nrow(X),
    series = length(unique(X$series)), r_level = stats::cor(X$pct_green, X[[gc]], method = "spearman"),
    r_weekly_anomaly = stats::cor(an(X$pct_green, X), an(X[[gc]], X), use = "complete.obs"),
    share_green_gsi_below_gu = mean(X[[gc]][X$pct_green >= 50] < 0.2), share_brown_gsi_above_120 = mean(X[[gc]][X$pct_green < 25] >= thr120)))) }))
for (i in seq_len(nrow(ST))) with(ST[i, ], log_msg(sprintf("  status %-3s %-20s n %6d: r level %.2f, r weekly anomaly %+.2f | %.0f %% of >= 50 %% green dates have GSI below GU; %.0f %% of < 25 %% green dates have LFM >= 120 %%",
  set, model, n, r_level, r_weekly_anomaly, 100 * share_green_gsi_below_gu, 100 * share_brown_gsi_above_120)))

# ---- write ------------------------------------------------------------------------------------------------
SCt <- do.call(rbind, Filter(Negate(is.null), SC)); rownames(SCt) <- NULL
utils::write.csv(SY[, setdiff(names(SY), "srow")], file.path(odir, "series_years.csv"), row.names = FALSE)
utils::write.csv(SCt, file.path(odir, "scores.csv"), row.names = FALSE); utils::write.csv(SW, file.path(odir, "sweep.csv"), row.names = FALSE)
utils::write.csv(PAR, file.path(odir, "params.csv"), row.names = FALSE); utils::write.csv(BIND, file.path(odir, "binding.csv"), row.names = FALSE)
utils::write.csv(ST, file.path(odir, "status.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(SCt, file.path(tdir, "npn_gsi_scores.csv"), row.names = FALSE); utils::write.csv(SW, file.path(tdir, "npn_gsi_sweep.csv"), row.names = FALSE)
utils::write.csv(PAR, file.path(tdir, "npn_gsi_params.csv"), row.names = FALSE); utils::write.csv(BIND, file.path(tdir, "npn_gsi_binding.csv"), row.names = FALSE)
utils::write.csv(ST, file.path(tdir, "npn_gsi_status.csv"), row.names = FALSE)
for (r in split(SCt, seq_len(nrow(SCt)))) if (r$set %in% c("all", "C3", "C4") && r$model != "v4 default" && r$target %in% c("greenup", "curing")) with(r, log_msg(sprintf("  %-34s %-8s %-3s MAE %5.1f d (typical %5.1f) | bias %+6.1f | never %3.0f %% | departure r %+.2f",
  model, target, set, mae, mae_typical, bias, 100 * never, dep_r)))
for (i in seq_len(nrow(PAR))) with(PAR[i, ], log_msg(sprintf("  params %-26s Tmin %5.1f..%5.1f C | VPD %4.2f..%4.2f kPa | day %5.2f..%5.2f h | precip %3d d, %5.1f mm | smooth %2d d | thr %.2f / %.2f",
  model, tmin_lo_C, tmin_hi_C, vpd_lo_kPa, vpd_hi_kPa, daylength_lo_h, daylength_hi_h, as.integer(precip_window_d), precip_hi_mm, as.integer(smoothing_d), thr_greenup, thr_curing)))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
