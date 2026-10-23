#!/usr/bin/env Rscript
# Step 58g (NPN grass curing study; laptop, lme4): 95 % intervals for the headline numbers by
# resampling clusters of nearby sites (single-linkage within 10 km, as 58e) with replacement, so
# neighbouring plots go in or out together. A resampled copy of a cluster gets new site labels so
# the mixed models treat it as a separate place. Statistics (relaxed set unless noted):
#   variance  53's partition on the stricter set (>= 12 readings): place share of all variation;
#             of the year-to-year part, the shares shared by a site's grasses and of the single plant
#             (green-up, 50 % curing, cured; all grasses)
#   weather   leave-one-year-out skill of the best single lead driver at the series-year scale (55f)
#   climate   days per SD of March-April temperature in the site-climate model of warm-season typical
#             green-up (58 Q4, species random effect)
#   c3c4      median gap C4 green-up - C3 cured over prairie sites with both (EPA 9.2-9.4; 58b)
# Output <run_dir>/npn_study/record/bootstrap.csv (statistic, estimate, lo, hi, B, valid) and
# report/tables/npn_bootstrap.csv.
#
# Usage: Rscript npn_study/scripts/58g_npn_bootstrap.R [profile] [B]   (default B = 200; after 51, 58, 58b)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
B <- if (length(args) >= 2) as.integer(args[2]) else 200L
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs lme4")

cfg <- load_npn_cfg(profile)
log_file <- log_init("58g_npn_bootstrap", profile)
ns <- cfg$npn_study %||% list()
base <- file.path(cfg$run_dir, "npn_study"); rdir <- file.path(base, "record")
t0 <- Sys.time(); set.seed(58)
log_msg("=== 58g_npn_bootstrap | profile: ", profile, " | B = ", B)
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
ctl <- lme4::lmerControl(calc.derivs = FALSE)

# ---- inputs ---------------------------------------------------------------------------------------
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id); cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$year %in% yrs & cm$pathway %in% c("C3", "C4") & cm$n_obs >= (ns$min_readings %||% 12), ]
u$region <- npn_site_regions(u$site_id); u$series <- paste(u$site_id, u$species_id)
cols <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
TY <- utils::read.csv(file.path(rdir, "series_typical.csv"), stringsAsFactors = FALSE, colClasses = c(region = "character"))
mt <- readRDS(file.path(base, "npn_metrics.rds")); N <- mt$normals; YM <- mt$site_years
TN <- merge(TY, N[, c("site_id", "t_marapr", "aridity", "p_jas_share")], by = "site_id")
TN <- TN[TN$pathway == "C4" & !is.na(TN$greenup) & stats::complete.cases(TN[, c("t_marapr", "aridity", "p_jas_share", "lat")]) & TN$region != "none", ]
TN$wt <- pmin(TN$years, 5)
gap_f <- file.path(rdir, "c3c4_gap_eco", "sites.csv")
GS <- if (file.exists(gap_f)) { g <- utils::read.csv(gap_f, stringsAsFactors = FALSE, colClasses = c(region = "character")); g[g$both %in% TRUE & g$region %in% c("9.2", "9.3", "9.4") & !is.na(g$gap_cured), ] } else NULL

# drivers (as 55f): site-year departures; lead set per target
drv <- setdiff(names(YM)[vapply(YM, is.numeric, TRUE)], c("site_id", "year"))
nyr <- ns$normals %||% list(start = 2002, end = 2025); YMn <- YM[YM$year >= nyr$start & YM$year <= nyr$end, ]
mu <- stats::aggregate(YMn[, drv], list(site_id = YMn$site_id), function(z) mean(z, na.rm = TRUE))
YA <- YM[, c("site_id", "year", drv)]; YA[, drv] <- YA[, drv] - mu[match(YA$site_id, mu$site_id), drv]
end_doy <- function(z) {
  e <- c(gdd5_mar31 = 90, gdd5_apr30 = 120, gdd5_may31 = 151, gdd10_apr30 = 120, gdd10_may31 = 151, gdd10_jun30 = 181, t_marapr = 120, t_mayjun = 181,
         t_julaug = 243, t_sepoct = 304, chill_days = 59, last_frz0 = 151, last_frz2 = 151, first_frz0 = 365, first_frz2 = 365, gdd5_at_last_frz2 = 151,
         p_son_prev = 0, p_djf = 59, p_mam = 151, p_jja = 243, p_son = 334, p_amj = 181, p_jas = 273, p_ann = 365, pet_ann = 365, vpd_mam = 151,
         vpd_jja = 243, events_amjjas = 273, dryspell_amjjas = 273, dryspell_jjas = 273, onset = 273)
  mo <- c(may1 = 121, jun1 = 152, jul1 = 182, aug1 = 213, sep1 = 244, oct1 = 274)
  v <- unname(e[z]); i <- is.na(v) & grepl("^ewma", z); v[i] <- mo[sub(".*_", "", z[i])]; v }
loyo <- function(X, dcols) {
  pred <- rep(NA_real_, nrow(X))
  for (yy in unique(X$year)) { tr <- X$year != yy; if (sum(tr) < 8) next
    r_tr <- vapply(dcols, function(z) { v <- suppressWarnings(stats::cor(X$y[tr], X[[z]][tr], method = "spearman", use = "complete.obs")); if (is.na(v)) 0 else abs(v) }, 0)
    b <- names(which.max(r_tr)); ok <- tr & !is.na(X[[b]]); if (sum(ok) < 6) next
    f <- stats::lm(X$y[ok] ~ X[[b]][ok]); xb <- X[[b]][!tr]; pred[!tr] <- ifelse(is.na(xb), 0, stats::coef(f)[1] + stats::coef(f)[2] * xb) }
  ok <- !is.na(pred); 1 - mean(abs(X$y[ok] - pred[ok])) / mean(abs(X$y[ok]))
}

# ---- clusters of sites within 10 km --------------------------------------------------------------
S <- unique(rbind(u[, c("site_id", "lat", "lon")], TY[, c("site_id", "lat", "lon")])); S <- stats::aggregate(cbind(lat, lon) ~ site_id, S, stats::median)
la <- S$lat * pi / 180; lo <- S$lon * pi / 180
cz <- outer(sin(la), sin(la)) + outer(cos(la), cos(la)) * cos(outer(lo, lo, "-")); cz[cz > 1] <- 1; cz[cz < -1] <- -1
clu <- stats::setNames(stats::cutree(stats::hclust(stats::as.dist(6371 * acos(cz)), method = "single"), h = 10), S$site_id)
K <- sort(unique(clu))
log_msg(sprintf("%d sites in %d clusters (10 km)", nrow(S), length(K)))

# ---- statistics on one (resampled) data set ------------------------------------------------------
# relabel(): a data frame with site_id -> copies with new site labels for each draw of its cluster
relabel <- function(d, draw) {
  if (is.null(draw)) return(d)
  do.call(rbind, lapply(seq_along(draw), function(j) { x <- d[clu[as.character(d$site_id)] == draw[j], , drop = FALSE]
    if (!nrow(x)) return(NULL); x$site_lab <- paste0(x$site_id, "_", j); x }))
}
stats_one <- function(draw = NULL) {
  out <- c()
  uu <- relabel(u, draw); if (is.null(draw)) uu$site_lab <- as.character(uu$site_id)
  for (tg in names(cols)) {
    x <- uu[!is.na(uu[[cols[[tg]]]]), ]; x$x <- x[[cols[[tg]]]]; x$ser <- paste(x$site_lab, x$species_id)
    med <- c(tapply(x$x, x$ser, cmed)); d <- cdiff(x$x, med[x$ser]); x <- x[abs(d) <= 90, ]; x$y <- med[x$ser] + cdiff(x$x, med[x$ser])
    x$site <- factor(x$site_lab); x$species <- factor(x$species_id); x$yearf <- factor(x$year); x$regionf <- factor(x$region)
    v <- tryCatch({ f <- suppressMessages(suppressWarnings(lme4::lmer(y ~ 1 + (1 | site) + (1 | species) + (1 | site:species) + (1 | yearf) + (1 | regionf:yearf) + (1 | site:yearf) + (1 | site:species:yearf),
                                                                      data = x, REML = TRUE, control = ctl)))
                    vc <- as.data.frame(lme4::VarCorr(f)); stats::setNames(vc$vcov, vc$grp) }, error = function(e) NULL)
    if (!is.null(v)) { pl <- c("site", "species", "site:species"); sh <- c("yearf", "regionf:yearf", "site:yearf"); w <- sum(v[setdiff(names(v), pl)])
      out[paste0("variance_place_", tg)] <- sum(v[pl]) / sum(v); out[paste0("variance_shared_", tg)] <- sum(v[sh]) / w; out[paste0("variance_plant_", tg)] <- v[["Residual"]] / w }
    sy <- stats::aggregate(y ~ ser + site_id + year, x, stats::median); sy$y <- sy$y - stats::ave(sy$y, sy$ser)
    ny <- table(sy$ser); sy <- sy[sy$ser %in% names(ny)[ny >= (ns$min_years %||% 3)], ]
    X <- merge(sy, YA, by = c("site_id", "year"))
    lead <- drv[!is.na(end_doy(drv)) & end_doy(drv) <= cmed(x$x) - (ns$lead_days %||% 14)]
    if (nrow(X) >= 30) out[paste0("weather_skill_", tg)] <- loyo(X, lead)
  }
  tn <- relabel(TN, draw); if (is.null(draw)) tn$site_lab <- as.character(tn$site_id)
  if (nrow(tn) >= 20) {
    tn$y <- cmed(tn$greenup) + cdiff(tn$greenup, cmed(tn$greenup)); for (z in c("t_marapr", "aridity", "p_jas_share", "lat")) tn[[paste0("z_", z)]] <- (tn[[z]] - mean(tn[[z]])) / stats::sd(tn[[z]])
    f <- tryCatch(suppressMessages(lme4::lmer(y ~ z_t_marapr + z_aridity + z_p_jas_share + z_lat + (1 | species_id), data = tn, weights = wt, REML = TRUE, control = ctl)), error = function(e) NULL)
    if (!is.null(f)) out["climate_c4_greenup_days_per_sd_t_marapr"] <- lme4::fixef(f)[["z_t_marapr"]]
  }
  if (!is.null(GS)) { g <- relabel(GS, draw); if (!is.null(g) && nrow(g)) out["c3c4_gap_prairie_days"] <- stats::median(g$gap_cured) }
  out
}

est <- stats_one(NULL)
log_msg("full data: ", paste(sprintf("%s %.2f", names(est), est), collapse = ", "))
R <- matrix(NA_real_, B, length(est), dimnames = list(NULL, names(est)))
for (b in seq_len(B)) {
  dr <- sample(K, length(K), replace = TRUE)
  s <- tryCatch(stats_one(dr), error = function(e) c())
  R[b, intersect(names(s), colnames(R))] <- s[intersect(names(s), colnames(R))]
  if (b %% 20 == 0) log_msg(sprintf("  %d / %d (%.1f min)", b, B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
O <- data.frame(statistic = names(est), estimate = unname(est), lo = apply(R, 2, stats::quantile, 0.025, na.rm = TRUE),
                hi = apply(R, 2, stats::quantile, 0.975, na.rm = TRUE), B = B, valid = colSums(!is.na(R)))
utils::write.csv(O, file.path(rdir, "bootstrap.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE); utils::write.csv(O, file.path(tdir, "npn_bootstrap.csv"), row.names = FALSE)
for (i in seq_len(nrow(O))) with(O[i, ], log_msg(sprintf("  %-44s %7.2f  [%7.2f, %7.2f]  (%d valid)", statistic, estimate, lo, hi, valid)))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(rdir, "bootstrap.csv"), log_file))
