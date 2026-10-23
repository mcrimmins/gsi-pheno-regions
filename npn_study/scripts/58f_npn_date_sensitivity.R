#!/usr/bin/env Rscript
# Step 58f (NPN grass curing study; laptop, lme4): are the results robust to how the curve dates
# are computed? Green-up and 50 % curing are recomputed for every usable plant-year from the
# readings (percent-green class midpoints, linear interpolation between readings) under variants:
#   own50     50 %, no limit on the gap between the two readings around the crossing (checks the
#             recomputation against the 02c dates used everywhere else)
#   th40, th60  40 % and 60 % instead of 50 %
#   gap21, gap14  50 %, crossing only if the readings around it are <= 21 / 14 days apart
#   firstcycle  50 %; curing = first fall below 50 % after green-up, not the fall after the peak
# The 02c dates (reference) and each variant go through the same steps:
#   agreement  plant-year dates vs the reference: n, share lost, median and mean absolute difference
#   typical    shift of the series typical dates (circular medians) vs the reference
#   variance   53's mixed model on the stricter set (>= 12 readings): place share of all variation;
#              of the year-to-year part, the shares shared by a site's grasses and of the single plant
#   weather    series-year medians; leave-one-year-out skill of the best single lead driver (as 55f)
# Cured dates are not recomputed (they come from the leaves-off and lowest classes in 02c).
# Output <run_dir>/npn_study/record/date_sensitivity.csv; report/tables/npn_date_sensitivity.csv.
#
# Usage: Rscript npn_study/scripts/58f_npn_date_sensitivity.R [profile]   (after 51)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs lme4")

cfg <- load_npn_cfg(profile)
log_file <- log_init("58f_npn_date_sensitivity", profile)
ns <- cfg$npn_study %||% list()
base <- file.path(cfg$run_dir, "npn_study"); rdir <- file.path(base, "record"); dir.create(rdir, showWarnings = FALSE)
t0 <- Sys.time()
log_msg("=== 58f_npn_date_sensitivity | profile: ", profile)
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id); cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$year %in% yrs & cm$pathway %in% c("C3", "C4"), ]
u$py <- paste(u$individual_id, u$year); u$region <- npn_site_regions(u$site_id); u$series <- paste(u$site_id, u$species_id)
go <- go[!is.na(go$pct_green), ]; go$py <- paste(go$individual_id, go$year); go <- go[go$py %in% u$py, ]
go <- go[order(go$py, go$doy), ]

# ---- recompute dates -------------------------------------------------------------------------------
cross <- function(d, p, th, maxgap, mode) {
  n <- length(p); up <- NA_real_; dn <- NA_real_
  if (n < 2) return(c(up, dn))
  iu <- which(p[-n] < th & p[-1] >= th & diff(d) <= maxgap)
  if (length(iu)) { i <- iu[1]; up <- d[i] + (th - p[i]) / (p[i + 1] - p[i]) * (d[i + 1] - d[i]) }
  start <- if (mode == "firstcycle") { if (length(iu)) iu[1] + 1 else 1 } else which.max(p)
  idn <- which(p[-n] >= th & p[-1] < th & diff(d) <= maxgap); idn <- idn[idn >= start]
  if (length(idn)) { i <- idn[1]; dn <- d[i] + (p[i] - th) / (p[i] - p[i + 1]) * (d[i + 1] - d[i]) }
  c(up, dn)
}
variants <- list(own50 = list(th = 50, gap = Inf, mode = "peak"), th40 = list(th = 40, gap = Inf, mode = "peak"), th60 = list(th = 60, gap = Inf, mode = "peak"),
                 gap21 = list(th = 50, gap = 21, mode = "peak"), gap14 = list(th = 50, gap = 14, mode = "peak"), firstcycle = list(th = 50, gap = Inf, mode = "firstcycle"))
idx <- split(seq_len(nrow(go)), go$py)
V <- list(reference = data.frame(py = u$py, greenup = u$g50_up_doy, curing = u$g50_down_doy))
for (vn in names(variants)) {
  v <- variants[[vn]]
  r <- t(vapply(idx, function(i) cross(go$doy[i], go$pct_green[i], v$th, v$gap, v$mode), c(0, 0)))
  V[[vn]] <- data.frame(py = names(idx), greenup = r[, 1], curing = r[, 2])
  log_msg(sprintf("  %-10s dates: green-up %d, curing %d plant-years", vn, sum(!is.na(r[, 1])), sum(!is.na(r[, 2]))))
}

# ---- drivers (as 55f): site-year departures, lead set per target -------------------------------
mt <- readRDS(file.path(base, "npn_metrics.rds")); YM <- mt$site_years
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

# ---- the same steps for the reference and every variant ------------------------------------------
fml <- y ~ 1 + (1 | site) + (1 | species) + (1 | site:species) + (1 | year) + (1 | region:year) + (1 | site:year) + (1 | site:species:year)
place <- c("site", "species", "site:species"); shared <- c("year", "region:year", "site:year")
ref_typ <- list(); OUT <- list()
for (vn in names(V)) for (tg in c("greenup", "curing")) {
  x <- merge(u[, c("py", "individual_id", "site_id", "species_id", "series", "year", "region", "pathway", "n_obs")], V[[vn]][, c("py", tg)], by = "py")
  names(x)[ncol(x)] <- "x"
  x$ref <- V$reference[[tg]][match(x$py, V$reference$py)]
  both <- !is.na(x$x) & !is.na(x$ref); dd <- abs(cdiff(x$x[both], x$ref[both]))
  x <- x[!is.na(x$x), ]
  typ <- c(tapply(x$x, x$series, cmed)); nyr_s <- tapply(x$year, x$series, function(z) length(unique(z)))
  if (vn == "reference") ref_typ[[tg]] <- typ
  sh <- cdiff(typ, ref_typ[[tg]][names(typ)])
  # variance partition (stricter set)
  s <- x[x$n_obs >= (ns$min_readings %||% 12), ]
  d <- cdiff(s$x, typ[s$series]); s <- s[abs(d) <= 90, ]; s$y <- typ[s$series] + cdiff(s$x, typ[s$series])
  s$site <- factor(s$site_id); s$species <- factor(s$species_id); s$yearf <- factor(s$year); s$region <- factor(s$region)
  vv <- tryCatch({ f <- suppressMessages(lme4::lmer(y ~ 1 + (1 | site) + (1 | species) + (1 | site:species) + (1 | yearf) + (1 | region:yearf) + (1 | site:yearf) + (1 | site:species:yearf),
                                                   data = s, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE)))
                   vc <- as.data.frame(lme4::VarCorr(f)); stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "plant", sub("yearf", "year", vc$grp))) }, error = function(e) NULL)
  pl <- if (is.null(vv)) NA else sum(vv[place]) / sum(vv); within <- if (is.null(vv)) NA else sum(vv[setdiff(names(vv), place)])
  # weather at the series-year scale, lead drivers
  sy <- stats::aggregate(y ~ series + site_id + year, s, stats::median); sy$y <- sy$y - stats::ave(sy$y, sy$series)
  ny <- table(sy$series); sy <- sy[sy$series %in% names(ny)[ny >= (ns$min_years %||% 3)], ]
  X <- merge(sy, YA, by = c("site_id", "year"))
  lead <- drv[!is.na(end_doy(drv)) & end_doy(drv) <= cmed(s$x) - (ns$lead_days %||% 14)]
  OUT[[length(OUT) + 1]] <- data.frame(variant = vn, target = tg, plant_years = nrow(x), lost = if (vn == "reference") 0 else 1 - nrow(x) / sum(!is.na(V$reference[[tg]])),
    median_abs_diff = if (vn == "reference") 0 else stats::median(dd), mean_abs_diff = if (vn == "reference") 0 else mean(dd),
    typical_shift_median = if (vn == "reference") 0 else stats::median(sh, na.rm = TRUE), typical_shift_iqr = if (vn == "reference") 0 else stats::IQR(sh, na.rm = TRUE),
    place_share = pl, shared_share = if (is.null(vv)) NA else sum(vv[shared]) / within, plant_share = if (is.null(vv)) NA else vv[["plant"]] / within,
    sd_year_to_year = sqrt(within), series_years = nrow(X), weather_skill = if (nrow(X) >= 30 && length(lead) >= 2) loyo(X, lead) else NA)
  with(OUT[[length(OUT)]], log_msg(sprintf("%-10s %-8s %5d plant-yrs, lost %3.0f %%, |diff| median %4.1f d | typical shift %+5.1f d (IQR %4.1f) | place %3.0f %%, shared %3.0f %%, plant %3.0f %% of yty (SD %4.1f d) | weather skill %+.2f",
    variant, target, plant_years, 100 * lost, median_abs_diff, typical_shift_median, typical_shift_iqr, 100 * place_share, 100 * shared_share, 100 * plant_share, sd_year_to_year, weather_skill)))
}
O <- do.call(rbind, OUT)
utils::write.csv(O, file.path(rdir, "date_sensitivity.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE); utils::write.csv(O, file.path(tdir, "npn_date_sensitivity.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(rdir, "date_sensitivity.csv"), log_file))
