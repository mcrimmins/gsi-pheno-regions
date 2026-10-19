#!/usr/bin/env Rscript
# Step 52 (NPN grass curing study): the model ladder, scored once per series-year (site x species
# x year) with leave-one-year-out cross-validation (every series loses the same year at once).
# Targets: green-up, curing (50 % down), cured. Crops excluded. A series needs >= 3 years with the
# target. Each date is unwrapped around its series' circular median; |departure| > 90 d counts as a
# season mismatch (reported, not scored), as in the satellite and ground checks.
#   L0  typical: the series' mean over the training years
#   L1  one driver (green-up: Mar-Apr mean temperature; curing and cured: Jun-Aug rain), slope per
#       pathway, on departures from the site's climatology
#   L2  process rule, threshold per series from its training years: green-up = day GDD5 (from
#       Jan 1) reaches the threshold; curing / cured = first day from typical green-up + 30 d on which
#       the 30-day P - PET memory falls below the threshold, or the first fall freeze (-2.2 C) if earlier
#   L3  forward-selected linear model on driver departures (<= 4 drivers, no intercept), per pathway
#   L4  random forest on the same departures plus site climate (aridity, warm-season rain share,
#       latitude), pathway and life cycle, pooled
# Two driver sets for curing / cured: "lead" uses only windows that end by the series' typical date
# minus 14 d (spring metrics, green-up-anchored rain, memory and heat at the cutoff), the
# operational case; "hindcast" adds the fixed summer and fall windows, which can run past the date.
# Inputs: <run_dir>/npn_study/npn_metrics.rds (51), npn_site_daily.rds (50).
# Output <run_dir>/npn_study/ladder/: predictions.csv (one row per series-year x target x set x
# level), scores.csv (MAE, skill vs L0, r of departures; all, by pathway, by region), selected.csv
# (L3 drivers per fold), importance.csv (L4).
#
# Usage: Rscript npn_study/scripts/52_npn_ladder.R [profile] [sensitivity]
#   sensitivity: a name under config npn_study: sensitivity (e.g. no1617); outputs then go to
#   <run_dir>/npn_study/sens/<name>/ladder/

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "full")

npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("52_npn_ladder", profile)
ns <- cfg$npn_study %||% list()
sens <- npn_sensitivity(ns, args)
sdir <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(sdir, sens$sub, "ladder")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
lead <- ns$lead_days %||% 14; min_years <- ns$min_years %||% 3; max_pred <- ns$max_predictors %||% 4
log_msg("=== 52_npn_ladder | profile: ", profile, if (nzchar(sens$name)) sprintf(" | sensitivity %s: excluding %s", sens$name, paste(sens$exclude_years, collapse = ", ")) else "")
set.seed(52)

M <- readRDS(file.path(sdir, "npn_metrics.rds"))
D <- readRDS(file.path(sdir, "npn_site_daily.rds"))
S <- M$series_years; S <- S[!(S$crop %in% "yes") & !(S$year %in% sens$exclude_years), ]
S <- cbind(S, M$anchor[match(paste(S$site_id, S$species_id, S$year), paste(M$anchor$site_id, M$anchor$species_id, M$anchor$year)), c("p_gu30", "p_gu60", "p_gu90", "dry_gu90", "events_gu90")])
YM <- M$site_years; N <- M$normals
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }

# daily fields for L2 and the lead-time drivers
dates <- D$dates; yr <- as.integer(format(dates, "%Y")); doy <- as.integer(format(dates, "%j"))
TM <- (D$tmin + D$tmax) / 2
PET <- hargreaves(D$tmin, D$tmax, ra_matrix(D$sites$lat, doy))
EW30 <- ewma_rows(D$ppt - PET, 30); EW60 <- ewma_rows(D$ppt - PET, 60)
G5 <- matrix(NA_real_, nrow(TM), ncol(TM))
for (Y in unique(yr)) { i <- which(yr == Y); G5[, i] <- t(apply(pmax(TM[, i, drop = FALSE] - 5, 0), 1, cumsum)) }
fz <- D$tmin <= -2.2
day_index <- function(Y, d) match(as.Date(sprintf("%d-01-01", Y)) + round(d) - 1, dates)

targets <- c("greenup", "curing", "cured")
spring <- c("gdd5_mar31", "gdd5_apr30", "gdd10_apr30", "t_marapr", "chill_days", "last_frz0", "last_frz2",
            "p_son_prev", "p_djf", "p_mam", "vpd_mam", "ewma30_may1", "ewma60_may1")
summer <- c("p_amj", "p_jja", "p_jas", "vpd_jja", "t_mayjun", "t_julaug", "t_sepoct", "events_amjjas",
            "dryspell_amjjas", "dryspell_jjas", "onset", "first_frz0", "first_frz2",
            paste0("ewma30_", c("jun1", "jul1", "aug1", "sep1")), paste0("ewma60_", c("jun1", "jul1", "aug1", "sep1")))
lead_vars <- c("ew30_cut", "ew60_cut", "p90_cut", "tx30_cut", "p_gu30", "p_gu60", "dry_gu90", "events_gu90")
sets <- list(greenup = list(all = spring),
             curing = list(lead = c(spring, lead_vars), hindcast = c(spring, lead_vars, summer)),
             cured = list(lead = c(spring, lead_vars), hindcast = c(spring, lead_vars, summer)))
l1_var <- c(greenup = "t_marapr", curing = "p_jja", cured = "p_jja")

pred_rows <- list(); sel_rows <- list(); imp_rows <- list()
for (tg in targets) {
  x <- S[!is.na(S[[tg]]), ]
  med <- c(tapply(x[[tg]], x$series, cmed))
  x$dep_raw <- cdiff(x[[tg]], med[x$series])
  x$y <- as.numeric(med[x$series]) + x$dep_raw                       # unwrapped around the series median
  x$mismatch <- abs(x$dep_raw) > 90
  ny <- tapply(!x$mismatch, x$series, sum); x <- x[x$series %in% names(ny)[ny >= min_years], ]
  # site-year drivers and the site climatology (mean over the site's PRISM years)
  m <- match(paste(x$site_id, x$year), paste(YM$site_id, YM$year))
  for (v in setdiff(unique(unlist(sets[[tg]])), lead_vars)) {
    clim <- tapply(YM[[v]], YM$site_id, mean, na.rm = TRUE)
    x[[paste0("d_", v)]] <- YM[[v]][m] - clim[as.character(x$site_id)]
  }
  # lead-time drivers at the cutoff (series typical date - lead): memory, rain and heat before it
  if (tg != "greenup") {
    cut_d <- med[x$series] - lead; srow <- match(x$site_id, D$sites$site_id)
    k <- day_index(x$year, cut_d)
    w <- function(M, n, f) vapply(seq_len(nrow(x)), function(i) if (is.na(k[i]) || k[i] <= n) NA_real_ else f(M[srow[i], (k[i] - n + 1):k[i]]), 0)
    x$ew30_cut <- EW30[cbind(srow, k)]; x$ew60_cut <- EW60[cbind(srow, k)]
    x$p90_cut <- w(D$ppt, 90, sum); x$tx30_cut <- w(D$tmax, 30, mean)
    # departures from the series' own mean of the same quantity (same calendar cutoff every year)
    for (v in lead_vars) x[[paste0("d_", v)]] <- x[[v]] - ave(x[[v]], x$series, FUN = function(z) mean(z, na.rm = TRUE))
  }
  x$aridity <- N$aridity[match(x$site_id, N$site_id)]; x$p_jas_share <- N$p_jas_share[match(x$site_id, N$site_id)]
  xs <- x[!x$mismatch, ]
  log_msg(sprintf("%s: %d series-years in %d series (%d mismatches dropped, %.1f %%)", tg, nrow(xs), length(unique(xs$series)),
                  sum(x$mismatch), 100 * mean(x$mismatch)))

  # L2 per-row quantities from the daily data
  srow <- match(xs$site_id, D$sites$site_id); k_obs <- day_index(xs$year, xs$y)
  if (tg == "greenup") {
    xs$q_obs <- G5[cbind(srow, k_obs)]
  } else {
    gu <- tapply(S$greenup, S$series, function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE))
    gu_s <- as.numeric(gu[xs$series]); gu_s[is.na(gu_s)] <- xs$gu_typical[is.na(gu_s)]
    gu_s[is.na(gu_s)] <- 100
    xs$start <- pmax(gu_s + 30, med[xs$series] - (ns$l2_window_days %||% 60))   # search typical +/- window
    xs$q_obs <- EW30[cbind(srow, k_obs)]
  }
  l2_date <- function(i, thr) {
    Y <- xs$year[i]; r <- srow[i]
    if (tg == "greenup") {
      j <- which(yr == Y); hit <- which(G5[r, j] >= thr)[1]
      return(if (is.na(hit)) NA_real_ else doy[j][hit])
    }
    a <- day_index(Y, xs$start[i]); if (is.na(a)) return(NA_real_)
    b <- min(day_index(Y, xs$base[i] + (ns$l2_window_days %||% 60)), length(dates), na.rm = TRUE); if (b <= a + 1) return(NA_real_); seg <- a:b
    e <- EW30[r, seg]; below <- e < thr
    # first downward crossing that then holds for 15 days; a season that starts below the
    # threshold and never rises above it keeps the typical date
    cr <- which(!below[-length(below)] & below[-1]) + 1
    cr <- cr[vapply(cr, function(j) all(below[j:min(j + 14, length(below))]), TRUE)]
    hit <- cr[1]
    fr <- which(fz[r, seg] & doy[seg] >= 200)[1]
    h <- suppressWarnings(min(hit, fr, na.rm = TRUE)); if (!is.finite(h)) return(NA_real_)
    as.numeric(dates[seg[h]] - as.Date(sprintf("%d-01-01", Y))) + 1
  }

  for (st in names(sets[[tg]])) {
    vars <- paste0("d_", sets[[tg]][[st]]); vars <- vars[vars %in% names(xs)]
    for (Y in sort(unique(xs$year))) {
      te <- which(xs$year == Y); tr <- which(xs$year != Y)
      base <- tapply(xs$y[tr], xs$series[tr], mean); nb <- tapply(xs$y[tr], xs$series[tr], length)
      te <- te[xs$series[te] %in% names(nb)[nb >= 2]]; if (!length(te)) next
      xs$base <- as.numeric(base[xs$series]); dep <- xs$y - xs$base
      out <- data.frame(target = tg, set = st, series = xs$series[te], site_id = xs$site_id[te], year = Y,
                        pathway = xs$pathway[te], region = xs$region[te], obs = xs$y[te], L0 = xs$base[te], L1 = NA_real_, L3 = NA_real_)
      for (pw in unique(xs$pathway[te])) {
        tp <- tr[xs$pathway[tr] %in% pw & !is.na(dep[tr])]; ep <- te[xs$pathway[te] %in% pw]; oe <- match(ep, te)
        X <- as.matrix(xs[, vars]); X[is.na(X)] <- 0
        # L1
        v1 <- paste0("d_", l1_var[[tg]]); if (!v1 %in% vars) v1 <- "d_ew30_cut"
        b1 <- if (length(tp) >= 20) sum(X[tp, v1] * dep[tp]) / sum(X[tp, v1]^2) else 0
        out$L1[oe] <- xs$base[ep] + pmax(pmin(b1 * X[ep, v1], 60), -60)
        # L3 forward selection
        sel <- character(); best <- 0; sst <- sum(dep[tp]^2)
        if (length(tp) >= 30) for (s in seq_len(max_pred)) {
          g <- vapply(setdiff(vars, sel), function(cn) 1 - sum(stats::lm.fit(X[tp, c(sel, cn), drop = FALSE], dep[tp])$residuals^2) / sst, 0)
          if (!length(g) || max(g) - best < 0.01) break
          sel <- c(sel, names(which.max(g))); best <- max(g)
        }
        if (length(sel)) {
          b <- stats::lm.fit(X[tp, sel, drop = FALSE], dep[tp])$coefficients; b[is.na(b)] <- 0
          out$L3[oe] <- xs$base[ep] + pmax(pmin(as.vector(X[ep, sel, drop = FALSE] %*% b), 60), -60)
          sel_rows[[length(sel_rows) + 1]] <- data.frame(target = tg, set = st, pathway = pw, year = Y, driver = sub("^d_", "", sel), b = b, r2_fit = best)
        } else out$L3[oe] <- xs$base[ep]
      }
      # L2 threshold rule: threshold = series training mean of the quantity at the observed date
      thr <- tapply(xs$q_obs[tr], xs$series[tr], function(z) mean(z, na.rm = TRUE))
      out$L2 <- vapply(te, function(i) l2_date(i, thr[[xs$series[i]]]), 0)
      out$L2 <- xs$base[te] + ifelse(is.na(out$L2), 0, pmax(pmin(cdiff(out$L2, xs$base[te]), 90), -90))
      # L4 random forest, pooled
      if (requireNamespace("ranger", quietly = TRUE) && length(tr) >= 60) {
        F <- data.frame(xs[, c(vars, "aridity", "p_jas_share", "lat")], pathway = factor(xs$pathway), life = factor(xs$life_cycle))
        F[vars][is.na(F[vars])] <- 0
        ok <- tr[!is.na(dep[tr]) & stats::complete.cases(F[tr, ])]
        rf <- ranger::ranger(x = F[ok, ], y = dep[ok], num.trees = 500, min.node.size = 10, importance = "permutation", seed = 52)
        Fe <- F[te, ]; Fe[is.na(Fe)] <- 0
        out$L4 <- xs$base[te] + stats::predict(rf, Fe)$predictions
        imp_rows[[length(imp_rows) + 1]] <- data.frame(target = tg, set = st, year = Y, driver = names(rf$variable.importance), imp = rf$variable.importance)
      } else out$L4 <- NA_real_
      pred_rows[[length(pred_rows) + 1]] <- out
    }
  }
}
P <- do.call(rbind, pred_rows)
utils::write.csv(P, file.path(odir, "predictions.csv"), row.names = FALSE)

# scores
lv <- c("L0", "L1", "L2", "L3", "L4")
score <- function(d, group, gval) {
  e0 <- abs(d$obs - d$L0); dep <- d$obs - d$L0
  do.call(rbind, lapply(lv, function(l) {
    e <- abs(d$obs - d[[l]]); ok <- !is.na(e)
    data.frame(target = d$target[1], set = d$set[1], group = group, value = gval, level = l, n = sum(ok),
               series = length(unique(d$series[ok])), mae = mean(e[ok]), mae_L0 = mean(e0[ok]),
               skill = 1 - mean(e[ok]) / mean(e0[ok]),
               r = if (l == "L0" || sum(ok) < 10) NA else stats::cor(d[[l]][ok] - d$L0[ok], dep[ok]))
  }))
}
SC <- list()
for (k in split(P, paste(P$target, P$set))) {
  SC[[length(SC) + 1]] <- score(k, "all", "all")
  for (g in unique(k$pathway)) { u <- k[k$pathway %in% g, ]; if (nrow(u) >= 30) SC[[length(SC) + 1]] <- score(u, "pathway", g) }
  for (g in setdiff(unique(stats::na.omit(k$region)), "none")) { u <- k[k$region %in% g, ]; if (nrow(u) >= 30) SC[[length(SC) + 1]] <- score(u, "region", g) }
}
SC <- do.call(rbind, SC)
utils::write.csv(SC, file.path(odir, "scores.csv"), row.names = FALSE)
if (length(sel_rows)) utils::write.csv(do.call(rbind, sel_rows), file.path(odir, "selected.csv"), row.names = FALSE)
if (length(imp_rows)) {
  I <- do.call(rbind, imp_rows); I <- stats::aggregate(imp ~ target + set + driver, I, mean)
  utils::write.csv(I[order(I$target, I$set, -I$imp), ], file.path(odir, "importance.csv"), row.names = FALSE)
}
for (i in which(SC$group == "all")) with(SC[i, ], log_msg(sprintf("%-8s %-8s %s: n %d | MAE %.1f d | skill %+.2f | r %s", target, set, level, n, mae, skill, ifelse(is.na(r), "-", sprintf("%.2f", r)))))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
