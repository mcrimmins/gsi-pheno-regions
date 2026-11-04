#!/usr/bin/env Rscript
# Step 11 (NPN x fire study; P720, minutes to ~1 h; xgboost): roadmap phase 3, machine learning. Does adding NPN
# status / intensity features improve out-of-sample prediction of weekly fire, where, and which phenophases?
# Fixed before results (Mike, 3 Nov 2026; roadmap doc claude/gsi-npn-fire-roadmap.md; config npn_fire: ml).
#
# Rows: Level II ecoregion x fuel class (timber TU+TL, shrub SH, herb GR+GS) x year x week, 2009-2024.
# Target: weekly fires in that fuel class, one model pair per fire record (config roadmap: outcomes).
# Weather-only features: week (and its sine / cosine), ecoregion, fuel class, ppet60 and its anomaly, 8-week
#   trailing tmean and its anomaly (fuel-class weight), and the training-years mean fire count for the
#   ecoregion x fuel x week (log). Anomalies are relative to training-year (2009-2017) weekly means.
# NPN features (added to the above): for each main plant group x pheno class, the weekly share-yes anomaly
#   (>= min_sites sites) and, for the percent classes, the mean-percent anomaly (>= min_sites percent sites);
#   missing where too few sites (xgboost handles missing values).
# Models: xgboost, Poisson objective, the same settings for both; rounds chosen by early stopping on the
#   validation years. Train 2009-2017, validate 2018-2020. Skill = 1 - deviance(weather + NPN) /
#   deviance(weather only) on the validation years, overall and by ecoregion and fuel class; 90 % interval
#   from a bootstrap over ecoregion-years. Importance: grouped permutation on the validation years (each plant
#   group x pheno class block shuffled across rows): increase in deviance.
# HOLD-OUT: 2021-2024 are never scored unless the script is run with `final` (roadmap phase 4, once): then the
#   models are refitted on 2009-2020 with the validation-chosen rounds and scored on 2021-2024.
# Outputs <run_dir>/npn_fire/ml/: skill.csv (outcome x scope), skill_eco.csv, importance.csv, rows.csv (counts);
#   final_*.csv in final mode.
# Usage: Rscript npn_fire/scripts/11_nf_ml.R [profile] [final]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && args[1] != "final") args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
final <- "final" %in% args
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
stub <- identical(Sys.getenv("NF_XGB_STUB"), "1")    # testing only: a Poisson GLM stands in for xgboost
if (!stub && !requireNamespace("xgboost", quietly = TRUE)) stop("needs xgboost: Rscript -e \"renv::install('xgboost')\"")

cfg <- load_nf_cfg(profile)
log_file <- log_init("11_nf_ml", profile)
rc <- cfg$npn_fire$roadmap; mc <- cfg$npn_fire$ml
outs <- unlist(rc$outcomes %||% list("fod_all"))
tr_y <- seq(mc$train_years[[1]] %||% 2009, mc$train_years[[2]] %||% 2017); va_y <- seq(mc$valid_years[[1]] %||% 2018, mc$valid_years[[2]] %||% 2020)
ho <- unlist(rc$holdout_years %||% list(2021, 2024)); te_y <- seq(ho[1], ho[2])
min_sites <- as.integer(mc$min_sites %||% 5); tw_w <- as.integer(mc$temp_weeks %||% 8); B <- as.integer(mc$bootstrap_B %||% 500)
groups <- unlist(mc$groups %||% list("deciduous_tree", "evergreen_tree", "conifer", "deciduous_shrub", "evergreen_shrub", "grass", "forb"))
pars <- mc$params %||% list(eta = 0.05, max_depth = 5, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 10)
max_rounds <- as.integer(mc$max_rounds %||% 2000); es <- as.integer(mc$early_stopping %||% 50); nthread <- as.integer(mc$nthread %||% 8)
out <- nf_out(cfg, "ml")
set.seed(mc$seed %||% 20261104)
t0 <- Sys.time()
log_msg("=== 11_nf_ml | profile: ", profile, " | train ", min(tr_y), "-", max(tr_y), " | validate ", min(va_y), "-", max(va_y),
        if (final) paste0(" | FINAL: scoring the hold-out ", min(te_y), "-", max(te_y)) else " | hold-out not scored",
        if (stub) " | STUB (GLM, testing only)" else "")
yrs <- if (final) c(tr_y, va_y, te_y) else c(tr_y, va_y)

FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
FW <- FW[FW$year %in% yrs, ]
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
if (!"tmean" %in% names(PW)) stop("prism/eco_weekly.csv has no tmean")
EW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "weekly", "eco_weekly"))); EW$eco_l2_code <- as.character(EW$eco_l2_code)
EW <- EW[EW$year %in% yrs & EW$group %in% groups, ]
fuel_cls <- list(timber = c("TU", "TL"), shrub = "SH", herb = c("GR", "GS"))
ecos <- sort(intersect(unique(FW$eco_l2_code), unique(PW$eco_l2_code)))

# --- base grid and weather features ---
G <- expand.grid(week = 1:52, year = yrs, fuel = names(fuel_cls), eco_l2_code = ecos, stringsAsFactors = FALSE)
G$key_w <- paste(G$eco_l2_code, G$fuel, G$year, G$week)
PWk <- paste(PW$eco_l2_code, PW$weight, PW$year, PW$week)
G$ppet60 <- PW$ppet60[match(paste(G$eco_l2_code, G$fuel, G$year, G$week), PWk)]
tm_at <- function(e, f, y, w) PW$tmean[match(paste(e, f, y, w), PWk)]
trail <- rowMeans(sapply(0:(tw_w - 1), function(k) { w <- G$week - k; y <- G$year - (w < 1); w[w < 1] <- w[w < 1] + 52; tm_at(G$eco_l2_code, G$fuel, y, w) }), na.rm = TRUE)
G$tmean8 <- trail
anom <- function(v, grp_key) { tr <- G$year %in% tr_y; m <- tapply(v[tr], grp_key[tr], mean, na.rm = TRUE); v - as.numeric(m[grp_key]) }
cw <- paste(G$eco_l2_code, G$fuel, G$week)
G$ppet60_anom <- anom(G$ppet60, cw); G$tmean8_anom <- anom(G$tmean8, cw)
G$wk_sin <- sin(2 * pi * G$week / 52); G$wk_cos <- cos(2 * pi * G$week / 52)
G$eco_id <- match(G$eco_l2_code, ecos); G$fuel_id <- match(G$fuel, names(fuel_cls))
for (e in seq_along(ecos)) G[[paste0("eco_", gsub("\\.", "_", ecos[e]))]] <- as.integer(G$eco_id == e)

# --- NPN features (anomalies vs training-year weekly means) ---
pct_cls <- c(2, 3, 4, 7, 12)
EW$k_ew <- paste(EW$eco_l2_code, EW$group, EW$pheno_class_id, EW$week)
npn_cols <- character(0); blocks <- list()
ekey <- paste(G$eco_l2_code, G$year, G$week)
for (g in groups) for (cl in 1:13) {
  d <- EW[EW$group == g & EW$pheno_class_id == cl, ]
  if (!nrow(d)) next
  cols <- character(0)
  for (var in c("share_yes", if (cl %in% pct_cls) "mean_pct")) {
    n <- if (var == "share_yes") d$sites else d$sites_pct
    v <- ifelse(n >= min_sites, d[[var]], NA_real_)
    trm <- tapply(v[d$year %in% tr_y], d$k_ew[d$year %in% tr_y], mean, na.rm = TRUE)
    a <- v - as.numeric(trm[d$k_ew])
    if (sum(is.finite(a)) < 200) next
    nm <- sprintf("npn_%s_c%02d_%s", g, cl, if (var == "share_yes") "share" else "pct")
    G[[nm]] <- a[match(ekey, paste(d$eco_l2_code, d$year, d$week))]
    cols <- c(cols, nm)
  }
  if (length(cols)) { npn_cols <- c(npn_cols, cols); blocks[[sprintf("%s_c%02d", g, cl)]] <- cols }
}
wx_cols <- c("week", "wk_sin", "wk_cos", "eco_id", "fuel_id", grep("^eco_[0-9]", names(G), value = TRUE), "ppet60", "ppet60_anom", "tmean8", "tmean8_anom", "clim_logfire")
log_msg(sprintf("grid rows: %d | NPN features: %d in %d blocks (plant group x pheno class) | ecoregions: %d", nrow(G), length(npn_cols), length(blocks), length(ecos)))

pdev <- function(y, mu) { mu <- pmax(mu, 1e-9); 2 * sum(ifelse(y > 0, y * log(y / mu), 0) - (y - mu)) }
fit_pred <- function(Xtr, ytr, Xva, yva, Xte = NULL, rounds = NULL) {
  if (stub) { D <- as.data.frame(Xtr); D$y <- ytr; keep <- vapply(D, function(z) all(is.finite(z)), TRUE); D <- D[, keep, drop = FALSE]
    m <- suppressWarnings(stats::glm(y ~ ., data = D, family = stats::poisson()))
    pr <- function(X) { X <- as.data.frame(X)[, setdiff(names(D), "y"), drop = FALSE]; X[!is.finite(as.matrix(X))] <- 0; stats::predict(m, X, type = "response") }
    return(list(rounds = 1, va = pr(Xva), te = if (is.null(Xte)) NULL else pr(Xte), predict = pr)) }
  p <- c(pars, list(objective = "count:poisson", nthread = nthread))
  dtr <- xgboost::xgb.DMatrix(Xtr, label = ytr, missing = NA); dva <- xgboost::xgb.DMatrix(Xva, label = yva, missing = NA)
  if (is.null(rounds)) {
    a <- list(params = p, data = dtr, nrounds = max_rounds, early_stopping_rounds = es, verbose = 0)
    if ("evals" %in% names(formals(xgboost::xgb.train))) a$evals <- list(val = dva) else a$watchlist <- list(val = dva)
    m0 <- do.call(xgboost::xgb.train, a)
    el <- m0$evaluation_log %||% attr(m0, "evaluation_log"); rounds <- if (!is.null(el)) which.min(el[[2]]) else max_rounds
  }
  m <- xgboost::xgb.train(params = p, data = dtr, nrounds = rounds, verbose = 0)
  pr <- function(X) stats::predict(m, xgboost::xgb.DMatrix(X, missing = NA))
  list(rounds = rounds, va = pr(Xva), te = if (is.null(Xte)) NULL else pr(Xte), predict = pr)
}

SK <- list(); SE <- list(); IMP <- list(); RW <- list(); FIN <- list()
for (o in outs) {
  fo <- FW[FW$outcome == o, ]
  if (!nrow(fo)) { log_msg("  ", o, ": no fire rows"); next }
  fo$fuel <- names(fuel_cls)[vapply(fo$fuel, function(f) { i <- which(vapply(fuel_cls, function(z) f %in% z, TRUE)); if (length(i)) i[1] else NA_integer_ }, 0L)]
  fo <- fo[!is.na(fo$fuel), ]
  s <- tapply(fo$value, paste(fo$eco_l2_code, fo$fuel, fo$year, fo$week), sum)
  D <- G; D$y <- round(as.numeric(s[D$key_w])); D$y[is.na(D$y)] <- 0
  D <- D[is.finite(D$ppet60) & is.finite(D$tmean8), ]
  tr <- D$year %in% tr_y
  cm <- tapply(D$y[tr], paste(D$eco_l2_code, D$fuel, D$week)[tr], mean)
  D$clim_logfire <- log1p(as.numeric(cm[paste(D$eco_l2_code, D$fuel, D$week)]))
  # drop ecoregion x fuel series with almost no fire in the training years
  tot <- tapply(D$y[tr], paste(D$eco_l2_code, D$fuel)[tr], sum); keep_s <- names(tot)[tot >= (mc$min_train_fires %||% 50)]
  D <- D[paste(D$eco_l2_code, D$fuel) %in% keep_s, ]
  tr <- D$year %in% tr_y; va <- D$year %in% va_y; te <- D$year %in% te_y
  RW[[o]] <- data.frame(outcome = o, rows_train = sum(tr), rows_valid = sum(va), series = length(keep_s), fires_train = sum(D$y[tr]), fires_valid = sum(D$y[va]))
  Xw <- as.matrix(D[, wx_cols]); Xn <- as.matrix(D[, c(wx_cols, npn_cols)])
  fw <- fit_pred(Xw[tr, ], D$y[tr], Xw[va, ], D$y[va])
  fn <- fit_pred(Xn[tr, ], D$y[tr], Xn[va, ], D$y[va])
  dv_w <- pdev(D$y[va], fw$va); dv_n <- pdev(D$y[va], fn$va); sk <- 1 - dv_n / dv_w
  # bootstrap over ecoregion-years in the validation years
  bk <- paste(D$eco_l2_code[va], D$year[va]); ub <- unique(bk); idx <- split(seq_along(bk), bk)
  dw_b <- vapply(idx, function(i) pdev(D$y[va][i], fw$va[i]), 0); dn_b <- vapply(idx, function(i) pdev(D$y[va][i], fn$va[i]), 0)
  bs <- vapply(seq_len(B), function(i) { s <- sample(length(ub), replace = TRUE); 1 - sum(dn_b[s]) / sum(dw_b[s]) }, 0)
  ci <- stats::quantile(bs, c(0.05, 0.95), na.rm = TRUE)
  SK[[o]] <- data.frame(outcome = o, rounds_wx = fw$rounds, rounds_npn = fn$rounds, dev_wx = round(dv_w, 1), dev_npn = round(dv_n, 1),
    skill = round(sk, 4), skill_lo = round(ci[[1]], 4), skill_hi = round(ci[[2]], 4))
  log_msg(sprintf("  %-18s validation skill of NPN over weather-only %+.4f [%+.4f, %+.4f] | rounds %d / %d | %d series", o, sk, ci[[1]], ci[[2]],
                  fw$rounds, fn$rounds, length(keep_s)))
  # by ecoregion x fuel
  sg <- paste(D$eco_l2_code[va], D$fuel[va])
  SE[[o]] <- do.call(rbind, lapply(split(seq_along(sg), sg), function(i) { a <- D[va, ][i[1], ]
    data.frame(outcome = o, eco_l2_code = a$eco_l2_code, fuel = a$fuel, fires_valid = sum(D$y[va][i]),
      npn_features_present = round(mean(rowSums(is.finite(Xn[va, npn_cols, drop = FALSE][i, , drop = FALSE])) > 0), 2),
      skill = round(1 - pdev(D$y[va][i], fn$va[i]) / pdev(D$y[va][i], fw$va[i]), 4)) }))
  # grouped permutation importance (NPN blocks) on the validation rows
  Xva <- Xn[va, ]
  IMP[[o]] <- do.call(rbind, lapply(names(blocks), function(b) { cols <- blocks[[b]]; Xp <- Xva; pi <- sample(nrow(Xp)); Xp[, cols] <- Xp[pi, cols]
    data.frame(outcome = o, block = b, features = length(cols), share_rows_with_data = round(mean(rowSums(is.finite(Xva[, cols, drop = FALSE])) > 0), 3),
               dev_increase = round(pdev(D$y[va], fn$predict(Xp)) - dv_n, 2)) }))
  if (final) {
    trf <- tr | va
    fw2 <- fit_pred(Xw[trf, ], D$y[trf], Xw[te, ], D$y[te], rounds = fw$rounds); fn2 <- fit_pred(Xn[trf, ], D$y[trf], Xn[te, ], D$y[te], rounds = fn$rounds)
    FIN[[o]] <- data.frame(outcome = o, dev_wx = pdev(D$y[te], fw2$va), dev_npn = pdev(D$y[te], fn2$va), skill = round(1 - pdev(D$y[te], fn2$va) / pdev(D$y[te], fw2$va), 4))
    log_msg(sprintf("  FINAL %-18s hold-out skill %+.4f", o, FIN[[o]]$skill))
  }
}
if (!length(SK)) stop("no outcome could be modelled")
utils::write.csv(do.call(rbind, SK), file.path(out, "skill.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, SE), file.path(out, "skill_eco.csv"), row.names = FALSE)
IM <- do.call(rbind, IMP); IM <- IM[order(IM$outcome, -IM$dev_increase), ]; utils::write.csv(IM, file.path(out, "importance.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, RW), file.path(out, "rows.csv"), row.names = FALSE)
if (final) utils::write.csv(do.call(rbind, FIN), file.path(out, "final_skill.csv"), row.names = FALSE)
for (o in names(IMP)) { x <- IM[IM$outcome == o, ][1:min(8, sum(IM$outcome == o)), ]
  log_msg("  top NPN blocks (", o, "): ", paste(sprintf("%s %+.1f", x$block, x$dev_increase), collapse = "; ")) }
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
