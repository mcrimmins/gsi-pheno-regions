#!/usr/bin/env Rscript
# Step 35 (evaluation, GSI curing forms): the GSI-per-region model test (34) showed that tuning
# the GSI ramps per region matches cell climatology for green-up but not for curing, and that
# the 28-day precipitation / KBDI moisture terms track year-to-year curing only weakly. This
# test changes the model's form, not the regions. Config `gsi_curing:` (anything not set there
# comes from `gsi_model:`); code R/gsi_model.R (+ R/gsi.R).
#   Forms = moisture option x crossing form, every combination:
#     moisture  kbdi (34's best), kbdi30 (30-day mean of KBDI), precip60, precip90 (60 / 90-day
#               precipitation sums): longer moisture memory
#     crossing  rel         green-up / curing at 20 % / 50 % of the season's amplitude (34)
#               rel_split   the two fractions fitted separately
#               abs_nfdrs4  absolute GSI levels as NFDRS4: green-up at GU (fitted), 50 % cured
#                           at (1 + GU) / 2 (NFDRS4 Cure())
#               abs_split   absolute green-up and curing levels fitted separately
#   Each form is fitted CONUS-wide and per region (A2 at k 7 and 13 by default) on odd years and
#   scored on even years, and the reverse, against the satellite green-up (sos20) and NDII7
#   curing (cure50) dates. "best" = per region and fold, the form with the lowest error in the
#   fitting years. Reference: each cell's median date in the fitting years (climatology).
# Needs the sample cache from 34 (eval/gsi_model/sample_<profile>.rds, with tmax): run 34 first.
# Missing dates (v2, Oct 10): where a fitted GSI gives no green-up or curing date, the cell's
#   typical date in the fitting years is used instead, in the fits and in the scores (config
#   `missing: fallback`; "penalty" = v1's fixed penalty_days). Scores also report coverage
#   (share of observed dates the GSI dates itself) and a paired comparison on the cell-years
#   the GSI does date: its error there and the typical date's error there (paired_*).
#   Each fit logs a line when it finishes ("fit i/n done"), so progress shows in the log.
# Outputs eval/gsi_curing/: scores.csv (model x form x partition x k), region_scores.csv (per
# region: best form, kbdi-rel reference, conus best, climatology), params.csv (fitted values).
# Parallel: forked workers on Linux (shared memory), sequential on Windows.
#
# Usage: Rscript scripts/35_gsi_curing.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/35_gsi_curing.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "gsi.R"))
source(here::here("R", "gsi_model.R"))

cfg <- load_cfg(profile)
log_file <- log_init("35_gsi_curing", profile)
notify_init(cfg)
gmc <- utils::modifyList(cfg$gsi_model, cfg$gsi_curing %||% list())
g <- cfg$features$gsi
out_dir <- file.path(cfg$run_dir, "eval", "gsi_curing"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
nw <- as.integer(Sys.getenv("GSI_WORKERS", gmc$workers %||% 1))
fork <- .Platform$OS.type == "unix" && nw > 1
pmap <- function(X, f) if (fork) parallel::mclapply(X, f, mc.cores = nw, mc.preschedule = FALSE) else lapply(X, f)
log_msg("=== 35_gsi_curing | profile: ", profile, " | ", if (fork) paste(nw, "forked workers") else "sequential")

# ---- sample and drivers (cache from 34) ---------------------------------------------------
cache <- file.path(cfg$run_dir, "eval", "gsi_model", sprintf("sample_%s.rds", profile))
if (!file.exists(cache)) stop("no ", cache, ": run scripts/34_gsi_model.R ", profile, " first")
S <- readRDS(cache)
if (is.null(S$drv$tmax) || is.null(S$cells$map_mm)) stop("sample cache lacks tmax / map_mm: rerun 34 (it extends the cache)")
years <- S$years; nY <- length(years); n <- nrow(S$cells)
S$drv$photo_s <- m_photoperiod(S$cells$lat, as.integer(format(S$dates, "%j")), g$photo_method)
for (w in c(28, 60, 90)) S$drv[[paste0("prec", w)]] <- m_roll_sum(S$drv$ppt, w)
S$drv$kbdi <- m_kbdi(S$drv$ppt, S$drv$tmax, S$cells$map_mm / 25.4, g$kbdi$interception_in)
S$drv$kbdi30 <- m_roll_mean(S$drv$kbdi, 30)
S$drv$ppt <- NULL; S$drv$tmax <- NULL; invisible(gc())
log_msg(sprintf("sample: %d cells x %d days (%s .. %s), years %d-%d", n, length(S$dates), min(S$dates), max(S$dates),
                min(years), max(years)))
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))

forms <- expand.grid(moisture = unlist(gmc$moisture), cross = unlist(gmc$cross), stringsAsFactors = FALSE)
forms$form <- paste(forms$moisture, forms$cross, sep = "-")
ref_form <- gmc$reference_form %||% "kbdi-rel"
log_msg(nrow(forms), " forms: ", paste(forms$form, collapse = ", "))

P <- eval_partitions(cfg, grid, unlist(gmc$k))
P <- P[intersect(unlist(gmc$partitions), names(P))]
if (!length(P)) stop("none of the partitions ", paste(unlist(gmc$partitions), collapse = ", "), " found")
units <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  reg <- P[[nm]][S$cells$row, kn]
  for (r in sort(unique(reg[!is.na(reg)]))) units[[length(units) + 1]] <- list(partition = nm, k = as.integer(sub("k", "", kn)), region = r, idx = which(reg == r))
}
# drivers for a set of cells; `need` = only those a form uses (keeps forked workers small:
# each fit copies its cells' series, ~0.11 GB per 1,000 cells per series on CONUS)
sub_drv <- function(idx, need = names(S$drv)) lapply(S$drv[need], function(m) m[idx, , drop = FALSE])
drv_need <- function(mo) c("tmin", "vpd_pa", "photo_s", if (!is.null(gm_moist[[mo]])) gm_moist[[mo]]$drv)
# observations for fitting, with each cell's typical date in the fitting years as the fallback
# where the GSI gives no date (v2, Oct 10: v1 counted a fixed 60 days, so fits in hard regions
# learned to skip curing dates)
fill_on <- !identical(gmc$missing, "penalty")
sub_obs <- function(idx, yi) {
  o <- list(sos = S$obs$sos[idx, yi, drop = FALSE], cure = S$obs$cure[idx, yi, drop = FALSE])
  if (fill_on) {
    rep_ <- function(m) matrix(apply(m, 1, stats::median, na.rm = TRUE), nrow(m), ncol(m))
    o$fill_sos <- rep_(o$sos); o$fill_cure <- rep_(o$cure)
  }
  o
}
fems <- g$param_sets$fems
start_for <- function(mo, cr) {
  p <- list(moisture = mo, cross = cr, tmin = unlist(fems$tmin), vpd = unlist(fems$vpd), photo_h = unlist(fems$photo_h))
  w <- c(precip = 28, precip60 = 60, precip90 = 90)
  if (mo %in% names(w)) p$moist <- c(fems$precip$lo, fems$precip$hi) * w[[mo]] / 28   # FEMS 28-day ramp, scaled
  if (mo %in% c("kbdi", "kbdi30")) p$moist <- c(g$kbdi$lo, g$kbdi$hi)
  p$x <- switch(cr, rel = NULL, rel_split = list(up = 0.2, down = 0.5), abs_nfdrs4 = list(gu = fems$greenup),
                abs_split = list(up = fems$greenup, down = (1 + fems$greenup) / 2))
  p
}

# ---- fits ---------------------------------------------------------------------------------
set.seed(gmc$seed + 1)
conus_idx <- sort(sample(n, min(gmc$conus_fit_cells, n)))
jobs <- list()
for (f in names(folds)) for (i in seq_len(nrow(forms))) {
  fm <- forms[i, ]
  base <- list(fold = f, form = fm$form, moisture = fm$moisture, cross = fm$cross)
  jobs[[length(jobs) + 1]] <- c(base, list(partition = "conus", k = 1L, region = 1L, idx = conus_idx))
  for (u in units) if (length(u$idx) >= gmc$min_cells) jobs[[length(jobs) + 1]] <- c(base, u)
}
# largest first, but no more than max_big CONUS-sized fits at once (memory): interleave them
# with the regional fits
jobs <- jobs[order(-vapply(jobs, function(j) length(j$idx), 0))]
big <- vapply(jobs, function(j) j$partition == "conus", TRUE); mb <- gmc$max_big %||% 8
if (sum(big) > mb && any(!big)) {
  sm <- jobs[!big]; bg <- jobs[big]; step <- max(1, floor(length(sm) / ceiling(length(bg) / mb)))
  out <- list(); while (length(bg) || length(sm)) {
    out <- c(out, bg[seq_len(min(mb, length(bg)))], sm[seq_len(min(step, length(sm)))])
    bg <- bg[-seq_len(min(mb, length(bg)))]; sm <- sm[-seq_len(min(step, length(sm)))]
  }
  jobs <- out
}
log_msg(length(units), " regions in ", length(P), " partition(s) x ", length(gmc$k), " k; ", length(jobs),
        " fits (maxit ", gmc$maxit, ", restarts ", gmc$restarts %||% 0, "); at most ", mb, " CONUS fits at once")
notify("35_gsi_curing started", sprintf("%d fits, %d workers", length(jobs), nw), priority = 2, tags = "hourglass")
nj <- length(jobs)
fits <- pmap(seq_along(jobs), function(ji) {
  j <- jobs[[ji]]; t1 <- Sys.time()
  tr <- folds[[j$fold]]
  ft <- tryCatch(gm_fit(sub_drv(j$idx, drv_need(j$moisture)), sub_obs(j$idx, tr), S$dates, years[tr], g, gmc, start_for(j$moisture, j$cross)),
                 error = function(e) list(error = conditionMessage(e)))
  log_msg(sprintf("fit %d/%d done: %s %s %s k%d r%d, %d cells, %s evaluations, %.1f min", ji, nj, j$fold, j$form,
                  j$partition, as.integer(j$k), as.integer(j$region), length(j$idx), ft$evals %||% "?",
                  as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  c(j[c("fold", "form", "moisture", "cross", "partition", "k", "region")], list(fit = ft, n = length(j$idx)))
})
fits <- lapply(fits, function(x) if (inherits(x, "try-error")) list(fit = list(error = as.character(x)), fold = "?",
                                    form = "?", partition = "?", k = NA, region = NA) else x)
for (b in Filter(function(x) !is.null(x$fit$error), fits))
  log_err("fit ", b$fold, " ", b$form, " ", b$partition, " k", b$k, " r", b$region, ": ", b$fit$error)
fits <- Filter(function(x) is.null(x$fit$error), fits)
key <- function(f, fm, p, k, r) paste(f, fm, p, k, r)
fit_of <- stats::setNames(lapply(fits, `[[`, "fit"), vapply(fits, function(x) key(x$fold, x$form, x$partition, x$k, x$region), ""))
params <- do.call(rbind, lapply(fits, function(x) data.frame(fold = x$fold, form = x$form, moisture = x$moisture, cross = x$cross,
  partition = x$partition, k = x$k, region = x$region, cells = x$n, loss = x$fit$loss, evals = x$fit$evals,
  t(gm_par_vec(x$fit$par)))))
utils::write.csv(params, file.path(out_dir, "params.csv"), row.names = FALSE)
log_msg(sprintf("fits done: %d (%.1f min); evaluations per fit: median %d, max %d; hit the limit: %.0f %%", nrow(params),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), as.integer(stats::median(params$evals)), max(params$evals),
                100 * mean(params$evals >= gmc$maxit)))
best_fit <- function(f, p, k, r) {
  cand <- Filter(Negate(is.null), lapply(forms$form, function(fm) fit_of[[key(f, fm, p, k, r)]]))
  if (!length(cand)) return(NULL)
  cand[[which.min(vapply(cand, `[[`, 0, "loss"))]]
}

# ---- predictions on held-out years and scores ---------------------------------------------
empty <- function() list(sos = matrix(NA_real_, n, nY), eos = matrix(NA_real_, n, nY))
pred_par <- function(par, idx) gm_predict(sub_drv(idx, drv_need(par$moisture)), par, gmc$smooth, S$dates, years, g)
cross_pred <- function(get_groups) {
  out <- empty()
  for (f in names(folds)) {
    te <- setdiff(seq_len(nY), folds[[f]])
    for (gr in get_groups(f)) {
      if (is.null(gr$par)) next
      pr <- pred_par(gr$par, gr$idx)
      out$sos[gr$idx, te] <- pr$sos[, te]; out$eos[gr$idx, te] <- pr$eos[, te]
    }
  }
  out
}
clim <- empty()
for (f in names(folds)) {
  tr <- folds[[f]]; te <- setdiff(seq_len(nY), tr)
  clim$sos[, te] <- apply(S$obs$sos[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
  clim$eos[, te] <- apply(S$obs$cure[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
}
obs_all <- list(sos = S$obs$sos, cure = S$obs$cure)
sc <- function(m, anom = TRUE) gm_score(m, obs_all, gmc$penalty_days, gmc$cap_days, anom, if (fill_on && anom) clim)
all_forms <- c(forms$form, "best")
getter <- function(fm, p, k) function(f, r) {
  ft <- if (fm == "best") best_fit(f, p, k, r) else fit_of[[key(f, fm, p, k, r)]]
  ft %||% (if (fm == "best") best_fit(f, "conus", 1, 1) else fit_of[[key(f, fm, "conus", 1, 1)]])   # small regions: CONUS fit
}
conus_m <- lapply(stats::setNames(all_forms, all_forms), function(fm) {
  gf <- getter(fm, "conus", 1)
  cross_pred(function(f) list(list(idx = seq_len(n), par = gf(f, 1)$par)))
})
reg_m <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]
  for (fm in all_forms) {
    gf <- getter(fm, nm, k)
    reg_m[[paste(nm, k, fm, sep = "|")]] <- cross_pred(function(f)
      lapply(sort(unique(reg[!is.na(reg)])), function(r) list(idx = which(reg == r), par = gf(f, r)$par)))
  }
}
form_cols <- function(fm) if (fm == "best") c(moisture = NA, cross = NA) else unlist(forms[forms$form == fm, c("moisture", "cross")])
scores <- do.call(rbind, c(
  list(data.frame(model = "climatology", form = NA, moisture = NA, cross = NA, partition = NA, k = NA, t(sc(clim, FALSE)))),
  lapply(all_forms, function(fm) data.frame(model = "conus", form = fm, t(form_cols(fm)), partition = NA, k = NA, t(sc(conus_m[[fm]])))),
  lapply(names(reg_m), function(nm) { p <- strsplit(nm, "|", fixed = TRUE)[[1]]
    data.frame(model = "regional", form = p[3], t(form_cols(p[3])), partition = p[1], k = as.integer(p[2]), t(sc(reg_m[[nm]]))) })))
utils::write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)
# which form wins per region and fold (fitting-year loss)
win <- do.call(rbind, lapply(names(folds), function(f) do.call(rbind, lapply(units, function(u) {
  b <- params[params$fold == f & params$partition == u$partition & params$k == u$k & params$region == u$region, ]
  if (!nrow(b)) return(NULL)
  data.frame(fold = f, partition = u$partition, k = u$k, region = u$region, form = b$form[which.min(b$loss)])
}))))
if (!is.null(win)) {
  utils::write.csv(win, file.path(out_dir, "best_form.csv"), row.names = FALSE)
  tw <- sort(table(win$form), decreasing = TRUE)
  log_msg("best form per region and fold: ", paste(sprintf("%s %d", names(tw), tw), collapse = ", "))
}

# per region: best form, reference form, CONUS best, climatology
rs <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]
  for (r in sort(unique(reg[!is.na(reg)]))) {
    i <- which(reg == r); ob <- list(sos = obs_all$sos[i, , drop = FALSE], cure = obs_all$cure[i, , drop = FALSE])
    subm <- function(m) list(sos = m$sos[i, , drop = FALSE], eos = m$eos[i, , drop = FALSE])
    mods <- list(regional_best = reg_m[[paste(nm, k, "best", sep = "|")]], regional_ref = reg_m[[paste(nm, k, ref_form, sep = "|")]],
                 conus_best = conus_m$best, climatology = clim)
    for (mn in names(mods)) if (!is.null(mods[[mn]]))
      rs[[length(rs) + 1]] <- data.frame(partition = nm, k = k, region = r, cells = length(i), fitted = length(i) >= gmc$min_cells,
                                         model = mn, t(gm_score(subm(mods[[mn]]), ob, gmc$penalty_days, gmc$cap_days, mn != "climatology",
                                                    if (fill_on && mn != "climatology") subm(clim))))
  }
}
utils::write.csv(do.call(rbind, rs), file.path(out_dir, "region_scores.csv"), row.names = FALSE)

for (i in seq_len(nrow(scores))) with(scores[i, ], log_msg(sprintf(
  "  %-11s %-20s %-10s %3s green-up MAE %5.1f d, r %5.2f | curing MAE %5.1f d, median %5.1f, <=30 d %3.0f%%, coverage %3.0f%%, r %5.2f",
  model, ifelse(is.na(form), "", form), ifelse(is.na(partition), "", partition), ifelse(is.na(k), "", k),
  sos_mae, sos_anom_r, cure_mae, cure_mdae, 100 * cure_within30, 100 * cure_coverage, cure_anom_r)))
msg <- sprintf("%d fits, %d cells, %.1f min", nrow(params), n, as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | outputs: ", out_dir, " | log: ", log_file)
notify("35_gsi_curing finished", msg, priority = 3, tags = "white_check_mark")
