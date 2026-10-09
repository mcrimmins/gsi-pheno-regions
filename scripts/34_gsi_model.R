#!/usr/bin/env Rscript
# Step 34 (evaluation, GSI per region): does tuning the GSI per region predict herbaceous
# green-up and curing better than the national defaults or one CONUS-wide tuning, and how
# many regions does that take? Config `gsi_model:`; code R/gsi_model.R (+ R/gsi.R).
#   1. Sample strict-herbaceous cells with satellite green-up in >= min_sat_years years;
#      read their daily PRISM series (tmin, vpdmax, ppt; 2001-2025) and their yearly
#      satellite sos20 (green-up) and cure50 (NDII7 curing). Cached in eval/gsi_model/.
#   2. Models, each scored on held-out years (two folds: fit on odd years, score on even, and
#      the reverse):
#        climatology    each cell's median satellite date in the fitting years (no weather)
#        default_fems   FEMS operational GSI (no precipitation ramp, 28-day smoothing)
#        default_nfdrs  NFDRS literature GSI with the 28-day precipitation ramp (21-day smoothing)
#        conus          one set of thresholds fitted to all regions
#        regional       thresholds fitted per region, for each partition x k in the config
#      Fitted thresholds: Tmin ramp, VPD ramp, start of the 1-hour daylength ramp and, for each
#      moisture option in the config (28-day precipitation sum, KBDI, none), that option's
#      ramp (FEMS smoothing). Each fit is a Nelder-Mead search with restarts. "best" picks,
#      per region and fold, the moisture option with the lowest error in the fitting years. GSI dates = the smoothed index's own sos20 / eos50
#      (relative crossings around its yearly peak, as for the satellite curves).
#   3. Scores: mean absolute error (days, each error capped at cap_days; a missing GSI date
#      where the satellite has one counts penalty_days), median absolute error, share within
#      30 days, coverage, and the correlation of yearly anomalies (how well each model tracks
#      early and late years within a cell; climatology has none). The fits minimize the capped
#      mean error.
# Outputs eval/gsi_model/: scores.csv (model x partition x k x moisture), region_scores.csv (per
# region), params.csv (fitted thresholds), sample_<profile>.rds (cache; an older cache without
# tmax is extended, not rebuilt).
# Parallel: forked workers on Linux (they share the sampled series), sequential on Windows.
#
# Usage: Rscript scripts/34_gsi_model.R [profile]
#   RStudio: gsi_profile <- "dev"; source("scripts/34_gsi_model.R")

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
log_file <- log_init("34_gsi_model", profile)
notify_init(cfg)
gmc <- cfg$gsi_model; g <- cfg$features$gsi
out_dir <- file.path(cfg$run_dir, "eval", "gsi_model"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
nw <- as.integer(Sys.getenv("GSI_WORKERS", gmc$workers %||% 1))
fork <- .Platform$OS.type == "unix" && nw > 1
pmap <- function(X, f) if (fork) parallel::mclapply(X, f, mc.cores = nw, mc.preschedule = FALSE) else lapply(X, f)
log_msg("=== 34_gsi_model | profile: ", profile, " | ", if (fork) paste(nw, "forked workers") else "sequential")

# ---- 1. sample, drivers, satellite targets (cached) --------------------------------------
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
pfile <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
sat_years <- cfg$years[file.exists(vapply(cfg$years, mfile, ""))]
years <- intersect(feature_years(cfg), sat_years)          # need Y-1 for the GSI spin-up
if (length(years) < 2) stop("need at least 2 years with satellite metrics and a previous PRISM year")
read_series <- function(v, cells, py) {
  parts <- pmap(py, function(y) {
    r <- terra::rast(pfile(v, y)); x <- as.matrix(terra::extract(r, cells)); colnames(x) <- names(r); x })
  if (any(vapply(parts, inherits, TRUE, "try-error"))) stop("reading ", v, " failed: ", parts[vapply(parts, inherits, TRUE, "try-error")][[1]])
  do.call(cbind, parts)
}
cache <- file.path(out_dir, sprintf("sample_%s.rds", profile))
if (file.exists(cache)) {
  S <- readRDS(cache); log_msg("sample from cache: ", nrow(S$cells), " cells")
} else {
  hm <- land_values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), grid)[, 1]
  sos_all <- vapply(years, function(Y) terra::values(terra::rast(mfile(Y))[["sos20"]], mat = FALSE)[land], numeric(length(land)))
  ok <- which(!is.na(hm) & hm > 0 & rowSums(!is.na(sos_all)) >= gmc$min_sat_years)
  set.seed(gmc$seed)
  pick <- sort(sample(ok, min(gmc$n_sample, length(ok))))
  cells <- data.frame(row = pick, cell = land[pick], lat = terra::yFromCell(grid, land[pick]))
  log_msg("sampled ", nrow(cells), " of ", length(ok), " eligible herbaceous cells")
  obs <- list(sos = sos_all[pick, , drop = FALSE],
              cure = vapply(years, function(Y) terra::values(terra::rast(mfile(Y))[["cure50"]], mat = FALSE)[land][pick], numeric(length(pick))))
  py <- seq(min(years) - 1, max(years))
  drv <- list(tmin = read_series("tmin", cells$cell, py), vpd_pa = read_series("vpdmax", cells$cell, py) * 100,
              ppt = read_series("ppt", cells$cell, py))
  dates <- as.Date(colnames(drv$tmin))
  if (any(diff(dates) != 1)) stop("daily series not contiguous")
  S <- list(cells = cells, obs = obs, drv = drv, dates = dates, years = years)
}
if (is.null(S$drv$tmax) || is.null(S$cells$map_mm)) {            # KBDI needs tmax and mean annual precipitation
  py <- seq(min(S$years) - 1, max(S$years))
  S$drv$tmax <- read_series("tmax", S$cells$cell, py)
  S$cells$map_mm <- terra::extract(terra::rast(file.path(static_dir(cfg), "kbdi_map.tif")), S$cells$cell)[, 1]
  saveRDS(S, cache)
  log_msg(sprintf("drivers: %d cells x %d days (%s .. %s), with tmax, cached", nrow(S$cells), length(S$dates), min(S$dates), max(S$dates)))
}
if (!file.exists(cache)) saveRDS(S, cache)
years <- S$years; nY <- length(years)
S$drv$photo_s <- m_photoperiod(S$cells$lat, as.integer(format(S$dates, "%j")), g$photo_method)
S$drv$prec28 <- m_roll_sum(S$drv$ppt, 28)
S$drv$kbdi <- m_kbdi(S$drv$ppt, S$drv$tmax, S$cells$map_mm / 25.4, g$kbdi$interception_in)
S$drv$ppt <- NULL; S$drv$tmax <- NULL; invisible(gc())
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))
moists <- unlist(gmc$moisture %||% "precip")

# regions of the sampled cells for each partition x k
P <- eval_partitions(cfg, grid, unlist(gmc$k))
P <- P[intersect(unlist(gmc$partitions), names(P))]
units <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  reg <- P[[nm]][S$cells$row, kn]
  for (r in sort(unique(reg[!is.na(reg)]))) units[[length(units) + 1]] <- list(partition = nm, k = as.integer(sub("k", "", kn)), region = r, idx = which(reg == r))
}
log_msg(length(units), " regions in ", length(P), " partitions x ", length(gmc$k), " k values; moisture: ",
        paste(moists, collapse = ", "), "; folds: ", paste(sprintf("%s (%d years)", names(folds), lengths(folds)), collapse = ", "))

sub_drv <- function(idx) lapply(S$drv, function(m) m[idx, , drop = FALSE])
sub_obs <- function(idx, yi) list(sos = S$obs$sos[idx, yi, drop = FALSE], cure = S$obs$cure[idx, yi, drop = FALSE])
fems <- g$param_sets$fems; nf <- g$param_sets$nfdrs
start_for <- function(mo) {
  p <- list(moisture = mo, tmin = unlist(fems$tmin), vpd = unlist(fems$vpd), photo_h = unlist(fems$photo_h))
  if (mo == "precip") p$moist <- c(fems$precip$lo, fems$precip$hi)
  if (mo == "kbdi") p$moist <- c(g$kbdi$lo, g$kbdi$hi)
  p
}

# ---- 2. fits ------------------------------------------------------------------------------
set.seed(gmc$seed + 1)
conus_idx <- sort(sample(nrow(S$cells), min(gmc$conus_fit_cells, nrow(S$cells))))
jobs <- list()
for (f in names(folds)) for (mo in moists) {
  jobs[[length(jobs) + 1]] <- list(fold = f, moisture = mo, partition = "conus", k = 1L, region = 1L, idx = conus_idx)
  for (u in units) if (length(u$idx) >= gmc$min_cells) jobs[[length(jobs) + 1]] <- c(list(fold = f, moisture = mo), u)
}
# largest jobs first so the long CONUS fits do not finish last
jobs <- jobs[order(-vapply(jobs, function(j) length(j$idx), 0))]
log_msg(length(jobs), " fits (maxit ", gmc$maxit, ", restarts ", gmc$restarts %||% 0, ")")
fits <- pmap(jobs, function(j) {
  tr <- folds[[j$fold]]
  ft <- tryCatch(gm_fit(sub_drv(j$idx), sub_obs(j$idx, tr), S$dates, years[tr], g, gmc, start_for(j$moisture)),
                 error = function(e) list(error = conditionMessage(e)))
  c(j[c("fold", "moisture", "partition", "k", "region")], list(fit = ft, n = length(j$idx)))
})
fits <- lapply(fits, function(x) if (inherits(x, "try-error")) list(fit = list(error = as.character(x)), fold = "?",
                                    moisture = "?", partition = "?", k = NA, region = NA) else x)
bad <- Filter(function(x) !is.null(x$fit$error), fits)
for (b in bad) log_err("fit ", b$fold, " ", b$moisture, " ", b$partition, " k", b$k, " r", b$region, ": ", b$fit$error)
fits <- Filter(function(x) is.null(x$fit$error), fits)
key <- function(f, mo, p, k, r) paste(f, mo, p, k, r)
fit_of <- setNames(lapply(fits, `[[`, "fit"), vapply(fits, function(x) key(x$fold, x$moisture, x$partition, x$k, x$region), ""))
params <- do.call(rbind, lapply(fits, function(x) data.frame(fold = x$fold, moisture = x$moisture, partition = x$partition,
  k = x$k, region = x$region, cells = x$n, loss = x$fit$loss, evals = x$fit$evals, t(gm_par_vec(x$fit$par)))))
utils::write.csv(params, file.path(out_dir, "params.csv"), row.names = FALSE)
log_msg(sprintf("fits done: %d (%.1f min); evaluations per fit: median %d, max %d", nrow(params),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), as.integer(stats::median(params$evals)), max(params$evals)))
# best moisture option per fold and unit (lowest loss in the fitting years)
best_fit <- function(f, p, k, r) {
  cand <- Filter(Negate(is.null), lapply(moists, function(mo) fit_of[[key(f, mo, p, k, r)]]))
  if (!length(cand)) return(NULL)
  cand[[which.min(vapply(cand, `[[`, 0, "loss"))]]
}

# ---- 3. predictions on held-out years and scores ----------------------------------------
n <- nrow(S$cells)
pred_par <- function(par, smooth, idx = seq_len(n)) gm_seasons(gm_gsi(sub_drv(idx), par, smooth), S$dates, years, g)
empty <- function() list(sos = matrix(NA_real_, n, nY), eos = matrix(NA_real_, n, nY))
models <- list()
models$default_fems <- pred_par(list(moisture = "none", tmin = unlist(fems$tmin), vpd = unlist(fems$vpd),
                                     photo_h = unlist(fems$photo_h)), fems$smooth)
models$default_nfdrs <- pred_par(list(moisture = "precip", tmin = unlist(nf$tmin), vpd = unlist(nf$vpd),
                                      photo_h = unlist(nf$photo_h), moist = c(nf$precip$lo, nf$precip$hi)), nf$smooth)
clim <- empty()
for (f in names(folds)) {
  tr <- folds[[f]]; te <- setdiff(seq_len(nY), tr)
  clim$sos[, te] <- apply(S$obs$sos[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
  clim$eos[, te] <- apply(S$obs$cure[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
}
models$climatology <- clim
cross <- function(get_groups) {               # each fold's parameters predict the other fold's years
  out <- empty()
  for (f in names(folds)) {
    te <- setdiff(seq_len(nY), folds[[f]])
    for (gr in get_groups(f)) {
      pr <- pred_par(gr$par, gmc$smooth, gr$idx)
      out$sos[gr$idx, te] <- pr$sos[, te]; out$eos[gr$idx, te] <- pr$eos[, te]
    }
  }
  out
}
conus_m <- list()
for (mo in c(moists, "best")) {
  get <- if (mo == "best") function(f) best_fit(f, "conus", 1, 1) else function(f) fit_of[[key(f, mo, "conus", 1, 1)]]
  conus_m[[mo]] <- cross(function(f) list(list(idx = seq_len(n), par = get(f)$par)))
}
reg_models <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]
  for (mo in c(moists, "best")) {
    getf <- function(f, r) {
      ft <- if (mo == "best") best_fit(f, nm, k, r) else fit_of[[key(f, mo, nm, k, r)]]
      ft %||% (if (mo == "best") best_fit(f, "conus", 1, 1) else fit_of[[key(f, mo, "conus", 1, 1)]])
    }
    reg_models[[paste(nm, k, mo, sep = "|")]] <- cross(function(f)
      lapply(sort(unique(reg[!is.na(reg)])), function(r) list(idx = which(reg == r), par = getf(f, r)$par)))
  }
}
obs_all <- list(sos = S$obs$sos, cure = S$obs$cure)
sc <- function(m, anom = TRUE) gm_score(m, obs_all, gmc$penalty_days, gmc$cap_days, anom)
scores <- do.call(rbind, c(
  lapply(names(models), function(nm) data.frame(model = nm, partition = NA, k = NA, moisture = NA, t(sc(models[[nm]], nm != "climatology")))),
  lapply(names(conus_m), function(mo) data.frame(model = "conus", partition = NA, k = NA, moisture = mo, t(sc(conus_m[[mo]])))),
  lapply(names(reg_models), function(nm) { p <- strsplit(nm, "|", fixed = TRUE)[[1]]
    data.frame(model = "regional", partition = p[1], k = as.integer(p[2]), moisture = p[3], t(sc(reg_models[[nm]]))) })))
utils::write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)

# per region: best-moisture regional fit vs best-moisture CONUS fit, defaults and climatology
rs <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]; rm_ <- reg_models[[paste(nm, k, "best", sep = "|")]]
  for (r in sort(unique(reg[!is.na(reg)]))) {
    i <- which(reg == r); ob <- list(sos = obs_all$sos[i, , drop = FALSE], cure = obs_all$cure[i, , drop = FALSE])
    subm <- function(m) list(sos = m$sos[i, , drop = FALSE], eos = m$eos[i, , drop = FALSE])
    for (mn in c("regional", "conus", "default_fems", "default_nfdrs", "climatology")) {
      m <- switch(mn, regional = rm_, conus = conus_m$best, models[[mn]])
      rs[[length(rs) + 1]] <- data.frame(partition = nm, k = k, region = r, cells = length(i),
                                         fitted = length(i) >= gmc$min_cells, model = mn,
                                         t(gm_score(subm(m), ob, gmc$penalty_days, gmc$cap_days, mn != "climatology")))
    }
  }
}
utils::write.csv(do.call(rbind, rs), file.path(out_dir, "region_scores.csv"), row.names = FALSE)

for (i in seq_len(nrow(scores))) with(scores[i, ], log_msg(sprintf(
  "  %-13s %-10s %3s %-7s green-up MAE %5.1f d, median %5.1f, <=30 d %3.0f%%, r %5.2f | curing MAE %5.1f d, median %5.1f, <=30 d %3.0f%%, r %5.2f",
  model, ifelse(is.na(partition), "", partition), ifelse(is.na(k), "", k), ifelse(is.na(moisture), "", moisture),
  sos_mae, sos_mdae, 100 * sos_within30, sos_anom_r, cure_mae, cure_mdae, 100 * cure_within30, cure_anom_r)))
msg <- sprintf("%d fits, %d cells, %.1f min", nrow(params), n, as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | outputs: ", out_dir, " | log: ", log_file)
notify("34_gsi_model finished", msg, priority = 3, tags = "white_check_mark")
