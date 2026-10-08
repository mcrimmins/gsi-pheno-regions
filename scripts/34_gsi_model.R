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
#      Fitted thresholds: Tmin ramp, VPD ramp, start of the 1-hour daylength ramp, 28-day
#      precipitation ramp (FEMS smoothing). GSI dates = the smoothed index's own sos20 / eos50
#      (relative crossings around its yearly peak, as for the satellite curves).
#   3. Scores: mean absolute error (days, each error capped at cap_days; a missing GSI date
#      where the satellite has one counts penalty_days), median absolute error, share within
#      30 days, coverage, and the correlation of yearly anomalies (how well each model tracks
#      early and late years within a cell; climatology has none). The fits minimize the capped
#      mean error.
# Outputs eval/gsi_model/: scores.csv (model x partition x k), region_scores.csv (per region),
# params.csv (fitted thresholds), sample.rds / drivers.rds (cache).
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
  rd <- function(v) {
    parts <- pmap(py, function(y) {
      r <- terra::rast(pfile(v, y)); x <- as.matrix(terra::extract(r, cells$cell))
      colnames(x) <- names(r); x })
    if (any(vapply(parts, inherits, TRUE, "try-error"))) stop("reading ", v, " failed: ", parts[vapply(parts, inherits, TRUE, "try-error")][[1]])
    do.call(cbind, parts)
  }
  drv <- list(tmin = rd("tmin"), vpd_pa = rd("vpdmax") * 100, ppt = rd("ppt"))
  dates <- as.Date(colnames(drv$tmin))
  if (any(diff(dates) != 1)) stop("daily series not contiguous")
  S <- list(cells = cells, obs = obs, drv = drv, dates = dates, years = years)
  saveRDS(S, cache)
  log_msg(sprintf("drivers: %d cells x %d days (%s .. %s) cached", nrow(cells), length(dates), min(dates), max(dates)))
}
S$drv$photo_s <- m_photoperiod(S$cells$lat, as.integer(format(S$dates, "%j")), g$photo_method)
years <- S$years; nY <- length(years)
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))

# regions of the sampled cells for each partition x k
P <- eval_partitions(cfg, grid, unlist(gmc$k))
P <- P[intersect(unlist(gmc$partitions), names(P))]
units <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  reg <- P[[nm]][S$cells$row, kn]
  for (r in sort(unique(reg[!is.na(reg)]))) units[[length(units) + 1]] <- list(partition = nm, k = as.integer(sub("k", "", kn)), region = r, idx = which(reg == r))
}
log_msg(length(units), " regions in ", length(P), " partitions x ", length(gmc$k), " k values; folds: ",
        paste(sprintf("%s (%d years)", names(folds), lengths(folds)), collapse = ", "))

sub_drv <- function(idx) lapply(S$drv, function(m) m[idx, , drop = FALSE])
sub_obs <- function(idx, yi) list(sos = S$obs$sos[idx, yi, drop = FALSE], cure = S$obs$cure[idx, yi, drop = FALSE])
fems <- g$param_sets$fems; nf <- g$param_sets$nfdrs
start <- list(tmin = unlist(fems$tmin), vpd = unlist(fems$vpd), photo_h = unlist(fems$photo_h),
              prec = c(fems$precip$lo, fems$precip$hi))

# ---- 2. fits ------------------------------------------------------------------------------
set.seed(gmc$seed + 1)
conus_idx <- sort(sample(nrow(S$cells), min(gmc$conus_fit_cells, nrow(S$cells))))
jobs <- list()
for (f in names(folds)) {
  jobs[[length(jobs) + 1]] <- list(fold = f, partition = "conus", k = 1L, region = 1L, idx = conus_idx)
  for (u in units) if (length(u$idx) >= gmc$min_cells) jobs[[length(jobs) + 1]] <- c(list(fold = f), u)
}
log_msg(length(jobs), " fits (maxit ", gmc$maxit, ")")
fits <- pmap(jobs, function(j) {
  tr <- folds[[j$fold]]
  ft <- tryCatch(gm_fit(sub_drv(j$idx), sub_obs(j$idx, tr), S$dates, years[tr], g, gmc, start),
                 error = function(e) list(error = conditionMessage(e)))
  c(j[c("fold", "partition", "k", "region")], list(fit = ft, n = length(j$idx)))
})
fits <- lapply(fits, function(x) if (inherits(x, "try-error")) list(fit = list(error = as.character(x)), fold = "?",
                                    partition = "?", k = NA, region = NA) else x)
bad <- Filter(function(x) !is.null(x$fit$error), fits)
for (b in bad) log_err("fit ", b$fold, " ", b$partition, " k", b$k, " r", b$region, ": ", b$fit$error)
fits <- Filter(function(x) is.null(x$fit$error), fits)
key <- function(f, p, k, r) paste(f, p, k, r)
fit_of <- setNames(lapply(fits, `[[`, "fit"), vapply(fits, function(x) key(x$fold, x$partition, x$k, x$region), ""))
params <- do.call(rbind, lapply(fits, function(x) data.frame(fold = x$fold, partition = x$partition, k = x$k,
  region = x$region, cells = x$n, loss = x$fit$loss, evals = x$fit$evals, t(gm_par_vec(x$fit$par)))))
utils::write.csv(params, file.path(out_dir, "params.csv"), row.names = FALSE)
log_msg(sprintf("fits done: %d (%.1f min)", nrow(params), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- 3. predictions on held-out years and scores ----------------------------------------
n <- nrow(S$cells)
pred_par <- function(par, smooth, idx = seq_len(n)) gm_seasons(gm_gsi(sub_drv(idx), par, smooth), S$dates, years, g)
empty <- function() list(sos = matrix(NA_real_, n, nY), eos = matrix(NA_real_, n, nY))
models <- list()
# defaults (no fitting; the same for both folds)
models$default_fems <- pred_par(list(tmin = unlist(fems$tmin), vpd = unlist(fems$vpd), photo_h = unlist(fems$photo_h), prec = NULL), fems$smooth)
models$default_nfdrs <- pred_par(list(tmin = unlist(nf$tmin), vpd = unlist(nf$vpd), photo_h = unlist(nf$photo_h),
                                      prec = c(nf$precip$lo, nf$precip$hi)), nf$smooth)
# climatology: median of the fitting years' satellite dates
clim <- empty()
for (f in names(folds)) {
  tr <- folds[[f]]; te <- setdiff(seq_len(nY), tr)
  ms <- apply(S$obs$sos[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
  mc <- apply(S$obs$cure[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
  clim$sos[, te] <- ms; clim$eos[, te] <- mc
}
models$climatology <- clim
# conus and regional fits: each fold's parameters predict the other fold's years
cross <- function(get_par) {
  out <- empty()
  for (f in names(folds)) {
    te <- setdiff(seq_len(nY), folds[[f]])
    groups <- get_par(f)                       # list of list(idx, par)
    for (gr in groups) {
      pr <- pred_par(gr$par, gmc$smooth, gr$idx)
      out$sos[gr$idx, te] <- pr$sos[, te]; out$eos[gr$idx, te] <- pr$eos[, te]
    }
  }
  out
}
models$conus <- cross(function(f) list(list(idx = seq_len(n), par = fit_of[[key(f, "conus", 1, 1)]]$par)))
reg_models <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]
  reg_models[[paste0("regional|", nm, "|", k)]] <- cross(function(f) {
    lapply(sort(unique(reg[!is.na(reg)])), function(r) {
      ft <- fit_of[[key(f, nm, k, r)]] %||% fit_of[[key(f, "conus", 1, 1)]]
      list(idx = which(reg == r), par = ft$par)
    })
  })
}
obs_all <- list(sos = S$obs$sos, cure = S$obs$cure)
sc <- function(m, anom = TRUE) gm_score(m, obs_all, gmc$penalty_days, gmc$cap_days, anom)
scores <- do.call(rbind, c(
  lapply(names(models), function(nm) data.frame(model = nm, partition = NA, k = NA, t(sc(models[[nm]], nm != "climatology")))),
  lapply(names(reg_models), function(nm) { p <- strsplit(nm, "|", fixed = TRUE)[[1]]
    data.frame(model = "regional", partition = p[2], k = as.integer(p[3]), t(sc(reg_models[[nm]]))) })))
utils::write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)

# per region (regional fit vs conus fit vs defaults vs climatology, on the region's cells)
rs <- list()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); reg <- P[[nm]][S$cells$row, kn]; rm_ <- reg_models[[paste0("regional|", nm, "|", k)]]
  for (r in sort(unique(reg[!is.na(reg)]))) {
    i <- which(reg == r); ob <- list(sos = obs_all$sos[i, , drop = FALSE], cure = obs_all$cure[i, , drop = FALSE])
    subm <- function(m) list(sos = m$sos[i, , drop = FALSE], eos = m$eos[i, , drop = FALSE])
    for (mn in c("regional", "conus", "default_fems", "default_nfdrs", "climatology")) {
      m <- if (mn == "regional") rm_ else models[[mn]]
      rs[[length(rs) + 1]] <- data.frame(partition = nm, k = k, region = r, cells = length(i),
                                         fitted = length(i) >= gmc$min_cells, model = mn,
                                         t(gm_score(subm(m), ob, gmc$penalty_days, gmc$cap_days, mn != "climatology")))
    }
  }
}
utils::write.csv(do.call(rbind, rs), file.path(out_dir, "region_scores.csv"), row.names = FALSE)

for (i in seq_len(nrow(scores))) with(scores[i, ], log_msg(sprintf(
  "  %-13s %-10s %3s  green-up MAE %5.1f d, median %5.1f, <=30 d %3.0f%%, r %5.2f | curing MAE %5.1f d, median %5.1f, <=30 d %3.0f%%, r %5.2f",
  model, ifelse(is.na(partition), "", partition), ifelse(is.na(k), "", k), sos_mae, sos_mdae, 100 * sos_within30,
  sos_anom_r, cure_mae, cure_mdae, 100 * cure_within30, cure_anom_r)))
msg <- sprintf("%d fits, %d cells, %.1f min", nrow(params), n, as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | outputs: ", out_dir, " | log: ", log_file)
notify("34_gsi_model finished", msg, priority = 3, tags = "white_check_mark")
