#!/usr/bin/env Rscript
# Step 67 (fire activity study; P720, ~20-40 min): the forecast test's phenology (scope decision 23c).
# Would phenology help if we only had curing dates a manager could know two weeks ahead? For every
# pilot region (not only the rules package's regions 3 and 11) the lead-time weather adjustment of
# 38 / 40 (R/cues.R adj_cues / adj_fit: typical curing date + linear adjustment of <= 3 cue
# departures, cues ending at the typical date minus rules: lead_days) is refit LEAVING EACH YEAR OUT:
# the year's predicted curing date comes from coefficients and cue climatology fitted on the other
# years, as the fire test does. Phenology cells, weights and typical dates as 65 (rules typical
# dates; herb cells, fuel share >= 0.35 in pheno_source regions, fuel weighted in
# herb_weighted_regions). Green-up stays at the typical date (the rules forecast curing only).
# Outputs <run_dir>/fire_study/forecast/:
#   forecast_weekly.csv   region x year x week: cured_share_fc (green from the typical green-up to the
#                         predicted curing date; neighbouring years typical), cured_share_typ_fc (all
#                         typical, same cells), curing_adj_med (weighted median predicted adjustment, d)
#   forecast_skill.csv    per region: held-out MAE of the typical vs the adjusted curing date against
#                         satellite eos50 (errors capped at 90 d), r of predicted vs observed departure
#   forecast_cues.csv     per region and held-out year: cues picked and coefficients
#   cues_<Y>_<profile>.rds  cue cache per year (delete to recompute)
# 66 reads forecast_weekly.csv when present (specs all_fc, human_unk_fc, px_fc).
#
# Usage: Rscript scripts/67_fire_forecast.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "fire_weekly.R"))
source(here::here("R", "cues.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("67_fire_forecast", profile)
notify_init(cfg)
fs <- cfg$fire_study; rc <- cfg$rules %||% list(); ca <- cfg$curing_adjust %||% list()
t0 <- Sys.time(); terra::terraOptions(progress = 0)
out_dir <- fire_dir(cfg, "forecast"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
regs <- as.integer(unlist(fs$pilot_regions))
yrs <- seq(as.integer(fs$fod$years$start), as.integer(fs$fod$years$end))
lead <- as.numeric(rc$lead_days %||% 14); max_adj <- as.numeric(rc$max_adj_days %||% 60)
jump <- as.numeric(ca$max_abs_anomaly %||% 90); cap <- as.numeric(ca$cap_days %||% 90)
nw <- as.integer(fs$week$n_weeks %||% 52)
sw <- step_workers(cfg, "fire_forecast")
log_msg("=== 67_fire_forecast | profile: ", profile, " | regions ", paste(regs, collapse = ", "), " | years ", min(yrs), "-", max(yrs),
        " | lead ", lead, " d | workers ", sw$n)

# ---- cells, weights, typical dates, satellite dates (as 65) ---------------------------------------------
fuel <- fire_fuel_share(cfg)
TD <- terra::rast(file.path(cfg$run_dir, "rules", "typical_dates.tif"))
tdv <- terra::values(TD[[c("greenup_doy", "curing_doy", "herb", "region")]])
region <- as.integer(tdv[, "region"])
ps <- fs$pheno_source %||% list()
cells_of <- function(r) {
  if (r %in% as.integer(unlist(ps$regions))) which(region == r & !is.na(fuel) & fuel >= as.numeric(ps$min_fuel_share %||% 0.35))
  else which(region == r & tdv[, "herb"] %in% 1)
}
CELLS <- lapply(stats::setNames(regs, regs), cells_of)
allc <- sort(unique(unlist(CELLS)))
TG <- tdv[allc, "greenup_doy"]; TC <- tdv[allc, "curing_doy"]
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
SOS <- EOS <- matrix(NA_real_, length(allc), length(yrs), dimnames = list(NULL, yrs))
for (j in seq_along(yrs)) if (file.exists(mfile(yrs[j]))) {
  m <- as.matrix(terra::extract(terra::rast(mfile(yrs[j]))[[c(fs$greenup_metric %||% "sos20", fs$curing_metric %||% "eos50")]], allc))
  SOS[, j] <- m[, 1]; EOS[, j] <- m[, 2]
}
log_msg(sprintf("%s phenology cells; satellite curing dates for %.0f %% of cell-years", format(length(allc), big.mark = ","), 100 * mean(!is.na(EOS))))

# ---- cues at the cutoff, per year (cached) ----------------------------------------------------------------
pfile <- function(v, y) prism_out_file(cfg, v, y)
cue_year <- function(j) {
  Y <- yrs[j]; f <- file.path(out_dir, sprintf("cues_%d_%s.rds", Y, profile))
  if (file.exists(f)) return(f)
  tryCatch({
    terra::terraOptions(memfrac = sw$memfrac, progress = 0)
    ys <- c(Y - 1, Y, Y + 1); ys <- ys[file.exists(pfile("ppt", ys)) & file.exists(pfile("tmax", ys))]
    if (!(Y %in% ys)) stop("no PRISM for ", Y)
    rd <- function(v) do.call(cbind, lapply(ys, function(y) as.matrix(terra::extract(terra::rast(pfile(v, y)), allc))))
    P <- rd("ppt"); TX <- rd("tmax")
    dates <- do.call(c, lapply(ys, function(y) as.Date(names(terra::rast(pfile("ppt", y))))))
    C <- adj_cues(P, TX, dates, Y, TC - lead, SOS[, j])
    saveRDS(C, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f); f
  }, error = function(e) paste("ERROR", Y, conditionMessage(e)))
}
need <- which(!file.exists(file.path(out_dir, sprintf("cues_%d_%s.rds", yrs, profile))))
log_msg(length(need), " year(s) of cues to compute")
if (length(need)) {
  res <- if (.Platform$OS.type == "unix" && sw$n > 1) parallel::mclapply(need, cue_year, mc.cores = min(sw$n, length(need)), mc.preschedule = FALSE) else lapply(need, cue_year)
  bad <- grep("^ERROR", unlist(res), value = TRUE)
  if (length(bad)) stop(paste(bad, collapse = "; "))
  log_msg(sprintf("cues done, %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
CU <- lapply(stats::setNames(adj_cue_names, adj_cue_names), function(x) matrix(NA_real_, length(allc), length(yrs)))
for (j in seq_along(yrs)) { C <- readRDS(file.path(out_dir, sprintf("cues_%d_%s.rds", yrs[j], profile))); for (cn in adj_cue_names) CU[[cn]][, j] <- C[, cn] }

# ---- leave-one-year-out adjustment per region ----------------------------------------------------------------
DEP <- EOS - TC                                         # observed departure from the typical curing date
ADJ <- matrix(NA_real_, length(allc), length(yrs))
SK <- list(); CO <- list()
for (r in regs) {
  i <- match(CELLS[[as.character(r)]], allc)
  if (length(i) < 25) { log_warn("region ", r, ": ", length(i), " cells, skipped"); next }
  for (j in seq_along(yrs)) {
    tr <- setdiff(seq_along(yrs), j)
    X <- lapply(CU, function(x) { m <- rowMeans(x[i, tr, drop = FALSE], na.rm = TRUE); x[i, , drop = FALSE] - m })
    dtr <- DEP[i, tr, drop = FALSE]; dtr[abs(dtr) > jump] <- NA
    fx <- adj_fit(dtr, lapply(X, function(x) x[, tr, drop = FALSE]), lapply(X, function(x) x[, j, drop = FALSE]),
                  ca$max_cues %||% 3, ca$min_gain %||% 0.01, max_adj)
    ADJ[i, j] <- fx$adj[, 1]
    if (length(fx$cues)) CO[[length(CO) + 1]] <- data.frame(region = r, held_out = yrs[j], cue = fx$cues, b = unname(fx$b), r2_fit = fx$r2_fit)
  }
  obs <- DEP[i, ]; a <- ADJ[i, ]; u <- !is.na(obs)
  nj <- u & abs(obs) <= jump
  SK[[length(SK) + 1]] <- data.frame(region = r, cells = length(i), cell_years = sum(u),
    mae_typical = mean(pmin(abs(obs[u]), cap)), mae_adjusted = mean(pmin(abs(obs[u] - a[u]), cap)),
    r_adj = if (sum(nj) > 30 && stats::sd(a[nj]) > 0) stats::cor(a[nj], obs[nj]) else NA_real_,
    sd_adj = stats::sd(a[u]))
  log_msg(sprintf("region %2d: %s cells | held-out curing MAE typical %.1f d, adjusted %.1f d | r %.2f | adjustment sd %.1f d",
                  r, format(length(i), big.mark = ","), SK[[length(SK)]]$mae_typical, SK[[length(SK)]]$mae_adjusted,
                  SK[[length(SK)]]$r_adj, SK[[length(SK)]]$sd_adj))
}
utils::write.csv(do.call(rbind, SK), file.path(out_dir, "forecast_skill.csv"), row.names = FALSE)
if (length(CO)) utils::write.csv(do.call(rbind, CO), file.path(out_dir, "forecast_cues.csv"), row.names = FALSE)

# ---- weekly forecast cured share ------------------------------------------------------------------------------------
mid <- fw_mid_doy(nw); OUT <- list()
for (r in regs) {
  i <- match(CELLS[[as.character(r)]], allc)
  if (length(i) < 25) next
  wt <- if (r %in% as.integer(unlist(fs$herb_weighted_regions))) fuel[allc[i]] else rep(1, length(i))
  tg <- TG[i]; tc <- TC[i]; okT <- !is.na(tg) & !is.na(tc)
  for (j in seq_along(yrs)) {
    Y <- yrs[j]
    lp <- as.numeric(as.Date(sprintf("%d-01-01", Y)) - as.Date(sprintf("%d-01-01", Y - 1)))
    lt <- as.numeric(as.Date(sprintf("%d-01-01", Y + 1)) - as.Date(sprintf("%d-01-01", Y)))
    a <- ADJ[i, j]; a[is.na(a)] <- 0
    G <- fw_green(mid, cbind(tg, tg, tg), cbind(tc, tc + a, tc), lp, lt)
    GT <- fw_green(mid, cbind(tg, tg, tg), cbind(tc, tc, tc), lp, lt)
    OUT[[length(OUT) + 1]] <- data.frame(region = r, year = Y, week = seq_len(nw),
      cured_share_fc = vapply(seq_len(nw), function(k) fw_wmean(as.numeric(!G[okT, k]), wt[okT]), 0),
      cured_share_typ_fc = vapply(seq_len(nw), function(k) fw_wmean(as.numeric(!GT[okT, k]), wt[okT]), 0),
      curing_adj_med = fw_wmedian(ADJ[i, j], wt))
  }
}
utils::write.csv(do.call(rbind, OUT), file.path(out_dir, "forecast_weekly.csv"), row.names = FALSE)
msg <- sprintf("%d regions x %d years, %.1f min", length(regs), length(yrs), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("67_fire_forecast finished", msg, priority = 3, tags = "white_check_mark")
