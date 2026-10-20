#!/usr/bin/env Rscript
# Step 65b (fire activity study; P720, ~15-30 min): a cured share a manager could have had each week
# (decided Oct 19: the "monitor curing" test). 65's cured share uses green-up and curing dates fitted to
# the whole season's smoothed curve, i.e. after the fact. Here each week sees only the 16-day
# composites already available that week (composite start + 16 d + latency_days <= the week's middle
# day), from Terra and Aqua (MOD13A2 / MYD13A2 NDVI on the grid, cutb/<product>/ndvi/; composites with
# fewer than cutb_curves: min_valid valid pixels or more than max_snow_frac snow are dropped).
# Per cell and week: g_now = mean NDVI of the composites starting in the 16 days before the cut,
# g_prev = the 16 days before that, each scaled to the cell's typical NDVI range (base and max = median
# of the cutb/metrics ndvi_base / ndvi_max of the OTHER years). Green when g_now >= curing_frac (0.5,
# as the curing date eos50) or when g_now >= greenup_frac (0.2, as sos20) and rising (g_now > g_prev);
# otherwise cured. No smoothing across the season, no future composites.
# Cells, weights: as 65 (herb cells; fuel share >= 0.35 in pheno_source regions; fuel-weighted in
# herb_weighted_regions). Weeks with no valid composite in a cell leave the cell out that week.
# Smoothed nowcast (Oct 19, *_rts): each week the weighted Whittaker smoother of the curves step (4-day
# nodes, cutb_curves lambda and envelope) is fitted to the composites available so far, and the state
# is read from the smoothed curve at the latest data node (g_prev 16 days earlier); no data in the last
# 32 days -> no state.
# Outputs <run_dir>/fire_study/realtime/: realtime_weekly.csv (region x year x week: cured_share_rt,
# cured_share_rt_typ = mean of the other years, coverage = share of cells with a state), realtime_check.csv
# (per region: correlation of the weekly real-time and after-the-fact cured shares, departures too).
# 66 reads realtime_weekly.csv when present (specs all_rt, human_unk_rt, px_rt).
#
# Usage: Rscript scripts/65b_fire_realtime.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "cutb.R"))
source(here::here("R", "cutb_curves.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "fire_weekly.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("65b_fire_realtime", profile)
notify_init(cfg)
fs <- cfg$fire_study; rt <- fs$realtime %||% list(); cc <- cfg$cutb_curves %||% list()
t0 <- Sys.time(); terra::terraOptions(progress = 0)
out_dir <- fire_dir(cfg, "realtime"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
regs <- as.integer(unlist(fs$pilot_regions))
yrs <- seq(as.integer(fs$fod$years$start), as.integer(fs$fod$years$end))
nw <- as.integer(fs$week$n_weeks %||% 52); mid <- fw_mid_doy(nw)
lat_d <- as.numeric(rt$latency_days %||% 3); gu <- as.numeric(rt$greenup_frac %||% 0.2); cu <- as.numeric(rt$curing_frac %||% 0.5)
min_valid <- as.numeric(cc$min_valid %||% 5); max_snow <- as.numeric(cc$max_snow_frac %||% 0.5)
full_valid <- as.numeric(cc$full_valid %||% 20); lam <- as.numeric(cc$lambda %||% 30)
env_iter <- as.integer(cc$envelope$iter %||% 1); env_below <- as.numeric(cc$envelope$below_weight %||% 0.5)
sw <- step_workers(cfg, "fire_realtime")
log_msg("=== 65b_fire_realtime | profile: ", profile, " | years ", min(yrs), "-", max(yrs), " | latency ", lat_d, " d, green-up ", gu,
        ", curing ", cu, " | workers ", sw$n)

fuel <- fire_fuel_share(cfg)
TD <- terra::rast(file.path(cfg$run_dir, "rules", "typical_dates.tif"))
tdv <- terra::values(TD[[c("herb", "region")]]); region <- as.integer(tdv[, "region"])
ps <- fs$pheno_source %||% list()
cells_of <- function(r) {
  if (r %in% as.integer(unlist(ps$regions))) which(region == r & !is.na(fuel) & fuel >= as.numeric(ps$min_fuel_share %||% 0.35))
  else which(region == r & tdv[, "herb"] %in% 1)
}
CELLS <- lapply(stats::setNames(regs, regs), cells_of)
allc <- sort(unique(unlist(CELLS)))

# typical NDVI range per cell, from the other years' metrics
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
my <- yrs[file.exists(mfile(yrs))]
MB <- MX <- matrix(NA_real_, length(allc), length(my), dimnames = list(NULL, my))
for (j in seq_along(my)) { m <- as.matrix(terra::extract(terra::rast(mfile(my[j]))[[c("ndvi_base", "ndvi_max")]], allc)); MB[, j] <- m[, 1]; MX[, j] <- m[, 2] }
log_msg(sprintf("%s cells; NDVI range from %d metric years", format(length(allc), big.mark = ","), length(my)))

year_state <- function(Y) {
  tryCatch({
    terra::terraOptions(memfrac = sw$memfrac, progress = 0)
    oth <- setdiff(colnames(MB), as.character(Y))
    base <- apply(MB[, oth, drop = FALSE], 1, stats::median, na.rm = TRUE)
    top <- apply(MX[, oth, drop = FALSE], 1, stats::median, na.rm = TRUE)
    amp <- top - base; amp[!is.finite(amp) | amp < 0.05] <- NA
    comp <- cc_composites(cfg, c(Y - 1, Y))
    if (is.null(comp)) stop("no composites for ", Y)
    comp <- comp[comp$start >= as.Date(sprintf("%d-10-01", Y - 1)) & comp$start <= as.Date(sprintf("%d-12-31", Y)), ]
    V <- NULL; WV <- NULL
    rd <- function(p, y, v, bands) as.matrix(terra::extract(terra::rast(cutb_out_file(cfg, p, v, y))[[bands]], allc))
    for (g in split(seq_len(nrow(comp)), paste(comp$product, comp$year))) {
      p <- comp$product[g[1]]; y <- comp$year[g[1]]; b <- comp$band[g]
      nd <- rd(p, y, "ndvi", b); nv <- rd(p, y, "n_valid", b); sn <- rd(p, y, "snow_frac", b)
      bad <- is.na(nv) | nv < min_valid | (!is.na(sn) & sn > max_snow) | is.na(nd)
      nd[bad] <- NA; wv <- pmin(nv / full_valid, 1); wv[bad] <- 0
      colnames(nd) <- format(comp$start[g]); V <- cbind(V, nd); WV <- cbind(WV, wv)
    }
    st <- as.Date(colnames(V))
    G <- (V - base) / amp                                     # relative greenness per composite
    S <- matrix(NA_real_, length(allc), nw)
    for (w in seq_len(nw)) {
      cut <- as.Date(sprintf("%d-01-01", Y)) + mid[w] - 1 - 16 - lat_d   # latest usable composite start
      a <- which(st > cut - 16 & st <= cut); b <- which(st > cut - 32 & st <= cut - 16)
      gn <- if (length(a)) rowMeans(G[, a, drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(allc))
      gp <- if (length(b)) rowMeans(G[, b, drop = FALSE], na.rm = TRUE) else rep(NA_real_, length(allc))
      gn[!is.finite(gn)] <- NA; gp[!is.finite(gp)] <- NA
      green <- gn >= cu | (gn >= gu & !is.na(gp) & gn > gp)
      S[, w] <- ifelse(is.na(gn), NA, ifelse(green, 0, 1))  # 1 = cured
    }
    # smoothed nowcast (Oct 19): each week, the weighted Whittaker smoother of the curves step (4-day
    # nodes, lambda and upper-envelope pass from cutb_curves) fitted to the composites available so far
    # only; state from the smoothed value at the latest data node and 16 days before it
    o <- order(st); st <- st[o]; V <- V[, o, drop = FALSE]; WV <- WV[, o, drop = FALSE]
    nodes <- seq(as.Date(sprintf("%d-10-01", Y - 1)), as.Date(sprintf("%d-12-31", Y)), by = 4)
    ni <- pmax(1L, pmin(length(nodes), round(as.numeric(st + 8 - nodes[1]) / 4) + 1L))
    YN <- matrix(0, length(allc), length(nodes)); WN <- YN
    for (k in seq_along(st)) {                                # composites sharing a node: weighted mean
      w0 <- WN[, ni[k]]; w1 <- WV[, k]; y1 <- V[, k]; y1[is.na(y1)] <- 0
      tot <- w0 + w1; YN[, ni[k]] <- ifelse(tot > 0, (YN[, ni[k]] * w0 + y1 * w1) / pmax(tot, 1e-9), 0); WN[, ni[k]] <- tot
    }
    S2 <- matrix(NA_real_, length(allc), nw)
    for (w in seq_len(nw)) {
      cut <- as.Date(sprintf("%d-01-01", Y)) + mid[w] - 1 - 16 - lat_d
      av <- which(st <= cut); if (!length(av)) next
      L <- max(ni[av]); if (L < 8) next
      Wc <- WN[, 1:L, drop = FALSE]; late <- ni[-av]; late <- late[late <= L]
      if (length(late)) Wc[, unique(late)] <- 0                # composites not yet available
      Z <- whit_envelope(YN[, 1:L, drop = FALSE], Wc, lam, env_iter, env_below)
      recent <- rowSums(Wc[, max(1, L - 7):L, drop = FALSE]) > 0   # data in the last 32 days
      gn <- (Z[, L] - base) / amp; gp <- (Z[, L - 4] - base) / amp
      green <- gn >= cu | (gn >= gu & gn > gp)
      S2[, w] <- ifelse(!recent | is.na(gn), NA, ifelse(green, 0, 1))
    }
    f <- file.path(out_dir, sprintf("state_%d_%s.rds", Y, profile)); saveRDS(list(raw = S, smooth = S2), f); f
  }, error = function(e) paste("ERROR", Y, conditionMessage(e)))
}
need <- yrs[!file.exists(file.path(out_dir, sprintf("state_%d_%s.rds", yrs, profile)))]
log_msg(length(need), " year(s) to compute")
if (length(need)) {
  res <- if (.Platform$OS.type == "unix" && sw$n > 1) parallel::mclapply(need, year_state, mc.cores = min(sw$n, length(need)), mc.preschedule = FALSE) else lapply(need, year_state)
  bad <- grep("^ERROR", unlist(res), value = TRUE)
  if (length(bad)) stop(paste(bad, collapse = "; "))
}

OUT <- list()
for (Y in yrs) {
  SS <- readRDS(file.path(out_dir, sprintf("state_%d_%s.rds", Y, profile)))
  if (!is.list(SS)) stop("old state cache for ", Y, ": delete state_*.rds and rerun")
  S <- SS$raw; S2 <- SS$smooth
  for (r in regs) {
    i <- match(CELLS[[as.character(r)]], allc); if (!length(i)) next
    wt <- if (r %in% as.integer(unlist(fs$herb_weighted_regions))) fuel[allc[i]] else rep(1, length(i))
    OUT[[length(OUT) + 1]] <- data.frame(region = r, year = Y, week = seq_len(nw),
      cured_share_rt = vapply(seq_len(nw), function(w) fw_wmean(S[i, w], wt), 0),
      coverage = colMeans(!is.na(S[i, , drop = FALSE])),
      cured_share_rts = vapply(seq_len(nw), function(w) fw_wmean(S2[i, w], wt), 0),
      coverage_s = colMeans(!is.na(S2[i, , drop = FALSE])))
  }
}
O <- do.call(rbind, OUT)
O$cured_share_rt_typ <- NA_real_; O$cured_share_rts_typ <- NA_real_
for (k in seq_len(nrow(O))) {                           # typical = mean of the other years, same region and week
  j <- O$region == O$region[k] & O$week == O$week[k] & O$year != O$year[k]
  O$cured_share_rt_typ[k] <- mean(O$cured_share_rt[j], na.rm = TRUE)
  O$cured_share_rts_typ[k] <- mean(O$cured_share_rts[j], na.rm = TRUE)
}
utils::write.csv(O, file.path(out_dir, "realtime_weekly.csv"), row.names = FALSE)

# check against 65's after-the-fact cured share
wf <- fire_dir(cfg, "weekly", "fire_weekly.csv")
if (file.exists(wf)) {
  W <- utils::read.csv(wf)
  m <- merge(O, W[, c("region", "year", "week", "cured_share", "cured_share_typ")], by = c("region", "year", "week"))
  CK <- do.call(rbind, lapply(regs, function(r) { x <- m[m$region == r, ]; u <- stats::complete.cases(x[, c("cured_share_rt", "cured_share")])
    v <- stats::complete.cases(x[, c("cured_share_rts", "cured_share")])
    data.frame(region = r, r_share = stats::cor(x$cured_share_rt[u], x$cured_share[u]),
               r_departure = stats::cor((x$cured_share_rt - x$cured_share_rt_typ)[u], (x$cured_share - x$cured_share_typ)[u]),
               mean_coverage = mean(x$coverage, na.rm = TRUE),
               r_share_smooth = stats::cor(x$cured_share_rts[v], x$cured_share[v]),
               r_departure_smooth = stats::cor((x$cured_share_rts - x$cured_share_rts_typ)[v], (x$cured_share - x$cured_share_typ)[v]),
               mean_coverage_smooth = mean(x$coverage_s, na.rm = TRUE)) }))
  utils::write.csv(CK, file.path(out_dir, "realtime_check.csv"), row.names = FALSE)
  for (k in seq_len(nrow(CK))) with(CK[k, ], log_msg(sprintf("region %2d: real-time vs after-the-fact cured share r %.2f (smoothed %.2f), departures r %.2f (smoothed %.2f), coverage %.0f %% (%.0f %%)",
                                                          region, r_share, r_share_smooth, r_departure, r_departure_smooth, 100 * mean_coverage, 100 * mean_coverage_smooth)))
}
msg <- sprintf("%d years, %.1f min", length(yrs), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("65b_fire_realtime finished", msg, priority = 3, tags = "white_check_mark")
