#!/usr/bin/env Rscript
# Step 65c (fire study monitoring test, Oct 20 2026; P720 or laptop, minutes): a weekly real-time cured
# share from daily VIIRS (05_viirs_daily.R), the near-daily version of 65b. The delay test (65b
# GSI_RT_LAG) found the growth signal needs greenness about a week old or less; here each week uses
# only days available by then: a composite of the composite_days (8) days ending latency_days (3)
# before the week's middle day.
# Per cell and week: g_now = composite NDVI of that window, g_prev = the window ending prev_offset_days
# (16) earlier, both scaled to the cell's own VIIRS range (range_quantiles 0.05 / 0.95 of its weekly
# composites in the OTHER years). Green when g_now >= curing_frac (0.5) or (g_now >= greenup_frac (0.2)
# and rising); else cured; the green / cured rule and fractions are 65b's (config fire_study: realtime).
# Two composites: vd = clear-pixel-weighted mean of the days with >= min_clear_px clear 500 m pixels;
# vdx = maximum of those days' NDVI (maximum-value compositing, robust to leftover cloud).
# Shares over the region's phenology cells (fuel-share weights in herb_weighted_regions, as 65 / 65b);
# typical = mean of the other years, same region and week.
# Outputs <run_dir>/fire_study/realtime/: realtime_weekly_vd.csv (region, year, week, cured_share_vd /
#   _typ, coverage_vd, cured_share_vdx / _typ, coverage_vdx), realtime_check_vd.csv (per region: r of
#   weekly shares and departures with 65's after-the-fact cured share and with 65b's 16-day real-time
#   share). 66 adds specs all_vd / human_unk_vd / px_vd and *_vdx (regions with VIIRS daily data only;
#   GSI_FIRE_REGIONS=3,10 restricts every spec to the same regions for comparison).
#
# Usage: Rscript scripts/65c_fire_viirs_daily.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "features.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "fire_weekly.R"))
source(here::here("R", "viirs_daily.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("65c_fire_viirs_daily", profile)
notify_init(cfg)
fs <- cfg$fire_study; vd <- vd_conf(cfg); rt <- fs$realtime %||% list()
t0 <- Sys.time()
out_dir <- fire_dir(cfg, "realtime"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
nw <- as.integer(fs$week$n_weeks %||% 52); mid <- fw_mid_doy(nw)
comp <- as.integer(vd$composite_days %||% 8); poff <- as.integer(vd$prev_offset_days %||% 16)
lat_d <- as.integer(vd$latency_days %||% 3); minpx <- as.numeric(vd$min_clear_px %||% 4)
rq <- as.numeric(unlist(vd$range_quantiles %||% list(0.05, 0.95)))
gu <- as.numeric(rt$greenup_frac %||% 0.2); cu <- as.numeric(rt$curing_frac %||% 0.5)
yrs <- seq(as.integer(vd$years$start), as.integer(vd$years$end))
log_msg("=== 65c_fire_viirs_daily | profile: ", profile, " | regions ", paste(vd$regions, collapse = ","), " | years ", min(yrs), "-", max(yrs),
        sprintf(" | composite %d d ending %d d before mid-week, prev %d d earlier, >= %g clear px/day", comp, lat_d, poff, minpx))

read_year <- function(r, Y) {
  f <- Sys.glob(sub("\\.rds$", "_c*.rds", vd_out_file(cfg, vd$product, r, Y)))
  if (!length(f)) return(NULL)
  L <- lapply(f, readRDS)
  cells <- L[[1]]$cells
  if (!all(vapply(L, function(x) identical(x$cells, cells), TRUE))) stop("region ", r, " ", Y, ": chunks have different cells")
  d <- do.call(c, lapply(L, `[[`, "dates")); o <- order(d); keep <- o[!duplicated(d[o])]
  list(cells = cells, dates = d[keep], ndvi = do.call(cbind, lapply(L, `[[`, "ndvi"))[, keep, drop = FALSE],
       n_clear = do.call(cbind, lapply(L, `[[`, "n_clear"))[, keep, drop = FALSE])
}

OUT <- list(); CK <- list()
for (r in vd$regions) {
  D <- lapply(stats::setNames(yrs, yrs), function(Y) read_year(r, Y))
  have <- yrs[!vapply(D, is.null, TRUE)]
  if (!length(have)) { log_warn("region ", r, ": no VIIRS daily files; skipped"); next }
  cells <- D[[as.character(have[1])]]$cells; nc <- length(cells)
  # one calendar-day matrix over all years (cells x days)
  all_d <- seq(as.Date(sprintf("%d-01-01", min(have))), as.Date(sprintf("%d-12-31", max(have))), by = 1)
  ND <- matrix(NA_real_, nc, length(all_d)); NC <- matrix(0, nc, length(all_d))
  for (Y in have) { x <- D[[as.character(Y)]]; j <- match(x$dates, all_d); ND[, j] <- x$ndvi; NC[, j] <- x$n_clear }
  rm(D); gc()
  ok <- NC >= minpx & is.finite(ND); ND[!ok] <- NA; NC[!ok] <- 0
  log_msg(sprintf("region %d: %s cells, %d years (%s), usable cell-days %.0f %%", r, format(nc, big.mark = ","), length(have),
                  paste(range(have), collapse = "-"), 100 * mean(ok)))
  # composites ending on day index e (window e - comp + 1 .. e)
  composite <- function(e) {
    j <- (e - comp + 1):e; j <- j[j >= 1]
    if (!length(j)) return(list(mean = rep(NA_real_, nc), max = rep(NA_real_, nc)))
    w <- NC[, j, drop = FALSE]; v <- ND[, j, drop = FALSE]; v0 <- v; v0[is.na(v0)] <- 0
    sw <- rowSums(w); mx <- do.call(pmax, c(lapply(seq_len(ncol(v)), function(i) v[, i]), na.rm = TRUE))
    list(mean = ifelse(sw > 0, rowSums(v0 * w) / sw, NA_real_), max = mx)
  }
  CN <- CP <- list(mean = array(NA_real_, c(nc, nw, length(yrs))), max = array(NA_real_, c(nc, nw, length(yrs))))
  for (k in seq_along(yrs)) {
    Y <- yrs[k]; if (!Y %in% have) next
    for (w in seq_len(nw)) {
      cut <- as.Date(sprintf("%d-01-01", Y)) + mid[w] - 1 - lat_d
      e <- match(cut, all_d); ep <- match(cut - poff, all_d)
      if (!is.na(e)) { a <- composite(e); CN$mean[, w, k] <- a$mean; CN$max[, w, k] <- a$max }
      if (!is.na(ep)) { b <- composite(ep); CP$mean[, w, k] <- b$mean; CP$max[, w, k] <- b$max }
    }
  }
  rm(ND, NC); gc()
  # weights as 65 / 65b
  wt <- if (r %in% as.integer(unlist(fs$herb_weighted_regions))) fire_fuel_share(cfg)[cells] else rep(1, nc)
  res <- list()
  for (st in c("mean", "max")) {
    S <- array(NA_real_, c(nc, nw, length(yrs)))
    for (k in seq_along(yrs)) {
      if (!yrs[k] %in% have) next
      oth <- setdiff(which(yrs %in% have), k)
      if (length(oth) < 2) next
      q <- row_quantiles(matrix(CN[[st]][, , oth], nc), rq)
      base <- q[, 1]; amp <- q[, 2] - q[, 1]; amp[!is.finite(amp) | amp < 0.05] <- NA
      gn <- (CN[[st]][, , k] - base) / amp; gp <- (CP[[st]][, , k] - base) / amp
      green <- gn >= cu | (gn >= gu & !is.na(gp) & gn > gp)
      S[, , k] <- ifelse(is.na(gn), NA, ifelse(green, 0, 1))
    }
    res[[st]] <- S
  }
  for (k in seq_along(yrs)) {
    if (!yrs[k] %in% have) next
    OUT[[length(OUT) + 1]] <- data.frame(region = r, year = yrs[k], week = seq_len(nw),
      cured_share_vd = vapply(seq_len(nw), function(w) fw_wmean(res$mean[, w, k], wt), 0),
      coverage_vd = colMeans(!is.na(res$mean[, , k, drop = TRUE])),
      cured_share_vdx = vapply(seq_len(nw), function(w) fw_wmean(res$max[, w, k], wt), 0),
      coverage_vdx = colMeans(!is.na(res$max[, , k, drop = TRUE])))
  }
}
if (!length(OUT)) stop("no VIIRS daily data for any region (run 05_viirs_daily.R)")
O <- do.call(rbind, OUT)
for (v in c("vd", "vdx")) {
  s <- paste0("cured_share_", v); O[[paste0(s, "_typ")]] <- NA_real_
  for (k in seq_len(nrow(O))) {
    j <- O$region == O$region[k] & O$week == O$week[k] & O$year != O$year[k]
    O[[paste0(s, "_typ")]][k] <- mean(O[[s]][j], na.rm = TRUE)
  }
}
utils::write.csv(O, file.path(out_dir, "realtime_weekly_vd.csv"), row.names = FALSE)

# check against 65's after-the-fact cured share and 65b's 16-day real-time share
wf <- fire_dir(cfg, "weekly", "fire_weekly.csv"); rf <- fire_dir(cfg, "realtime", "realtime_weekly.csv")
if (file.exists(wf)) {
  W <- utils::read.csv(wf)[, c("region", "year", "week", "cured_share", "cured_share_typ")]
  m <- merge(O, W, by = c("region", "year", "week"))
  if (file.exists(rf)) m <- merge(m, utils::read.csv(rf)[, c("region", "year", "week", "cured_share_rt", "cured_share_rt_typ")], by = c("region", "year", "week"), all.x = TRUE)
  rr <- function(a, b) { u <- stats::complete.cases(a, b); if (sum(u) > 10) stats::cor(a[u], b[u]) else NA_real_ }
  CK <- do.call(rbind, lapply(unique(m$region), function(r) { x <- m[m$region == r, ]
    dep <- function(v) x[[paste0("cured_share_", v)]] - x[[paste0("cured_share_", v, "_typ")]]
    data.frame(region = r, r_share_vd = rr(x$cured_share_vd, x$cured_share), r_dep_vd = rr(dep("vd"), x$cured_share - x$cured_share_typ),
               r_share_vdx = rr(x$cured_share_vdx, x$cured_share), r_dep_vdx = rr(dep("vdx"), x$cured_share - x$cured_share_typ),
               r_dep_rt = if ("cured_share_rt" %in% names(x)) rr(dep("rt"), x$cured_share - x$cured_share_typ) else NA_real_,
               coverage_vd = mean(x$coverage_vd, na.rm = TRUE), coverage_vdx = mean(x$coverage_vdx, na.rm = TRUE)) }))
  utils::write.csv(CK, file.path(out_dir, "realtime_check_vd.csv"), row.names = FALSE)
  for (k in seq_len(nrow(CK))) with(CK[k, ], log_msg(sprintf("region %2d: daily vs after-the-fact share r %.2f (max comp %.2f); departures r %.2f (max %.2f; 16-day real-time %.2f); coverage %.0f %%",
                                                         region, r_share_vd, r_share_vdx, r_dep_vd, r_dep_vdx, r_dep_rt, 100 * coverage_vd)))
}
msg <- sprintf("%d region-years, %.1f min", nrow(unique(O[, c("region", "year")])), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("65c_fire_viirs_daily finished", msg, priority = 3, tags = "white_check_mark")
