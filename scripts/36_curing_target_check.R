#!/usr/bin/env Rscript
# Step 36 (evaluation, target check): is the satellite curing date a good target where the GSI
# fails? The GSI model tests (34, 35) score the GSI against satellite green-up (NDVI sos20) and
# curing (NDII7 cure50). In the monsoon Southwest, the low deserts and the California annual
# grasslands the GSI does worse than a cell's typical date. Before blaming the GSI, compare the
# satellite dates with ground observations of the same events (config `target_check: pairs`):
# NPN grass green-up / curing / cured dates and Globe-LFMC herbaceous 50 % decline dates.
#   Level 1, site medians (always): each site's median ground date vs its 4 km cell's median
#     satellite date (summaries/cutb_metrics_median.tif). Per working region and CONUS: sites,
#     median bias (satellite - ground), mean absolute difference, share within 30 days,
#     rank correlation across sites (does the satellite put late sites late?).
#   Level 2, year by year (if per-year satellite metrics are present, cutb/metrics/
#     metrics_<Y>.tif): each site-year vs the cell's satellite date that year. Same scores,
#     plus the anomaly correlation (site-years minus the site's own mean): does the satellite
#     follow early and late years on the ground? Also the ground's own year-to-year spread.
# Dates are compared on a circular year (a satellite date of -20 = Dec 12 of the year before),
# so a difference is at most half a year; |difference| > 90 days is counted as a season
# mismatch (e.g. spring vs monsoon).
# Outputs <run_dir>/eval/target_check/: sites.csv, site_years.csv (level 2), summary.csv,
# figs/target_check_sites.png. Run on the laptop (study data live there); level 2 needs the
# per-year metrics copied from the P720 (~25 files).
#
# Usage: Rscript scripts/36_curing_target_check.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/36_curing_target_check.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("36_curing_target_check", profile)
tc <- cfg$target_check
out_dir <- file.path(cfg$run_dir, "eval", "target_check"); fig_dir <- file.path(out_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
pairs <- do.call(rbind, lapply(tc$pairs, as.data.frame))
log_msg("=== 36_curing_target_check | profile: ", profile, " | ", nrow(pairs), " comparisons")

# ---- ground observations, cells, regions ---------------------------------------------------
ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
E <- obs_events(ddir, unique(pairs$set))
if (!nrow(E)) stop("no ground observations for ", paste(unique(pairs$set), collapse = ", "), " in ", ddir)
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
E$cell <- terra::cellFromXY(grid, as.matrix(E[, c("lon", "lat")]))
E$row <- match(E$cell, land)
E <- E[!is.na(E$row), ]
rs <- list(partition = tc$partition, k = tc$k)
D <- region_display(cfg, grid, list(rs))[[1]]
if (is.null(D)) stop("region set ", tc$partition, " k", tc$k, " not found under clusters/")
E$region <- D$num[E$row]
hm_f <- file.path(static_dir(cfg), "herb_mask.tif")
herb <- if (file.exists(hm_f)) { v <- land_values(terra::rast(hm_f), grid)[, 1]; !is.na(v) & v > 0 } else rep(NA, length(land))
E$herb_cell <- herb[E$row]
log_msg(sprintf("ground: %d site-years at %d sites (%s); region set %s k%d", nrow(E), length(unique(E$site)),
                paste(sprintf("%s %d", names(table(E$set)), table(E$set)), collapse = ", "), tc$partition, as.integer(tc$k)))

# circular difference (days), in (-182.5, 182.5]
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
score <- function(g, s, site = NULL) {          # g ground, s satellite (same rows)
  d <- cdiff(s, g); ok <- !is.na(d); d <- d[ok]; g <- g[ok]
  if (length(d) < 3) return(c(n = length(d), bias = NA, mad = NA, within30 = NA, mismatch90 = NA, r = NA))
  sa <- g + d                                  # satellite date aligned to the ground's year
  c(n = length(d), bias = stats::median(d), mad = mean(abs(d)), within30 = mean(abs(d) <= 30),
    mismatch90 = mean(abs(d) > 90),
    r = if (stats::sd(g) > 0 && stats::sd(sa) > 0) stats::cor(g, sa, method = "spearman") else NA)
}

# ---- level 1: site medians vs median satellite dates ----------------------------------------
med <- terra::rast(summary_file(cfg, "cutb_metrics", "median"))
miss <- setdiff(unique(pairs$metric), names(med))
if (length(miss)) stop("cutb_metrics_median lacks ", paste(miss, collapse = ", "))
M <- land_values(med, grid, unique(pairs$metric))
S <- stats::aggregate(cbind(doy, lon, lat) ~ set + site + row + region + herb_cell, E, stats::median)
ny <- stats::aggregate(year ~ set + site, E, function(y) length(unique(y)))
S$n_years <- ny$year[match(paste(S$set, S$site), paste(ny$set, ny$site))]
summ <- data.frame(); sites_out <- data.frame()
regs <- c(list(CONUS = NULL), stats::setNames(as.list(sort(unique(S$region))), sort(unique(S$region))))
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i, ]; s <- S[S$set == p$set, ]; s$sat <- M[s$row, p$metric]
  s$diff <- cdiff(s$sat, s$doy)
  sites_out <- rbind(sites_out, data.frame(set = p$set, metric = p$metric, s[, c("site", "lon", "lat", "region", "herb_cell", "n_years", "doy", "sat", "diff")]))
  for (rn in names(regs)) {
    u <- if (rn == "CONUS") rep(TRUE, nrow(s)) else s$region %in% regs[[rn]]
    for (dom in c("all", "herb")) {
      v <- u & (dom == "all" | s$herb_cell %in% TRUE)
      summ <- rbind(summ, data.frame(level = "site_medians", set = p$set, metric = p$metric, region = rn, domain = dom,
                                     t(score(s$doy[v], s$sat[v])), anom_r = NA, ground_spread = NA))
    }
  }
}
utils::write.csv(sites_out, file.path(out_dir, "sites.csv"), row.names = FALSE)

# ---- level 2: site-years vs that year's satellite date --------------------------------------
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
yrs <- sort(unique(E$year)); have <- yrs[file.exists(vapply(yrs, mfile, ""))]
if (!length(have)) {
  log_warn("no per-year satellite metrics under ", file.path(cfg$run_dir, "cutb", "metrics"),
           ": level 2 (year by year) skipped. Copy cutb/metrics/metrics_<year>.tif from the P720 to add it.")
} else {
  log_msg("level 2: per-year satellite metrics for ", length(have), " years (", min(have), "-", max(have), ")")
  Y <- E[E$year %in% have, ]
  Y$key <- paste(Y$row, Y$year)
  sat_y <- list()
  for (yy in have) {
    r <- terra::rast(mfile(yy)); rows <- unique(Y$row[Y$year == yy])
    v <- terra::extract(r[[unique(pairs$metric)]], land[rows])
    sat_y[[as.character(yy)]] <- data.frame(row = rows, year = yy, v, check.names = FALSE)
  }
  SY <- do.call(rbind, sat_y)
  yr_out <- data.frame()
  for (i in seq_len(nrow(pairs))) {
    p <- pairs[i, ]; y <- Y[Y$set == p$set, ]
    y$sat <- SY[[p$metric]][match(y$key, paste(SY$row, SY$year))]
    y$diff <- cdiff(y$sat, y$doy)
    # anomalies: sites with >= min_years years of both
    ok <- !is.na(y$diff)
    cnt <- stats::ave(as.numeric(ok), y$site, FUN = sum)
    a <- ok & cnt >= (tc$min_years %||% 3)
    y$g_anom <- NA_real_; y$s_anom <- NA_real_
    if (any(a)) {
      sa <- y$doy + y$diff                                       # satellite aligned to ground's year
      y$g_anom[a] <- y$doy[a] - stats::ave(y$doy[a], y$site[a])
      y$s_anom[a] <- sa[a] - stats::ave(sa[a], y$site[a])
    }
    yr_out <- rbind(yr_out, data.frame(set = p$set, metric = p$metric,
                                       y[, c("site", "year", "region", "herb_cell", "doy", "sat", "diff", "g_anom", "s_anom")]))
    for (rn in names(regs)) {
      u <- if (rn == "CONUS") rep(TRUE, nrow(y)) else y$region %in% regs[[rn]]
      for (dom in c("all", "herb")) {
        v <- u & (dom == "all" | y$herb_cell %in% TRUE)
        w <- v & !is.na(y$g_anom) & !is.na(y$s_anom)
        ar <- if (sum(w) >= 10 && stats::sd(y$g_anom[w]) > 0 && stats::sd(y$s_anom[w]) > 0) stats::cor(y$g_anom[w], y$s_anom[w]) else NA
        summ <- rbind(summ, data.frame(level = "site_years", set = p$set, metric = p$metric, region = rn, domain = dom,
                                       t(score(y$doy[v], y$sat[v])), anom_r = ar,
                                       ground_spread = if (sum(w)) mean(abs(y$g_anom[w])) else NA))
      }
    }
  }
  utils::write.csv(yr_out, file.path(out_dir, "site_years.csv"), row.names = FALSE)
}
utils::write.csv(summ, file.path(out_dir, "summary.csv"), row.names = FALSE)

# ---- log and figure ---------------------------------------------------------------------------
for (lv in unique(summ$level)) for (i in seq_len(nrow(pairs))) {
  d <- summ[summ$level == lv & summ$set == pairs$set[i] & summ$metric == pairs$metric[i] & summ$domain == "all" & summ$n >= 3, ]
  log_msg(sprintf("%-12s %-18s vs %-7s %s", lv, pairs$set[i], pairs$metric[i],
                  paste(sprintf("%s: n %d, bias %+.0f, |d| %.0f, r %.2f%s", d$region, as.integer(d$n), d$bias, d$mad, d$r,
                                ifelse(is.na(d$anom_r), "", sprintf(", anom r %.2f", d$anom_r))), collapse = " | ")))
}
nc <- min(3, nrow(pairs)); nr <- ceiling(nrow(pairs) / nc)
png(file.path(fig_dir, "target_check_sites.png"), 560 * nc, 520 * nr, res = 130)
par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
for (i in seq_len(nrow(pairs))) {
  s <- sites_out[sites_out$set == pairs$set[i] & sites_out$metric == pairs$metric[i] & !is.na(sites_out$diff), ]
  sa <- s$doy + s$diff; lim <- range(c(s$doy, sa, 0, 365), na.rm = TRUE)
  plot(s$doy, sa, col = D$pal[s$region], pch = 16, cex = 0.8, xlim = lim, ylim = lim,
       xlab = paste("ground:", pairs$set[i], "(day of year)"), ylab = paste("satellite:", pairs$metric[i]),
       main = pairs$label[i] %||% paste(pairs$set[i], "vs", pairs$metric[i]), cex.main = 0.8)
  abline(0, 1); abline(30, 1, lty = 3); abline(-30, 1, lty = 3)
}
dev.off()
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
