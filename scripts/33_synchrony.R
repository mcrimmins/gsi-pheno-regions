#!/usr/bin/env Rscript
# Step 33 (evaluation, year-to-year synchrony): do cells in the same region green up and cure
# early or late in the same years? A region is a good domain for a phenology model when its
# cells move together from year to year, not only when their average seasons match.
# For each satellite season metric (config `synchrony: metrics`; per-year files from
# 06_cutb_curves.R) each cell's yearly anomaly is its value minus its own median across years
# (cells with fewer than `min_years` valid years are dropped). Anomalies larger than
# `max_abs_anomaly` days are dropped: they are jumps to a different season (e.g. a spring
# instead of a monsoon peak), not earlier or later timing of the same season. For every partition in config
# `evaluation: partitions` (+ the geographic baseline) and k:
#   sync_r2   share of the cells' anomaly variance explained by their region's mean anomaly in
#             the same year: 1 - sum((a - region-year mean)^2) / sum(a^2). 0 = regions share
#             nothing year to year; 1 = every cell follows its region exactly. The CONUS-wide
#             yearly mean (one region) is reported as k = 1.
# per region, for the region sets in config `region_profiles: sets`: sync_r2 within the region
# and the median correlation of a cell's anomalies with its region's mean anomaly series; and a
# raster of that per-cell correlation (cells that do not follow their region).
# Domains: strict herbaceous cells (default) and all cells.
# Outputs <run_dir>/eval/sync/: sync_r2.csv, sync_regions_<set>.csv, cellcorr_<set>.tif.
# Needs cutb/metrics/metrics_<Y>.tif (on the P720) and the cluster rasters. Minutes.
#
# Usage: Rscript scripts/33_synchrony.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/33_synchrony.R")

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
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("33_synchrony", profile)
sc <- cfg$synchrony; ev <- cfg$evaluation
ks <- seq(ev$k[[1]], ev$k[[2]])
out_dir <- file.path(cfg$run_dir, "eval", "sync"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
mfiles <- file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", cfg$years))
yrs <- cfg$years[file.exists(mfiles)]; mfiles <- mfiles[file.exists(mfiles)]
if (length(yrs) < 3) stop("need per-year Cut B metrics (cutb/metrics/metrics_<Y>.tif) for at least 3 years")
mets <- unlist(sc$metrics)
log_msg("=== 33_synchrony | profile: ", profile, " | ", length(yrs), " years (", min(yrs), "-", max(yrs),
        ") | metrics: ", paste(mets, collapse = ", "))

# anomalies: cells x years per metric
A <- list()
for (mt in mets) {
  X <- vapply(mfiles, function(f) terra::values(terra::rast(f)[[mt]], mat = FALSE)[land], numeric(length(land)))
  nv <- rowSums(!is.na(X))
  med <- apply(X, 1, stats::median, na.rm = TRUE)
  a <- X - med; a[nv < sc$min_years, ] <- NA
  big <- !is.na(a) & abs(a) > sc$max_abs_anomaly      # season jumped (e.g. spring vs monsoon peak)
  a[big] <- NA
  A[[mt]] <- a
  log_msg(sprintf("  %-14s %.1f %% of cell-years dropped as jumps > %d days", mt,
                  100 * sum(big) / max(1, sum(!is.na(X))), sc$max_abs_anomaly))
  log_msg(sprintf("  %-14s cells with >= %d years: %d; anomaly SD (median cell) %.1f days", mt, sc$min_years,
                  sum(nv >= sc$min_years), stats::median(apply(a, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)))
}
hm <- land_values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), grid)[, 1]
domains <- list(herb = !is.na(hm) & hm > 0, all = rep(TRUE, length(land)))

# share of anomaly variance explained by region-year means
sync_r2 <- function(a, g) {
  ok <- !is.na(g); a <- a[ok, , drop = FALSE]; g <- g[ok]
  tot <- sum(a^2, na.rm = TRUE); if (!is.finite(tot) || tot == 0) return(NA_real_)
  res <- 0
  for (j in seq_len(ncol(a))) {
    x <- a[, j]; u <- !is.na(x); if (!any(u)) next
    mu <- tapply(x[u], g[u], mean)
    res <- res + sum((x[u] - mu[as.character(g[u])])^2)
  }
  1 - res / tot
}

P <- eval_partitions(cfg, grid, ks)
rows <- list()
for (dn in names(domains)) {
  u <- domains[[dn]]
  for (mt in mets) {
    a <- A[[mt]][u, , drop = FALSE]
    rows[[length(rows) + 1]] <- data.frame(partition = "conus", k = 1L, domain = dn, metric = mt,
                                           sync_r2 = sync_r2(a, rep(1L, nrow(a))))
    for (nm in names(P)) for (kn in colnames(P[[nm]]))
      rows[[length(rows) + 1]] <- data.frame(partition = nm, k = as.integer(sub("k", "", kn)), domain = dn,
                                             metric = mt, sync_r2 = sync_r2(a, P[[nm]][u, kn]))
  }
}
R2 <- do.call(rbind, rows)
utils::write.csv(R2, file.path(out_dir, "sync_r2.csv"), row.names = FALSE)
for (k in intersect(unlist(ev$report_k), ks)) for (mt in mets) {
  d <- R2[R2$domain == "herb" & R2$metric == mt & R2$k %in% c(1, k), ]
  log_msg(sprintf("  herb k=%2d %-14s %s", k, mt, paste(sprintf("%s %.2f", d$partition, d$sync_r2), collapse = "  ")))
}

# per region, for the profiled region sets
D <- Filter(Negate(is.null), region_display(cfg, grid, cfg$region_profiles$sets))
for (d in D) {
  g <- d$num; k <- max(g, na.rm = TRUE); u <- domains[[sc$region_domain %||% "herb"]]
  reg <- data.frame()
  cc <- matrix(NA_real_, length(land), length(mets), dimnames = list(NULL, mets))
  for (mt in mets) {
    a <- A[[mt]]
    for (j in seq_len(k)) {
      i <- which(g == j & u); if (length(i) < 10) next
      aa <- a[i, , drop = FALSE]
      ser <- colMeans(aa, na.rm = TRUE)
      r <- apply(aa, 1, function(x) { ok <- !is.na(x); if (sum(ok) < sc$min_years) NA else suppressWarnings(stats::cor(x[ok], ser[ok])) })
      cc[i, mt] <- r
      reg <- rbind(reg, data.frame(region = j, metric = mt, cells = length(i),
                                   sync_r2 = sync_r2(aa, rep(1L, length(i))),
                                   median_cell_cor = stats::median(r, na.rm = TRUE)))
    }
  }
  utils::write.csv(reg, file.path(out_dir, sprintf("sync_regions_%s.csv", d$id)), row.names = FALSE)
  write_feature_matrix({ m <- matrix(NA_real_, terra::ncell(grid), length(mets)); m[land, ] <- cc; m },
                       grid, paste0("cor_", mets), file.path(out_dir, sprintf("cellcorr_%s.tif", d$id)))
  log_msg(d$id, ": per-region synchrony written (", k, " regions)")
}
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s",
                as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
