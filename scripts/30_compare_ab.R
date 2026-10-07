#!/usr/bin/env Rscript
# Step 30 (evaluation, Cut A vs Cut B): how well does each partition work as a set of domains
# for the observed (satellite) seasonal curves? For every partition in config `evaluation:
# partitions` (Cut A and Cut B k-means rasters) plus a geographic baseline (k-means on location
# only) and every k:
#   r2_curves   share of the Cut B curve variance explained by the regions: between-region
#               sum of squares / total, on the Cut B PCA features (NDVI and NDII7 curves, each
#               index weighted equally; from 07_cutb_regions.R)
#   r2_<metric> the same for single phenology metrics (across-year medians from 06): green-up,
#               peak, end of season, season length, curing, amplitudes
# Each for all cells with curves and for the strict herbaceous mask. Cut B partitions are the
# reference: they are fitted to these curves, so they show how much a partition of that size
# can explain at most; Cut A is independent of the curves. The geographic baseline shows how
# much any compact regions of the same number explain.
# Outputs <run_dir>/eval/ab/: r2_curves.csv, r2_metrics.csv, figs/r2_curves.png,
#   figs/r2_metrics.png. Runs in a few minutes (one process).
#
# Usage: Rscript scripts/30_compare_ab.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/30_compare_ab.R")

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

cfg <- load_cfg(profile)
log_file <- log_init("30_compare_ab", profile)
ev <- cfg$evaluation
ks <- seq(ev$k[[1]], ev$k[[2]])
out_dir <- file.path(cfg$run_dir, "eval", "ab"); fig_dir <- file.path(out_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
log_msg("=== 30_compare_ab | profile: ", profile, " | k ", min(ks), "-", max(ks))

grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
ff <- file.path(cfg$run_dir, "clusters", ev$curve_features)
if (!file.exists(ff)) stop("missing ", ff, " (run 07_cutb_regions.R)")
Y <- land_values(terra::rast(ff), grid)
mf <- summary_file(cfg, "cutb_metrics", "median")
M <- land_values(terra::rast(mf), grid, intersect(unlist(ev$metrics), names(terra::rast(mf))))
hm <- land_values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), grid)[, 1]
domains <- list(all = stats::complete.cases(Y), herb = stats::complete.cases(Y) & !is.na(hm) & hm > 0)
log_msg("cells with curves: ", sum(domains$all), "; strict herbaceous mask: ", sum(domains$herb))
P <- eval_partitions(cfg, grid, ks)
log_msg("partitions: ", paste(names(P), collapse = ", "))

rc <- data.frame(); rm_ <- data.frame()
for (nm in names(P)) for (kn in colnames(P[[nm]])) {
  k <- as.integer(sub("k", "", kn)); g <- P[[nm]][, kn]
  for (dn in names(domains)) {
    u <- domains[[dn]]
    rc <- rbind(rc, data.frame(partition = nm, k = k, domain = dn, r2 = r2_groups(Y[u, , drop = FALSE], g[u])))
    for (mn in colnames(M)) rm_ <- rbind(rm_, data.frame(partition = nm, k = k, domain = dn, metric = mn,
                                                         r2 = r2_groups(M[u, mn, drop = FALSE], g[u])))
  }
}
utils::write.csv(rc, file.path(out_dir, "r2_curves.csv"), row.names = FALSE)
utils::write.csv(rm_, file.path(out_dir, "r2_metrics.csv"), row.names = FALSE)
for (k in intersect(unlist(ev$report_k), ks)) {
  d <- rc[rc$k == k, ]
  log_msg(sprintf("  k=%2d  curve R^2 (all / herb mask): %s", k,
                  paste(sprintf("%s %.2f/%.2f", unique(d$partition),
                                d$r2[d$domain == "all"], d$r2[d$domain == "herb"]), collapse = "; ")))
}

# figures
pal <- c("#2a78d6", "#1a7f37", "#eb6834", "#c0392b", "#7b4fb3", "#888888")
pn <- names(P); pcol <- setNames(pal[seq_along(pn)], pn); pcol["geo"] <- "#888888"
plty <- setNames(ifelse(grepl("^B_", pn), 2, 1), pn); plty["geo"] <- 3
png(file.path(fig_dir, "r2_curves.png"), 1800, 750, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), las = 1)
for (dn in names(domains)) {
  d <- rc[rc$domain == dn, ]
  plot(NA, xlim = range(ks), ylim = c(0, 1), xlab = "number of regions (k)", ylab = "share of curve variance explained",
       main = if (dn == "all") "Satellite curves, all cells" else "Satellite curves, strict herbaceous mask")
  for (nm in pn) { e <- d[d$partition == nm, ]; lines(e$k, e$r2, col = pcol[nm], lty = plty[nm], lwd = 2) }
  if (dn == "all") legend("bottomright", pn, col = pcol, lty = plty, lwd = 2, bty = "n", cex = 0.8)
}
dev.off()
mets <- colnames(M); nc <- min(4, length(mets)); nr <- ceiling(length(mets) / nc)
png(file.path(fig_dir, "r2_metrics.png"), 520 * nc, 430 * nr + 60, res = 150)
par(mfrow = c(nr, nc), mar = c(4, 4, 2.2, 0.5), las = 1)
for (mn in mets) {
  d <- rm_[rm_$metric == mn & rm_$domain == "herb", ]
  plot(NA, xlim = range(ks), ylim = c(0, 1), xlab = "k", ylab = "R^2", main = paste(mn, "(herb mask)"), cex.main = 0.9)
  for (nm in pn) { e <- d[d$partition == nm, ]; lines(e$k, e$r2, col = pcol[nm], lty = plty[nm], lwd = 2) }
}
plot.new(); legend("center", pn, col = pcol, lty = plty, lwd = 2, bty = "n")
dev.off()

msg <- sprintf("%d partitions x %d k, %.1f min", length(P), length(ks), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | outputs: ", out_dir, " | log: ", log_file)
