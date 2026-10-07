#!/usr/bin/env Rscript
# Step 07 (Cut B regions): cluster the across-year median seasonal curves from 06 and compare
# the result with the Cut A partitions. Config `cutb_regions:` and `clustering:` (stability
# settings, seed, check points).
# Features: summaries/cutb_ndvi_curve_median.tif and cutb_ndii_curve_median.tif (46 values
#   each). With equally spaced, already smoothed curves, functional PCA reduces to PCA of the
#   discretized curves: PCA per index (keep pca_var), each index scaled to total variance 1.
# Variants:
#   raw    curve level and shape (productivity and timing)
#   shape  each cell's curve minus its mean, divided by max(amplitude, min_amp_scale):
#          timing / shape only, so a sparse and a dense grassland with the same calendar match
# Methods: spatial-block stability over k (as 20_), k-means for every k, two-stage Ward.
# Comparison: ARI of every Cut B k against every Cut A k (all cells and the strict herbaceous
#   mask), best-matching A k per B k; check-site regions; cluster-mean curves.
# Outputs in <run_dir>/clusters/b_<variant>/: b_kmeans.tif (k04..), b_ward.tif, b_features.tif,
#   stability.csv, kmeans_fit.csv, ari_vs_A.csv, check_points.csv, curves_k<K>.csv, figs/*.png
#
# Usage: Rscript scripts/07_cutb_regions.R [profile] [variant[,variant]]
#   RStudio: gsi_profile <- "dev"; source("scripts/07_cutb_regions.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
variant_arg <- if (length(args) >= 2) args[2] else if (exists("gsi_variant", envir = globalenv())) {
  base::get("gsi_variant", envir = globalenv())
} else NULL

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))

cfg <- load_cfg(profile)
log_file <- log_init("07_cutb_regions", profile)
notify_init(cfg)
cb <- cfg$cutb_regions; cc <- cfg$clustering
variants <- if (is.null(variant_arg)) unlist(cb$default_variants) else strsplit(variant_arg, ",")[[1]]
ks <- seq(cb$k_range[[1]], cb$k_range[[2]])
terra::terraOptions(progress = 0)
t_all <- Sys.time()
log_msg("=== 07_cutb_regions | profile: ", profile, " | variants: ", paste(variants, collapse = ", "),
        " | k ", min(ks), "-", max(ks))

grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land_all <- which(!is.na(terra::values(grid, mat = FALSE)))
cf <- function(ix) summary_file(cfg, sprintf("cutb_%s_curve", ix), "median")
for (ix in c("ndvi", "ndii")) if (!file.exists(cf(ix))) stop("missing ", cf(ix), " (run 06_cutb_curves.R)")
CV <- list(ndvi = land_values(terra::rast(cf("ndvi")), grid), ndii = land_values(terra::rast(cf("ndii")), grid))
ok <- stats::complete.cases(CV$ndvi, CV$ndii)
log_msg(length(land_all), " land cells; ", sum(ok), " with complete curves (", sum(!ok), " left unlabeled)")
xy <- terra::xyFromCell(grid, land_all)
doy <- as.integer(sub("doy", "", colnames(CV$ndvi)))
hm_file <- file.path(static_dir(cfg), "herb_mask.tif")
hmask <- if (file.exists(hm_file)) { v <- land_values(terra::rast(hm_file), grid)[ok, 1]; !is.na(v) & v > 0 } else rep(FALSE, sum(ok))
log_msg("strict herbaceous mask: ", sum(hmask), " of ", sum(ok), " cells")
full_ok <- function(L) { M <- matrix(NA_integer_, length(land_all), ncol(L), dimnames = list(NULL, colnames(L))); M[ok, ] <- L; M }
lab_rast <- function(v) { r <- terra::rast(grid); m <- rep(NA_real_, terra::ncell(grid)); m[land_all[ok]] <- v; terra::values(r) <- m; r }
states <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL
# Cut A partitions to compare against
A <- list()
for (nm in names(cb$compare_to)) {
  f <- file.path(clusters_dir(cfg, cb$compare_to[[nm]]), "a1_kmeans.tif")
  if (file.exists(f)) A[[nm]] <- land_values(terra::rast(f), grid)[ok, , drop = FALSE] else log_warn("no ", f, "; skipping ", nm)
}
cp <- cc$check_points
cp_df <- if (length(cp)) data.frame(name = vapply(cp, `[[`, "", "name"), lon = vapply(cp, function(p) as.numeric(p$lon), 0),
                                    lat = vapply(cp, function(p) as.numeric(p$lat), 0)) else NULL
if (!is.null(cp_df)) {
  cp_df$row <- match(terra::cellFromXY(grid, as.matrix(cp_df[, c("lon", "lat")])), land_all[ok])
  cp_df <- cp_df[!is.na(cp_df$row), , drop = FALSE]
}

run_variant <- function(vname) {
  v <- cb$variants[[vname]]
  out_dir <- clusters_dir(cfg, paste0("b_", vname)); fig_dir <- file.path(out_dir, "figs")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  log_msg("--- variant ", vname, " -> ", out_dir)
  X <- lapply(CV, function(m) m[ok, , drop = FALSE])
  if (isTRUE(v$shape)) {
    X <- lapply(X, function(m) {
      amp <- pmax(apply(m, 1, max) - apply(m, 1, min), v$min_amp_scale)
      (m - rowMeans(m)) / amp
    })
  }
  parts <- lapply(names(X), function(ix) {
    r <- block_reduce(X[[ix]], 1, cb$pca_var)
    colnames(r$scores) <- paste0(ix, "_PC", seq_len(ncol(r$scores)))
    log_msg(sprintf("  %s: 46 values -> %d PCs (%.1f %% variance; PC1 %.1f %%)", ix, r$keep,
                    100 * sum(r$eigen[seq_len(r$keep)]) / sum(r$eigen), 100 * r$eigen[1] / sum(r$eigen)))
    r
  })
  Z <- do.call(cbind, lapply(parts, `[[`, "scores"))
  M <- matrix(NA_real_, terra::ncell(grid), ncol(Z)); M[land_all[ok], ] <- Z
  write_feature_matrix(M, grid, colnames(Z), file.path(out_dir, "b_features.tif"))

  # stability and k-means
  t1 <- Sys.time()
  stab <- stability_k(Z, xy[ok, ], ks, cc$stability, cc$kmeans_nstart, cc$seed)
  utils::write.csv(stab, file.path(out_dir, "stability.csv"), row.names = FALSE)
  stab_m <- stats::aggregate(ari ~ k, stab, mean)
  log_msg(sprintf("  stability (%.1f min): %s", as.numeric(difftime(Sys.time(), t1, units = "mins")),
                  paste(sprintf("k%d=%.2f", stab_m$k, stab_m$ari), collapse = " ")))
  set.seed(cc$seed)
  KM <- matrix(NA_integer_, nrow(Z), length(ks), dimnames = list(NULL, sprintf("k%02d", ks)))
  cols <- list(); tss <- wtss(Z); fit <- data.frame()
  for (k in ks) {
    km <- kmeans_fit(Z, k, cc$kmeans_nstart)
    oc <- order_and_colors(km$centers)
    KM[, sprintf("k%02d", k)] <- match(seq_len(k), oc$order)[km$cluster]
    cols[[as.character(k)]] <- oc$col[oc$order]
    fit <- rbind(fit, data.frame(k = k, between_ss_share = 1 - km$tot.withinss / tss, min_size = min(km$size)))
  }
  utils::write.csv(merge(fit, stab_m), file.path(out_dir, "kmeans_fit.csv"), row.names = FALSE)
  write_labels(full_ok(KM), grid, file.path(out_dir, "b_kmeans.tif"))
  # two-stage Ward (nested)
  set.seed(cc$seed + 2)
  micro <- stats::kmeans(Z, cc$ward$n_micro, iter.max = 50, algorithm = "MacQueen")
  hc <- ward_two_stage(micro$centers, micro$size)
  WD <- sapply(ks, function(k) stats::cutree(hc, k)[micro$cluster]); colnames(WD) <- sprintf("cut%02d", ks)
  write_labels(full_ok(WD), grid, file.path(out_dir, "b_ward.tif"))
  ward_ari <- vapply(ks, function(k) ari(KM[, sprintf("k%02d", k)], WD[, sprintf("cut%02d", k)]), 0)
  log_msg("  Ward vs k-means ARI: ", paste(sprintf("k%d=%.2f", ks, ward_ari), collapse = " "))

  # comparison with Cut A: ARI for every (B k, A k), all cells and the strict herb mask
  cmp <- data.frame()
  for (nm in names(A)) for (kb in ks) for (ka in colnames(A[[nm]])) {
    a <- A[[nm]][, ka]; b <- KM[, sprintf("k%02d", kb)]; use <- !is.na(a)
    cmp <- rbind(cmp, data.frame(cut_a = nm, k_a = as.integer(sub("k", "", ka)), k_b = kb,
                                 ari_all = ari(a[use], b[use]), ari_herb = ari(a[use & hmask], b[use & hmask])))
  }
  if (nrow(cmp)) {
    utils::write.csv(cmp, file.path(out_dir, "ari_vs_A.csv"), row.names = FALSE)
    for (nm in names(A)) {
      d <- cmp[cmp$cut_a == nm & cmp$k_a == cmp$k_b, ]
      log_msg(sprintf("  ARI B vs %s, same k (all / herb mask): %s", nm,
                      paste(sprintf("k%d=%.2f/%.2f", d$k_b, d$ari_all, d$ari_herb), collapse = " ")))
      for (kk in intersect(c(7, 10, 13), ks)) {
        d2 <- cmp[cmp$cut_a == nm & cmp$k_b == kk, ]; b2 <- d2[which.max(d2$ari_all), ]
        log_msg(sprintf("    B k=%d best matches %s k=%d (ARI %.2f; herb mask %.2f)", kk, nm, b2$k_a, b2$ari_all, b2$ari_herb))
      }
    }
  }
  # check sites
  if (!is.null(cp_df) && nrow(cp_df)) {
    cpt <- data.frame(site = cp_df$name, KM[cp_df$row, , drop = FALSE], check.names = FALSE)
    utils::write.csv(cpt, file.path(out_dir, "check_points.csv"), row.names = FALSE)
    ksh <- intersect(c(4, 6, 7, 8, 10, 13, 16), ks)
    for (i in seq_len(nrow(cpt))) log_msg(sprintf("  %-46s %s", cpt$site[i],
      paste(sprintf("k%02d=%2d", ksh, unlist(cpt[i, sprintf("k%02d", ksh)])), collapse = " ")))
  }

  # figures: stability, maps, cluster-mean curves
  png(file.path(fig_dir, "b_stability.png"), 1500, 700, res = 150)
  par(mar = c(4.5, 4.5, 2.5, 1), las = 1)
  plot(stab_m$k, stab_m$ari, type = "b", pch = 19, col = "#2a78d6", ylim = c(0.4, 1), xlab = "k",
       ylab = "stability (ARI between half-samples)", main = sprintf("Cut B (%s) k-means stability", vname))
  dev.off()
  for (k in intersect(unlist(cb$curve_fig_k), ks)) {
    kn <- sprintf("k%02d", k); kc <- cols[[as.character(k)]]; lab <- KM[, kn]
    png(file.path(fig_dir, sprintf("b_kmeans_%s.png", kn)), 2000, 1300, res = 150)
    terra::plot(lab_rast(lab), col = kc, breaks = seq(0.5, k + 0.5), legend = FALSE, axes = FALSE,
                mar = c(2.5, 0.5, 2, 0.5), main = sprintf("Cut B (%s) k-means, k = %d", vname, k))
    if (!is.null(states)) lines(states$x, states$y, col = "#ffffffaa", lwd = 0.5)
    legend("bottom", legend = seq_len(k), fill = kc, ncol = min(k, 13), bty = "n", cex = 0.8, border = NA,
           x.intersp = 0.4, inset = c(0, -0.06), xpd = NA)
    if (!is.null(cp_df)) { points(cp_df$lon, cp_df$lat, pch = 21, bg = "white", cex = 0.9)
      text(cp_df$lon, cp_df$lat, lab[cp_df$row], pos = 3, cex = 0.65, font = 2) }
    dev.off()
    # cluster-mean raw curves (always in index units, also for the shape variant)
    mc <- do.call(rbind, lapply(seq_len(k), function(j) {
      i <- lab == j
      data.frame(cluster = j, n_cells = sum(i), herb_share = mean(hmask[i]), doy = doy,
                 ndvi = colMeans(CV$ndvi[ok, , drop = FALSE][i, , drop = FALSE]),
                 ndii = colMeans(CV$ndii[ok, , drop = FALSE][i, , drop = FALSE]))
    }))
    utils::write.csv(mc, file.path(out_dir, sprintf("curves_%s.csv", kn)), row.names = FALSE)
    nr <- ceiling(k / 5)
    png(file.path(fig_dir, sprintf("b_curves_%s.png", kn)), 2200, 420 * nr + 120, res = 150)
    par(mfrow = c(nr, 5), mar = c(2.2, 2.5, 2, 0.5), oma = c(0, 0, 2, 0))
    for (j in seq_len(k)) {
      d <- mc[mc$cluster == j, ]
      plot(d$doy, d$ndvi, type = "l", lwd = 2.5, col = kc[j], ylim = c(-0.2, 0.85), xaxt = "n", xlab = "", ylab = "",
           main = sprintf("%d  (%.0f%% cells, herb %.0f%%)", j, 100 * d$n_cells[1] / nrow(Z), 100 * d$herb_share[1]), cex.main = 0.9)
      lines(d$doy, d$ndii, lwd = 2, lty = 2, col = kc[j])
      axis(1, at = c(1, 91, 182, 274), labels = c("J", "A", "J", "O"))
    }
    mtext(sprintf("Cut B (%s), k = %d: cluster-mean NDVI (solid) and NDII7 (dashed)", vname, k), outer = TRUE, cex = 0.9)
    dev.off()
  }
  if (nrow(cmp)) {   # ARI heat maps vs each Cut A
    png(file.path(fig_dir, "b_ari_vs_A.png"), 900 * length(A), 800, res = 150)
    par(mfrow = c(1, length(A)), mar = c(4.5, 4.5, 2.5, 1))
    for (nm in names(A)) {
      d <- cmp[cmp$cut_a == nm, ]; ka <- sort(unique(d$k_a))
      Mx <- matrix(NA, length(ks), length(ka)); Mx[cbind(match(d$k_b, ks), match(d$k_a, ka))] <- d$ari_all
      image(ks, ka, Mx, col = grDevices::hcl.colors(20, "YlGnBu", rev = TRUE), zlim = c(0, max(cmp$ari_all)),
            xlab = "Cut B k", ylab = paste(nm, "k"), main = paste("ARI: Cut B vs", nm))
      abline(0, 1, col = "#ffffff", lty = 2)
    }
    dev.off()
  }
  saveRDS(list(variant = vname, rotation = lapply(parts, `[[`, "rotation"), hc = hc, ok = ok, doy = doy),
          file.path(out_dir, "b_run.rds"))
  sprintf("%s: %.1f min", vname, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

msgs <- vapply(variants, function(vn) tryCatch(run_variant(vn), error = function(e) {
  log_err("variant ", vn, " failed: ", conditionMessage(e)); paste(vn, "FAILED") }), "")
msg <- sprintf("%s | total %.1f min", paste(msgs, collapse = "; "), as.numeric(difftime(Sys.time(), t_all, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("07_cutb_regions finished", msg, priority = 3, tags = "white_check_mark")
