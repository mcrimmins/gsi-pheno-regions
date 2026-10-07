#!/usr/bin/env Rscript
# Step 20 (Cut A, run A1): cluster CONUS land cells on Block 1 across-year medians +
# static (elevation, latitude). Config `clustering:`.
#
# Features: the Block 1 medians listed in clustering: a1: block1 (one freeze threshold,
#   0 C; ffs is dropped as fff - lsf - 1; t_ann / t_range dropped as linear in the seasonal
#   means), log10 of p_ann / ai_ann / ai_warm. Standardize, PCA within the block (keep
#   pca_var), scale the block to total variance 1. Static (elev, lat) standardized and scaled
#   to static_weight. A sensitivity run without static (sensitivity_static_weight).
# Methods:
#   stability   spatial-block half-samples, k-means on each, ARI on a common evaluation set,
#               for every k in k_range -> picks k_detail (or use config k_detail)
#   kmeans      all cells, every k in k_range (one band per k, ids ordered by centroid PC1)
#   gmm         mclust (one model) fit on a subsample at each k_detail, all cells
#               predicted: class + max membership probability
#   ward        k-means micro-clusters (n_micro), Ward on their centroids; cuts at
#               ward: cuts + k_detail
# Outputs in <run_dir>/clusters/a1/:
#   a1_features.tif (reduced features), a1_kmeans.tif (k04..), a1_kmeans_nostatic.tif,
#   a1_gmm_k<K>.tif (class, maxprob), a1_ward.tif (cut<K> bands), stability.csv,
#   kmeans_fit.csv, agreement.csv, pca.csv, profiles_k<K>.csv, figs/*.png, a1_run.rds
# Not resume-safe per piece (whole run takes minutes); rerun overwrites.
#
# Usage: Rscript scripts/20_cluster_a1.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/20_cluster_a1.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

suppressPackageStartupMessages(library(mclust))   # Mclust needs it attached
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))

cfg <- load_cfg(profile)
log_file <- log_init("20_cluster_a1", profile)
notify_init(cfg)
cc <- cfg$clustering; a1 <- cc$a1
out_dir <- clusters_dir(cfg, "a1"); fig_dir <- file.path(out_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
ks <- seq(cc$k_range[[1]], cc$k_range[[2]])
log_msg("=== 20_cluster_a1 | profile: ", profile, " | k ", min(ks), "-", max(ks))

# ---- Features ----------------------------------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
med <- terra::rast(summary_file(cfg, "block1", "median"))
b1_names <- unlist(a1$block1)
miss <- setdiff(b1_names, names(med))
if (length(miss)) stop("Block 1 median lacks: ", paste(miss, collapse = ", "))
X1 <- land_values(med, grid, b1_names)
lat_file <- file.path(static_dir(cfg), "lat.tif")
lat_r <- if (file.exists(lat_file)) terra::rast(lat_file) else terra::init(grid, "y")
static_r <- c(terra::rast(file.path(static_dir(cfg), "elev.tif")), lat_r)
names(static_r) <- c("elev", "lat")
XS <- land_values(static_r, grid, unlist(a1$static))
xy <- terra::xyFromCell(grid, which(!is.na(terra::values(grid, mat = FALSE))))
ok <- stats::complete.cases(X1, XS)
log_msg(nrow(X1), " land cells; ", sum(!ok), " dropped for missing values (left unlabeled)")

rb1 <- block_reduce(X1[ok, ], 1, cc$pca_var, unlist(a1$log))
rst <- block_reduce(XS[ok, , drop = FALSE], a1$static_weight, 1)   # keep all static axes
Z <- cbind(rb1$scores, rst$scores)
colnames(Z) <- c(paste0("b1_PC", seq_len(rb1$keep)), paste0("st_PC", seq_len(rst$keep)))
Z0 <- rb1$scores                                                       # sensitivity: no static
if (a1$sensitivity_static_weight > 0) {
  Z0 <- cbind(rb1$scores, block_reduce(XS[ok, , drop = FALSE], a1$sensitivity_static_weight, 1)$scores)
}
log_msg(sprintf("Block 1: %d features -> %d PCs (%.1f %% variance); static: %d axes, weight %.2f",
                ncol(X1), rb1$keep, 100 * sum(rb1$eigen[1:rb1$keep]) / sum(rb1$eigen),
                rst$keep, a1$static_weight))
pca_tab <- data.frame(feature = rownames(rb1$rotation), round(rb1$rotation, 3))
utils::write.csv(rbind(pca_tab, data.frame(feature = "variance_share",
  t(round(rb1$eigen[1:rb1$keep] / sum(rb1$eigen), 3)) |> `colnames<-`(colnames(rb1$rotation)))),
  file.path(out_dir, "pca.csv"), row.names = FALSE)
full_ok <- function(L) { M <- matrix(NA_integer_, nrow(X1), ncol(L), dimnames = list(NULL, colnames(L))); M[ok, ] <- L; M }
write_feature_matrix({ M <- matrix(NA_real_, terra::ncell(grid), ncol(Z))
  M[which(!is.na(terra::values(grid, mat = FALSE)))[ok], ] <- Z; M },
  grid, colnames(Z), file.path(out_dir, "a1_features.tif"))

# ---- Stability over k --------------------------------------------------------------------
t1 <- Sys.time()
stab <- stability_k(Z, xy[ok, ], ks, cc$stability, cc$kmeans_nstart, cc$seed)
utils::write.csv(stab, file.path(out_dir, "stability.csv"), row.names = FALSE)
stab_m <- stats::aggregate(ari ~ k, stab, function(v) c(mean = mean(v), q10 = unname(stats::quantile(v, .1)),
                                                        q90 = unname(stats::quantile(v, .9))))
stab_m <- data.frame(k = stab_m$k, stab_m$ari)
k_detail <- if (is.null(cc$k_detail)) pick_k(stab, cc$n_detail) else unlist(cc$k_detail)
log_msg(sprintf("stability: %d boots x %d k in %.1f min; mean ARI %s",
                cc$stability$n_boot, length(ks), as.numeric(difftime(Sys.time(), t1, units = "mins")),
                paste(sprintf("k%d=%.2f", stab_m$k, stab_m$mean), collapse = " ")))
log_msg("k_detail: ", paste(k_detail, collapse = ", "))

# ---- k-means, every k ---------------------------------------------------------------------
set.seed(cc$seed)
tss <- sum(scale(Z, scale = FALSE)^2)
KM <- matrix(NA_integer_, nrow(Z), length(ks), dimnames = list(NULL, sprintf("k%02d", ks)))
KM0 <- KM; cents <- list(); cols <- list(); fit <- data.frame()
for (k in ks) {
  km <- kmeans_fit(Z, k, cc$kmeans_nstart)
  oc <- order_and_colors(km$centers)
  relab <- match(seq_len(k), oc$order)                 # old id -> new id
  KM[, sprintf("k%02d", k)] <- relab[km$cluster]
  cents[[as.character(k)]] <- km$centers[oc$order, , drop = FALSE]
  cols[[as.character(k)]] <- oc$col[oc$order]
  km0 <- kmeans_fit(Z0, k, cc$kmeans_nstart)
  KM0[, sprintf("k%02d", k)] <- match_labels(KM[, sprintf("k%02d", k)], km0$cluster)
  fit <- rbind(fit, data.frame(k = k, between_ss_share = 1 - km$tot.withinss / tss,
                               min_size = min(km$size), ari_vs_nostatic = ari(km$cluster, km0$cluster)))
}
utils::write.csv(merge(fit, stab_m), file.path(out_dir, "kmeans_fit.csv"), row.names = FALSE)
write_labels(full_ok(KM), grid, file.path(out_dir, "a1_kmeans.tif"))
write_labels(full_ok(KM0), grid, file.path(out_dir, "a1_kmeans_nostatic.tif"))
log_msg("k-means: ", length(ks), " k values done")

# ---- Gaussian mixture at k_detail ---------------------------------------------------------
gm <- list(); agree <- data.frame()
set.seed(cc$seed + 1)
fit_i <- sample(nrow(Z), min(cc$gmm$n_fit, nrow(Z)))
for (k in k_detail) {
  t2 <- Sys.time()
  m <- mclust::Mclust(Z[fit_i, ], G = k, modelNames = cc$gmm$model, verbose = FALSE,
                      initialization = list(subset = sample(length(fit_i), cc$gmm$n_init_subset)))
  if (is.null(m)) { log_warn("GMM k=", k, " failed to fit"); next }
  cl <- integer(nrow(Z)); pmax_ <- numeric(nrow(Z))
  for (s in seq(1, nrow(Z), by = 50000)) {
    i <- s:min(s + 49999, nrow(Z))
    p <- stats::predict(m, Z[i, ])
    cl[i] <- p$classification; pmax_[i] <- apply(p$z, 1, max)
  }
  kmk <- KM[, sprintf("k%02d", k)]
  gm[[as.character(k)]] <- list(class = match_labels(kmk, cl), pmax = pmax_, bic = m$bic)
  r <- terra::rast(grid, nlyrs = 2); M <- matrix(NA_real_, terra::ncell(grid), 2)
  M[which(!is.na(terra::values(grid, mat = FALSE)))[ok], ] <- cbind(gm[[as.character(k)]]$class, pmax_)
  write_feature_matrix(M, grid, c("class", "maxprob"), file.path(out_dir, sprintf("a1_gmm_k%02d.tif", k)))
  agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "gmm", ari = ari(kmk, cl)))
  log_msg(sprintf("GMM k=%d: %.1f min; median max-prob %.2f; %.1f %% cells < 0.6; ARI vs k-means %.2f",
                  k, as.numeric(difftime(Sys.time(), t2, units = "mins")), stats::median(pmax_),
                  100 * mean(pmax_ < 0.6), ari(kmk, cl)))
}

# ---- Two-stage Ward ------------------------------------------------------------------------
t3 <- Sys.time()
set.seed(cc$seed + 2)
micro <- stats::kmeans(Z, cc$ward$n_micro, iter.max = 50, algorithm = "MacQueen")
hc <- ward_two_stage(micro$centers, micro$size)
cuts <- sort(unique(c(unlist(cc$ward$cuts), k_detail)))
WD <- sapply(cuts, function(k) stats::cutree(hc, k)[micro$cluster])
colnames(WD) <- sprintf("cut%02d", cuts)
for (k in k_detail) {
  kmk <- KM[, sprintf("k%02d", k)]; w <- WD[, sprintf("cut%02d", k)]
  WD[, sprintf("cut%02d", k)] <- match_labels(kmk, w)
  agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "ward", ari = ari(kmk, w)))
  if (!is.null(gm[[as.character(k)]]))
    agree <- rbind(agree, data.frame(k = k, a = "gmm", b = "ward", ari = ari(gm[[as.character(k)]]$class, w)))
  agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "kmeans_nostatic",
                                   ari = ari(kmk, KM0[, sprintf("k%02d", k)])))
}
write_labels(full_ok(WD), grid, file.path(out_dir, "a1_ward.tif"))
utils::write.csv(agree, file.path(out_dir, "agreement.csv"), row.names = FALSE)
log_msg(sprintf("Ward: %d micro-clusters, %.1f min", cc$ward$n_micro,
                as.numeric(difftime(Sys.time(), t3, units = "mins"))))

# ---- Cluster profiles (original units) -----------------------------------------------------
herb_r <- terra::rast(file.path(static_dir(cfg), "herb_share.tif"))   # older builds lack open_herb_share
herb <- land_values(herb_r, grid, intersect(c("herb_share", "crop_share", "open_herb_share"),
                                            names(herb_r)))[ok, , drop = FALSE]
hm_file <- file.path(static_dir(cfg), "herb_mask.tif")
hmask <- if (file.exists(hm_file)) !is.na(land_values(terra::rast(hm_file), grid)[ok, 1]) &
  land_values(terra::rast(hm_file), grid)[ok, 1] > 0 else rep(NA, sum(ok))
allf <- cbind(land_values(med, grid)[ok, ], XS[ok, , drop = FALSE], herb, herb_mask = hmask)
for (k in k_detail) {
  kmk <- KM[, sprintf("k%02d", k)]
  prof <- do.call(rbind, lapply(seq_len(k), function(j) {
    v <- allf[kmk == j, , drop = FALSE]
    data.frame(cluster = j, n_cells = nrow(v), share = nrow(v) / nrow(allf),
               lon = stats::median(xy[ok, ][kmk == j, 1]), lat_c = stats::median(xy[ok, ][kmk == j, 2]),
               t(apply(v, 2, stats::median, na.rm = TRUE)),
               herb_mask_share = mean(v[, "herb_mask"], na.rm = TRUE), color = cols[[as.character(k)]][j])
  }))
  prof$herb_mask <- NULL
  utils::write.csv(prof, file.path(out_dir, sprintf("profiles_k%02d.csv", k)), row.names = FALSE)
}

saveRDS(list(rb1 = rb1[c("rotation", "center", "scale", "eigen", "keep", "factor", "log_cols")],
             rst = rst[c("center", "scale", "factor")], k_detail = k_detail, cents = cents,
             cols = cols, hc = hc, micro_centers = micro$centers, micro_size = micro$size,
             gmm_bic = lapply(gm, `[[`, "bic"), ok = ok, cfg = cc),
        file.path(out_dir, "a1_run.rds"))

# ---- Figures --------------------------------------------------------------------------------
land_all <- which(!is.na(terra::values(grid, mat = FALSE)))
lab_rast <- function(v) { r <- terra::rast(grid); m <- rep(NA_real_, terra::ncell(grid)); m[land_all[ok]] <- v; terra::values(r) <- m; r }
states <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL
asp <- 1 / cos(mean(as.vector(terra::ext(grid))[3:4]) * pi / 180)
map_cat <- function(v, col, title) {
  terra::plot(lab_rast(v), col = col, breaks = seq(0.5, length(col) + 0.5), legend = FALSE,
              axes = FALSE, mar = c(0.5, 0.5, 2, 0.5), main = title, cex.main = 1.1)
  if (!is.null(states)) lines(states$x, states$y, col = "#ffffffaa", lwd = 0.5)
}
png(file.path(fig_dir, "a1_stability.png"), 1500, 700, res = 150)
par(mar = c(4.5, 4.5, 2.5, 4.5), las = 1)
plot(stab_m$k, stab_m$mean, type = "n", ylim = range(c(stab_m$q10, stab_m$q90, 0.5, 1)),
     xlab = "k", ylab = "stability (ARI between half-samples)", main = "A1 k-means: stability and fit by k")
polygon(c(stab_m$k, rev(stab_m$k)), c(stab_m$q10, rev(stab_m$q90)), col = "#2a78d633", border = NA)
lines(stab_m$k, stab_m$mean, lwd = 2, col = "#2a78d6"); points(stab_m$k, stab_m$mean, pch = 19, col = "#2a78d6", cex = 0.7)
abline(v = k_detail, lty = 3, col = "#52514e")
par(new = TRUE); plot(fit$k, fit$between_ss_share, type = "l", col = "#eb6834", lwd = 2, axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = "#eb6834"); mtext("between-cluster share of variance", 4, 3, col = "#eb6834", las = 0)
dev.off()
for (k in k_detail) {
  png(file.path(fig_dir, sprintf("a1_kmeans_k%02d.png", k)), 2000, 1250, res = 150)
  map_cat(KM[, sprintf("k%02d", k)], cols[[as.character(k)]], sprintf("A1 k-means, k = %d (ids ordered by centroid PC1)", k))
  e <- terra::ext(grid)
  legend("bottomleft", legend = seq_len(k), fill = cols[[as.character(k)]], ncol = ceiling(k / 2),
         bty = "n", cex = 0.8, border = NA, x.intersp = 0.5)
  dev.off()
  png(file.path(fig_dir, sprintf("a1_methods_k%02d.png", k)), 2400, 1500, res = 150)
  par(mfrow = c(2, 2))
  map_cat(KM[, sprintf("k%02d", k)], cols[[as.character(k)]], sprintf("k-means, k = %d", k))
  if (!is.null(gm[[as.character(k)]])) {
    ck <- cols[[as.character(k)]]; ck <- c(ck, grDevices::gray.colors(max(0, max(gm[[as.character(k)]]$class) - k)))
    map_cat(gm[[as.character(k)]]$class, ck,
            sprintf("Gaussian mixture (ARI vs k-means %.2f)", agree$ari[agree$k == k & agree$b == "gmm"]))
  } else plot.new()
  map_cat(WD[, sprintf("cut%02d", k)], c(cols[[as.character(k)]], grDevices::gray.colors(5)),
          sprintf("Ward on micro-clusters (ARI vs k-means %.2f)", agree$ari[agree$k == k & agree$b == "ward"]))
  map_cat(KM0[, sprintf("k%02d", k)], c(cols[[as.character(k)]], grDevices::gray.colors(5)),
          sprintf("k-means without static (ARI %.2f)", agree$ari[agree$k == k & agree$b == "kmeans_nostatic"]))
  dev.off()
  if (!is.null(gm[[as.character(k)]])) {
    png(file.path(fig_dir, sprintf("a1_gmm_maxprob_k%02d.png", k)), 2000, 1250, res = 150)
    br <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0001)
    terra::plot(lab_rast(gm[[as.character(k)]]$pmax), breaks = br, col = grDevices::hcl.colors(7, "Inferno"),
                axes = FALSE, mar = c(0.5, 0.5, 2, 5), main = sprintf("Gaussian mixture k = %d: max membership probability (low = transition zone)", k))
    if (!is.null(states)) lines(states$x, states$y, col = "#ffffff88", lwd = 0.5)
    dev.off()
  }
}
png(file.path(fig_dir, "a1_ward_cuts.png"), 2400, 1500, res = 150)
par(mfrow = c(2, 2))
for (k in intersect(c(2, 4, 6, 8), cuts)) {
  w <- WD[, sprintf("cut%02d", k)]
  cc_ <- order_and_colors(do.call(rbind, lapply(seq_len(k), function(j) colMeans(Z[w == j, , drop = FALSE]))))
  map_cat(match(w, cc_$order), cc_$col[cc_$order], sprintf("Ward tree cut at %d", k))
}
dev.off()

msg <- sprintf("k_detail %s; %.1f min", paste(k_detail, collapse = ","),
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | outputs: ", out_dir, " | log: ", log_file)
notify("20_cluster_a1 finished", msg, priority = 3, tags = "white_check_mark")
