#!/usr/bin/env Rscript
# Step 20 (Cut A, run A1): cluster CONUS land cells on Block 1 across-year medians +
# static (elevation, latitude). Config `clustering:`.
#
# Features: the Block 1 medians in clustering: a1: block1 (freeze threshold 0 C only; ffs
#   dropped as fff - lsf - 1; t_ann / t_range dropped as linear in the seasonal means),
#   log10 of p_ann / ai_ann / ai_warm, standardized, PCA (keep pca_var). Block 1 has total
#   variance 1; static (elev, lat) is scaled to static_weight. A k-means run without static
#   is the sensitivity comparison.
# Variants (config a1: variants; output clusters/a1 for `base`, clusters/a1_<name> else):
#   groups       Block 1 split into feature groups (temperature / moisture amount / moisture
#                seasonality) that share the block's variance equally, so temperature no
#                longer dominates by feature count. null = one group (the first run).
#   cell_weight  weight cells by a herb_share.tif layer (floor + (1 - floor) * share) in
#                every fit: weighted k-means, weight-proportional samples for stability,
#                GMM and the Ward micro-clusters. Every cell is still labeled.
# Methods:
#   stability   spatial-block half-samples, k-means on each, ARI on a common evaluation set,
#               for every k in k_range -> picks k_detail (or use config k_detail)
#   kmeans      all cells, every k in k_range (one band per k, ids ordered by centroid PC1)
#   gmm         mclust (one model) fit on a subsample at each k_detail, all cells
#               predicted: class + max membership probability
#   ward        k-means micro-clusters (n_micro), Ward on their centroids; cuts at
#               ward: cuts + k_detail
# Outputs in <run_dir>/clusters/<a1 | a1_variant>/:
#   a1_features.tif (reduced features), a1_kmeans.tif (k04..), a1_kmeans_nostatic.tif,
#   a1_gmm_k<K>.tif (class, maxprob), a1_ward.tif (cut<K> bands), stability.csv,
#   kmeans_fit.csv, agreement.csv (incl. ARI vs the base run), pca.csv, profiles_k<K>.csv,
#   check_points.csv (region of each check site per k), figs/*.png, a1_run.rds
# Not resume-safe per piece (each variant takes ~20-30 min on CONUS); rerun overwrites.
#
# Usage: Rscript scripts/20_cluster_a1.R [profile] [variant[,variant...]]
#   RStudio: gsi_profile <- "dev"; gsi_variant <- "groups"; source("scripts/20_cluster_a1.R")
#   Variants default to config a1: default_variants.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
variant_arg <- if (length(args) >= 2) args[2] else if (exists("gsi_variant", envir = globalenv())) {
  base::get("gsi_variant", envir = globalenv())
} else NULL

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
variants <- if (is.null(variant_arg)) unlist(a1$default_variants) else strsplit(variant_arg, ",")[[1]]
bad <- setdiff(variants, names(a1$variants))
if (length(bad)) stop("unknown variant(s): ", paste(bad, collapse = ", "))
terra::terraOptions(progress = 0)
t_all <- Sys.time()
ks <- seq(cc$k_range[[1]], cc$k_range[[2]])
log_msg("=== 20_cluster_a1 | profile: ", profile, " | variants: ", paste(variants, collapse = ", "),
        " | k ", min(ks), "-", max(ks))

# ---- Inputs (shared by all variants) -------------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land_all <- which(!is.na(terra::values(grid, mat = FALSE)))
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
xy <- terra::xyFromCell(grid, land_all)
herb_r <- terra::rast(file.path(static_dir(cfg), "herb_share.tif"))
herb_layers <- intersect(c("herb_share", "crop_share", "open_herb_share"), names(herb_r))
HV <- land_values(herb_r, grid, herb_layers)
# Older herb_share.tif builds lack open_herb_share: rebuild it as 03_herb_mask.R does,
# (broad classes) / (1 - water), from nlcd_shares.tif.
ns_file <- file.path(static_dir(cfg), "nlcd_shares.tif")
if (!"open_herb_share" %in% colnames(HV) && file.exists(ns_file)) {
  ns <- terra::rast(ns_file)
  broad <- unlist(cfg$herb_mask$broad_classes)
  if (all(c(broad, "water") %in% names(ns))) {
    S <- land_values(ns, grid, c(broad, "water"))
    lnd <- 1 - S[, "water"]; lnd[lnd <= 0.01] <- NA
    HV <- cbind(HV, open_herb_share = rowSums(S[, broad, drop = FALSE]) / lnd)
    log_warn("herb_share.tif has no open_herb_share; rebuilt from nlcd_shares.tif (rerun 03_herb_mask.R to store it)")
  }
}
ok <- stats::complete.cases(X1, XS)
log_msg(nrow(X1), " land cells; ", sum(!ok), " dropped for missing values (left unlabeled)")
hm_file <- file.path(static_dir(cfg), "herb_mask.tif")
hmask <- if (file.exists(hm_file)) { v <- land_values(terra::rast(hm_file), grid)[ok, 1]; !is.na(v) & v > 0 } else
  rep(NA, sum(ok))
allf <- cbind(land_values(med, grid)[ok, ], XS[ok, , drop = FALSE], HV[ok, , drop = FALSE], herb_mask = hmask)
full_ok <- function(L) { M <- matrix(NA_integer_, nrow(X1), ncol(L), dimnames = list(NULL, colnames(L))); M[ok, ] <- L; M }
lab_rast <- function(v) { r <- terra::rast(grid); m <- rep(NA_real_, terra::ncell(grid)); m[land_all[ok]] <- v; terra::values(r) <- m; r }
# check sites -> row in the ok-cell matrices
cp <- cc$check_points
cp_df <- if (length(cp)) data.frame(name = vapply(cp, `[[`, "", "name"),
                                    lon = vapply(cp, function(p) as.numeric(p$lon), 0),
                                    lat = vapply(cp, function(p) as.numeric(p$lat), 0)) else NULL
if (!is.null(cp_df)) {
  cell <- terra::cellFromXY(grid, as.matrix(cp_df[, c("lon", "lat")]))
  cp_df$row <- match(cell, land_all[ok])
  cp_df <- cp_df[!is.na(cp_df$row), , drop = FALSE]
}
`%||%` <- function(a, b) if (is.null(a)) b else a
# Extra feature sources for variants with more blocks (A2): "<source>:<band>" names.
feature_sources <- list(b1med = summary_file(cfg, "block1", "median"), b1iqr = summary_file(cfg, "block1", "iqr"),
                        b2med = summary_file(cfg, "block2", "median"), b2iqr = summary_file(cfg, "block2", "iqr"),
                        b2var = summary_file(cfg, "block2", "var"))
src_values <- function(names_) {
  sp <- do.call(rbind, strsplit(names_, ":", fixed = TRUE))
  bad <- setdiff(unique(sp[, 1]), names(feature_sources))
  if (length(bad)) stop("unknown feature source(s): ", paste(bad, collapse = ", "))
  X <- do.call(cbind, lapply(unique(sp[, 1]), function(s) {
    f <- feature_sources[[s]]
    if (!file.exists(f)) stop("missing ", f, " (run the block's feature script first)")
    r <- terra::rast(f); b <- sp[sp[, 1] == s, 2]
    miss <- setdiff(b, names(r)); if (length(miss)) stop(basename(f), " lacks: ", paste(miss, collapse = ", "))
    v <- land_values(r, grid, b); colnames(v) <- paste0(s, ":", b); v
  }))
  X[, names_, drop = FALSE]
}
states <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL

run_variant <- function(vname) {
  v <- a1$variants[[vname]]
  out_name <- v$out %||% (if (vname == "base") "a1" else paste0("a1_", vname))
  out_dir <- clusters_dir(cfg, out_name)
  base_km_file <- file.path(clusters_dir(cfg, v$compare_to %||% "a1"), "a1_kmeans.tif")
  fig_dir <- file.path(out_dir, "figs")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  log_msg("--- variant ", vname, " -> ", out_dir)

  # ---- Features and weights ---------------------------------------------------------------
  rb1 <- groups_reduce(X1[ok, ], v$groups, 1, cc$pca_var, unlist(a1$log))
  rst <- block_reduce(XS[ok, , drop = FALSE], a1$static_weight, 1)   # keep all static axes
  colnames(rst$scores) <- paste0("static_PC", seq_len(ncol(rst$scores)))
  Z <- cbind(rb1$scores, rst$scores)
  Z0 <- rb1$scores                                                     # sensitivity: no static
  if (!is.null(v$block2)) {                                            # A2: Block 2 as its own block
    b2 <- v$block2
    nm2 <- unique(unlist(b2$groups))
    X2 <- src_values(nm2)[ok, , drop = FALSE]
    for (cn in intersect(unlist(b2$log1p), colnames(X2))) X2[, cn] <- log10(pmax(X2[, cn], 0) + 1)
    nna <- colSums(is.na(X2))
    for (cn in names(nna)[nna > 0]) X2[is.na(X2[, cn]), cn] <- stats::median(X2[, cn], na.rm = TRUE)
    if (any(nna > 0)) log_warn("  block 2: NA filled with column medians: ", paste(sprintf("%s (%d)", names(nna)[nna > 0], nna[nna > 0]), collapse = ", "))
    rb2 <- groups_reduce(X2, b2$groups, b2$weight %||% 1, cc$pca_var)
    colnames(rb2$scores) <- paste0("b2_", colnames(rb2$scores))
    names(rb2$parts) <- paste0("b2_", names(rb2$parts))
    Z <- cbind(rb1$scores, rb2$scores, rst$scores); Z0 <- cbind(Z0, rb2$scores)
    rb1$parts <- c(rb1$parts, rb2$parts)                               # report all groups below
  }
  if (a1$sensitivity_static_weight > 0) {
    Z0 <- cbind(Z0, block_reduce(XS[ok, , drop = FALSE], a1$sensitivity_static_weight, 1)$scores)
  }
  for (g in names(rb1$parts)) {
    p <- rb1$parts[[g]]
    log_msg(sprintf("  group %-12s %2d features -> %d PCs (%.1f %% of its variance)", g,
                    nrow(p$rotation), p$keep, 100 * sum(p$eigen[seq_len(p$keep)]) / sum(p$eigen)))
  }
  w <- NULL
  if (!is.null(v$cell_weight)) {
    lyr <- v$cell_weight$layer
    if (!lyr %in% colnames(HV)) stop("cell_weight layer '", lyr, "' not in herb_share.tif (", paste(colnames(HV), collapse = ", "), ")")
    w <- cell_weights(HV[ok, lyr], v$cell_weight$floor)
    w[!is.finite(w)] <- v$cell_weight$floor
    log_msg(sprintf("  cell weights: %s, floor %.2f; mean %.2f; effective share of cells with weight > 0.5: %.1f %%",
                    lyr, v$cell_weight$floor, mean(w), 100 * sum(w[w > 0.5]) / sum(w)))
  }
  pad <- function(d, n) { for (j in seq_len(n)) { cn <- paste0("PC", j); if (!cn %in% names(d)) d[[cn]] <- NA }
                          d[, c("group", "feature", paste0("PC", seq_len(n)))] }
  npc <- max(vapply(rb1$parts, function(p) p$keep, 0))
  pca_tab <- do.call(rbind, lapply(names(rb1$parts), function(g) {
    p <- rb1$parts[[g]]
    rbind(data.frame(group = g, feature = rownames(p$rotation), round(p$rotation, 3), check.names = FALSE),
          data.frame(group = g, feature = "variance_share",
                     t(round(p$eigen[seq_len(p$keep)] / sum(p$eigen), 3)) |> `colnames<-`(colnames(p$rotation)),
                     check.names = FALSE)) |> pad(npc)
  }))
  utils::write.csv(pca_tab, file.path(out_dir, "pca.csv"), row.names = FALSE)   # NA = PC not kept in that group
  M <- matrix(NA_real_, terra::ncell(grid), ncol(Z)); M[land_all[ok], ] <- Z
  write_feature_matrix(M, grid, colnames(Z), file.path(out_dir, "a1_features.tif"))

  # ---- Stability over k ------------------------------------------------------------------
  t1 <- Sys.time()
  stab <- stability_k(Z, xy[ok, ], ks, cc$stability, cc$kmeans_nstart, cc$seed, w)
  utils::write.csv(stab, file.path(out_dir, "stability.csv"), row.names = FALSE)
  stab_m <- stats::aggregate(ari ~ k, stab, function(a) c(mean = mean(a), q10 = unname(stats::quantile(a, .1)),
                                                         q90 = unname(stats::quantile(a, .9))))
  stab_m <- data.frame(k = stab_m$k, stab_m$ari)
  k_detail <- if (is.null(cc$k_detail)) pick_k(stab, cc$n_detail) else unlist(cc$k_detail)
  log_msg(sprintf("  stability: %d boots x %d k in %.1f min; mean ARI %s",
                  cc$stability$n_boot, length(ks), as.numeric(difftime(Sys.time(), t1, units = "mins")),
                  paste(sprintf("k%d=%.2f", stab_m$k, stab_m$mean), collapse = " ")))
  log_msg("  k_detail: ", paste(k_detail, collapse = ", "))

  # ---- k-means, every k -------------------------------------------------------------------
  set.seed(cc$seed)
  tss <- wtss(Z, w)
  KM <- matrix(NA_integer_, nrow(Z), length(ks), dimnames = list(NULL, sprintf("k%02d", ks)))
  KM0 <- KM; cents <- list(); cols <- list(); fit <- data.frame()
  for (k in ks) {
    km <- wkmeans(Z, w, k, cc$kmeans_nstart, cc$weighted_sample)
    oc <- order_and_colors(km$centers)
    relab <- match(seq_len(k), oc$order)                 # old id -> new id
    KM[, sprintf("k%02d", k)] <- relab[km$cluster]
    cents[[as.character(k)]] <- km$centers[oc$order, , drop = FALSE]
    cols[[as.character(k)]] <- oc$col[oc$order]
    km0 <- wkmeans(Z0, w, k, cc$kmeans_nstart, cc$weighted_sample)
    KM0[, sprintf("k%02d", k)] <- match_labels(KM[, sprintf("k%02d", k)], km0$cluster)
    fit <- rbind(fit, data.frame(k = k, between_ss_share = 1 - km$tot.withinss / tss,
                                 min_size = min(km$size), ari_vs_nostatic = ari(km$cluster, km0$cluster)))
  }
  utils::write.csv(merge(fit, stab_m), file.path(out_dir, "kmeans_fit.csv"), row.names = FALSE)
  write_labels(full_ok(KM), grid, file.path(out_dir, "a1_kmeans.tif"))
  write_labels(full_ok(KM0), grid, file.path(out_dir, "a1_kmeans_nostatic.tif"))
  log_msg("  k-means: ", length(ks), " k values done")

  # ---- Gaussian mixture at k_detail -------------------------------------------------------
  gm <- list(); agree <- data.frame()
  set.seed(cc$seed + 1)
  fit_i <- sample(nrow(Z), min(cc$gmm$n_fit, nrow(Z)), prob = w)
  for (k in k_detail) {
    t2 <- Sys.time()
    m <- NULL
    for (mn in unlist(cc$gmm$model)) {   # first model in the list that fits (full covariance can be singular)
      m <- tryCatch(mclust::Mclust(Z[fit_i, ], G = k, modelNames = mn, verbose = FALSE,
                                   initialization = list(subset = sample(length(fit_i), cc$gmm$n_init_subset))),
                    error = function(e) NULL)
      if (!is.null(m)) break
    }
    if (is.null(m)) { log_warn("  GMM k=", k, " failed to fit (", paste(unlist(cc$gmm$model), collapse = ", "), ")"); next }
    cl <- integer(nrow(Z)); pmax_ <- numeric(nrow(Z))
    for (s in seq(1, nrow(Z), by = 50000)) {
      i <- s:min(s + 49999, nrow(Z))
      p <- stats::predict(m, Z[i, ])
      cl[i] <- p$classification; pmax_[i] <- apply(p$z, 1, max)
    }
    kmk <- KM[, sprintf("k%02d", k)]
    gm[[as.character(k)]] <- list(class = match_labels(kmk, cl), pmax = pmax_, bic = m$bic)
    M <- matrix(NA_real_, terra::ncell(grid), 2)
    M[land_all[ok], ] <- cbind(gm[[as.character(k)]]$class, pmax_)
    write_feature_matrix(M, grid, c("class", "maxprob"), file.path(out_dir, sprintf("a1_gmm_k%02d.tif", k)))
    agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "gmm", ari = ari(kmk, cl)))
    log_msg(sprintf("  GMM k=%d (%s): %.1f min; median max-prob %.2f; %.1f %% cells < 0.6; ARI vs k-means %.2f",
                    k, m$modelName, as.numeric(difftime(Sys.time(), t2, units = "mins")), stats::median(pmax_),
                    100 * mean(pmax_ < 0.6), ari(kmk, cl)))
  }

  # ---- Two-stage Ward --------------------------------------------------------------------
  t3 <- Sys.time()
  set.seed(cc$seed + 2)
  if (is.null(w)) {
    micro <- stats::kmeans(Z, cc$ward$n_micro, iter.max = 50, algorithm = "MacQueen")
    mc <- micro$centers; mcl <- micro$cluster; msz <- micro$size
  } else {   # micro-clusters fit on a weighted sample; sizes = summed cell weights
    idx <- sample.int(nrow(Z), min(cc$weighted_sample, nrow(Z)), prob = w)
    mc <- stats::kmeans(Z[idx, ], cc$ward$n_micro, iter.max = 50, algorithm = "MacQueen")$centers
    mcl <- assign_nearest(Z, mc)
    msz <- as.vector(rowsum(w, factor(mcl, levels = seq_len(nrow(mc)))))
    keep <- msz > 0
    if (!all(keep)) { mc <- mc[keep, , drop = FALSE]; mcl <- match(mcl, which(keep)); msz <- msz[keep] }
  }
  hc <- ward_two_stage(mc, msz)
  cuts <- sort(unique(c(unlist(cc$ward$cuts), k_detail)))
  WD <- sapply(cuts, function(k) stats::cutree(hc, k)[mcl])
  colnames(WD) <- sprintf("cut%02d", cuts)
  base_km <- if (out_name != basename(dirname(base_km_file)) && file.exists(base_km_file)) land_values(terra::rast(base_km_file), grid)[ok, , drop = FALSE] else NULL
  for (k in k_detail) {
    kn <- sprintf("k%02d", k); kmk <- KM[, kn]; wk <- WD[, sprintf("cut%02d", k)]
    WD[, sprintf("cut%02d", k)] <- match_labels(kmk, wk)
    agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "ward", ari = ari(kmk, wk)))
    if (!is.null(gm[[as.character(k)]]))
      agree <- rbind(agree, data.frame(k = k, a = "gmm", b = "ward", ari = ari(gm[[as.character(k)]]$class, wk)))
    agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "kmeans_nostatic", ari = ari(kmk, KM0[, kn])))
    if (!is.null(base_km) && kn %in% colnames(base_km))
      agree <- rbind(agree, data.frame(k = k, a = "kmeans", b = "kmeans_base_run", ari = ari(kmk, base_km[, kn])))
  }
  write_labels(full_ok(WD), grid, file.path(out_dir, "a1_ward.tif"))
  utils::write.csv(agree, file.path(out_dir, "agreement.csv"), row.names = FALSE)
  log_msg(sprintf("  Ward: %d micro-clusters, %.1f min", nrow(mc), as.numeric(difftime(Sys.time(), t3, units = "mins"))))
  if (any(agree$b == "kmeans_base_run"))
    log_msg("  ARI vs comparison run (", basename(dirname(base_km_file)), "): ", paste(sprintf("k%d=%.2f", agree$k[agree$b == "kmeans_base_run"],
                                                agree$ari[agree$b == "kmeans_base_run"]), collapse = " "))

  # ---- Cluster profiles (original units) and check sites ---------------------------------
  for (k in k_detail) {
    kmk <- KM[, sprintf("k%02d", k)]
    prof <- do.call(rbind, lapply(seq_len(k), function(j) {
      a <- allf[kmk == j, , drop = FALSE]
      data.frame(cluster = j, n_cells = nrow(a), share = nrow(a) / nrow(allf),
                 weight_share = if (is.null(w)) NA else sum(w[kmk == j]) / sum(w),
                 lon = stats::median(xy[ok, ][kmk == j, 1]), lat_c = stats::median(xy[ok, ][kmk == j, 2]),
                 t(apply(a[, colnames(a) != "herb_mask", drop = FALSE], 2, stats::median, na.rm = TRUE)),
                 herb_mask_share = mean(a[, "herb_mask"], na.rm = TRUE), color = cols[[as.character(k)]][j])
    }))
    utils::write.csv(prof, file.path(out_dir, sprintf("profiles_k%02d.csv", k)), row.names = FALSE)
  }
  if (!is.null(cp_df) && nrow(cp_df)) {
    kshow <- sort(unique(c(4, 6, 7, 8, 10, 12, k_detail)))
    kshow <- kshow[kshow %in% ks]
    cpt <- data.frame(site = cp_df$name, KM[cp_df$row, sprintf("k%02d", kshow), drop = FALSE], check.names = FALSE)
    utils::write.csv(cpt, file.path(out_dir, "check_points.csv"), row.names = FALSE)
    for (i in seq_len(nrow(cpt)))
      log_msg(sprintf("  %-46s %s", cpt$site[i], paste(sprintf("%s=%2d", names(cpt)[-1], unlist(cpt[i, -1])), collapse = " ")))
  }

  saveRDS(list(variant = vname, groups = v$groups, cell_weight = v$cell_weight,
               parts = lapply(rb1$parts, `[`, c("rotation", "center", "scale", "eigen", "keep", "factor", "log_cols")),
               rst = rst[c("center", "scale", "factor")], k_detail = k_detail, cents = cents,
               cols = cols, hc = hc, micro_centers = mc, micro_size = msz,
               gmm_bic = lapply(gm, `[[`, "bic"), ok = ok, cfg = cc),
          file.path(out_dir, "a1_run.rds"))

  # ---- Figures ---------------------------------------------------------------------------
  vlab <- if (vname == "base") "A1" else if (!is.null(v$block2)) paste0("A2 (", vname, ")") else paste0("A1 (", vname, ")")
  map_cat <- function(lbl, col, title) {
    terra::plot(lab_rast(lbl), col = col, breaks = seq(0.5, length(col) + 0.5), legend = FALSE,
                axes = FALSE, mar = c(2.5, 0.5, 2, 0.5), main = title, cex.main = 1.1)
    if (!is.null(states)) lines(states$x, states$y, col = "#ffffffaa", lwd = 0.5)
  }
  cat_legend <- function(col) {
    k <- length(col)
    legend("bottom", legend = seq_len(k), fill = col, ncol = min(k, 12), bty = "n", cex = 0.8,
           border = NA, x.intersp = 0.4, inset = c(0, -0.06), xpd = NA)
  }
  extra_cols <- function(col, n) c(col, grDevices::gray.colors(max(0, n - length(col))))
  png(file.path(fig_dir, "a1_stability.png"), 1500, 700, res = 150)
  par(mar = c(4.5, 4.5, 2.5, 4.5), las = 1)
  plot(stab_m$k, stab_m$mean, type = "n", ylim = range(c(stab_m$q10, stab_m$q90, 0.5, 1)),
       xlab = "k", ylab = "stability (ARI between half-samples)", main = paste(vlab, "k-means: stability and fit by k"))
  polygon(c(stab_m$k, rev(stab_m$k)), c(stab_m$q10, rev(stab_m$q90)), col = "#2a78d633", border = NA)
  lines(stab_m$k, stab_m$mean, lwd = 2, col = "#2a78d6"); points(stab_m$k, stab_m$mean, pch = 19, col = "#2a78d6", cex = 0.7)
  abline(v = k_detail, lty = 3, col = "#52514e")
  par(new = TRUE); plot(fit$k, fit$between_ss_share, type = "l", col = "#eb6834", lwd = 2, axes = FALSE, xlab = "", ylab = "")
  axis(4, col.axis = "#eb6834"); mtext("between-cluster share of variance", 4, 3, col = "#eb6834", las = 0)
  dev.off()
  for (k in sort(unique(c(k_detail, intersect(c(7, 10, 12), ks))))) {
    kc <- cols[[as.character(k)]]
    png(file.path(fig_dir, sprintf("a1_kmeans_k%02d.png", k)), 2000, 1300, res = 150)
    map_cat(KM[, sprintf("k%02d", k)], kc, sprintf("%s k-means, k = %d (ids ordered by centroid PC1)", vlab, k))
    cat_legend(kc)
    if (!is.null(cp_df)) {
      points(cp_df$lon, cp_df$lat, pch = 21, bg = "white", col = "black", cex = 0.9)
      text(cp_df$lon, cp_df$lat, KM[cp_df$row, sprintf("k%02d", k)], pos = 3, cex = 0.65, font = 2)
    }
    dev.off()
  }
  for (k in k_detail) {
    kc <- cols[[as.character(k)]]
    png(file.path(fig_dir, sprintf("a1_methods_k%02d.png", k)), 2400, 1500, res = 150)
    par(mfrow = c(2, 2))
    map_cat(KM[, sprintf("k%02d", k)], kc, sprintf("k-means, k = %d", k))
    if (!is.null(gm[[as.character(k)]])) {
      map_cat(gm[[as.character(k)]]$class, extra_cols(kc, max(gm[[as.character(k)]]$class)),
              sprintf("Gaussian mixture (ARI vs k-means %.2f)", agree$ari[agree$k == k & agree$b == "gmm"]))
    } else plot.new()
    map_cat(WD[, sprintf("cut%02d", k)], extra_cols(kc, max(WD[, sprintf("cut%02d", k)])),
            sprintf("Ward on micro-clusters (ARI vs k-means %.2f)", agree$ari[agree$k == k & agree$b == "ward"]))
    map_cat(KM0[, sprintf("k%02d", k)], extra_cols(kc, max(KM0[, sprintf("k%02d", k)])),
            sprintf("k-means without static (ARI %.2f)", agree$ari[agree$k == k & agree$b == "kmeans_nostatic"]))
    dev.off()
    if (!is.null(gm[[as.character(k)]])) {
      png(file.path(fig_dir, sprintf("a1_gmm_maxprob_k%02d.png", k)), 2000, 1250, res = 150)
      br <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0001)
      terra::plot(lab_rast(gm[[as.character(k)]]$pmax), breaks = br, col = grDevices::hcl.colors(7, "Inferno"),
                  axes = FALSE, mar = c(0.5, 0.5, 2, 5),
                  main = sprintf("Gaussian mixture k = %d: max membership probability (low = transition zone)", k))
      if (!is.null(states)) lines(states$x, states$y, col = "#ffffff88", lwd = 0.5)
      dev.off()
    }
  }
  png(file.path(fig_dir, "a1_ward_cuts.png"), 2400, 1500, res = 150)
  par(mfrow = c(2, 2))
  for (k in intersect(c(2, 4, 6, 8), cuts)) {
    wk <- WD[, sprintf("cut%02d", k)]
    cc_ <- order_and_colors(do.call(rbind, lapply(seq_len(k), function(j) colMeans(Z[wk == j, , drop = FALSE]))))
    map_cat(match(wk, cc_$order), cc_$col[cc_$order], sprintf("Ward tree cut at %d", k))
  }
  dev.off()
  if (!is.null(w)) {
    png(file.path(fig_dir, "a1_cell_weight.png"), 2000, 1250, res = 150)
    terra::plot(lab_rast(w), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0001), col = grDevices::hcl.colors(5, "YlGn", rev = TRUE),
                axes = FALSE, mar = c(0.5, 0.5, 2, 5), main = sprintf("Cell weight: %.1f + %.1f x %s",
                v$cell_weight$floor, 1 - v$cell_weight$floor, v$cell_weight$layer))
    if (!is.null(states)) lines(states$x, states$y, col = "#55555588", lwd = 0.5)
    dev.off()
  }
  msg <- sprintf("%s: k_detail %s; %.1f min", vname, paste(k_detail, collapse = ","),
                 as.numeric(difftime(Sys.time(), t0, units = "mins")))
  log_msg("  done | ", msg)
  msg
}

msgs <- vapply(variants, function(vn) tryCatch(run_variant(vn), error = function(e) {
  log_err("variant ", vn, " failed: ", conditionMessage(e)); paste(vn, "FAILED")
}), "")
msg <- sprintf("%s | total %.1f min", paste(msgs, collapse = "; "),
               as.numeric(difftime(Sys.time(), t_all, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("20_cluster_a1 finished", msg, priority = 3, tags = "white_check_mark")
