# Cut A clustering helpers: feature matrix, block PCA/weighting, stability over k,
# k-means / Gaussian mixture / two-stage Ward, label matching and agreement.
# All functions work on a cells x features matrix of land cells (rows in grid order).

clusters_dir <- function(cfg, run) file.path(cfg$run_dir, "clusters", run)

# Land-cell values of selected bands of a raster (rows = land cells of `grid`, grid order).
land_values <- function(r, grid, bands = names(r)) {
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  v <- terra::values(r[[bands]], mat = TRUE)[land, , drop = FALSE]
  colnames(v) <- bands
  v
}

# Standardize, PCA, keep `pca_var` of the variance, scale so the block's total variance is
# `weight`. Returns scores plus everything needed to describe or reapply the transform.
block_reduce <- function(X, weight = 1, pca_var = 0.95, log_cols = character()) {
  for (v in intersect(log_cols, colnames(X))) X[, v] <- log10(pmax(X[, v], 1e-3))
  ctr <- colMeans(X); sds <- apply(X, 2, stats::sd)
  Z <- sweep(sweep(X, 2, ctr), 2, sds, "/")
  if (ncol(Z) == 1) {
    S <- Z; rot <- matrix(1, 1, 1, dimnames = list(colnames(X), "PC1")); ev <- 1; keep <- 1
  } else {
    p <- stats::prcomp(Z, center = FALSE, scale. = FALSE)
    ev <- p$sdev^2
    keep <- which(cumsum(ev) / sum(ev) >= pca_var)[1]
    S <- p$x[, seq_len(keep), drop = FALSE]; rot <- p$rotation[, seq_len(keep), drop = FALSE]
  }
  f <- sqrt(weight / sum(ev[seq_len(keep)]))
  list(scores = S * f, rotation = rot, center = ctr, scale = sds, eigen = ev, keep = keep,
       factor = f, log_cols = intersect(log_cols, colnames(X)))
}

# Nearest-centroid assignment, in chunks (squared Euclidean).
assign_nearest <- function(X, centers, chunk = 50000) {
  out <- integer(nrow(X)); c2 <- rowSums(centers^2)
  for (s in seq(1, nrow(X), by = chunk)) {
    i <- s:min(s + chunk - 1, nrow(X))
    d <- -2 * X[i, , drop = FALSE] %*% t(centers) + matrix(c2, length(i), nrow(centers), byrow = TRUE)
    out[i] <- max.col(-d, ties.method = "first")
  }
  out
}

kmeans_fit <- function(X, k, nstart = 5, iter.max = 100) {
  km <- suppressWarnings(stats::kmeans(X, k, nstart = nstart, iter.max = iter.max,
                                       algorithm = "Hartigan-Wong"))
  if (km$ifault == 4) {   # Hartigan-Wong can stop on quick-transfer steps; fall back to Lloyd
    km <- stats::kmeans(X, km$centers, iter.max = 200, algorithm = "Lloyd")
  }
  km
}

# Adjusted Rand index between two labelings.
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  c2 <- function(x) sum(x * (x - 1) / 2)
  s <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / c2(n)
  (s - e) / ((sa + sb) / 2 - e)
}

# Stability over k by spatial-block subsampling: two fits on independent half-samples of
# spatial blocks, both applied to a common evaluation set, compared by ARI.
# With cell weights w, fits and the evaluation set are drawn in proportion to w.
stability_k <- function(X, xy, ks, st, nstart, seed, w = NULL) {
  set.seed(seed)
  blk <- paste(floor(xy[, 1] / st$block_deg), floor(xy[, 2] / st$block_deg))
  ub <- unique(blk)
  ev <- sample(nrow(X), min(st$n_eval, nrow(X)), prob = w)
  draw <- function() {
    b <- sample(ub, round(length(ub) * st$frac))
    pool <- which(blk %in% b)
    pool[sample.int(length(pool), min(st$n_fit, length(pool)), prob = if (is.null(w)) NULL else w[pool])]
  }
  res <- expand.grid(boot = seq_len(st$n_boot), k = ks)
  res$ari <- NA_real_
  for (b in seq_len(st$n_boot)) {
    i1 <- draw(); i2 <- draw()
    for (k in ks) {
      k1 <- kmeans_fit(X[i1, ], k, nstart); k2 <- kmeans_fit(X[i2, ], k, nstart)
      res$ari[res$boot == b & res$k == k] <-
        ari(assign_nearest(X[ev, ], k1$centers), assign_nearest(X[ev, ], k2$centers))
    }
  }
  res
}

# Pick n k values at local maxima of the mean stability curve (ties broken by stability),
# falling back to the highest-stability k values.
pick_k <- function(stab, n = 3, min_k = 5) {
  m <- stats::aggregate(ari ~ k, stab, mean)
  m <- m[order(m$k), ]
  loc <- which(m$k >= min_k & m$ari >= c(-Inf, head(m$ari, -1)) & m$ari >= c(tail(m$ari, -1), -Inf))
  cand <- m[loc, ]
  if (nrow(cand) < n) cand <- rbind(cand, m[m$k >= min_k & !(m$k %in% cand$k), ])
  sort(head(cand[order(-cand$ari), "k"], n))
}

# Two-stage Ward: micro-cluster centroids with sizes; initial dissimilarity
# d_ij^2 = 2 n_i n_j / (n_i + n_j) |c_i - c_j|^2 makes the Lance-Williams Ward update exact
# for merges of whole micro-clusters.
ward_two_stage <- function(centers, sizes) {
  D <- as.matrix(stats::dist(centers))
  w <- sqrt(2 * outer(sizes, sizes) / outer(sizes, sizes, "+"))
  stats::hclust(stats::as.dist(D * w), method = "ward.D2", members = sizes)
}

# Relabel `b` to best match reference labels `a` (greedy on the contingency table).
match_labels <- function(a, b) {
  tab <- table(factor(b), factor(a)); map <- integer(0)
  rl <- rownames(tab); cl <- colnames(tab); used_r <- used_c <- character()
  while (length(used_r) < length(rl) && length(used_c) < length(cl)) {
    t2 <- tab[setdiff(rl, used_r), setdiff(cl, used_c), drop = FALSE]
    ix <- which(t2 == max(t2), arr.ind = TRUE)[1, ]
    r <- rownames(t2)[ix[1]]; c <- colnames(t2)[ix[2]]
    map[r] <- as.integer(c); used_r <- c(used_r, r); used_c <- c(used_c, c)
  }
  extra <- setdiff(rl, used_r)
  if (length(extra)) map[extra] <- max(as.integer(cl)) + seq_along(extra)
  unname(map[as.character(b)])
}

# Order cluster ids so similar clusters get nearby ids (by the first PC of the centroids)
# and colors from PC1/PC2 angle of the centroids, so similar regimes look alike on maps.
order_and_colors <- function(centers) {
  p <- stats::prcomp(centers, center = TRUE, scale. = FALSE)
  s <- p$x[, 1:2, drop = FALSE]
  ord <- order(s[, 1])
  h <- (atan2(s[, 2], s[, 1]) * 180 / pi) %% 360
  r <- sqrt(rowSums(s^2)); r <- r / max(r)
  col <- grDevices::hcl(h, 35 + 45 * r, 72 - 22 * r)
  list(order = ord, col = col)
}

# Write integer labels (land cells, grid order) as bands of one raster.
write_labels <- function(L, grid, out, datatype = "INT2U") {
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  r <- terra::rast(grid, nlyrs = ncol(L))
  m <- matrix(NA_real_, terra::ncell(grid), ncol(L)); m[land, ] <- L
  terra::values(r) <- m; names(r) <- colnames(L)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(out, ".tmp.tif")
  terra::writeRaster(r, tmp, overwrite = TRUE, datatype = datatype,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES", "INTERLEAVE=BAND"))
  unlink(out); file.rename(tmp, out)
  invisible(out)
}

# ---- Variants: feature groups and cell weights ---------------------------------------------

# Reduce Block 1 by feature groups: each group standardized + PCA, scaled so all groups share
# the block's total variance equally (block total = weight). groups = NULL: one group.
groups_reduce <- function(X, groups, weight = 1, pca_var = 0.99, log_cols = character()) {
  if (is.null(groups)) groups <- list(all = colnames(X))
  gw <- weight / length(groups)
  parts <- lapply(names(groups), function(g) {
    r <- block_reduce(X[, unlist(groups[[g]]), drop = FALSE], gw, pca_var, log_cols)
    colnames(r$scores) <- paste0(g, "_PC", seq_len(ncol(r$scores)))
    r
  })
  names(parts) <- names(groups)
  list(scores = do.call(cbind, lapply(parts, `[[`, "scores")), parts = parts)
}

# Cell weights from a static layer: floor + (1 - floor) * share, so no cell has zero weight.
cell_weights <- function(v, floor = 0.1) floor + (1 - floor) * pmin(pmax(v, 0), 1)

# k-means with optional cell weights. Unweighted: kmeans_fit on all cells. Weighted: start
# from k-means on a weight-proportional sample, then weighted Lloyd iterations on all cells
# (centroids = weighted means), which minimizes the weighted within-cluster sum of squares.
wkmeans <- function(X, w, k, nstart = 5, n_sample = 150000, max_iter = 50) {
  if (is.null(w)) {
    km <- kmeans_fit(X, k, nstart)
    return(list(centers = km$centers, cluster = km$cluster, size = km$size,
                tot.withinss = km$tot.withinss, iter = km$iter))
  }
  idx <- sample.int(nrow(X), min(n_sample, nrow(X)), prob = w)
  C <- kmeans_fit(X[idx, , drop = FALSE], k, nstart)$centers
  for (it in seq_len(max_iter)) {
    cl <- assign_nearest(X, C)
    sw <- as.vector(rowsum(w, factor(cl, levels = seq_len(k))))
    Cn <- rowsum(X * w, factor(cl, levels = seq_len(k))) / sw
    empty <- !is.finite(Cn[, 1])
    if (any(empty)) Cn[empty, ] <- C[empty, ]
    moved <- max(abs(Cn - C)); C <- Cn
    if (moved < 1e-6) break
  }
  cl <- assign_nearest(X, C)
  wss <- sum(w * rowSums((X - C[cl, , drop = FALSE])^2))
  list(centers = C, cluster = cl, size = tabulate(cl, k), tot.withinss = wss, iter = it)
}

# Weighted total sum of squares (about the weighted mean).
wtss <- function(X, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(X))
  m <- colSums(X * w) / sum(w)
  sum(w * rowSums(sweep(X, 2, m)^2))
}
