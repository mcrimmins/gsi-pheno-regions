# Evaluation helpers: partitions (k-means label rasters with bands k04, k05, ...), a purely
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a
# geographic baseline, and explained-variance scores. Used by 30_compare_ab.R and
# 31_obs_regions.R. Config `evaluation:`.

# Label matrix (land cells x k bands) for each configured partition that exists.
eval_partitions <- function(cfg, grid, ks) {
  ev <- cfg$evaluation
  P <- list()
  for (nm in names(ev$partitions)) {
    f <- file.path(cfg$run_dir, "clusters", ev$partitions[[nm]])
    if (!file.exists(f)) { log_warn("partition ", nm, ": no ", f, " (skipped)"); next }
    r <- terra::rast(f)
    kb <- intersect(sprintf("k%02d", ks), names(r))
    if (!length(kb)) { log_warn("partition ", nm, ": none of k ", paste(range(ks), collapse = "-"), " in ", f); next }
    P[[nm]] <- land_values(r, grid, kb)
  }
  if (isTRUE(ev$geo_baseline)) P[["geo"]] <- geo_partition(cfg, grid, ks)
  P
}

# Geographic baseline: k-means on location only (lon scaled by cos(lat)), fitted on a sample,
# every land cell assigned to the nearest centre. Cached in clusters/geo/geo_kmeans.tif.
geo_partition <- function(cfg, grid, ks) {
  f <- file.path(cfg$run_dir, "clusters", "geo", "geo_kmeans.tif")
  if (file.exists(f)) {
    r <- terra::rast(f)
    if (all(sprintf("k%02d", ks) %in% names(r))) return(land_values(r, grid, sprintf("k%02d", ks)))
  }
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  xy <- terra::xyFromCell(grid, land)
  X <- cbind(xy[, 1] * cos(mean(xy[, 2]) * pi / 180), xy[, 2])
  set.seed(cfg$clustering$seed + 7)
  i <- sample(nrow(X), min(60000, nrow(X)))
  L <- sapply(ks, function(k) assign_nearest(X, kmeans_fit(X[i, ], k, nstart = 5)$centers))
  colnames(L) <- sprintf("k%02d", ks)
  write_labels(L, grid, f)
  L
}

# Share of variance explained by groups g: between-group SS / total SS, summed over the
# columns of Y (multivariate R^2). Rows with NA in g or Y are dropped.
r2_groups <- function(Y, g) {
  Y <- as.matrix(Y)
  use <- !is.na(g) & stats::complete.cases(Y)
  Y <- Y[use, , drop = FALSE]; g <- g[use]
  if (nrow(Y) < 3 || length(unique(g)) < 2) return(NA_real_)
  tss <- sum(scale(Y, scale = FALSE)^2)
  gm <- rowsum(Y, g) / as.vector(table(g)[rownames(rowsum(Y, g))])
  wss <- sum((Y - gm[as.character(g), , drop = FALSE])^2)
  1 - wss / tss
}

# Univariate R^2, adjusted R^2 (penalizes the number of groups) and counts.
r2_adj <- function(y, g) {
  use <- !is.na(y) & !is.na(g); y <- y[use]; g <- g[use]
  n <- length(y); m <- length(unique(g))
  if (n < 5 || m < 2 || n <= m) return(c(n = n, groups = m, r2 = NA, adj_r2 = NA))
  r2 <- r2_groups(matrix(y), g)
  c(n = n, groups = m, r2 = r2, adj_r2 = 1 - (1 - r2) * (n - 1) / (n - m))
}

# ---- Region display labels (report maps and region profiles) ------------------------------
# Qualitative palette for up to 20 regions (distinct hues, two lightness levels).
region_pal <- c("#2a78d6", "#eb6834", "#1a7f37", "#c0392b", "#7b4fb3", "#d4a020", "#17a2b8",
                "#8c564b", "#e377c2", "#6b8e23", "#9ec5f4", "#f5b386", "#98df8a", "#f19c99",
                "#c5b0d5", "#f0d58c", "#9edae5", "#c49c94", "#f7b6d2", "#4d4d4d")

# Raw k-means labels (land cells) of one region set list(partition, k); NULL if missing.
region_set_labels <- function(cfg, grid, s) {
  f <- file.path(cfg$run_dir, "clusters", cfg$evaluation$partitions[[s$partition]] %||% "none")
  if (!file.exists(f)) return(NULL)
  r <- terra::rast(f); kn <- sprintf("k%02d", s$k)
  if (!kn %in% names(r)) return(NULL)
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  terra::values(r[[kn]], mat = FALSE)[land]
}
region_set_id <- function(s) sprintf("%s_k%02d", s$partition, as.integer(s$k))

# Display labels for region sets, shared by the report maps and the profiles so numbers and
# colours agree everywhere. The reference set's regions are ordered south to north (mean
# latitude); every other set's regions are matched to them by overlap. Returns per set:
# col (colour slot per land cell, index into region_pal), num (region number 1..k per land
# cell: colour slots in increasing order), pal (colour per region number).
region_display <- function(cfg, grid, sets, ref = cfg$report$region_maps[[1]]) {
  land <- which(!is.na(terra::values(grid, mat = FALSE)))
  lat <- terra::yFromCell(grid, land)
  r0 <- region_set_labels(cfg, grid, ref)
  if (is.null(r0)) stop("reference region set ", region_set_id(ref), " not found")
  m <- tapply(lat, r0, mean); r0 <- match(r0, as.integer(names(m))[order(m)])
  out <- lapply(sets, function(s) {
    l <- region_set_labels(cfg, grid, s)
    if (is.null(l)) return(NULL)
    col <- if (identical(s$partition, ref$partition) && identical(as.integer(s$k), as.integer(ref$k))) r0 else match_labels(r0, l)
    slots <- sort(unique(col[!is.na(col)]))
    list(id = region_set_id(s), set = s, col = col, num = match(col, slots),
         pal = rep_len(region_pal, max(slots))[slots])
  })
  names(out) <- vapply(sets, region_set_id, "")
  out
}
