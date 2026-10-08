#!/usr/bin/env Rscript
# Step 32 (evaluation, region profiles): what each candidate region looks like, to help choose
# the working regions. For every region set in config `region_profiles: sets`:
#   features.csv  climate and site features per region (all cells): 25 %, median, 75 % across
#                 cells (config `region_profiles: features`)
#   metrics.csv   satellite season metrics per region (cutb_metrics_median; green-up, peak,
#                 end, curing, amplitudes), for all cells and the strict herbaceous cells
#   curves.csv    median NDVI and NDII7 seasonal curves per region (25 %, median, 75 % across
#                 cells of the per-cell median curves), all cells and herbaceous cells
#   obs.csv       ground-observation dates per region (site medians from 31_obs_regions.R)
#   regions.csv   region key: number, colour, cells, share of land, herbaceous cells, label
#                 point (the region cell nearest its mean location)
# plus crosswalk_<fine>__<coarse>.csv for each config `region_profiles: crosswalk` pair (share
# of each fine region's cells in each coarse region).
# Region numbers and colours match the project-page maps (R/evaluate.R region_display()).
# Outputs <run_dir>/eval/profiles/<partition>_k<KK>/. Needs the Block 1/2 summaries, Cut B
# summaries, static layers and cluster rasters; obs.csv needs eval/obs/sites.csv (31), so run
# it on the laptop. A few minutes.
#
# Usage: Rscript scripts/32_region_profiles.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/32_region_profiles.R")

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
log_file <- log_init("32_region_profiles", profile)
rp <- cfg$region_profiles
out_root <- file.path(cfg$run_dir, "eval", "profiles")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
log_msg("=== 32_region_profiles | profile: ", profile, " | sets: ",
        paste(vapply(rp$sets, region_set_id, ""), collapse = ", "))

grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
xy <- terra::xyFromCell(grid, land)
D <- region_display(cfg, grid, rp$sets)
D <- Filter(Negate(is.null), D)
if (!length(D)) stop("none of the region sets found under ", file.path(cfg$run_dir, "clusters"))

# ---- Inputs (land cells) ------------------------------------------------------------------
sources <- list(b1med = summary_file(cfg, "block1", "median"), b2med = summary_file(cfg, "block2", "median"),
                b2var = summary_file(cfg, "block2", "var"), static = file.path(static_dir(cfg), "elev.tif"),
                herb = file.path(static_dir(cfg), "herb_share.tif"))
feat <- do.call(rbind, lapply(rp$features, as.data.frame))
sp <- do.call(rbind, strsplit(feat$name, ":", fixed = TRUE))
X <- matrix(NA_real_, length(land), nrow(feat), dimnames = list(NULL, feat$name))
for (s in unique(sp[, 1])) {
  f <- sources[[s]]
  if (is.null(f) || !file.exists(f)) { log_warn("feature source ", s, ": missing ", f %||% "(unknown)"); next }
  r <- terra::rast(f)
  if (s == "static") names(r)[1] <- "elev"
  b <- sp[sp[, 1] == s, 2]; ok <- b %in% names(r)
  if (any(!ok)) log_warn("feature source ", s, ": no band(s) ", paste(b[!ok], collapse = ", "))
  if (any(ok)) X[, feat$name[sp[, 1] == s][ok]] <- land_values(r, grid, b[ok])
}
met <- do.call(rbind, lapply(rp$metrics, as.data.frame))
mf <- summary_file(cfg, "cutb_metrics", "median")
M <- if (file.exists(mf)) land_values(terra::rast(mf), grid, intersect(met$name, names(terra::rast(mf)))) else NULL
if (is.null(M)) log_warn("missing ", mf, " (season metrics skipped)")
curve <- list()
for (ix in c("ndvi", "ndii")) {
  f <- summary_file(cfg, sprintf("cutb_%s_curve", ix), "median")
  if (file.exists(f)) curve[[ix]] <- land_values(terra::rast(f), grid) else log_warn("missing ", f)
}
hm <- land_values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), grid)[, 1]
herb <- !is.na(hm) & hm > 0
sites_f <- file.path(cfg$run_dir, "eval", "obs", "sites.csv")
S <- if (file.exists(sites_f)) utils::read.csv(sites_f, stringsAsFactors = FALSE) else NULL
if (is.null(S)) log_warn("no ", sites_f, " (run 31_obs_regions.R; obs.csv skipped)") else
  S <- S[S$set %in% unlist(rp$obs_sets), ]
log_msg(sprintf("inputs: %d land cells, %d features, %d metrics, curves: %s, %d herbaceous cells, %d sites",
                length(land), sum(colSums(!is.na(X)) > 0), if (is.null(M)) 0 else ncol(M),
                paste(names(curve), collapse = "+"), sum(herb), if (is.null(S)) 0 else nrow(S)))

# 25 %, median, 75 % of each column of V for the cells in idx (rows), with the count
qtab <- function(V, idx) {
  V <- V[idx, , drop = FALSE]
  q <- apply(V, 2, function(v) { v <- v[!is.na(v)]
    if (length(v) < 3) return(c(NA, NA, NA, length(v)))
    c(stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE), length(v)) })
  data.frame(var = colnames(V), q25 = q[1, ], med = q[2, ], q75 = q[3, ], n = q[4, ], row.names = NULL)
}

# ---- Per set ------------------------------------------------------------------------------
for (d in D) {
  od <- file.path(out_root, d$id); dir.create(od, showWarnings = FALSE)
  g <- d$num; k <- max(g, na.rm = TRUE)
  # region key with a label point: region cell nearest the region's mean location
  key <- do.call(rbind, lapply(seq_len(k), function(j) {
    i <- which(g == j); m <- colMeans(xy[i, , drop = FALSE])
    c0 <- i[which.min((xy[i, 1] - m[1])^2 * cos(m[2] * pi / 180)^2 + (xy[i, 2] - m[2])^2)]
    data.frame(region = j, colour = d$pal[j], cells = length(i), share = length(i) / length(g),
               herb_cells = sum(herb[i]), lon = xy[c0, 1], lat = xy[c0, 2],
               mean_lon = m[1], mean_lat = m[2])
  }))
  utils::write.csv(key, file.path(od, "regions.csv"), row.names = FALSE)
  # features (all cells)
  ft <- do.call(rbind, lapply(seq_len(k), function(j) cbind(region = j, qtab(X, which(g == j)))))
  ft$label <- feat$label[match(ft$var, feat$name)]; ft$unit <- feat$unit[match(ft$var, feat$name)]
  utils::write.csv(ft, file.path(od, "features.csv"), row.names = FALSE)
  # season metrics and curves, all cells and herbaceous cells
  mt <- data.frame(); cv <- data.frame()
  for (j in seq_len(k)) for (dom in c("all", "herb")) {
    idx <- which(g == j & (dom == "all" | herb))
    if (dom == "herb" && length(idx) < 3) next
    if (!is.null(M)) mt <- rbind(mt, cbind(region = j, domain = dom, qtab(M, idx)))
    for (ix in names(curve)) {
      q <- qtab(curve[[ix]], idx)
      cv <- rbind(cv, data.frame(region = j, domain = dom, index = ix,
                                 doy = as.integer(sub("doy", "", q$var)), q[, c("q25", "med", "q75", "n")]))
    }
  }
  if (nrow(mt)) {
    mt$label <- met$label[match(mt$var, met$name)]; mt$unit <- met$unit[match(mt$var, met$name)]
    utils::write.csv(mt, file.path(od, "metrics.csv"), row.names = FALSE)
  }
  if (nrow(cv)) utils::write.csv(cv, file.path(od, "curves.csv"), row.names = FALSE)
  # ground observations
  if (!is.null(S) && nrow(S)) {
    S$region <- g[S$row]
    ob <- do.call(rbind, lapply(split(S, list(S$set, S$region), drop = TRUE), function(s)
      data.frame(set = s$set[1], region = s$region[1], n_sites = nrow(s),
                 q25 = stats::quantile(s$doy, 0.25, names = FALSE), med = stats::median(s$doy),
                 q75 = stats::quantile(s$doy, 0.75, names = FALSE))))
    ob <- ob[order(ob$set, ob$region), ]
    utils::write.csv(ob, file.path(od, "obs.csv"), row.names = FALSE)
    utils::write.csv(S[!is.na(S$region), c("set", "site", "doy", "lon", "lat", "n_years", "region")],
                     file.path(od, "obs_sites.csv"), row.names = FALSE)
  }
  log_msg(sprintf("%s: %d regions, cells %s; herbaceous cells per region %s", d$id, k,
                  paste(key$cells, collapse = "/"), paste(key$herb_cells, collapse = "/")))
}

# ---- Crosswalks ---------------------------------------------------------------------------
for (p in rp$crosswalk %||% list()) {
  a <- D[[p[[1]]]]; b <- D[[p[[2]]]]
  if (is.null(a) || is.null(b)) { log_warn("crosswalk ", p[[1]], " -> ", p[[2]], ": set missing"); next }
  tab <- table(fine = a$num, coarse = b$num)
  cw <- as.data.frame(tab, stringsAsFactors = FALSE); cw <- cw[cw$Freq > 0, ]
  cw$fine <- as.integer(cw$fine); cw$coarse <- as.integer(cw$coarse)
  cw$share_of_fine <- cw$Freq / as.vector(rowSums(tab))[cw$fine]
  names(cw)[names(cw) == "Freq"] <- "cells"
  cw <- cw[order(cw$fine, -cw$share_of_fine), ]
  f <- file.path(out_root, sprintf("crosswalk_%s__%s.csv", p[[1]], p[[2]]))
  utils::write.csv(cw, f, row.names = FALSE)
  main <- cw[!duplicated(cw$fine), ]
  log_msg(sprintf("crosswalk %s -> %s: %s", p[[1]], p[[2]],
                  paste(sprintf("%d->%d (%.0f%%)", main$fine, main$coarse, 100 * main$share_of_fine), collapse = ", ")))
}

log_msg(sprintf("=== done | %d sets, %.1f min | outputs: %s | log: %s", length(D),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), out_root, log_file))
