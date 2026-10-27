#!/usr/bin/env Rscript
# Step 73 (fuel phenology by ecoregion x fuel type; laptop, ~5 min): does the ecoregion x fuel split
# earn its keep? Groupings of the same 4 km cells are scored on how well they explain satellite and
# ground green-up / curing timing. Decided Oct 26 2026 (Mike): run this comparison, with the season
# split of two-season units as one of the groupings, before more units work.
#
# Groupings (every land cell gets a label; fuel = the cell's dominant group from 71):
#   eco2               EPA / CEC Level II ecoregion
#   eco2_fuel          Level II x fuel group (the 72 units, wall to wall)
#   eco2_fuel_season   ... each unit split in two (early / late season) where its cells' typical
#                      (green-up, curing) form two clusters whose centres are > season_split_days apart
#                      and the smaller holds >= season_min_share (2-means; units >= season_min_cells cells)
#   eco2_fuel_elev     ... x 500 m elevation bands (PRISM DEM)
#   eco3 / eco3_fuel   EPA Level III ecoregion (downloaded once; config ecofuel: ecoregions_l3) [x fuel]
#   a2 / a2_fuel       the climate regions (A2, 13 regions) [x fuel]
# Scores
#   satellite (cells of the 72 climatology units at >= 50 % purity, so all groupings use the same cells):
#     r2 / adj_r2      share of the variance of the cells' typical dates (median over 2001-2025) explained
#     cv_r2            the same predicting each cell from its group's mean in cells outside its
#                      block_deg x block_deg block (groups with no cells outside: CONUS mean); penalises
#                      groupings that only fit nearby cells
#     sync             share of cells' year-to-year departures shared by their group that year (as 72)
#   ground (NPN record sites; site date = median over its species series; labels from the site's cell):
#     cv_r2 / cv_mae   each site predicted by the median of the other sites in its group, leaving out
#                      sites within cluster_km (none left: median of all other sites); own_group = share
#                      of sites predicted from their own group
# Outputs <run_dir>/ecofuel/compare/: satellite.csv, ground.csv, groupings.csv (groups, cells, median
# cells per group), season_split.csv (units split, centres, shares); ecofuel/eco_l3.tif + lookup.
# Usage: Rscript scripts/73_ecofuel_compare.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "ecofuel.R"))
if (!requireNamespace("sf", quietly = TRUE)) stop("needs the sf package")

cfg <- load_cfg(profile)
log_file <- log_init("73_ecofuel_compare", profile)
ec <- cfg$ecofuel; cl <- ec$climatology; cp <- ec$compare
odir <- ecofuel_dir(cfg, "compare"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
mets <- unlist(cp$metrics %||% list("sos20", "eos50"))
yrs <- seq(as.integer(cl$years$start), as.integer(cl$years$end))
t0 <- Sys.time()
log_msg("=== 73_ecofuel_compare | profile: ", profile, " | metrics ", paste(mets, collapse = ", "))

mask_file <- file.path(static_dir(cfg), "grid_mask.tif")
mask <- terra::rast(mask_file)
land <- which(!is.na(terra::values(mask, mat = FALSE)))
lv <- function(f, band = 1) terra::values(terra::rast(f)[[band]], mat = FALSE)[land]

# ---- labels ---------------------------------------------------------------------------------------
L2 <- utils::read.csv(ecofuel_dir(cfg, "eco_l2_lookup.csv"), stringsAsFactors = FALSE, colClasses = c(code = "character"))
eco2 <- L2$code[match(lv(ecofuel_dir(cfg, "eco_l2.tif"), "eco"), L2$eco)]          # codes, e.g. "9.4"
dom_f <- ecofuel_dir(cfg, "inventory", "dominant.tif")
if (!file.exists(dom_f)) stop("no ", dom_f, " (run 71)")
DG <- utils::read.csv(ecofuel_dir(cfg, "inventory", "dominant_groups.csv"), stringsAsFactors = FALSE)
fuel <- DG$group[match(lv(dom_f, "dominant"), DG$dominant)]
l3_file <- ecofuel_dir(cfg, "eco_l3.tif"); l3_lookup <- ecofuel_dir(cfg, "eco_l3_lookup.csv")
if (!file.exists(l3_file) || !file.exists(l3_lookup)) {
  G <- terra::rast(ecofuel_dir(cfg, "fuel_groups.tif"))
  ecofuel_eco_raster(cfg, ec$ecoregions_l3, mask_file, l3_file, l3_lookup, water = terra::values(G[["NB_water"]], mat = FALSE))
}
L3 <- utils::read.csv(l3_lookup, stringsAsFactors = FALSE, colClasses = c(code = "character"))
eco3 <- L3$code[match(lv(l3_file, "eco"), L3$eco)]
elev_f <- file.path(static_dir(cfg), "elev.tif")
elev <- if (file.exists(elev_f)) lv(elev_f) else { log_warn("no static/elev.tif: elevation grouping skipped"); NULL }
a2 <- region_set_labels(cfg, mask, list(partition = cp$a2$partition %||% "A2_groups", k = as.integer(cp$a2$k %||% 13)))
if (is.null(a2)) log_warn("A2 k", cp$a2$k, " labels not found under clusters/: A2 groupings skipped")

# typical dates per cell (median over the years; summaries from 06)
med_f <- file.path(cfg$run_dir, "summaries", "cutb_metrics_median.tif")
if (!file.exists(med_f)) stop("no ", med_f)
TY <- terra::values(terra::rast(med_f)[[mets]])[land, , drop = FALSE]

# season split of eco2 x fuel units (2-means on the cells' typical dates)
lab <- function(...) { x <- paste(..., sep = "|"); x[Reduce(`|`, lapply(list(...), is.na))] <- NA; x }
ef <- lab(eco2, fuel)
season <- rep("all", length(land)); SPL <- list()
min_cells <- as.integer(cp$season_min_cells %||% 200)
set.seed(cfg$clustering$seed %||% 1)
for (u in names(which(table(ef) >= min_cells))) {
  i <- which(ef == u & stats::complete.cases(TY))
  if (length(i) < min_cells) next
  km <- tryCatch(stats::kmeans(TY[i, , drop = FALSE], centers = 2, nstart = 5), error = function(e) NULL)
  if (is.null(km)) next
  d <- sqrt(sum((km$centers[1, ] - km$centers[2, ])^2)); sh <- min(table(km$cluster)) / length(i)
  if (d > as.numeric(cp$season_split_days %||% 60) && sh >= as.numeric(cp$season_min_share %||% 0.1)) {
    early <- which.min(km$centers[, 1])
    season[i] <- ifelse(km$cluster == early, "early", "late")
    SPL[[u]] <- data.frame(unit = u, cells = length(i), early_share = mean(km$cluster == early), centre_distance = d,
                           early_sos20 = km$centers[early, 1], late_sos20 = km$centers[3 - early, 1],
                           early_eos50 = km$centers[early, 2], late_eos50 = km$centers[3 - early, 2])
  }
}
SPL <- do.call(rbind, SPL)
if (!is.null(SPL)) { SPL$unit <- sub("\\|", " ", SPL$unit); utils::write.csv(SPL, file.path(odir, "season_split.csv"), row.names = FALSE) }
log_msg(sprintf("season split: %d units split in two: %s", if (is.null(SPL)) 0 else nrow(SPL),
                if (is.null(SPL)) "-" else paste(sprintf("%s (%.0f%% early)", SPL$unit, 100 * SPL$early_share), collapse = ", ")))

GR <- list(eco2 = lab(eco2), eco2_fuel = ef, eco2_fuel_season = lab(eco2, fuel, season))
if (!is.null(elev)) GR$eco2_fuel_elev <- lab(eco2, fuel, floor(elev / as.numeric(cp$elev_band_m %||% 500)))
GR$eco3 <- lab(eco3); GR$eco3_fuel <- lab(eco3, fuel)
if (!is.null(a2)) { GR$a2 <- lab(a2); GR$a2_fuel <- lab(a2, fuel) }

# ---- satellite: the 72 climatology unit cells ------------------------------------------------------
UM <- ecofuel_unit_cells(cfg, as.numeric(cl$purity))
UT <- utils::read.csv(ecofuel_dir(cfg, "climatology", "units.csv"), stringsAsFactors = FALSE)
UM <- UM[UM$unit %in% UT$unit[UT$eligible], ]
dom_rows <- sort(unique(UM$row))
log_msg(sprintf("satellite domain: %d cells in %d climatology units", length(dom_rows), length(unique(UM$unit))))
xy <- terra::xyFromCell(mask, land[dom_rows])
bd <- as.numeric(cp$block_deg %||% 2)
block <- paste(floor(xy[, 1] / bd), floor(xy[, 2] / bd))
cv_r2 <- function(y, g, block) {
  ok <- !is.na(y) & !is.na(g); y <- y[ok]; g <- g[ok]; b <- block[ok]
  pred <- rep(NA_real_, length(y)); mu <- mean(y)
  for (bb in unique(b)) {
    te <- b == bb; tr <- !te
    s <- rowsum(y[tr], g[tr]); n <- rowsum(rep(1, sum(tr)), g[tr])
    m <- (s / n)[, 1]
    p <- m[g[te]]; p[is.na(p)] <- mu
    pred[te] <- p
  }
  1 - sum((y - pred)^2) / sum((y - mu)^2)
}
# per-year departures for sync (as 72)
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
have <- yrs[file.exists(vapply(yrs, mfile, ""))]
DEP <- setNames(lapply(mets, function(m) matrix(NA_real_, length(dom_rows), length(have))), mets)
cells_dom <- land[dom_rows]
for (k in seq_along(have)) {
  v <- terra::values(terra::rast(mfile(have[k]))[[mets]])[cells_dom, , drop = FALSE]
  for (m in mets) DEP[[m]][, k] <- v[, m]
}
for (m in mets) {
  ty <- apply(DEP[[m]], 1, stats::median, na.rm = TRUE)
  d <- DEP[[m]] - ty; d[abs(d) > as.numeric(cl$max_jump %||% 90)] <- NA; DEP[[m]] <- d
}
sync_g <- function(D, g) {
  ok_r <- !is.na(g); D <- D[ok_r, , drop = FALSE]; g <- g[ok_r]
  tot <- 0; within <- 0; allv <- D[!is.na(D)]; mu <- mean(allv)
  for (j in seq_len(ncol(D))) {
    z <- D[, j]; o <- !is.na(z); if (sum(o) < 2) next
    gm <- (rowsum(z[o], g[o]) / rowsum(rep(1, sum(o)), g[o]))[, 1]
    within <- within + sum((z[o] - gm[g[o]])^2); tot <- tot + sum((z[o] - mu)^2)
  }
  if (tot > 0) 1 - within / tot else NA_real_
}
SAT <- do.call(rbind, lapply(names(GR), function(nm) {
  g <- GR[[nm]][dom_rows]
  do.call(rbind, lapply(mets, function(m) {
    y <- TY[dom_rows, m]; ra <- r2_adj(y, g)
    data.frame(grouping = nm, metric = m, cells = unname(ra["n"]), groups = unname(ra["groups"]), r2 = unname(ra["r2"]),
               adj_r2 = unname(ra["adj_r2"]), cv_r2 = cv_r2(y, g, block), sync = sync_g(DEP[[m]], g))
  }))
}))
utils::write.csv(SAT, file.path(odir, "satellite.csv"), row.names = FALSE)

GS <- do.call(rbind, lapply(names(GR), function(nm) {
  g <- GR[[nm]]; t <- table(g); td <- table(g[dom_rows])
  data.frame(grouping = nm, groups_all = length(t), groups_domain = length(td), median_cells_domain = stats::median(as.numeric(td)))
}))
utils::write.csv(GS, file.path(odir, "groupings.csv"), row.names = FALSE)

# ---- ground: NPN record sites -------------------------------------------------------------------
rec_f <- file.path(cfg$run_dir, "npn_study", "record", "series_typical.csv")
if (file.exists(rec_f)) {
  ST <- utils::read.csv(rec_f, stringsAsFactors = FALSE)
  S <- unique(ST[, c("site_id", "lat", "lon")]); S <- S[!duplicated(S$site_id), ]
  S$row <- match(terra::cellFromXY(mask, cbind(S$lon, S$lat)), land)
  P <- sf::st_as_sf(S, coords = c("lon", "lat"), crs = 4326)
  DM <- matrix(as.numeric(sf::st_distance(P)), nrow(S)) / 1000          # great-circle km
  near <- DM <= as.numeric(cp$cluster_km %||% 10)
  GD <- list()
  for (tg in intersect(c("greenup", "curing", "cured"), names(ST))) {
    a <- ST[!is.na(ST[[tg]]), c("site_id", tg)]
    sv <- tapply(a[[tg]], a$site_id, stats::median)
    k <- match(names(sv), S$site_id); y <- as.numeric(sv)
    for (nm in names(GR)) {
      g <- GR[[nm]][S$row[k]]
      pred <- own <- rep(NA, length(y))
      for (i in seq_along(y)) {
        keep <- !near[k[i], k]; keep[i] <- FALSE
        same <- keep & !is.na(g) & g == g[i]
        if (!is.na(g[i]) && any(same)) { pred[i] <- stats::median(y[same]); own[i] <- TRUE } else { pred[i] <- stats::median(y[keep]); own[i] <- FALSE }
      }
      GD[[length(GD) + 1]] <- data.frame(grouping = nm, target = tg, sites = length(y), groups_with_sites = length(unique(g[!is.na(g)])),
        cv_r2 = if (length(y) >= 5) 1 - sum((y - pred)^2) / sum((y - mean(y))^2) else NA_real_, cv_mae = mean(abs(y - pred)), own_group = mean(own))
    }
  }
  GD <- do.call(rbind, GD)
  utils::write.csv(GD, file.path(odir, "ground.csv"), row.names = FALSE)
} else { GD <- NULL; log_warn("no NPN record (", rec_f, "): ground scores skipped (run on the laptop)") }

# ---- log --------------------------------------------------------------------------------------------
for (nm in names(GR)) {
  z <- SAT[SAT$grouping == nm, ]; g <- GS[GS$grouping == nm, ]
  gd <- if (is.null(GD)) NULL else GD[GD$grouping == nm, ]
  gtxt <- if (is.null(gd)) "" else paste(sprintf("%s cv %.2f (own %.0f%%)", gd$target, gd$cv_r2, 100 * gd$own_group), collapse = " | ")
  log_msg(sprintf("  %-17s groups %4d | sat %s | ground %s", nm, g$groups_domain,
                  paste(sprintf("%s r2 %.2f cv %.2f sync %.2f", z$metric, z$r2, z$cv_r2, z$sync), collapse = " | "), gtxt))
}
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
