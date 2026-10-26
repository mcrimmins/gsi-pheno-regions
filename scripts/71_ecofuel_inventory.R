#!/usr/bin/env Rscript
# Step 71 (fuel phenology by ecoregion x fuel type; laptop or P720, < 1 min): inventory of the
# candidate units, ecoregion x fuel group, to settle the open questions (fuel grouping, minimum unit
# size, whether shrub / timber groups and fire are worth carrying) BEFORE any climatology is built.
# Nothing here is a decision: it reports areas, purity and data counts.
#
# Inputs (read only):
#   ecofuel/fuel_groups.tif, fbfm_shares.tif, eco_l<L>.tif, eco_l<L>_lookup.csv   (70)
#   static/herb_mask.tif                                                         (03)
#   npn_study/data/npn_site_ecoregions.csv, npn_site_fuelmodels.csv              (NPN study, in git)
#   <run_dir>/npn_study/record/series_typical.csv  (optional: sites in the usable record; laptop)
#   fire_study/fod/fod_fires.parquet|rds           (optional: FPA-FOD fires on the grid; 60)
#
# Definitions
#   share of land  = group share of mapped pixels / (1 - water share)   (as 03_herb_mask.R)
#   area           = group share x mapped share x cell area (km2)
#   dominant group = the group with the largest share of land in the cell (all 10 land groups,
#                    nonburnable included); pure_p50 / pure_p75 = cells where one group has >= 0.5 /
#                    >= 0.75 of the land (config ecofuel: purity). Cells less than half mapped by
#                    LANDFIRE (borders, coast; config ecofuel: min_mapped) count toward area only.
#   NPN sites      = site fuel group at the point if burnable, else the most common burnable model
#                    within 90 m (59b's rule, as 59c); counted in the site's Level II ecoregion (59)
#   fires          = FPA-FOD wildfires >= 10 ac, debris burning removed (the fire study's primary
#                    outcome), by the 4 km cell: fires_dom = in cells where the group is dominant;
#                    fires_wt = sum over fires of the group's share of land in the fire's cell. The
#                    4 km share is used, not the 30 m fuel model at the fire point: LF2023 maps fuels
#                    after the fire, and FPA-FOD points are often imprecise.
#
# Outputs <run_dir>/ecofuel/inventory/
#   units.csv       ecoregion x group: land_km2, share_of_eco, cells_dom, cells_p50, cells_p75,
#                   herb_cells_dom, npn_sites, npn_sites_record, fires_dom, fires_wt, fires_wt_per_yr
#   ecoregions.csv  per ecoregion: cells, land_km2, share burnable / agriculture / urban, dominant
#                   group mix, NPN sites, fires
#   models.csv      ecoregion x FBFM40 model: km2 (for grouping decisions, e.g. GR1-2 vs GR3+)
#   purity.csv      CONUS cells by dominant group and max share class (for the purity question)
#   dominant.tif    per cell: dominant (index into dominant_groups.csv), dominant_share (for 79's maps)
# Level of ecoregion from config ecofuel: ecoregions: level (2).
#
# Usage: Rscript scripts/71_ecofuel_inventory.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "ecofuel.R"))

cfg <- load_cfg(profile)
log_file <- log_init("71_ecofuel_inventory", profile)
ec <- cfg$ecofuel %||% list()
lev <- as.integer(ec$ecoregions$level %||% 2)
pur <- as.numeric(unlist(ec$purity %||% list(0.5, 0.75)))
odir <- ecofuel_dir(cfg, "inventory"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 71_ecofuel_inventory | profile: ", profile, " | Level ", lev, " | purity ", paste(pur, collapse = "/"))

need <- c(ecofuel_dir(cfg, "fuel_groups.tif"), ecofuel_dir(cfg, "fbfm_shares.tif"),
          ecofuel_dir(cfg, sprintf("eco_l%d.tif", lev)), ecofuel_dir(cfg, sprintf("eco_l%d_lookup.csv", lev)))
if (!all(file.exists(need))) stop("missing 70 outputs: ", paste(basename(need[!file.exists(need)]), collapse = ", "),
                                  " (run 70, or copy ecofuel/ from the P720 with 91 ecofuel)")
mask <- grid_mask(cfg)
land <- which(!is.na(terra::values(mask, mat = FALSE)))
G <- terra::values(terra::rast(need[1]))[land, , drop = FALSE]
S <- terra::values(terra::rast(need[2]))[land, , drop = FALSE]
EC <- terra::values(terra::rast(need[3]))[land, , drop = FALSE]
L <- utils::read.csv(need[4], stringsAsFactors = FALSE, colClasses = c(code = "character"))
ca <- terra::values(terra::cellSize(mask, unit = "km"), mat = FALSE)[land]
hm_f <- file.path(static_dir(cfg), "herb_mask.tif")
herb <- if (file.exists(hm_f)) { h <- terra::values(terra::rast(hm_f), mat = FALSE)[land]; !is.na(h) & h > 0 } else {
  log_warn("no static/herb_mask.tif: herb_cells_dom left NA"); rep(NA, length(land)) }

grp <- setdiff(ecofuel_groups(), "NB_water")              # the 10 land groups
FT <- fbfm_table()
water <- G[, "NB_water"]; water[is.na(water)] <- 0
mapped <- G[, "mapped"]; mapped[is.na(mapped)] <- 0
lw <- 1 - water
SL <- G[, grp, drop = FALSE] / ifelse(lw > 0.01, lw, NA)      # share of land
AREA <- G[, grp, drop = FALSE] * mapped * ca                  # km2 of each group in the cell
eco <- EC[, "eco"]
min_mapped <- as.numeric(ec$min_mapped %||% 0.5)
ok <- !is.na(eco) & rowSums(!is.na(SL)) > 0
log_msg(sprintf("cells with < %.0f %% of their area mapped by LANDFIRE (border, coast): %d; their area counts, but not dominance / purity",
                100 * min_mapped, sum(ok & mapped < min_mapped)))
SLd <- SL; SLd[mapped < min_mapped, ] <- NA                   # shares used for dominance / purity
okd <- ok & rowSums(!is.na(SLd)) > 0
dom_i <- rep(NA_integer_, nrow(SL)); dom_i[okd] <- max.col(replace(SLd[okd, , drop = FALSE], is.na(SLd[okd, , drop = FALSE]), -1), ties.method = "first")
dom <- grp[dom_i]
dom_share <- rep(NA_real_, nrow(SL)); dom_share[okd] <- SLd[cbind(which(okd), dom_i[okd])]
log_msg(sprintf("land cells %d | with a mapped fuel and an ecoregion %d | median dominant share %.2f",
                length(land), sum(ok), stats::median(dom_share, na.rm = TRUE)))

# ---- NPN sites (read only) ----------------------------------------------------------------------
npn <- NULL
fe <- here::here("npn_study", "data", "npn_site_ecoregions.csv"); ff <- here::here("npn_study", "data", "npn_site_fuelmodels.csv")
if (lev == 2 && file.exists(fe) && file.exists(ff)) {
  SE <- utils::read.csv(fe, stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
  SF <- utils::read.csv(ff, stringsAsFactors = FALSE)
  npn <- merge(SE[, c("site_id", "eco_l2_code")], SF[, c("site_id", "fbfm_group", "fbfm_3x3_group", "fbfm", "fbfm_3x3")], by = "site_id")
  burn <- !is.na(npn$fbfm_group) & npn$fbfm_group != "NB"
  npn$group <- ifelse(burn, npn$fbfm_group, npn$fbfm_3x3_group)
  npn$model <- ifelse(burn, npn$fbfm, npn$fbfm_3x3)
  npn$group[is.na(npn$group)] <- "none"
  rec_f <- file.path(cfg$run_dir, "npn_study", "record", "series_typical.csv")
  npn$record <- if (file.exists(rec_f)) npn$site_id %in% utils::read.csv(rec_f)$site_id else NA
  log_msg(sprintf("NPN sites: %d (usable record: %s) | groups: %s", nrow(npn),
                  if (file.exists(rec_f)) sum(npn$record) else "file not here",
                  paste(names(table(npn$group)), table(npn$group), collapse = ", ")))
} else log_warn("NPN site tables not found (or ecoregion level != 2): site counts left NA")

# ---- fires (read only) --------------------------------------------------------------------------
fires <- NULL
fstem <- fire_dir(cfg, "fod", "fod_fires")
if (!is.na(fire_table_file(fstem))) {
  F0 <- fire_read_table(fstem)
  yrs <- ec$fire_years %||% cfg$fire_study$fod$years %||% list(start = 2001, end = 2024)
  yr <- as.integer(format(as.Date(F0$date), "%Y"))
  keep <- F0$c_plus %in% TRUE & !(F0$debris %in% TRUE) & yr >= yrs$start & yr <= yrs$end
  li <- match(F0$cell[keep], land)                       # row in the land-cell tables
  fires <- data.frame(row = li[!is.na(li)])
  n_years <- yrs$end - yrs$start + 1
  log_msg(sprintf("FPA-FOD: %d C+ wildfires (no debris) %d-%d on grid land cells", nrow(fires), yrs$start, yrs$end))
} else log_warn("no fire_study/fod/fod_fires table here: fire counts left NA")

# ---- units: ecoregion x group --------------------------------------------------------------------
rows <- list()
for (e in L$eco) {
  ie <- which(ok & eco == e)
  eco_land <- sum(AREA[ie, ], na.rm = TRUE)
  for (g in grp) {
    j <- match(g, grp)
    a <- sum(AREA[ie, j], na.rm = TRUE)
    sl <- SLd[ie, j]
    r <- data.frame(eco = e, code = L$code[e], name = L$name[e], group = g,
                    land_km2 = round(a), share_of_eco = if (eco_land > 0) a / eco_land else NA,
                    cells_dom = sum(dom[ie] == g, na.rm = TRUE), cells_p50 = sum(sl >= pur[1], na.rm = TRUE),
                    cells_p75 = sum(sl >= pur[length(pur)], na.rm = TRUE),
                    herb_cells_dom = if (all(is.na(herb))) NA else sum(dom[ie] == g & herb[ie], na.rm = TRUE),
                    npn_sites = NA, npn_sites_record = NA, fires_dom = NA, fires_wt = NA, fires_wt_per_yr = NA)
    if (!is.null(npn)) {
      s <- npn$eco_l2_code == L$code[e] & npn$group == g
      r$npn_sites <- sum(s); if (!all(is.na(npn$record))) r$npn_sites_record <- sum(s & npn$record)
    }
    if (!is.null(fires)) {
      fr <- fires$row[fires$row %in% ie]
      r$fires_dom <- sum(dom[fr] == g, na.rm = TRUE)
      r$fires_wt <- round(sum(SL[fr, j], na.rm = TRUE), 1)
      r$fires_wt_per_yr <- round(r$fires_wt / n_years, 2)
    }
    rows[[length(rows) + 1]] <- r
  }
}
U <- do.call(rbind, rows)
utils::write.csv(U, file.path(odir, "units.csv"), row.names = FALSE)

# ---- per ecoregion ------------------------------------------------------------------------------
burn <- ecofuel_burnable()
ER <- do.call(rbind, lapply(L$eco, function(e) {
  u <- U[U$eco == e, ]; ie <- which(ok & eco == e)
  tot <- sum(u$land_km2)
  data.frame(eco = e, code = L$code[e], name = L$name[e], l1_name = L$l1_name[e], cells = length(ie),
             land_km2 = tot, share_burnable = sum(u$land_km2[u$group %in% burn]) / max(tot, 1),
             share_ag = sum(u$land_km2[u$group == "NB_ag"]) / max(tot, 1),
             share_urban = sum(u$land_km2[u$group == "NB_urban"]) / max(tot, 1),
             top_group = u$group[which.max(u$land_km2)],
             burnable_mix = paste(sprintf("%s %.0f%%", u$group[u$group %in% burn], 100 * u$share_of_eco[u$group %in% burn]), collapse = ", "),
             median_dom_share = stats::median(dom_share[ie], na.rm = TRUE),
             npn_sites = if (is.null(npn)) NA else sum(npn$eco_l2_code == L$code[e]),
             npn_sites_nonburnable = if (is.null(npn)) NA else sum(npn$eco_l2_code == L$code[e] & npn$group == "none"),
             fires = if (is.null(fires)) NA else sum(fires$row %in% ie))
}))
ER <- ER[order(-ER$land_km2), ]
utils::write.csv(ER, file.path(odir, "ecoregions.csv"), row.names = FALSE)

# ---- per ecoregion x model ----------------------------------------------------------------------
MA <- S[, FT$model, drop = FALSE] * S[, "mapped"] * ca
M <- do.call(rbind, lapply(L$eco, function(e) {
  ie <- which(ok & eco == e)
  data.frame(eco = e, code = L$code[e], name = L$name[e], model = FT$model, group = FT$group,
             km2 = round(colSums(MA[ie, , drop = FALSE], na.rm = TRUE)),
             npn_sites = if (is.null(npn)) NA else vapply(FT$model, function(m) sum(npn$eco_l2_code == L$code[e] & npn$model %in% m), 0))
}))
M <- M[M$km2 > 0 | (!is.na(M$npn_sites) & M$npn_sites > 0), ]
utils::write.csv(M, file.path(odir, "models.csv"), row.names = FALSE)

# ---- purity --------------------------------------------------------------------------------------
br <- c(0, 0.25, 0.5, 0.75, 1.0001)
cls <- cut(dom_share, br, right = FALSE, labels = c("<0.25", "0.25-0.5", "0.5-0.75", ">=0.75"))
P <- as.data.frame(table(dominant = factor(dom, levels = grp), max_share = cls), stringsAsFactors = FALSE)
names(P)[3] <- "cells"
utils::write.csv(P, file.path(odir, "purity.csv"), row.names = FALSE)

# ---- per-cell dominant group and its share (for the maps in 79) -----------------------------------
D <- terra::rast(mask, nlyrs = 2); names(D) <- c("dominant", "dominant_share")
dv <- matrix(NA_real_, terra::ncell(mask), 2); dv[land, 1] <- dom_i; dv[land, 2] <- dom_share
terra::values(D) <- dv
tmp <- file.path(odir, "dominant.tmp.tif")
terra::writeRaster(D, tmp, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
unlink(file.path(odir, "dominant.tif")); file.rename(tmp, file.path(odir, "dominant.tif"))
utils::write.csv(data.frame(dominant = seq_along(grp), group = grp), file.path(odir, "dominant_groups.csv"), row.names = FALSE)

# ---- log ---------------------------------------------------------------------------------------
for (i in seq_len(nrow(ER))) with(ER[i, ], log_msg(sprintf("  %-5s %-34s %7.0f km2 | burnable %3.0f %% ag %3.0f %% | NPN %3s | fires %6s | %s",
  code, substr(name, 1, 34), land_km2, 100 * share_burnable, 100 * share_ag,
  ifelse(is.na(npn_sites), "-", npn_sites), ifelse(is.na(fires), "-", fires), burnable_mix)))
ub <- U[U$group %in% burn, ]
for (t in c(50, 200, 1000)) {
  z <- ub[ub$cells_p50 >= t, ]
  log_msg(sprintf("burnable units with >= %d cells at >= %.0f %% purity: %d (GR %d, GS %d, SH %d, TU %d, TL %d)", t, 100 * pur[1], nrow(z),
                  sum(z$group == "GR"), sum(z$group == "GS"), sum(z$group == "SH"), sum(z$group == "TU"), sum(z$group == "TL")))
}
if (!is.null(npn)) {
  z <- ub[ub$group %in% c("GR", "GS"), ]
  log_msg(sprintf("GR/GS units: %d with any NPN site, %d with >= 3; sites in GR/GS units: %d", sum(z$npn_sites > 0), sum(z$npn_sites >= 3), sum(z$npn_sites)))
}
if (!is.null(fires)) log_msg(sprintf("fires by dominant group: %s", paste(sprintf("%s %d", grp, vapply(grp, function(g) sum(U$fires_dom[U$group == g]), 0)), collapse = ", ")))
log_msg("=== done | ", odir, " | log: ", log_file)
