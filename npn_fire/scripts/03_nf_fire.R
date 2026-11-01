#!/usr/bin/env Rscript
# Step 03 (NPN x fire study; P720, minutes; arrow, terra): fire activity by EPA Level II ecoregion x
# year x week x LANDFIRE fuel group, for every fire outcome we have. Decided 1 Nov 2026 (Mike):
# LANDFIRE fuel groups (ecofuel 70) as the fuel layer; FPA-FOD all sizes and C+; VIIRS from the fire
# study's CONUS-wide link screen (63 flags every CONUS detection).
#
# Inputs, read-only (parent pipeline data outputs; no parent code is sourced):
#   fire_study/fod/fod_fires            (60) FPA-FOD 7th ed., one row per fire on a grid land cell
#   fire_study/mtbs/mtbs_fires, _cells  (61) MTBS perimeters, fire x grid cell burned km2
#   fire_study/link/<SENSOR>/link_<SENSOR>_<Y>  (63) VIIRS detections with keep = wildfire & !crop
#   ecofuel/fuel_groups.tif             (70) shares of GR GS SH TU TL SB NB_* + mapped per 4 km cell
#   ecofuel/eco_l2.tif + eco_l2_lookup.csv (70) Level II ecoregion per cell (majority of sub-cells)
# Fuel weights per cell: share of each burnable group (GR GS SH TU TL SB) and of NB (ag, urban,
# barren, snow) over mapped non-water land, so they sum to 1; cells mapped < min_mapped -> "unmapped".
# Each fire is split over its cell's weights (FPA-FOD, VIIRS: the detection's cell; MTBS: burned km2
# per cell x that cell's weights). Matching fuel groups to plant groups is left to 05.
# Outcomes (value column; week = fixed 7-day blocks, days 365-366 in week 52):
#   fod_all          FPA-FOD wildfires, all sizes, debris / open burning removed, by discovery week
#   fod_all_lightning, fod_all_human, fod_all_unknown   the same by cause (cause coding changes after 2019)
#   fod_cplus        the same, fires >= min_size_ac (10 ac)
#   fod_cplus_acres  final acres of the C+ fires, at their discovery week
#   mtbs_fires       MTBS wildfires (Wildfire, Wildland Fire Use) by ignition week, split by burned area
#   mtbs_km2         MTBS wildfire burned km2 by ignition week
#   viirs_px         unique 375 m pixel-days of 63 `keep` detections (S-NPP; NOAA-20 fills whole UTC
#                    days without S-NPP data), by solar date
# Outputs <run_dir>/npn_fire/fire/:
#   fire_weekly.parquet  outcome, eco_l2_code, year, week, fuel, value (long; zero rows not stored)
#   fire_inventory.csv   outcome x ecoregion x fuel: total, per year, years and weeks with fire
#   fire_year.csv        outcome x year: CONUS total and share assigned to an ecoregion / mapped fuel
#   fire_filters.csv     records in / removed at each screen, per source
#   cell_weights.csv.gz  cell, eco_l2_code, weights (for 05's local checks)
#   viirs_fill.csv       outage days filled from NOAA-20, by year
#
# Usage: Rscript npn_fire/scripts/03_nf_fire.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("terra", quietly = TRUE)) stop("needs the terra package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("03_nf_fire", profile)
fc <- cfg$npn_fire$fire
fsd <- fc$fire_study_dir %||% file.path(cfg$run_dir, "fire_study")
efd <- fc$ecofuel_dir %||% file.path(cfg$run_dir, "ecofuel")
yrs <- seq(fc$years$start, fc$years$end)
out <- nf_out(cfg, "fire")
t0 <- Sys.time()
log_msg("=== 03_nf_fire | profile: ", profile, " | fire_study: ", fsd, " | ecofuel: ", efd)

# ---- cell weights and ecoregions ---------------------------------------------------------------------
G <- terra::rast(file.path(efd, "fuel_groups.tif"))
E <- terra::rast(file.path(efd, "eco_l2.tif"))
gm <- fc$grid_mask %||% file.path(cfg$run_dir, "static", "grid_mask.tif")
if (file.exists(gm)) {
  M <- terra::rast(gm)
  if (!isTRUE(all.equal(dim(M)[1:2], dim(G)[1:2])) || !isTRUE(all.equal(as.vector(terra::ext(M)), as.vector(terra::ext(G)), tolerance = 1e-6)))
    stop("fuel_groups.tif is not on the fire study's grid (", gm, "): cell numbers would not match")
  log_msg("grid check: fuel_groups.tif matches ", basename(gm), " (", terra::ncell(M), " cells)")
} else log_warn("no grid mask at ", gm, ": cell alignment not checked")
LK <- utils::read.csv(file.path(efd, "eco_l2_lookup.csv"), stringsAsFactors = FALSE, colClasses = c(code = "character"))
GV <- terra::values(G); ev <- terra::values(E[["eco"]], mat = FALSE)
fg <- unlist(fc$fuel_groups); nbg <- unlist(fc$nonburnable)
miss <- setdiff(c(fg, nbg, "NB_water", "mapped"), colnames(GV))
if (length(miss)) stop("fuel_groups.tif lacks band(s) ", paste(miss, collapse = ", "))
z <- function(v) { v[is.na(v)] <- 0; v }
B <- sapply(fg, function(g) z(GV[, g])); NB <- rowSums(sapply(nbg, function(g) z(GV[, g])))
tot <- rowSums(B) + NB
mapped <- z(GV[, "mapped"])
W <- cbind(B, NB = NB) / ifelse(tot > 0, tot, NA)
unm <- mapped < (fc$min_mapped %||% 0.5) | !is.finite(rowSums(W))
W[unm, ] <- 0
W <- cbind(W, unmapped = as.numeric(unm))
fuels <- colnames(W)
eco_code <- LK$code[match(ev, LK$eco)]
log_msg(sprintf("cells: %d with an ecoregion | %d unmapped (mapped < %s) | fuels: %s", sum(!is.na(eco_code)),
                sum(unm & !is.na(eco_code)), fc$min_mapped %||% 0.5, paste(fuels, collapse = " ")))
cw <- data.frame(cell = which(!is.na(eco_code)), eco_l2_code = eco_code[!is.na(eco_code)], round(W[!is.na(eco_code), , drop = FALSE], 4))
con <- gzfile(file.path(out, "cell_weights.csv.gz"), "w"); utils::write.csv(cw, con, row.names = FALSE); close(con)

FILT <- list()
filt <- function(source, step, n_in, n_out) FILT[[length(FILT) + 1]] <<- data.frame(source = source, step = step, records = n_out, removed = n_in - n_out)

# Spread events over fuels: cell, year, week, value -> long table outcome x eco x year x week x fuel
spread <- function(outcome, cell, year, week, value = 1) {
  ok <- !is.na(cell) & !is.na(eco_code[cell]) & !is.na(week)
  if (!any(ok)) return(NULL)
  cell <- cell[ok]; year <- year[ok]; week <- week[ok]; value <- rep_len(value, length(ok))[ok]
  key <- paste(eco_code[cell], year, week, sep = "|")
  S <- rowsum(W[cell, , drop = FALSE] * value, key, reorder = FALSE)
  k <- do.call(rbind, strsplit(rownames(S), "|", fixed = TRUE))
  L <- data.frame(outcome = outcome, eco_l2_code = rep(k[, 1], ncol(S)), year = as.integer(rep(k[, 2], ncol(S))),
                  week = as.integer(rep(k[, 3], ncol(S))), fuel = rep(colnames(S), each = nrow(S)), value = as.vector(S),
                  stringsAsFactors = FALSE)
  L[L$value > 0, ]
}
ALL <- list(); YR <- list()
year_tot <- function(outcome, year, value, assigned) YR[[length(YR) + 1]] <<- data.frame(outcome = outcome,
  year = sort(unique(year)), total = as.numeric(tapply(value, year, sum)[as.character(sort(unique(year)))]),
  in_ecoregion = as.numeric(tapply(value * assigned, year, sum)[as.character(sort(unique(year)))]))

# ---- FPA-FOD -----------------------------------------------------------------------------------------
F <- nf_read_parent(file.path(fsd, "fod", "fod_fires"))
n0 <- nrow(F); F <- F[F$year %in% yrs, ]; filt("fod", "study years", n0, nrow(F))
n0 <- nrow(F); F <- F[!F$debris %in% TRUE, ]; filt("fod", "debris / open burning removed", n0, nrow(F))
F$week <- nf_week(F$date)
asg <- !is.na(eco_code[F$cell])
filt("fod", "on a cell with a Level II ecoregion", nrow(F), sum(asg))
ALL$fod_all <- spread("fod_all", F$cell, F$year, F$week)
for (cz in c("lightning", "human", "unknown")) { i <- F$cause %in% cz
  ALL[[paste0("fod_all_", cz)]] <- spread(paste0("fod_all_", cz), F$cell[i], F$year[i], F$week[i]) }
cp <- !is.na(F$size_ac) & F$size_ac >= (fc$min_size_ac %||% 10)
ALL$fod_cplus <- spread("fod_cplus", F$cell[cp], F$year[cp], F$week[cp])
ALL$fod_cplus_acres <- spread("fod_cplus_acres", F$cell[cp], F$year[cp], F$week[cp], F$size_ac[cp])
year_tot("fod_all", F$year, rep(1, nrow(F)), asg); year_tot("fod_cplus", F$year[cp], rep(1, sum(cp)), asg[cp])
log_msg(sprintf("FPA-FOD: %d wildfires (no debris) %d-%d, %d C+ | %.1f %% on a Level II cell",
                nrow(F), min(yrs), max(yrs), sum(cp), 100 * mean(asg)))
rm(F); invisible(gc())

# ---- MTBS --------------------------------------------------------------------------------------------
MF <- nf_read_parent(file.path(fsd, "mtbs", "mtbs_fires"))
MC <- nf_read_parent(file.path(fsd, "mtbs", "mtbs_cells"))
n0 <- nrow(MF); MF <- MF[MF$year %in% yrs, ]; filt("mtbs", "study years", n0, nrow(MF))
n0 <- nrow(MF); MF <- MF[MF$type %in% unlist(fc$mtbs_wildfire_types), ]; filt("mtbs", "wildfire types only", n0, nrow(MF))
MF$week <- nf_week(MF$ig_date)
MC <- MC[MC$event_id %in% MF$event_id & MC$km2 > 0, ]
i <- match(MC$event_id, MF$event_id)
area <- tapply(MC$km2, MC$event_id, sum)
MC$share <- MC$km2 / as.numeric(area[MC$event_id])
ALL$mtbs_km2 <- spread("mtbs_km2", MC$cell, MF$year[i], MF$week[i], MC$km2)
ALL$mtbs_fires <- spread("mtbs_fires", MC$cell, MF$year[i], MF$week[i], MC$share)
asg_m <- !is.na(eco_code[MC$cell])
year_tot("mtbs_km2", MF$year[i], MC$km2, asg_m)
log_msg(sprintf("MTBS: %d wildfires %d-%d, %.0f km2 on %d fire x cell rows | %.1f %% of area on a Level II cell",
                nrow(MF), min(MF$year), max(MF$year), sum(MC$km2), nrow(MC), 100 * sum(MC$km2[asg_m]) / sum(MC$km2)))
rm(MF, MC); invisible(gc())

# ---- VIIRS pixel-days --------------------------------------------------------------------------------
vc <- fc$viirs; vy <- seq(vc$years$start, vc$years$end)
read_link <- function(sensor, Y) {
  stem <- file.path(fsd, "link", sensor, sprintf("link_%s_%d", sensor, Y))
  if (!any(file.exists(paste0(stem, c(".parquet", ".rds"))))) return(NULL)
  d <- nf_read_parent(stem)
  d[, intersect(c("lat", "lon", "date_utc", "date_solar", "cell", "keep"), names(d))]
}
FILL <- list(); VP <- list()
for (Y in vy) {
  P <- read_link(vc$primary, Y)
  if (is.null(P)) { log_warn("no ", vc$primary, " link table for ", Y); next }
  days <- seq(as.Date(sprintf("%d-01-01", Y)), as.Date(sprintf("%d-12-31", Y)), by = "day")
  if (Y == 2012) days <- days[days >= as.Date("2012-01-20")]                 # S-NPP record starts
  seen <- unique(as.Date(P$date_utc))
  out_days <- days[!days %in% seen]
  n_fill <- 0
  if (length(out_days)) {
    Q <- read_link(vc$fill, Y)
    if (!is.null(Q)) { Q <- Q[as.Date(Q$date_utc) %in% out_days, ]; n_fill <- length(unique(as.Date(Q$date_utc))); P <- rbind(P, Q) }
  }
  FILL[[length(FILL) + 1]] <- data.frame(year = Y, outage_days = length(out_days), filled_days = n_fill)
  n_all <- nrow(P); P <- P[P$keep %in% TRUE, ]
  filt("viirs", sprintf("%d: keep (linked wildfire, not crop)", Y), n_all, nrow(P))
  key <- nf_pixel_key(P$lon, P$lat, as.Date(P$date_solar))
  P <- P[!duplicated(key), ]
  P$year <- as.integer(format(as.Date(P$date_solar), "%Y")); P$week <- nf_week(P$date_solar)
  P <- P[P$year == Y, ]                                                       # solar date can cross Dec 31
  VP[[length(VP) + 1]] <- spread("viirs_px", P$cell, P$year, P$week)
  year_tot("viirs_px", P$year, rep(1, nrow(P)), !is.na(eco_code[P$cell]))
  log_msg(sprintf("VIIRS %d: %d kept pixel-days | %d outage days, %d filled from %s", Y, nrow(P), length(out_days), n_fill, vc$fill))
}
ALL$viirs_px <- do.call(rbind, VP)
utils::write.csv(do.call(rbind, FILL), file.path(out, "viirs_fill.csv"), row.names = FALSE)

# ---- write ------------------------------------------------------------------------------------------
FW <- do.call(rbind, ALL[!vapply(ALL, is.null, TRUE)])
rownames(FW) <- NULL
invisible(nf_write_table(FW, file.path(out, "fire_weekly")))
utils::write.csv(do.call(rbind, FILT), file.path(out, "fire_filters.csv"), row.names = FALSE)
YT <- do.call(rbind, YR); YT$share_in_ecoregion <- round(YT$in_ecoregion / YT$total, 4)
utils::write.csv(YT, file.path(out, "fire_year.csv"), row.names = FALSE)

key <- paste(FW$outcome, FW$eco_l2_code, FW$fuel, sep = "|")
k <- do.call(rbind, strsplit(unique(key), "|", fixed = TRUE))
INV <- data.frame(outcome = k[, 1], eco_l2_code = k[, 2], fuel = k[, 3], stringsAsFactors = FALSE)
uk <- unique(key)
INV$total <- as.numeric(tapply(FW$value, key, sum)[uk])
INV$years_with <- as.integer(tapply(FW$year, key, function(z) length(unique(z)))[uk])
INV$weeks_with <- as.integer(tapply(paste(FW$year, FW$week), key, function(z) length(unique(z)))[uk])
ny <- ifelse(INV$outcome == "viirs_px", length(vy), length(yrs))
INV$per_year <- round(INV$total / ny, 2)
INV$eco_l2_name <- LK$name[match(INV$eco_l2_code, LK$code)]
INV$mtbs_incomplete_last_year <- INV$outcome %in% c("mtbs_fires", "mtbs_km2") & max(yrs) > (fc$mtbs_complete_end %||% 2023)
INV <- INV[order(INV$outcome, INV$eco_l2_code, -INV$total), ]
utils::write.csv(INV, file.path(out, "fire_inventory.csv"), row.names = FALSE)

for (o in c("fod_all", "fod_cplus", "mtbs_fires", "viirs_px")) {
  s <- INV[INV$outcome == o, ]
  if (!nrow(s)) next
  bf <- tapply(s$total, s$fuel, sum); bf <- round(100 * bf / sum(bf))
  log_msg(sprintf("%-10s %10.0f in %2d ecoregions | by fuel %%: %s", o, sum(s$total), length(unique(s$eco_l2_code)),
                  paste(names(bf), bf, sep = " ", collapse = ", ")))
}
log_msg(sprintf("=== done | %d rows in fire_weekly | %.1f min | %s | log: %s", nrow(FW),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
