#!/usr/bin/env Rscript
# Step 61 (fire activity study; P720): MTBS burned area boundaries -> large-fire perimeters on the
# 4 km grid, with ignition dates, region, herbaceous fuel share, and a completeness check against
# FPA-FOD (60).
#   1. Download mtbs_perimeter_data.zip once (config fire_study: mtbs: url) into
#      <data_root>/_shared/raw/fire/mtbs/ and unzip. A zip or .shp placed there by hand is used as
#      is. MTBS v12.0 (April 2025) covers 1984-2024; per the MTBS release schedule (checked Oct 17
#      2026) 2023 and 2024 are still being mapped, so recent years may be incomplete (step 5).
#   2. Field check (config mtbs: fields; names as in the MTBS burned area boundaries dataset).
#   3. Filters (logged): ignition year in mtbs: years -> not AK / HI / PR (state = first two
#      characters of Event_ID) -> perimeter touches a grid land cell. Incid_Type is kept as a
#      column and a flag: wildfire = Incid_Type in mtbs: wildfire_types (Wildfire, Wildland Fire
#      Use); prescribed fires and "Unknown" stay in the table but are not wildfires.
#   4. Perimeters -> grid: each perimeter (made valid, projected to the grid CRS) is overlaid on
#      the 4 km grid with exact cell fractions (terra::extract exact = TRUE); burned area per cell =
#      covered fraction x cell area. Per fire: burned area on the grid, the region holding most of
#      it, its share on strict herbaceous cells, and its fuel-share-weighted area (fuel share =
#      grass + shrub + pasture + herbaceous wetland of land, config fire_study: fuel_classes).
#   5. Checks: (a) MTBS polygon area vs BurnBndAc; (b) completeness by year: share of FPA-FOD
#      wildfires >= check_min_acres whose MTBS_ID is in this table (a drop in the last years =
#      MTBS not finished mapping them); (c) ignition date vs FPA-FOD discovery date for linked
#      fires.
# Outputs <run_dir>/fire_study/mtbs/:
#   mtbs_fires.parquet   one row per fire: event_id, irwin_id, name, type, wildfire, ig_date, year,
#                        doy, acres (BurnBndAc), state, lon / lat (centroid), area_km2 (on the grid),
#                        region (most burned area), region_share, herb_mask_km2, fuel_km2, pilot,
#                        fod_id / fod_discovery (linked FPA-FOD record, if any)
#   mtbs_cells.parquet   one row per fire x grid cell: event_id, cell, frac, km2, region, fuel_share
#   mtbs_fields.csv, mtbs_values.csv (Incid_Type, Map_Prog, Asmnt_Type counts), mtbs_filters.csv,
#   mtbs_completeness.csv (year x FPA-FOD fires >= check_min_acres: total, with MTBS_ID, found
#   here), mtbs_region_year.csv (region x year: wildfires, burned km2, fuel-weighted km2)
#
# Usage: Rscript scripts/61_fire_mtbs.R [profile] [force]
#   Needs 60's fod_fires table for the completeness check and link (skipped with a warning if
#   missing). RStudio: gsi_profile <- "full"; source("scripts/61_fire_mtbs.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
force <- "force" %in% args

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("61_fire_mtbs", profile)
notify_init(cfg)
fs <- cfg$fire_study; mc <- fs$mtbs
out_dir <- fire_dir(cfg, "mtbs"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stem_f <- file.path(out_dir, "mtbs_fires"); stem_c <- file.path(out_dir, "mtbs_cells")
terra::terraOptions(progress = 0)
t0 <- Sys.time()
yrs <- seq(as.integer(mc$years$start), as.integer(mc$years$end))
log_msg("=== 61_fire_mtbs | profile: ", profile, " | years ", min(yrs), "-", max(yrs))

if (!force && !is.na(fire_table_file(stem_f)) && !is.na(fire_table_file(stem_c))) {
  log_msg("tables exist (", fire_table_file(stem_f), "); nothing to do (add 'force' to rebuild)")
} else {
notify("61 MTBS start", sprintf("years %d-%d", min(yrs), max(yrs)), priority = 2)

# ---- 1. download and unzip -----------------------------------------------------------------------
cache <- fire_cache_dir(cfg, "mtbs"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)
find_shp <- function() list.files(cache, "\\.shp$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
shp <- find_shp()
if (!length(shp)) {
  zf <- file.path(cache, basename(mc$url))
  if (!file.exists(zf)) {
    log_msg("downloading ", mc$url)
    t1 <- Sys.time()
    fire_download(mc$url, zf, ua = cfg$prism$user_agent %||% "gsi-pheno-regions (R httr2)")
    log_msg(sprintf("  %.0f MB in %s", file.size(zf) / 2^20, fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))))
  } else log_msg("zip on disk: ", zf)
  utils::unzip(zf, exdir = cache)
  shp <- find_shp()
}
if (length(shp) > 1 && any(grepl("perim", basename(shp), ignore.case = TRUE))) shp <- shp[grepl("perim", basename(shp), ignore.case = TRUE)]
if (length(shp) != 1) stop("expected one perimeter .shp under ", cache, ", found: ", paste(basename(shp), collapse = ", "))
log_msg("shapefile ", basename(shp), " (", format(file.mtime(shp), "%Y-%m-%d"), ")")

# ---- 2. read and field check ---------------------------------------------------------------------
t1 <- Sys.time()
V <- terra::vect(shp)
need <- unlist(mc$fields)
# Field names change case between MTBS releases (v12.0 file: all lower case); match ignoring case
# and rename to the names used here.
nm <- names(V); hit <- match(tolower(nm), tolower(need))
if (any(!is.na(hit) & nm != need[hit])) {
  log_msg("  renaming ", sum(!is.na(hit) & nm != need[hit]), " fields to the expected case (", paste(utils::head(nm[!is.na(hit) & nm != need[hit]], 4), collapse = ", "), ", ...)")
  nm[!is.na(hit)] <- need[hit[!is.na(hit)]]; names(V) <- nm
}
have <- names(V)
FLD <- data.frame(field = union(need, have), required = union(need, have) %in% need, present = union(need, have) %in% have)
utils::write.csv(FLD, file.path(out_dir, "mtbs_fields.csv"), row.names = FALSE)
log_msg(sprintf("read %s perimeters in %s | CRS %s | fields: %s", format(nrow(V), big.mark = ","),
                fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs"))), terra::crs(V, describe = TRUE)$name,
                paste(have, collapse = " ")))
miss <- setdiff(need, have)
if (length(miss)) stop("required MTBS fields missing: ", paste(miss, collapse = ", "), " (see mtbs_fields.csv)")

A <- as.data.frame(V)
vc <- function(col) { t <- sort(table(A[[col]], useNA = "ifany"), decreasing = TRUE)
  data.frame(field = col, value = ifelse(is.na(names(t)), "<NA>", names(t)), n = as.integer(t)) }
VAL <- do.call(rbind, lapply(intersect(c("Incid_Type", "Map_Prog", "Asmnt_Type"), have), vc))
utils::write.csv(VAL, file.path(out_dir, "mtbs_values.csv"), row.names = FALSE)
for (col in unique(VAL$field)) {
  v <- VAL[VAL$field == col, ]
  log_msg(sprintf("  %s: %s", col, paste(sprintf("%s %s", v$value, format(v$n, big.mark = ",", trim = TRUE)), collapse = " | ")))
}
wtypes <- unlist(mc$wildfire_types)
bad <- setdiff(wtypes, VAL$value[VAL$field == "Incid_Type"])
if (length(bad)) log_warn("wildfire_types not found in Incid_Type: ", paste(bad, collapse = ", "))

# ---- 3. filters ----------------------------------------------------------------------------------
ig <- fire_parse_date(A$Ig_Date)
log_msg(sprintf("  Ig_Date class %s, parsed for %.2f %%, range %s to %s", class(A$Ig_Date)[1], 100 * mean(!is.na(ig)),
                format(min(ig, na.rm = TRUE)), format(max(ig, na.rm = TRUE))))
if (mean(!is.na(ig)) < 0.99) stop("Ig_Date parsed for only ", round(100 * mean(!is.na(ig)), 1), " % of perimeters; examples: ",
                                  paste(utils::head(A$Ig_Date[is.na(ig)], 5), collapse = ", "))
yy <- as.integer(format(ig, "%Y"))
log_msg("  perimeters per ignition year (all states, all types), last 6: ",
        paste(sprintf("%d %d", tail(sort(unique(yy)), 6), as.integer(table(yy)[as.character(tail(sort(unique(yy)), 6))])), collapse = ", "))
state <- toupper(substr(A$Event_ID, 1, 2))
FL <- list(data.frame(step = "all perimeters", records = nrow(A), removed = 0L))
step <- function(name, keep) {
  FL[[length(FL) + 1]] <<- data.frame(step = name, records = sum(keep), removed = length(keep) - sum(keep))
  log_msg(sprintf("  filter %-26s %8s left, %7s removed", name, format(sum(keep), big.mark = ","), format(length(keep) - sum(keep), big.mark = ",")))
  keep
}
keep <- step("study years", !is.na(yy) & yy %in% yrs)
V <- V[keep]; A <- A[keep, ]; ig <- ig[keep]; state <- state[keep]
keep <- step("not AK / HI / PR", !state %in% toupper(unlist(mc$exclude_states)))
V <- V[keep]; A <- A[keep, ]; ig <- ig[keep]; state <- state[keep]

# ---- 4. perimeters on the grid -------------------------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
if (!terra::same.crs(V, grid)) V <- terra::project(V, terra::crs(grid))
nv <- sum(!terra::is.valid(V))
if (nv) { log_msg("  making ", nv, " invalid perimeters valid"); V <- terra::makeValid(V) }
gv <- terra::values(grid, mat = FALSE)
csz <- terra::values(terra::cellSize(grid, unit = "km"), mat = FALSE)
reg <- fire_region_values(cfg, grid, fs$partition %||% "A2_groups", fs$k %||% 13)
log_msg("regions from ", reg$src)
fuel <- fire_fuel_share(cfg)
hmv <- terra::values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), mat = FALSE)
herb <- !is.na(hmv) & hmv > 0

chunk <- 500L; idx <- split(seq_len(nrow(V)), ceiling(seq_len(nrow(V)) / chunk))
t1 <- Sys.time(); CL <- vector("list", length(idx))
for (j in seq_along(idx)) {
  i <- idx[[j]]
  e <- terra::extract(grid, V[i], cells = TRUE, exact = TRUE, ID = TRUE)
  e <- e[!is.na(e$cell) & e$fraction > 0, c("ID", "cell", "fraction")]
  e$row <- i[e$ID]
  CL[[j]] <- e[, c("row", "cell", "fraction")]
  if (j %% 5 == 0 || j == length(idx))
    log_msg(sprintf("  overlay %d / %d perimeters (%s)", max(i), nrow(V), fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))))
}
C <- do.call(rbind, CL); rm(CL)
if (!nrow(C)) stop("no perimeter overlaps a grid land cell")
C <- C[!is.na(gv[C$cell]), ]                                   # land cells only
C$km2 <- C$fraction * csz[C$cell]
C$region <- as.integer(reg$v[C$cell]); C$fuel_share <- fuel[C$cell]; C$herb <- herb[C$cell]
on_grid <- seq_len(nrow(V)) %in% C$row
keep <- step("on a grid land cell", on_grid)
FLT <- do.call(rbind, FL)
utils::write.csv(FLT, file.path(out_dir, "mtbs_filters.csv"), row.names = FALSE)

# per-fire summaries from the cells
area <- tapply(C$km2, C$row, sum)
herb_km2 <- tapply(C$km2 * C$herb, C$row, sum)
fuel_km2 <- tapply(C$km2 * ifelse(is.na(C$fuel_share), 0, C$fuel_share), C$row, sum)
RA <- stats::aggregate(km2 ~ row + region, C, sum)
RA <- RA[order(RA$row, -RA$km2), ]; RA <- RA[!duplicated(RA$row), ]
rk <- as.character(seq_len(nrow(V)))
F1 <- data.frame(
  event_id = A$Event_ID, irwin_id = A$irwinID, name = A$Incid_Name, type = A$Incid_Type,
  wildfire = A$Incid_Type %in% wtypes, ig_date = ig, year = as.integer(format(ig, "%Y")),
  doy = as.integer(format(ig, "%j")), acres = as.numeric(A$BurnBndAc), state = state,
  lon = suppressWarnings(as.numeric(A$BurnBndLon)), lat = suppressWarnings(as.numeric(A$BurnBndLat)),
  map_prog = A$Map_Prog, asmnt_type = A$Asmnt_Type,
  area_km2 = as.numeric(area[rk]), herb_mask_km2 = as.numeric(herb_km2[rk]), fuel_km2 = as.numeric(fuel_km2[rk]),
  region = RA$region[match(seq_len(nrow(V)), RA$row)], stringsAsFactors = FALSE)
F1$region_share <- RA$km2[match(seq_len(nrow(V)), RA$row)] / F1$area_km2
F1 <- F1[keep, ]
F1$pilot <- F1$region %in% as.integer(unlist(fs$pilot_regions))
C$event_id <- A$Event_ID[C$row]
C <- C[, c("event_id", "cell", "fraction", "km2", "region", "fuel_share")]

# (a) polygon area vs reported acres
ratio <- F1$area_km2 / (F1$acres * 0.00404686)
log_msg(sprintf("  area on the grid / BurnBndAc: median %.3f, 5-95 %% %.2f-%.2f (fires partly off the grid or in water lower it)",
                stats::median(ratio, na.rm = TRUE), stats::quantile(ratio, 0.05, na.rm = TRUE), stats::quantile(ratio, 0.95, na.rm = TRUE)))
log_msg(sprintf("  types kept: %s | wildfires (%s): %s", paste(sprintf("%s %d", names(table(F1$type)), as.integer(table(F1$type))), collapse = ", "),
                paste(wtypes, collapse = " + "), format(sum(F1$wildfire), big.mark = ",")))

# ---- 5. FPA-FOD link and completeness -------------------------------------------------------------
fod_stem <- fire_dir(cfg, "fod", "fod_fires")
F1$fod_id <- NA_integer_; F1$fod_discovery <- as.Date(NA)
if (is.na(fire_table_file(fod_stem))) {
  log_warn("no FPA-FOD table (run 60 first): link and completeness check skipped")
} else {
  FOD <- fire_read_table(fod_stem)
  FOD$date <- as.Date(FOD$date)
  lk <- FOD[!is.na(FOD$mtbs_id) & nzchar(FOD$mtbs_id), c("fod_id", "mtbs_id", "date", "size_ac", "year", "debris")]
  lk <- lk[order(-lk$size_ac), ]; lk1 <- lk[!duplicated(lk$mtbs_id), ]
  m <- match(F1$event_id, lk1$mtbs_id)
  F1$fod_id <- lk1$fod_id[m]; F1$fod_discovery <- lk1$date[m]
  w <- F1$wildfire
  log_msg(sprintf("  MTBS wildfires with an FPA-FOD record: %.1f %% (%d of %d)", 100 * mean(!is.na(m[w])), sum(!is.na(m[w])), sum(w)))
  dd <- as.numeric(F1$ig_date - F1$fod_discovery)
  if (any(!is.na(dd))) log_msg(sprintf("  ignition - discovery date (linked): median %+.0f d, |d| <= 1 d %.0f %%, |d| > 7 d %.1f %%",
                                       stats::median(dd, na.rm = TRUE), 100 * mean(abs(dd) <= 1, na.rm = TRUE), 100 * mean(abs(dd) > 7, na.rm = TRUE)))
  big <- FOD[FOD$size_ac >= (mc$check_min_acres %||% 1000) & !FOD$debris, ]
  CP <- do.call(rbind, lapply(yrs, function(y) {
    b <- big[big$year == y, ]
    data.frame(year = y, fod_large = nrow(b), with_mtbs_id = sum(!is.na(b$mtbs_id) & nzchar(b$mtbs_id)),
               found_here = sum(b$mtbs_id %in% F1$event_id), mtbs_wildfires = sum(F1$wildfire & F1$year == y))
  }))
  CP$share_found <- round(CP$found_here / CP$fod_large, 3)
  utils::write.csv(CP, file.path(out_dir, "mtbs_completeness.csv"), row.names = FALSE)
  typical <- stats::median(CP$share_found[CP$year <= max(yrs) - 3], na.rm = TRUE)
  log_msg(sprintf("  completeness: FPA-FOD wildfires >= %g ac found in MTBS, median %.2f through %d; last years: %s",
                  mc$check_min_acres %||% 1000, typical, max(yrs) - 3,
                  paste(sprintf("%d %.2f (%d MTBS wildfires)", tail(CP$year, 4), tail(CP$share_found, 4), tail(CP$mtbs_wildfires, 4)), collapse = ", ")))
  low <- CP$year[which(CP$share_found < typical - 0.15)]
  if (length(low)) log_warn("  years with clearly lower MTBS coverage (likely not finished mapping): ", paste(low, collapse = ", "))
}

f_out <- fire_write_table(F1, stem_f); c_out <- fire_write_table(C, stem_c)
log_msg(sprintf("wrote %s fires -> %s; %s fire-cells -> %s", format(nrow(F1), big.mark = ","), basename(f_out),
                format(nrow(C), big.mark = ","), basename(c_out)))

# ---- 6. region x year summary ----------------------------------------------------------------------
W <- F1[F1$wildfire, ]
Cw <- C[C$event_id %in% W$event_id, ]; Cw$year <- W$year[match(Cw$event_id, W$event_id)]
Cw$fuel_km2 <- Cw$km2 * ifelse(is.na(Cw$fuel_share), 0, Cw$fuel_share)
RY <- if (nrow(Cw)) stats::aggregate(cbind(km2, fuel_km2) ~ region + year, Cw, sum) else
  data.frame(region = integer(), year = integer(), km2 = numeric(), fuel_km2 = numeric())
NF <- if (nrow(W)) stats::aggregate(list(wildfires = rep(1L, nrow(W))), W[, c("region", "year")], sum) else NULL
if (!is.null(NF)) RY <- merge(RY, NF, all = TRUE)
RY <- RY[order(RY$region, RY$year), ]
RY$km2 <- round(RY$km2, 1); RY$fuel_km2 <- round(RY$fuel_km2, 1)
utils::write.csv(RY, file.path(out_dir, "mtbs_region_year.csv"), row.names = FALSE)
P <- RY[RY$region %in% as.integer(unlist(fs$pilot_regions)), ]
if (nrow(P)) {
  S <- stats::aggregate(cbind(wildfires, km2, fuel_km2) ~ region, P, sum, na.action = stats::na.pass)
  log_msg("pilot regions, MTBS wildfires ", min(yrs), "-", max(yrs), " (by region holding most of each fire's area for counts; area by cell):")
  for (k in seq_len(nrow(S))) log_msg(sprintf("  region %2d: %5d fires | %8.0f km2 burned | %5.0f %% of it fuel-weighted",
                                              S$region[k], S$wildfires[k], S$km2[k], 100 * S$fuel_km2[k] / S$km2[k]))
}
msg <- sprintf("%s MTBS fires on the grid (%s wildfires), %s fire-cells; %s", format(nrow(F1), big.mark = ","),
               format(sum(F1$wildfire), big.mark = ","), format(nrow(C), big.mark = ","),
               fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("61 MTBS done", msg, priority = 3, tags = "white_check_mark")
}
