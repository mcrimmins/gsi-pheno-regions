#!/usr/bin/env Rscript
# Step 60 (fire activity study; P720): FPA-FOD 7th edition fire records -> one table of CONUS
# wildfires with cause, size, discovery date, region and herbaceous cover at the point.
#   1. Download the GeoPackage once (config fire_study: fod: url) into
#      <data_root>/_shared/raw/fire/fpa_fod/ and unzip it. A zip or .gpkg placed there by hand
#      is used as is.
#   2. Field check: every field the study needs must exist in the Fires layer (names verified
#      Oct 17 2026 against the 7th edition metadata; unchanged from the 6th). All field names
#      and the cause / size-class value counts are written out, and the configured
#      drop_general_causes values must occur in the data.
#   3. Filters (each logged with the records it removes): study years (config fod: years) ->
#      valid discovery date -> not AK / HI / PR -> point on a land cell of the analysis grid.
#      Nothing else is dropped here; the screens used later are FLAGS:
#        cause     lightning (NWCG_CAUSE_CLASSIFICATION "Natural"), human, unknown
#        c_plus    FIRE_SIZE >= min_size_ac (10 ac, size class C and up)
#        debris    NWCG_GENERAL_CAUSE in drop_general_causes ("Debris and open burning")
#        rx_name   fire name looks like a prescribed burn (FPA-FOD excludes planned ignitions
#                  except escapes that needed suppression, Short 2014; this counts what slips in)
#        pilot     region in fire_study: pilot_regions
#   4. Per fire: grid cell, region (A2 k13 map number), strict herbaceous mask and NLCD shares
#      of the 4 km cell (herb_share, open_herb_share, crop_share).
# Outputs <run_dir>/fire_study/fod/:
#   fod_fires.parquet (.rds without arrow)  one row per fire
#   fod_fields.csv          fields in the layer, type, required / present
#   fod_values.csv          value counts: cause class, general cause, size class, source type (all years)
#   fod_filters.csv         records left after each filter
#   fod_region_summary.csv  per region x cause: all fires, C+ fires, C+ without debris, debris
#                           share, rx-like names, median herb share of the cell
#   fod_region_year.csv     per region x year x cause: C+ wildfires without debris burning
# Discovery week is NOT set here (65 decides it); the table keeps date and day of year.
#
# Usage: Rscript scripts/60_fire_fod.R [profile] [force]
#   force = rebuild the table even if it exists (the download is never repeated).
#   RStudio: gsi_profile <- "full"; source("scripts/60_fire_fod.R")

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
log_file <- log_init("60_fire_fod", profile)
notify_init(cfg)
fs <- cfg$fire_study; fc <- fs$fod
out_dir <- fire_dir(cfg, "fod"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stem <- file.path(out_dir, "fod_fires")
terra::terraOptions(progress = 0)
t0 <- Sys.time()
yrs <- seq(as.integer(fc$years$start), as.integer(fc$years$end))
log_msg("=== 60_fire_fod | profile: ", profile, " | ", fc$edition, " | years ", min(yrs), "-", max(yrs))
if (!requireNamespace("arrow", quietly = TRUE)) log_warn("arrow not installed: table written as .rds (renv::install('arrow') for Parquet)")

if (!force && !is.na(fire_table_file(stem))) {
  log_msg("table exists (", fire_table_file(stem), "); nothing to do (add 'force' to rebuild)")
} else {
notify("60 FPA-FOD start", sprintf("%s, years %d-%d", fc$edition, min(yrs), max(yrs)), priority = 2)

# ---- 1. download and unzip -----------------------------------------------------------------------
cache <- fire_cache_dir(cfg, "fpa_fod"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)
gpkg <- list.files(cache, "\\.gpkg$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(gpkg)) {
  zf <- file.path(cache, basename(fc$url))
  if (!file.exists(zf)) {
    log_msg("downloading ", fc$url)
    t1 <- Sys.time()
    fire_download(fc$url, zf, ua = cfg$prism$user_agent %||% "gsi-pheno-regions (R httr2)")
    log_msg(sprintf("  %.0f MB in %s", file.size(zf) / 2^20, fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))))
  } else log_msg("zip on disk: ", zf)
  utils::unzip(zf, exdir = cache)
  gpkg <- list.files(cache, "\\.gpkg$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}
if (length(gpkg) != 1) stop("expected one .gpkg under ", cache, ", found ", length(gpkg))
if (!grepl(fc$edition, basename(gpkg), fixed = TRUE)) log_warn("GeoPackage name ", basename(gpkg), " does not contain ", fc$edition)
lay <- terra::vector_layers(gpkg)
layer <- lay[tolower(lay) == tolower(fc$layer %||% "Fires")]
if (length(layer) != 1) stop("no layer '", fc$layer, "' in ", basename(gpkg), " (layers: ", paste(lay, collapse = ", "), ")")
log_msg("GeoPackage ", basename(gpkg), " | layers: ", paste(lay, collapse = ", "), " | using ", layer)

# ---- 2. field check ------------------------------------------------------------------------------
head5 <- terra::vect(gpkg, layer = layer, query = sprintf('SELECT * FROM "%s" LIMIT 5', layer), what = "attributes")
have <- names(head5)
need <- unlist(fc$fields)
FLD <- data.frame(field = union(need, have), required = union(need, have) %in% need, present = union(need, have) %in% have,
                  type = vapply(union(need, have), function(f) if (f %in% have) class(head5[[f]])[1] else NA_character_, ""))
utils::write.csv(FLD, file.path(out_dir, "fod_fields.csv"), row.names = FALSE)
log_msg(sprintf("fields in %s: %d | %s", layer, length(have), paste(have, collapse = " ")))
miss <- setdiff(need, have)
if (length(miss)) stop("required FPA-FOD fields missing: ", paste(miss, collapse = ", "), " (see fod_fields.csv)")
log_msg("all ", length(need), " required fields present")

# ---- 3. read, value counts, filters --------------------------------------------------------------
t1 <- Sys.time()
q <- sprintf('SELECT %s FROM "%s"', paste(sprintf('"%s"', need), collapse = ", "), layer)
F0 <- as.data.frame(terra::vect(gpkg, layer = layer, query = q, what = "attributes"))
log_msg(sprintf("read %s records (%d-%d) in %s", format(nrow(F0), big.mark = ","), min(F0$FIRE_YEAR, na.rm = TRUE),
                max(F0$FIRE_YEAR, na.rm = TRUE), fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))))

vc <- function(col) { t <- sort(table(F0[[col]], useNA = "ifany"), decreasing = TRUE)
  data.frame(field = col, value = ifelse(is.na(names(t)), "<NA>", names(t)), n = as.integer(t)) }
VAL <- do.call(rbind, lapply(c("NWCG_CAUSE_CLASSIFICATION", "NWCG_GENERAL_CAUSE", "FIRE_SIZE_CLASS", "SOURCE_SYSTEM_TYPE"), vc))
utils::write.csv(VAL, file.path(out_dir, "fod_values.csv"), row.names = FALSE)
for (col in unique(VAL$field)) {
  v <- VAL[VAL$field == col, ]
  log_msg(sprintf("  %s: %s", col, paste(sprintf("%s %s", v$value, format(v$n, big.mark = ",", trim = TRUE)), collapse = " | ")))
}
drop_gc <- unlist(fc$drop_general_causes)
bad <- setdiff(drop_gc, VAL$value[VAL$field == "NWCG_GENERAL_CAUSE"])
if (length(bad)) stop("drop_general_causes not found in NWCG_GENERAL_CAUSE: ", paste(bad, collapse = ", "))
cls <- VAL$value[VAL$field == "NWCG_CAUSE_CLASSIFICATION"]
for (x in c(fc$lightning_class, fc$human_class)) if (!x %in% cls) stop("cause class '", x, "' not in NWCG_CAUSE_CLASSIFICATION")

FL <- list(); step <- function(name, keep) {
  FL[[length(FL) + 1]] <<- data.frame(step = name, records = sum(keep), removed = length(keep) - sum(keep))
  log_msg(sprintf("  filter %-28s %10s left, %9s removed", name, format(sum(keep), big.mark = ","), format(length(keep) - sum(keep), big.mark = ",")))
  keep
}
FL[[1]] <- data.frame(step = "all records", records = nrow(F0), removed = 0L)
F0 <- F0[step("study years", F0$FIRE_YEAR %in% yrs), ]

# discovery date: FIRE_YEAR + DISCOVERY_DOY is the reference; DISCOVERY_DATE must agree
d_doy <- as.Date(sprintf("%d-01-01", as.integer(F0$FIRE_YEAR))) + as.integer(F0$DISCOVERY_DOY) - 1L
d_txt <- fire_parse_date(F0$DISCOVERY_DATE)
both <- !is.na(d_doy) & !is.na(d_txt)
log_msg(sprintf("  DISCOVERY_DATE parsed for %.2f %% (class %s); agrees with FIRE_YEAR + DISCOVERY_DOY for %.3f %% of those",
                100 * mean(!is.na(d_txt)), class(F0$DISCOVERY_DATE)[1], 100 * mean(d_doy[both] == d_txt[both])))
if (any(both & d_doy != d_txt)) log_warn(sprintf("  %d records where DISCOVERY_DATE and the DOY differ: DOY used", sum(both & d_doy != d_txt)))
date <- ifelse(is.na(d_doy), d_txt, d_doy); class(date) <- "Date"
keep <- step("valid discovery date", !is.na(date) & as.integer(format(date, "%Y")) == F0$FIRE_YEAR)
F0 <- F0[keep, ]; date <- date[keep]
keep <- step("not AK / HI / PR", !toupper(F0$STATE) %in% toupper(unlist(fc$exclude_states)))
F0 <- F0[keep, ]; date <- date[keep]

# point -> analysis grid cell (PRISM 4 km land cells; dev = the dev box)
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
xy <- cbind(as.numeric(F0$LONGITUDE), as.numeric(F0$LATITUDE))
log_msg("grid CRS: ", terra::crs(grid, describe = TRUE)$name, " (fire points are NAD83 lon/lat)")
if (!terra::is.lonlat(grid)) xy <- terra::crds(terra::project(terra::vect(xy, crs = "EPSG:4269"), terra::crs(grid)))
cell <- terra::cellFromXY(grid, xy)
gv <- terra::values(grid, mat = FALSE)
keep <- step("on a grid land cell", !is.na(cell) & !is.na(gv[cell]))
F0 <- F0[keep, ]; date <- date[keep]; cell <- cell[keep]
FLT <- do.call(rbind, FL)
utils::write.csv(FLT, file.path(out_dir, "fod_filters.csv"), row.names = FALSE)

# ---- 4. flags and cell attributes ----------------------------------------------------------------
cause <- ifelse(F0$NWCG_CAUSE_CLASSIFICATION %in% fc$lightning_class, "lightning",
         ifelse(F0$NWCG_CAUSE_CLASSIFICATION %in% fc$human_class, "human", "unknown"))
size <- as.numeric(F0$FIRE_SIZE)
c_plus <- !is.na(size) & size >= fc$min_size_ac
cls_cplus <- toupper(F0$FIRE_SIZE_CLASS) %in% c("C", "D", "E", "F", "G")
if (any(c_plus != cls_cplus)) log_warn(sprintf("  FIRE_SIZE >= %g and size class C+ disagree for %d fires (FIRE_SIZE used)",
                                              fc$min_size_ac, sum(c_plus != cls_cplus)))
rx <- grepl(fc$rx_name_pattern, toupper(F0$FIRE_NAME), perl = TRUE)

reg <- fire_region_values(cfg, grid, fs$partition %||% "A2_groups", fs$k %||% 13)
log_msg("regions from ", reg$src)
hs <- terra::rast(file.path(static_dir(cfg), "herb_share.tif"))
hm <- terra::rast(file.path(static_dir(cfg), "herb_mask.tif"))
uc <- sort(unique(cell))
H <- terra::extract(hs[[c("herb_share", "open_herb_share", "crop_share")]], uc)
hmv <- terra::extract(hm, uc)[, 1]
ix <- match(cell, uc)
ll <- terra::xyFromCell(grid, cell)

FD <- data.frame(
  fod_id = F0$FOD_ID, fpa_id = F0$FPA_ID, year = as.integer(F0$FIRE_YEAR), date = date,
  doy = as.integer(format(date, "%j")), cause = cause, cause_class = F0$NWCG_CAUSE_CLASSIFICATION,
  general_cause = F0$NWCG_GENERAL_CAUSE, size_ac = size, size_class = F0$FIRE_SIZE_CLASS, c_plus = c_plus,
  debris = F0$NWCG_GENERAL_CAUSE %in% drop_gc, rx_name = rx,
  cont_date = fire_parse_date(F0$CONT_DATE), mtbs_id = F0$MTBS_ID, ics209_id = F0$ICS_209_PLUS_INCIDENT_JOIN_ID,
  lon = as.numeric(F0$LONGITUDE), lat = as.numeric(F0$LATITUDE), state = F0$STATE, fips = F0$FIPS_CODE,
  owner = F0$OWNER_DESCR, source_type = F0$SOURCE_SYSTEM_TYPE, agency = F0$NWCG_REPORTING_AGENCY,
  cell = as.integer(cell), cell_lon = ll[, 1], cell_lat = ll[, 2],
  region = as.integer(reg$v[cell]), herb_mask = as.integer(!is.na(hmv[ix]) & hmv[ix] > 0),
  herb_share = H$herb_share[ix], open_herb_share = H$open_herb_share[ix], crop_share = H$crop_share[ix],
  stringsAsFactors = FALSE)
FD$pilot <- FD$region %in% as.integer(unlist(fs$pilot_regions))
rm(F0); invisible(gc())
f_out <- fire_write_table(FD, stem)
log_msg(sprintf("wrote %s fires -> %s", format(nrow(FD), big.mark = ","), f_out))

# ---- 5. summaries --------------------------------------------------------------------------------
agg <- function(d) data.frame(fires = nrow(d), c_plus = sum(d$c_plus), c_plus_wildfire = sum(d$c_plus & !d$debris),
  debris_share_c_plus = if (sum(d$c_plus)) round(mean(d$debris[d$c_plus]), 3) else NA,
  rx_name = sum(d$rx_name), c_plus_herb_mask = sum(d$c_plus & !d$debris & d$herb_mask == 1),
  median_herb_share = round(stats::median(d$herb_share, na.rm = TRUE), 3))
grp <- split(FD, list(FD$region, FD$cause), drop = TRUE)
RS <- do.call(rbind, lapply(names(grp), function(n) cbind(region = grp[[n]]$region[1], cause = grp[[n]]$cause[1], agg(grp[[n]]))))
RS <- RS[order(RS$region, RS$cause), ]
RS$pilot <- RS$region %in% as.integer(unlist(fs$pilot_regions))
utils::write.csv(RS, file.path(out_dir, "fod_region_summary.csv"), row.names = FALSE)
W <- FD[FD$c_plus & !FD$debris, ]
RY <- if (nrow(W)) stats::aggregate(list(c_plus_wildfire = rep(1L, nrow(W))), W[, c("region", "year", "cause")], sum) else
  data.frame(region = integer(), year = integer(), cause = character(), c_plus_wildfire = integer())
utils::write.csv(RY[order(RY$region, RY$year, RY$cause), ], file.path(out_dir, "fod_region_year.csv"), row.names = FALSE)

P <- RS[RS$pilot, ]
log_msg("pilot regions, C+ wildfires (no debris burning) by cause, ", min(yrs), "-", max(yrs), ":")
for (r in sort(unique(P$region))) {
  p <- P[P$region == r, ]
  log_msg(sprintf("  region %2d: %s | debris share of C+ %s | rx-like names %d | C+ on herb-mask cells %d",
                  r, paste(sprintf("%s %d", p$cause, p$c_plus_wildfire), collapse = ", "),
                  paste(sprintf("%s %.2f", p$cause, p$debris_share_c_plus), collapse = "/"), sum(p$rx_name), sum(p$c_plus_herb_mask)))
}
msg <- sprintf("%s fires on the grid, %s C+ wildfires (%s in pilot regions); %s",
               format(nrow(FD), big.mark = ","), format(sum(FD$c_plus & !FD$debris), big.mark = ","),
               format(sum(FD$c_plus & !FD$debris & FD$pilot), big.mark = ","), fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("60 FPA-FOD done", msg, priority = 3, tags = "white_check_mark")
}
