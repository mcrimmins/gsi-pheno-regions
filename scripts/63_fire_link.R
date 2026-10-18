#!/usr/bin/env Rscript
# Step 63 (fire activity study; P720): link VIIRS active-fire detections (62) to reported fires,
# and flag the screens of decision 12 (non-wildfire burning). Nothing is dropped here; 65 uses the
# flags (`keep` is the primary screen).
#   1. MTBS perimeters (all types), buffered by mtbs_buffer_m (VIIRS 375 m pixel + geolocation):
#      a detection inside a perimeter whose ignition date is within [ig - 1, ig + mtbs_max_days]
#      links to it (nearest ignition date if several). MTBS type is kept: prescribed fires and
#      "Other" are not wildfires.
#   2. Otherwise FPA-FOD fires (all sizes, debris burning included so it can be flagged): a
#      detection within radius r = fod_base_km + the fire's equivalent-circle radius
#      (sqrt(acres / pi)) of the fire's point, and with solar date in [discovery - 1,
#      end + 1], links to the nearest such fire. end = containment date when given (within a
#      year), else discovery + the 90th-percentile duration of the fire's size class
#      (fod_default_days, from the 7th edition: A-B 1, C 2, D 6, E 15, F 47, G 84 days).
#   3. Cropland: crop share of the 4 km cell >= crop_max.
#   Flags per detection: link (mtbs / fod / none), mtbs_id, mtbs_type, fod_id, dist_km,
#   wildfire (linked to an MTBS wildfire type, or to an FPA-FOD record that is not debris
#   burning), crop, keep = wildfire & !crop (decision 12: drop cropland and unlinked detections).
#   The burn-season county-week screen (sensitivity) needs county polygons and is not built yet.
# Outputs <run_dir>/fire_study/link/:
#   <SENSOR>/link_<SENSOR>_<year>.parquet   the 62 detections plus the flags above
#   link_summary.csv   sensor x year x region: detections, linked to MTBS (wildfire / other), to
#                      FPA-FOD (wildfire / debris), unlinked, crop, kept; share removed by each
#                      screen (pilot regions and all CONUS)
#
# Usage: Rscript scripts/63_fire_link.R [profile] [sensors] [years] [force]
#   Needs 60 (fod_fires), 61 (MTBS shapefile in the cache, mtbs_fires) and 62 (FIRMS tables).
#   Forked workers by sensor-year on Linux (config workers_by_step: fire_link), sequential on Windows.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
force <- "force" %in% args
rest <- setdiff(args[-1], "force")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("63_fire_link", profile)
notify_init(cfg)
fs <- cfg$fire_study; lc <- fs$link %||% list(); fc <- fs$firms
sensors <- c(fc$primary, fc$check)
if (length(rest) && grepl("VIIRS", rest[1])) { sensors <- strsplit(rest[1], ",")[[1]]; rest <- rest[-1] }
yrs <- seq(as.integer(fc$years$start), as.integer(fc$years$end))
if (length(rest)) { r <- as.integer(strsplit(rest[1], "-")[[1]]); yrs <- seq(r[1], r[length(r)]) }
out_dir <- fire_dir(cfg, "link"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
buf_m <- as.numeric(lc$mtbs_buffer_m %||% 750); mtbs_max <- as.integer(lc$mtbs_max_days %||% 180)
base_km <- as.numeric(lc$fod_base_km %||% 1.5); crop_max <- as.numeric(lc$crop_max %||% 0.5)
dflt <- unlist(lc$fod_default_days %||% list(A = 1, B = 1, C = 2, D = 6, E = 15, F = 47, G = 84))
log_msg("=== 63_fire_link | profile: ", profile, " | sensors ", paste(sensors, collapse = ", "), " | years ", min(yrs), "-", max(yrs),
        sprintf(" | MTBS buffer %g m, max %d d | FPA-FOD radius %g km + fire radius | crop >= %g", buf_m, mtbs_max, base_km, crop_max))

# ---- inputs shared by all years ---------------------------------------------------------------------
FOD <- fire_read_table(fire_dir(cfg, "fod", "fod_fires"))
FOD$date <- as.Date(FOD$date); FOD$cont_date <- as.Date(FOD$cont_date)
dur <- as.numeric(FOD$cont_date - FOD$date)
dd <- dflt[toupper(FOD$size_class)]; dd[is.na(dd)] <- dflt[["A"]]
FOD$end <- FOD$date + ifelse(!is.na(dur) & dur >= 0 & dur <= 365, dur, dd)
FOD$r_km <- base_km + sqrt(pmax(FOD$size_ac, 0) * 0.00404686 / pi)
FOD <- FOD[, c("fod_id", "date", "end", "lat", "lon", "r_km", "debris", "year", "size_ac")]
log_msg(sprintf("FPA-FOD: %s fires; search radius median %.2f km, 99th pct %.1f km; window median %d d",
                format(nrow(FOD), big.mark = ","), stats::median(FOD$r_km), stats::quantile(FOD$r_km, 0.99),
                as.integer(stats::median(FOD$end - FOD$date))))

MF <- fire_read_table(fire_dir(cfg, "mtbs", "mtbs_fires"))
MF$ig_date <- as.Date(MF$ig_date)
shp <- list.files(fire_cache_dir(cfg, "mtbs"), "\\.shp$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(shp) > 1 && any(grepl("perim", basename(shp), ignore.case = TRUE))) shp <- shp[grepl("perim", basename(shp), ignore.case = TRUE)]
if (length(shp) != 1) stop("expected one MTBS perimeter shapefile under ", fire_cache_dir(cfg, "mtbs"), " (run 61)")
P <- terra::vect(shp)
nm <- names(P); i <- match(tolower(nm), "event_id"); names(P)[which(!is.na(i))] <- "Event_ID"
P <- P[P$Event_ID %in% MF$event_id]
P <- terra::makeValid(terra::project(P, "EPSG:5070"))
P <- terra::buffer(P, buf_m)
P$ig_date <- MF$ig_date[match(P$Event_ID, MF$event_id)]
P$type <- MF$type[match(P$Event_ID, MF$event_id)]
P$wildfire <- MF$wildfire[match(P$Event_ID, MF$event_id)]
pf <- file.path(out_dir, sprintf("mtbs_buffered_%s.rds", profile))      # wrapped for forked workers
saveRDS(terra::wrap(P), pf); rm(P)
log_msg(sprintf("MTBS: %s perimeters buffered by %g m", format(nrow(MF), big.mark = ","), buf_m))

# ---- one sensor-year ----------------------------------------------------------------------------------
link_year <- function(sen, y) {
  inp <- fire_dir(cfg, "firms", sen, sprintf("firms_%s_%d", sen, y))
  of <- file.path(out_dir, sen, sprintf("link_%s_%d", sen, y))
  if (is.na(fire_table_file(inp))) return(list(msg = sprintf("%s %d: no 62 table", sen, y)))
  if (!force && !is.na(fire_table_file(of))) return(list(msg = sprintf("%s %d: exists", sen, y), sum = NULL))
  t1 <- Sys.time()
  X <- fire_read_table(inp); n <- nrow(X)
  X$date_solar <- as.Date(X$date_solar)
  X$link <- "none"; X$mtbs_id <- NA_character_; X$mtbs_type <- NA_character_; X$fod_id <- NA_integer_; X$dist_km <- NA_real_
  # 1. MTBS: points in buffered perimeters with a date match
  P <- terra::unwrap(readRDS(pf))
  P <- P[!is.na(P$ig_date) & P$ig_date >= as.Date(sprintf("%d-01-01", y)) - mtbs_max - 1 & P$ig_date <= as.Date(sprintf("%d-12-31", y)) + 1]
  if (nrow(P) && n) {
    pts <- terra::project(terra::vect(cbind(X$lon, X$lat), crs = "EPSG:4269"), "EPSG:5070")
    rel <- terra::relate(pts, P, "intersects", pairs = TRUE)
    if (length(rel) && nrow(rel)) {
      rel <- as.data.frame(rel); names(rel) <- c("pt", "poly")
      lag <- as.numeric(X$date_solar[rel$pt] - P$ig_date[rel$poly])
      rel <- rel[lag >= -1 & lag <= mtbs_max, , drop = FALSE]; lag <- lag[lag >= -1 & lag <= mtbs_max]
      if (nrow(rel)) {
        rel <- rel[order(rel$pt, abs(lag)), ]; rel <- rel[!duplicated(rel$pt), ]
        X$link[rel$pt] <- "mtbs"; X$mtbs_id[rel$pt] <- P$Event_ID[rel$poly]; X$mtbs_type[rel$pt] <- P$type[rel$poly]
      }
    }
  }
  # 2. FPA-FOD: candidate fires by grid cell neighbourhood, then exact distance and date window
  F1 <- FOD[FOD$end >= as.Date(sprintf("%d-01-01", y)) - 1 & FOD$date <= as.Date(sprintf("%d-12-31", y)) + 1, ]
  todo <- which(X$link == "none")
  if (nrow(F1) && length(todo)) {
    deg_lat <- F1$r_km / 111.2; deg_lon <- F1$r_km / (111.2 * cos(F1$lat * pi / 180))
    key <- function(lon, lat) paste(floor(lon * 10), floor(lat * 10))    # 0.1-degree buckets
    # expand each fire to the 0.1-degree buckets its search circle touches
    nb <- pmax(ceiling(pmax(deg_lat, deg_lon) * 10), 1)
    E <- do.call(rbind, lapply(sort(unique(nb)), function(k) {
      j <- which(nb == k); off <- expand.grid(dx = -k:k, dy = -k:k)
      data.frame(f = rep(j, each = nrow(off)),
                 key = paste(rep(floor(F1$lon[j] * 10), each = nrow(off)) + rep(off$dx, length(j)),
                             rep(floor(F1$lat[j] * 10), each = nrow(off)) + rep(off$dy, length(j))))
    }))
    cand <- split(E$f, E$key)
    xk <- key(X$lon[todo], X$lat[todo])
    for (kk in intersect(unique(xk), names(cand))) {
      ip <- todo[xk == kk]; jf <- cand[[kk]]
      dlat <- outer(X$lat[ip], F1$lat[jf], "-") * 111.2
      dlon <- outer(X$lon[ip], F1$lon[jf], "-") * 111.2 * cos(X$lat[ip] * pi / 180)
      dkm <- sqrt(dlat^2 + dlon^2)
      ok <- dkm <= matrix(F1$r_km[jf], length(ip), length(jf), byrow = TRUE) &
        outer(X$date_solar[ip], F1$date[jf] - 1, ">=") & outer(X$date_solar[ip], F1$end[jf] + 1, "<=")
      dkm[!ok] <- Inf
      best <- max.col(-dkm, ties.method = "first"); bd <- dkm[cbind(seq_along(ip), best)]
      hit <- is.finite(bd)
      if (any(hit)) {
        X$link[ip[hit]] <- "fod"; X$fod_id[ip[hit]] <- F1$fod_id[jf[best[hit]]]; X$dist_km[ip[hit]] <- round(bd[hit], 2)
      }
    }
  }
  deb <- FOD$debris[match(X$fod_id, FOD$fod_id)]
  X$wildfire <- (X$link == "mtbs" & X$mtbs_type %in% unlist(fs$mtbs$wildfire_types)) | (X$link == "fod" & !deb %in% TRUE)
  X$crop <- !is.na(X$crop_share) & X$crop_share >= crop_max
  X$keep <- X$wildfire & !X$crop
  fire_write_table(X, of)
  grp <- function(sub, region) data.frame(sensor = sen, year = y, region = region, detections = nrow(sub),
    mtbs_wildfire = sum(sub$link == "mtbs" & sub$wildfire), mtbs_other = sum(sub$link == "mtbs" & !sub$wildfire),
    fod_wildfire = sum(sub$link == "fod" & sub$wildfire), fod_debris = sum(sub$link == "fod" & !sub$wildfire),
    unlinked = sum(sub$link == "none"), crop = sum(sub$crop), crop_wildfire = sum(sub$crop & sub$wildfire), kept = sum(sub$keep))
  S <- rbind(grp(X, "all"), grp(X[X$pilot, ], "pilot"),
             do.call(rbind, lapply(sort(unique(X$region[X$pilot])), function(r) grp(X[X$region %in% r, ], as.character(r)))))
  list(msg = sprintf("%s %d: %s detections | MTBS %.0f %% | FPA-FOD %.0f %% | unlinked %.0f %% | crop %.0f %% | kept %.0f %% (%s)",
                     sen, y, format(n, big.mark = ","), 100 * mean(X$link == "mtbs"), 100 * mean(X$link == "fod"),
                     100 * mean(X$link == "none"), 100 * mean(X$crop), 100 * mean(X$keep),
                     fmt_dur(as.numeric(difftime(Sys.time(), t1, units = "secs")))), sum = S)
}

jobs <- expand.grid(y = yrs, sen = sensors, stringsAsFactors = FALSE)
w <- step_workers(cfg, "fire_link")$n
log_msg(sprintf("%d sensor-years, %d workers", nrow(jobs), w))
run1 <- function(i) tryCatch(link_year(jobs$sen[i], jobs$y[i]), error = function(e) list(msg = sprintf("%s %d: ERROR %s", jobs$sen[i], jobs$y[i], conditionMessage(e))))
res <- if (.Platform$OS.type != "windows" && w > 1) parallel::mclapply(seq_len(nrow(jobs)), run1, mc.cores = w, mc.preschedule = FALSE) else lapply(seq_len(nrow(jobs)), run1)
for (r in res) if (grepl("ERROR", r$msg)) log_err(r$msg) else log_msg("  ", r$msg)
S <- do.call(rbind, lapply(res, `[[`, "sum"))
if (!is.null(S) && nrow(S)) {
  f <- file.path(out_dir, "link_summary.csv")
  old <- if (file.exists(f)) utils::read.csv(f, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[!paste(old$sensor, old$year) %in% paste(S$sensor, S$year), ]
  S <- rbind(old, S); S <- S[order(S$sensor, S$year, S$region), ]
  utils::write.csv(S, f, row.names = FALSE)
  P <- S[S$region == "pilot" & S$sensor == fc$primary, ]
  if (nrow(P)) {
    tot <- colSums(P[, c("detections", "mtbs_wildfire", "mtbs_other", "fod_wildfire", "fod_debris", "unlinked", "crop", "kept")])
    log_msg(sprintf("pilot regions, %s %d-%d: %s detections | MTBS wildfire %.0f %% | MTBS prescribed / other %.0f %% | FPA-FOD wildfire %.0f %% | FPA-FOD debris %.0f %% | unlinked %.0f %% | crop cells %.0f %% | kept %.0f %%",
                    sub("VIIRS_", "", fc$primary), min(P$year), max(P$year), format(tot[["detections"]], big.mark = ","),
                    100 * tot[["mtbs_wildfire"]] / tot[["detections"]], 100 * tot[["mtbs_other"]] / tot[["detections"]],
                    100 * tot[["fod_wildfire"]] / tot[["detections"]], 100 * tot[["fod_debris"]] / tot[["detections"]],
                    100 * tot[["unlinked"]] / tot[["detections"]], 100 * tot[["crop"]] / tot[["detections"]], 100 * tot[["kept"]] / tot[["detections"]]))
  }
}
unlink(pf)
msg <- sprintf("%d sensor-years linked; %s", sum(!vapply(res, function(r) is.null(r$sum), TRUE)), fmt_dur(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("63 link done", msg, priority = 3, tags = "white_check_mark")
