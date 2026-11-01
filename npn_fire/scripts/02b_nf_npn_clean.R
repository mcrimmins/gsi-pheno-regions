#!/usr/bin/env Rscript
# Step 02b (NPN x fire study; P720, minutes per year; arrow, sf): clean the NPN status records from
# 01 into one row per plant x phenophase x year with onset / end dates and intensity-threshold
# dates, plus site, site-year and ecoregion effort tables.
# Inputs: _shared/raw/npn/status/<Y>/*.parquet (01), lookups (01: phenophases with pheno_class_id),
# npn_fire/data/nf_species_groups.csv (02a), EPA / CEC Level II polygons (config ecoregions).
#
# Per year (resumable: <out>/clean/years/<Y>_*.parquet / .csv are skipped if present):
#   - plant phenophase classes only (pheno_class_id 1-13); status -1 (uncertain) dropped; repeat
#     visits on one day merged (yes if any observer said yes; intensity = highest)
#   - per individual x phenophase x year (`plant_years`):
#       n_obs, n_yes, first_obs_doy, last_obs_doy
#       onset_doy = first yes; onset_gap = days since the previous "no" (NA: no earlier no that year)
#       end_doy = last yes; end_gap = days to the next "no" (NA: no later no that year)
#       n_cycles = runs of yes separated by a no; last_onset_doy / last_onset_gap = start of the last run
#       percent-class phenophases (intensity like "25-49%"; a "no" counts as 0 %): n_pct, peak_pct,
#       peak_doy, up50_doy / up50_gap (first reading >= threshold; gap from the previous reading
#       below it), down50_doy / down50_gap (first reading below the threshold after the peak; gap
#       from the previous reading at or above it). Threshold: config clean: pct_threshold (50).
#     Gaps are stored, not filtered: 05+ apply the 14-day rule (7 / 30 sensitivity) when reading.
#   - per site x year (`site_years`): n_dates, n_plants, n_species, n_observers, n_records, share of
#     records from NEON (dataset_id 6 / 16 or partner group NEON), main partner-group-free source
#   - sites: median lat / lon / elevation, state, first / last year
# Then: Level II ecoregion per site (point in polygon; sites outside every polygon take the nearest,
# flagged) and effort by ecoregion x year x group x pheno class (sites, plants, plant-years with an
# onset with gap <= 14 d, NEON share).
# Outputs <run_dir>/npn_fire/clean/: plant_years.parquet, site_years.csv, sites.csv (coordinates:
# never published), effort_eco.csv, phenophase_map.csv, intensity_map.csv, clean_summary.csv.
#
# Usage: Rscript npn_fire/scripts/02b_nf_npn_clean.R [profile] [years] [force]
# Needs: arrow, sf. Run 02a (laptop) first and git pull its data/nf_species_groups.csv.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && !grepl("^[0-9]|^force$", args[1])) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
rest <- if (length(args) && identical(args[1], profile)) args[-1] else args
force <- "force" %in% rest; rest <- setdiff(rest, "force")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("sf", quietly = TRUE)) stop("needs the sf package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("02b_nf_npn_clean", profile)
nf <- cfg$npn_fire; cc <- nf$clean
years <- nf_parse_years(rest[grepl("^[0-9]", rest)][1]) %||% seq(nf$years$start, nf$years$end)
thr <- cc$pct_threshold %||% 50; gap_ok <- cc$gap_days %||% 14
out <- nf_out(cfg, "clean"); ydir <- nf_out(cfg, "clean", "years")
log_msg("=== 02b_nf_npn_clean | profile: ", profile, " | years ", paste(range(years), collapse = "-"),
        " | threshold ", thr, " % | force: ", force)

# ---- lookups -------------------------------------------------------------------------------------
lk <- file.path(cfg$raw_dir, "npn", "lookups")
PH <- utils::read.csv(file.path(lk, "phenophases.csv"), stringsAsFactors = FALSE)
PC <- utils::read.csv(file.path(lk, "pheno_classes.csv"), stringsAsFactors = FALSE)
gf <- nf_path("data", "nf_species_groups.csv")
if (!file.exists(gf)) stop("no ", gf, ": run 02a on the laptop, commit, and git pull here")
SG <- utils::read.csv(gf, stringsAsFactors = FALSE)
plant_classes <- unlist(cc$pheno_classes %||% 1:13)

# ---- per-year cleaning ------------------------------------------------------------------------------
read_year <- function(Y) {
  fs <- list.files(file.path(cfg$raw_dir, "npn", "status", Y), "^npn_status_.*\\.(parquet|csv\\.gz)$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  keep <- c("observation_id", "site_id", "latitude", "longitude", "elevation_in_meters", "state", "species_id",
            "individual_id", "phenophase_id", "observation_date", "phenophase_status", "intensity_category_id",
            "intensity_value", "observedby_person_id", "partner_group", "dataset_id")
  do.call(rbind, lapply(fs, function(f) {
    d <- nf_read_table(f)
    for (k in setdiff(keep, names(d))) d[[k]] <- NA
    d <- d[, keep]
    d$intensity_value <- as.character(d$intensity_value); d$partner_group <- as.character(d$partner_group)
    d$state <- as.character(d$state); d$observation_date <- as.character(d$observation_date)
    d
  }))
}

# Lagged values within key runs (data sorted by key, then date)
lag1 <- function(x, key) { y <- c(NA, x[-length(x)]); y[c(TRUE, key[-1] != key[-length(key)])] <- NA; y }
lead1 <- function(x, key) { y <- c(x[-1], NA); y[c(key[-1] != key[-length(key)], TRUE)] <- NA; y }
# Row index of the first / last TRUE of `mask` within each key (keys with no TRUE give no row)
first_in_key <- function(mask, key) { mask <- mask %in% TRUE; w <- which(mask); w[!duplicated(key[w])] }
last_in_key <- function(mask, key) { mask <- mask %in% TRUE; w <- which(mask); w[!duplicated(key[w], fromLast = TRUE)] }

clean_year <- function(Y) {
  f_py <- file.path(ydir, sprintf("plant_years_%d", Y)); f_sy <- file.path(ydir, sprintf("site_years_%d.csv", Y))
  f_im <- file.path(ydir, sprintf("intensity_%d.csv", Y))
  if (!force && !is.na(nf_table_path(f_py)) && file.exists(f_sy)) { log_msg(Y, ": on disk, skipped"); return(invisible()) }
  t0 <- Sys.time()
  d <- read_year(Y)
  if (is.null(d) || !nrow(d)) { log_warn(Y, ": no status files"); return(invisible()) }
  n0 <- nrow(d)
  d$pheno_class_id <- PH$pheno_class_id[match(d$phenophase_id, PH$phenophase_id)]
  d <- d[d$pheno_class_id %in% plant_classes & d$species_id %in% SG$species_id, ]
  n1 <- nrow(d)
  d <- d[d$phenophase_status %in% c(0, 1), ]
  n2 <- nrow(d)
  d$date <- as.Date(d$observation_date); d$doy <- as.integer(format(d$date, "%j"))
  d$neon <- d$dataset_id %in% c(6, 16) | grepl("NEON|National Ecological Observatory", d$partner_group %||% "")

  # site-years (before merging repeat visits)
  sk <- d$site_id
  SY <- data.frame(year = Y, site_id = sort(unique(sk)))
  agg_n <- function(v) tapply(v, sk, function(z) length(unique(z[!is.na(z)])))[as.character(SY$site_id)]
  SY$n_records <- as.integer(table(sk)[as.character(SY$site_id)])
  SY$n_dates <- as.integer(agg_n(d$date)); SY$n_plants <- as.integer(agg_n(d$individual_id))
  SY$n_species <- as.integer(agg_n(d$species_id)); SY$n_observers <- as.integer(agg_n(d$observedby_person_id))
  SY$neon_share <- round(as.numeric(tapply(d$neon, sk, mean)[as.character(SY$site_id)]), 3)

  # intensity value -> class type and midpoint (percent classes only)
  iv <- unique(d$intensity_value[!is.na(d$intensity_value)])
  IM <- data.frame(intensity_value = iv, type = nf_intensity_type(iv), pct_mid = nf_pct_mid(iv), stringsAsFactors = FALSE)
  utils::write.csv(IM, f_im, row.names = FALSE)
  pct_ph <- unique(d$phenophase_id[d$intensity_value %in% IM$intensity_value[IM$type == "percent"]])
  d$pct <- IM$pct_mid[match(d$intensity_value, IM$intensity_value)]

  # merge repeat visits: one row per individual x phenophase x date
  d$key <- paste(d$individual_id, d$phenophase_id, sep = "_")
  d <- d[order(d$key, d$date, -d$phenophase_status, -ifelse(is.na(d$pct), -1, d$pct)), ]
  d <- d[!duplicated(paste(d$key, d$date)), ]        # first row = yes if any yes, highest intensity
  is_pct <- d$phenophase_id %in% pct_ph
  d$pct[is_pct & d$phenophase_status == 0] <- 0      # "no" = 0 % for percent classes
  k <- d$key; s <- d$phenophase_status; dt <- as.integer(d$date)
  ps <- lag1(s, k); pd <- lag1(dt, k); ns <- lead1(s, k); nd <- lead1(dt, k)
  first_row <- !duplicated(k); last_row <- !duplicated(k, fromLast = TRUE)
  ki <- match(k, unique(k)); nk <- max(ki)
  yes <- s == 1
  # first / last yes per key; start of the last yes-run
  fy <- first_in_key(yes, k); ly <- last_in_key(yes, k)
  cyc_start <- yes & (is.na(ps) | ps == 0)
  lc <- last_in_key(cyc_start, k)
  P <- data.frame(year = Y, key = unique(k), stringsAsFactors = FALSE)
  P$individual_id <- d$individual_id[first_row]; P$site_id <- d$site_id[first_row]
  P$species_id <- d$species_id[first_row]; P$phenophase_id <- d$phenophase_id[first_row]
  P$pheno_class_id <- d$pheno_class_id[first_row]
  P$n_obs <- tabulate(ki, nk); P$n_yes <- tabulate(ki[yes], nk)
  P$first_obs_doy <- d$doy[first_row]; P$last_obs_doy <- d$doy[last_row]
  P$onset_doy <- NA_integer_; P$onset_gap <- NA_integer_; P$end_doy <- NA_integer_; P$end_gap <- NA_integer_
  P$onset_doy[ki[fy]] <- d$doy[fy]; P$onset_gap[ki[fy]] <- ifelse(ps[fy] %in% 0, dt[fy] - pd[fy], NA)
  P$end_doy[ki[ly]] <- d$doy[ly]; P$end_gap[ki[ly]] <- ifelse(ns[ly] %in% 0, nd[ly] - dt[ly], NA)
  P$n_cycles <- tabulate(ki[cyc_start], nk)
  P$last_onset_doy <- NA_integer_; P$last_onset_gap <- NA_integer_
  P$last_onset_doy[ki[lc]] <- d$doy[lc]; P$last_onset_gap[ki[lc]] <- ifelse(ps[lc] %in% 0, dt[lc] - pd[lc], NA)

  # percent-class thresholds
  P$pct_class <- P$phenophase_id %in% pct_ph
  P$n_pct <- 0L; P$peak_pct <- NA_real_; P$peak_doy <- NA_integer_
  P$up50_doy <- NA_integer_; P$up50_gap <- NA_integer_; P$down50_doy <- NA_integer_; P$down50_gap <- NA_integer_
  r <- which(is_pct & !is.na(d$pct))
  if (length(r)) {
    kr <- ki[r]; v <- d$pct[r]; dd <- dt[r]; doy <- d$doy[r]
    pv <- lag1(v, kr); pdd <- lag1(dd, kr)
    P$n_pct[unique(kr)] <- tabulate(kr, nk)[unique(kr)]
    mx <- as.numeric(tapply(v, kr, max)); ukr <- sort(unique(kr))
    P$peak_pct[ukr] <- mx
    at_peak <- v == mx[match(kr, ukr)]
    fp <- first_in_key(at_peak, kr)
    P$peak_doy[kr[fp]] <- doy[fp]
    up <- v >= thr
    fu <- first_in_key(up, kr)
    P$up50_doy[kr[fu]] <- doy[fu]; P$up50_gap[kr[fu]] <- ifelse(!is.na(pv[fu]) & pv[fu] < thr, dd[fu] - pdd[fu], NA)
    peak_date <- rep(NA_integer_, nk); peak_date[kr[fp]] <- dd[fp]
    after <- dd > peak_date[kr] & v < thr & P$peak_pct[kr] >= thr
    fd <- first_in_key(after, kr)
    P$down50_doy[kr[fd]] <- doy[fd]; P$down50_gap[kr[fd]] <- ifelse(!is.na(pv[fd]) & pv[fd] >= thr, dd[fd] - pdd[fd], NA)
  }
  P$key <- NULL
  SG_i <- match(P$species_id, SG$species_id); P$group <- SG$group[SG_i]
  invisible(nf_write_table(P, f_py))
  utils::write.csv(SY, f_sy, row.names = FALSE)

  # sites seen this year
  S <- stats::aggregate(cbind(latitude, longitude, elevation_in_meters) ~ site_id, data = d, FUN = stats::median, na.action = na.pass)
  S$state <- d$state[match(S$site_id, d$site_id)]; S$year <- Y
  utils::write.csv(S, file.path(ydir, sprintf("sites_%d.csv", Y)), row.names = FALSE)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  log_msg(sprintf("%d: %9d records -> %9d plant phenophase -> %9d yes/no | %7d plant-years (%d with onset gap <= %d d) | %4d sites | %s s",
                  Y, n0, n1, n2, nrow(P), sum(P$onset_gap <= gap_ok, na.rm = TRUE), gap_ok, nrow(SY), secs))
}
for (Y in years) clean_year(Y)

# ---- combine ---------------------------------------------------------------------------------------
all_years <- seq(nf$years$start, nf$years$end)
pyf <- stats::na.omit(vapply(all_years, function(Y) nf_table_path(file.path(ydir, sprintf("plant_years_%d", Y))), ""))
if (!length(pyf)) stop("no yearly tables in ", ydir)
P <- do.call(rbind, lapply(pyf, nf_read_table))
P$group <- SG$group[match(P$species_id, SG$species_id)]     # always the current 02a groups (yearly files may be older)
SY <- do.call(rbind, lapply(list.files(ydir, "^site_years_\\d{4}\\.csv$", full.names = TRUE), utils::read.csv))
S <- do.call(rbind, lapply(list.files(ydir, "^sites_\\d{4}\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE))
IM <- unique(do.call(rbind, lapply(list.files(ydir, "^intensity_\\d{4}\\.csv$", full.names = TRUE), utils::read.csv, stringsAsFactors = FALSE)))
Sx <- stats::aggregate(cbind(latitude, longitude, elevation_in_meters) ~ site_id, data = S, FUN = stats::median, na.action = na.pass)
Sx$state <- S$state[match(Sx$site_id, S$site_id)]
Sx$first_year <- as.integer(tapply(S$year, S$site_id, min)[as.character(Sx$site_id)])
Sx$last_year <- as.integer(tapply(S$year, S$site_id, max)[as.character(Sx$site_id)])
Sx$n_years <- as.integer(tapply(S$year, S$site_id, function(z) length(unique(z)))[as.character(Sx$site_id)])

# ---- Level II ecoregions --------------------------------------------------------------------------
eco_shp <- function() {
  ec <- cfg$ecoregions; edir <- nf_raw(cfg, "ecoregions", paste0("l", ec$level))
  shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
  if (!length(shp)) {
    alt <- list.files(file.path(cfg$run_dir, "npn_study", "ecoregions"), "\\.shp$", recursive = TRUE, full.names = TRUE)
    if (length(alt)) shp <- alt
  }
  if (!length(shp)) {
    zf <- file.path(edir, basename(ec$url))
    utils::download.file(ec$url, zf, mode = "wb"); utils::unzip(zf, exdir = edir)
    shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
  }
  shp[grepl("l2|level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
}
shp <- eco_shp(); E <- sf::st_read(shp, quiet = TRUE)
g <- attr(E, "sf_column"); names(E)[names(E) != g] <- toupper(names(E)[names(E) != g])
E <- sf::st_make_valid(E); E <- E[!E$NA_L2CODE %in% c("0.0", "0"), ]
ok <- is.finite(Sx$latitude) & is.finite(Sx$longitude)
Pt <- sf::st_transform(sf::st_as_sf(Sx[ok, ], coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), sf::st_crs(E))
hit <- sf::st_intersects(Pt, E); i <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, 0L)
near <- is.na(i); if (any(near)) i[near] <- sf::st_nearest_feature(Pt[near, ], E)
dist <- rep(0, length(i)); if (any(near)) dist[near] <- as.numeric(sf::st_distance(Pt[near, ], E[i[near], ], by_element = TRUE)) / 1000
Sx$eco_l2_code <- NA_character_; Sx$eco_l2_name <- NA_character_; Sx$eco_l1_code <- NA_character_; Sx$eco_l1_name <- NA_character_
Sx$eco_nearest <- NA; Sx$eco_dist_km <- NA_real_
Sx$eco_l2_code[ok] <- as.character(E$NA_L2CODE[i]); Sx$eco_l2_name[ok] <- as.character(E$NA_L2NAME[i])
Sx$eco_l1_code[ok] <- as.character(E$NA_L1CODE[i]); Sx$eco_l1_name[ok] <- as.character(E$NA_L1NAME[i])
Sx$eco_nearest[ok] <- near; Sx$eco_dist_km[ok] <- round(dist, 1)
far <- Sx$eco_dist_km > (cc$eco_max_km %||% 25) & !is.na(Sx$eco_dist_km)
if (any(far)) log_warn(sum(far), " sites more than ", cc$eco_max_km %||% 25, " km outside every Level II polygon (kept, flagged)")
log_msg(sprintf("sites: %d | %d Level II ecoregions | %d by nearest polygon", nrow(Sx), length(unique(stats::na.omit(Sx$eco_l2_code))), sum(near)))

# ---- write combined tables and effort ---------------------------------------------------------------
P$eco_l2_code <- Sx$eco_l2_code[match(P$site_id, Sx$site_id)]
SY$eco_l2_code <- Sx$eco_l2_code[match(SY$site_id, Sx$site_id)]
invisible(nf_write_table(P, file.path(out, "plant_years")))
utils::write.csv(SY, file.path(out, "site_years.csv"), row.names = FALSE)
utils::write.csv(Sx, file.path(out, "sites.csv"), row.names = FALSE)
utils::write.csv(IM, file.path(out, "intensity_map.csv"), row.names = FALSE)
PM <- PH[PH$pheno_class_id %in% plant_classes, ]
PM$pheno_class_name <- PC$name[match(PM$pheno_class_id, PC$id)]
PM$plant_years <- as.integer(table(factor(P$phenophase_id, levels = PM$phenophase_id)))
PM$percent_class <- PM$phenophase_id %in% unique(P$phenophase_id[P$pct_class])
utils::write.csv(PM, file.path(out, "phenophase_map.csv"), row.names = FALSE)

ok_on <- !is.na(P$onset_gap) & P$onset_gap <= gap_ok
P$neon <- SY$neon_share[match(paste(P$site_id, P$year), paste(SY$site_id, SY$year))] > 0.5
grp <- list(eco = P$eco_l2_code, year = P$year, group = P$group, pheno_class_id = P$pheno_class_id)
EF <- stats::aggregate(list(plant_years = rep(1L, nrow(P))), grp, sum)
f2 <- function(v, fun) stats::aggregate(list(x = v), grp, fun)$x
EF$plants <- f2(P$individual_id, function(z) length(unique(z)))
EF$sites <- f2(P$site_id, function(z) length(unique(z)))
EF$onset_ok <- f2(ok_on, sum)
EF$sites_onset_ok <- stats::aggregate(list(x = ifelse(ok_on, P$site_id, NA)), grp, function(z) length(unique(z[!is.na(z)])))$x
EF$neon_plant_share <- round(f2(P$neon %in% TRUE, mean), 3)
names(EF)[1] <- "eco_l2_code"
EF$eco_l2_name <- Sx$eco_l2_name[match(EF$eco_l2_code, Sx$eco_l2_code)]
EF$pheno_class_name <- PC$name[match(EF$pheno_class_id, PC$id)]
utils::write.csv(EF[order(EF$eco_l2_code, EF$group, EF$pheno_class_id, EF$year), ], file.path(out, "effort_eco.csv"), row.names = FALSE)

SM <- stats::aggregate(list(plant_years = rep(1L, nrow(P))), list(year = P$year), sum)
SM$plants <- as.integer(tapply(P$individual_id, P$year, function(z) length(unique(z))))
SM$sites <- as.integer(tapply(P$site_id, P$year, function(z) length(unique(z))))
SM$onset_ok <- as.integer(tapply(ok_on, P$year, sum))
SM$pct_class_plant_years <- as.integer(tapply(P$pct_class, P$year, sum))
SM$down50_ok <- as.integer(tapply(!is.na(P$down50_gap) & P$down50_gap <= gap_ok, P$year, sum))
utils::write.csv(SM, file.path(out, "clean_summary.csv"), row.names = FALSE)
for (j in seq_len(nrow(SM))) with(SM[j, ], log_msg(sprintf("  %d: %8d plant-years | %6d plants | %5d sites | onset ok %7d | down50 ok %6d",
                                                       year, plant_years, plants, sites, onset_ok, down50_ok)))
tg <- table(P$group[ok_on]); log_msg("onset ok by group: ", paste(names(tg), tg, collapse = ", "))
log_msg("=== done | ", out, " | log: ", log_file)
