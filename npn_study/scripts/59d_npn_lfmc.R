#!/usr/bin/env Rscript
# Step 59d (NPN grass curing study; laptop, sf, < 2 min): where do NPN grass percent-green records
# and measured herbaceous live fuel moisture (LFMC) meet? LFMC: Globe-LFMC 2.0, U.S. records (the
# GSI study's lfmc_us_all.rds, which carries the National Fuel Moisture Database samples), herbaceous
# class (grasses and forbs), records flagged out of range / elevation / NFMD outlier dropped.
#   1. Inventory: herbaceous LFMC samples, sites, years and species; sampling interval; by EPA
#      Level II ecoregion next to U.S. grass fuel area (59c) and NPN plant-years.
#   2. Co-location: distance from every NPN study site (relaxed set, 58) to the nearest herbaceous
#      LFMC site; counts within 1 / 5 / 10 / 25 / 50 km; years in common.
#   3. Pairs: for NPN-LFMC site pairs within npn_study: lfmc: max_km (default 10) and the same
#      year, each LFMC sample is matched with the NPN readings at the paired site within +/- max_days
#      (default 7): the mean percent green of each NPN species (and of all its grasses) on the
#      nearest reading date. Same species / same genus flagged. Where there are enough pairs
#      (>= 20): Spearman r of percent green vs LFMC, overall and by pair; and by site-year, the LFMC
#      fall through 120 % (NFDRS herbaceous transition) after the peak vs the NPN 50 % curing date.
# Inputs: lfmc_us_all.rds, npn_green_obs.rds, npn_curing_metrics.rds (GSI_OBS_DIR or config
# study_dir); <run_dir>/npn_study/record/series_typical.csv (58); ecoregions (59);
# report/tables/npn_representativeness_ecoregion.csv (59c, optional).
# Outputs <run_dir>/npn_study/lfmc/: lfmc_sites.csv, colocation.csv, pairs.csv, pair_summary.csv,
# inventory_ecoregion.csv; report/tables/npn_lfmc_inventory.csv, npn_lfmc_colocation.csv,
# npn_lfmc_pairs.csv (summary only; no site names anywhere).
#
# Usage: Rscript npn_study/scripts/59d_npn_lfmc.R [profile]   (after 58, 59, 59c)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("sf", quietly = TRUE)) stop("needs sf")

cfg <- load_npn_cfg(profile)
log_file <- log_init("59d_npn_lfmc", profile)
ns <- cfg$npn_study %||% list(); lp <- ns$lfmc %||% list()
max_km <- as.numeric(lp$max_km %||% 10); max_days <- as.numeric(lp$max_days %||% 7)
base <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(base, "lfmc"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
log_msg(sprintf("=== 59d_npn_lfmc | profile: %s | pairs within %g km and +/- %g days", profile, max_km, max_days))
hav <- function(a1, o1, a2, o2) { r <- pi / 180; 6371 * 2 * asin(sqrt(sin((a2 - a1) * r / 2)^2 + cos(a1 * r) * cos(a2 * r) * sin((o2 - o1) * r / 2)^2)) }

L <- readRDS(file.path(npn_obs_dir(cfg), "lfmc_us_all.rds"))
n0 <- nrow(L)
L <- L[L$fuel_class %in% "herb" & !(L$flag_range %in% TRUE) & !(L$flag_elev %in% TRUE) & !(L$flag_nfmd_outlier %in% TRUE) & !is.na(L$lfmc), ]
L$grass <- L$ftype_orig %in% "Grass"; L$genus <- sub(" .*", "", L$species)
log_msg(sprintf("LFMC: %d U.S. records; herbaceous after flags %d (%d grass, %d forb) at %d sites, %d-%d", n0, nrow(L), sum(L$grass), sum(!L$grass),
                length(unique(L$site_key)), min(L$year), max(L$year)))

# ---- sites and ecoregions -------------------------------------------------------------------------
LS <- do.call(rbind, lapply(split(L, L$site_key), function(d) data.frame(site_key = d$site_key[1], lat = d$lat[1], lon = d$lon[1],
  samples = nrow(d), grass_samples = sum(d$grass), years = length(unique(d$year)), first_year = min(d$year), last_year = max(d$year),
  species = length(unique(d$species)),
  interval_days = { z <- unlist(lapply(split(d$date, d$year), function(x) diff(sort(unique(x))))); if (length(z)) stats::median(as.numeric(z)) else NA })))
LS$lfmc_id <- sprintf("H%03d", seq_len(nrow(LS)))
edir <- file.path(base, "ecoregions")
shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
shp <- shp[grepl("l2|level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
LS$region <- "none"
if (length(shp) && !is.na(shp)) {
  sf::sf_use_s2(FALSE)
  E <- sf::st_read(shp, quiet = TRUE); gc <- attr(E, "sf_column"); names(E)[names(E) != gc] <- toupper(names(E)[names(E) != gc]); E <- sf::st_make_valid(E)
  P <- sf::st_transform(sf::st_as_sf(LS, coords = c("lon", "lat"), crs = 4326), sf::st_crs(E))
  i <- vapply(sf::st_intersects(P, E), function(h) if (length(h)) h[1] else NA_integer_, 0L)
  if (any(is.na(i))) i[is.na(i)] <- sf::st_nearest_feature(P[is.na(i), ], E)
  LS$region <- as.character(E$NA_L2CODE[i])
} else log_warn("no ecoregion shapefile (run 59): LFMC sites get region 'none'")
L$region <- LS$region[match(L$site_key, LS$site_key)]; L$lfmc_id <- LS$lfmc_id[match(L$site_key, LS$site_key)]
utils::write.csv(LS[, setdiff(names(LS), "site_key")], file.path(odir, "lfmc_sites.csv"), row.names = FALSE)

# ---- the NPN record -----------------------------------------------------------------------------
TY <- utils::read.csv(file.path(base, "record", "series_typical.csv"), stringsAsFactors = FALSE, colClasses = c(region = "character"))
NS <- stats::aggregate(plant_years ~ site_id + lat + lon + region, TY, sum)
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
nyears <- tapply(cm$year[cm$site_id %in% NS$site_id], cm$site_id[cm$site_id %in% NS$site_id], function(z) sort(unique(z)))

# ---- 1. inventory by ecoregion --------------------------------------------------------------------
rg <- sort(unique(c(LS$region, NS$region))); rg <- rg[rg != "none"]
IV <- data.frame(region = rg,
  lfmc_sites = vapply(rg, function(k) sum(LS$region == k), 0), lfmc_samples = vapply(rg, function(k) sum(L$region == k), 0),
  lfmc_grass_samples = vapply(rg, function(k) sum(L$region == k & L$grass), 0),
  npn_sites = vapply(rg, function(k) sum(NS$region == k), 0), npn_plant_years = vapply(rg, function(k) sum(NS$plant_years[NS$region == k]), 0))
IV$lfmc_sample_share <- IV$lfmc_samples / sum(IV$lfmc_samples); IV$npn_plant_year_share <- IV$npn_plant_years / sum(IV$npn_plant_years)
rf <- npn_path("report", "tables", "npn_representativeness_ecoregion.csv")
if (file.exists(rf)) { R <- utils::read.csv(rf, stringsAsFactors = FALSE, colClasses = c(region = "character"))
  IV$grass_fuel_area_share <- R$grass_fuel_area_share[match(IV$region, R$region)]; IV$name <- R$name[match(IV$region, R$region)] } else { IV$grass_fuel_area_share <- NA; IV$name <- NA }
IV <- IV[order(-ifelse(is.na(IV$grass_fuel_area_share), 0, IV$grass_fuel_area_share)), ]

# ---- 2. co-location ----------------------------------------------------------------------------------
CO <- do.call(rbind, lapply(seq_len(nrow(NS)), function(i) {
  dk <- hav(NS$lat[i], NS$lon[i], LS$lat, LS$lon); j <- which.min(dk)
  ly <- sort(unique(L$year[L$site_key == LS$site_key[j]])); ny <- nyears[[as.character(NS$site_id[i])]]
  data.frame(site_id = NS$site_id[i], region = NS$region[i], plant_years = NS$plant_years[i], nearest_km = dk[j], nearest_lfmc = LS$lfmc_id[j],
             years_in_common = length(intersect(ly, ny)), lfmc_within_max = sum(dk <= max_km)) }))
cut_km <- c(1, 5, 10, 25, 50)
CK <- data.frame(within_km = cut_km, npn_sites = vapply(cut_km, function(k) sum(CO$nearest_km <= k), 0),
                 npn_sites_years_in_common = vapply(cut_km, function(k) sum(CO$nearest_km <= k & CO$years_in_common > 0), 0),
                 plant_years = vapply(cut_km, function(k) sum(CO$plant_years[CO$nearest_km <= k]), 0))
CK$share_of_sites <- CK$npn_sites / nrow(NS)
for (i in seq_len(nrow(CK))) with(CK[i, ], log_msg(sprintf("  NPN sites within %2g km of herbaceous LFMC: %3d of %d (%.0f %%), %d with years in common", within_km, npn_sites, nrow(NS), 100 * share_of_sites, npn_sites_years_in_common)))

# ---- 3. pairs ---------------------------------------------------------------------------------------------
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds")); go <- go[!is.na(go$pct_green) & go$site_id %in% NS$site_id, ]
go$date <- as.Date(go$observation_date)
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
go$sp_name <- paste(go$genus, go$species); go$pathway <- sp$pathway[match(go$species_id, sp$species_id)]
PAIRS <- list()
near <- CO[CO$nearest_km <= max_km, ]
for (i in seq_len(nrow(near))) {
  s <- near$site_id[i]; nlat <- NS$lat[NS$site_id == s]; nlon <- NS$lon[NS$site_id == s]
  ks <- LS$site_key[hav(nlat, nlon, LS$lat, LS$lon) <= max_km]
  g <- go[go$site_id == s, ]
  for (k in ks) {
    l <- L[L$site_key == k & L$year %in% unique(g$year), ]; if (!nrow(l)) next
    dkm <- hav(nlat, nlon, l$lat[1], l$lon[1])
    for (j in seq_len(nrow(l))) {
      gg <- g[abs(as.numeric(g$date - l$date[j])) <= max_days, ]; if (!nrow(gg)) next
      # nearest reading date per species, then mean percent green over that species' plants
      by_sp <- do.call(rbind, lapply(split(gg, gg$species_id), function(z) { z <- z[abs(as.numeric(z$date - l$date[j])) == min(abs(as.numeric(z$date - l$date[j]))), ]
        data.frame(species_id = z$species_id[1], sp_name = z$sp_name[1], genus = z$genus[1], pathway = z$pathway[1], pct_green = mean(z$pct_green), plants = length(unique(z$individual_id)),
                   lag_days = as.numeric(z$date[1] - l$date[j])) }))
      PAIRS[[length(PAIRS) + 1]] <- cbind(data.frame(site_id = s, lfmc_id = l$lfmc_id[j], km = dkm, date = l$date[j], year = l$year[j], lfmc = l$lfmc[j],
        lfmc_species = l$species[j], lfmc_grass = l$grass[j]), by_sp,
        same_species = by_sp$sp_name == l$species[j], same_genus = by_sp$genus == l$genus[j])
    }
  }
}
PA <- if (length(PAIRS)) do.call(rbind, PAIRS) else NULL
PS <- data.frame(stat = character(), value = numeric())
add <- function(k, v) PS <<- rbind(PS, data.frame(stat = k, value = v))
add("npn_sites", nrow(NS)); add("lfmc_herb_sites", nrow(LS)); add("lfmc_herb_samples", nrow(L)); add("lfmc_grass_samples", sum(L$grass))
add("npn_sites_within_max_km", nrow(near)); add("max_km", max_km); add("max_days", max_days)
if (!is.null(PA) && nrow(PA)) {
  utils::write.csv(PA, file.path(odir, "pairs.csv"), row.names = FALSE)
  add("pairs", nrow(PA)); add("pairs_same_species", sum(PA$same_species)); add("pairs_same_genus", sum(PA$same_genus))
  add("pair_sites", length(unique(PA$site_id))); add("lfmc_samples_paired", length(unique(paste(PA$lfmc_id, PA$date, PA$lfmc_species))))
  sm <- function(X, tag) if (nrow(X) >= 20) { r <- suppressWarnings(stats::cor(X$pct_green, X$lfmc, method = "spearman")); add(paste0("spearman_", tag), r); add(paste0("n_", tag), nrow(X)) }
  sm(PA, "all"); sm(PA[PA$same_species, ], "same_species"); sm(PA[PA$same_genus, ], "same_genus"); sm(PA[PA$lfmc_grass, ], "lfmc_grass")
  # site-year timing: LFMC falls through 120 % after its peak vs NPN 50 % curing
  cu <- cm[cm$site_id %in% PA$site_id, c("site_id", "year", "species_id", "g50_down_doy")]
  TT <- do.call(rbind, lapply(split(PA, paste(PA$site_id, PA$lfmc_id, PA$year, PA$lfmc_species)), function(z) {
    l <- L[L$lfmc_id == z$lfmc_id[1] & L$year == z$year[1] & L$species == z$lfmc_species[1], ]; l <- l[order(l$date), ]
    if (nrow(l) < 4) return(NULL)
    d <- as.integer(format(l$date, "%j")); v <- l$lfmc; ip <- which.max(v); x <- NA
    for (q in ip:(length(v) - 1)) if (v[q] >= 120 && v[q + 1] < 120) { x <- d[q] + (v[q] - 120) / (v[q] - v[q + 1]) * (d[q + 1] - d[q]); break }
    c50 <- stats::median(cu$g50_down_doy[cu$site_id == z$site_id[1] & cu$year == z$year[1]], na.rm = TRUE)
    data.frame(site_id = z$site_id[1], lfmc_id = z$lfmc_id[1], year = z$year[1], lfmc_species = z$lfmc_species[1], lfmc_120_doy = x, npn_curing50_doy = c50) }))
  if (!is.null(TT) && nrow(TT)) {
    utils::write.csv(TT, file.path(odir, "timing.csv"), row.names = FALSE)
    ok <- !is.na(TT$lfmc_120_doy) & !is.na(TT$npn_curing50_doy); add("timing_site_years", sum(ok))
    if (sum(ok) >= 5) { add("timing_r", stats::cor(TT$lfmc_120_doy[ok], TT$npn_curing50_doy[ok])); add("timing_median_diff_days", stats::median(TT$npn_curing50_doy[ok] - TT$lfmc_120_doy[ok])) }
  }
  log_msg(sprintf("pairs: %d (LFMC samples %d; NPN sites %d); same species %d, same genus %d", nrow(PA), PS$value[PS$stat == "lfmc_samples_paired"], length(unique(PA$site_id)),
                  sum(PA$same_species), sum(PA$same_genus)))
  for (k in grep("^spearman_|^timing_", PS$stat, value = TRUE)) log_msg(sprintf("  %-26s %.2f", k, PS$value[PS$stat == k]))
} else { add("pairs", 0); log_msg("pairs: none within ", max_km, " km and ", max_days, " days") }

utils::write.csv(CO, file.path(odir, "colocation.csv"), row.names = FALSE)
utils::write.csv(IV, file.path(odir, "inventory_ecoregion.csv"), row.names = FALSE)
utils::write.csv(PS, file.path(odir, "pair_summary.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(IV, file.path(tdir, "npn_lfmc_inventory.csv"), row.names = FALSE)
utils::write.csv(CK, file.path(tdir, "npn_lfmc_colocation.csv"), row.names = FALSE)
utils::write.csv(PS, file.path(tdir, "npn_lfmc_pairs.csv"), row.names = FALSE)
for (i in seq_len(min(8, nrow(IV)))) with(IV[i, ], log_msg(sprintf("  %-5s grass fuel %4.1f %% | NPN plant-years %4.1f %% | herbaceous LFMC %4.1f %% of samples (%d sites)", region,
  100 * ifelse(is.na(grass_fuel_area_share), NA, grass_fuel_area_share), 100 * npn_plant_year_share, 100 * lfmc_sample_share, lfmc_sites)))
log_msg("=== done | ", odir, " | log: ", log_file)
