# Ground observations from the GSI study tables (GrowingSeasonIndex/study/data/derived/*.rds):
# one row per site-year with set, site, lon, lat, year, doy. Same sets and filters as
# scripts/31_obs_regions.R (which builds the same table inline); used by 36_curing_target_check.R.
#   grass_greenup / grass_curing   NPN graminoid % green rising / falling through 50 %
#   grass_cured                    NPN graminoid cured date (< 5 % green)
#   herb_onset / herb_end / woody_onset   NPN phenometrics
#   lfmc_<herb|woody>_peak / _decline     Globe-LFMC peak date and 50 % decline date
#   inat_flower                    iNaturalist flowering onset

obs_study_dir <- function(cfg) {
  os <- if (.Platform$OS.type == "windows") "windows" else "linux"
  study <- Sys.getenv("GSI_STUDY_DIR", "")
  if (!nzchar(study)) study <- cfg$evaluation$study_dir[[os]]
  file.path(path.expand(study), "data", "derived")
}

obs_events <- function(ddir, sets = NULL) {
  rd <- function(f) { p <- file.path(ddir, f); if (file.exists(p)) readRDS(p) else { log_warn("missing ", p); NULL } }
  rows <- list()
  add <- function(set, site, lon, lat, year, doy) {
    if (!is.null(sets) && !set %in% sets) return(invisible())
    d <- data.frame(set = set, site = as.character(site), lon = as.numeric(lon), lat = as.numeric(lat),
                    year = as.integer(year), doy = as.numeric(doy))
    d <- d[stats::complete.cases(d), ]
    if (nrow(d)) rows[[length(rows) + 1]] <<- d
  }
  cm <- rd("npn_curing_metrics.rds")
  if (!is.null(cm)) {
    cm <- cm[cm$usable %in% TRUE, ]
    id <- paste0("npn", cm$site_id, "|", cm$species_id)
    add("grass_greenup", id, cm$lon, cm$lat, cm$year, cm$g50_up_doy)
    add("grass_curing", id, cm$lon, cm$lat, cm$year, cm$g50_down_doy)
    add("grass_cured", id, cm$lon, cm$lat, cm$year, cm$cured_doy)
  }
  ph <- rd("npn_phenometrics.rds")
  if (!is.null(ph)) {
    ph <- ph[ph$in_conus %in% TRUE & !(ph$flag_conflict %in% TRUE), ]
    pick <- function(role, fc) ph[ph$role == role & ph$fuel_class == fc, ]
    a <- pick("onset", "herb"); add("herb_onset", paste0("npn", a$site_id, "|", a$species_id), a$longitude, a$latitude, a$event_year, a$event_doy)
    a <- pick("leaves", "herb"); add("herb_end", paste0("npn", a$site_id, "|", a$species_id), a$longitude, a$latitude, a$event_year, a$event_doy)
    a <- pick("onset", "woody"); add("woody_onset", paste0("npn", a$site_id, "|", a$species_id), a$longitude, a$latitude, a$event_year, a$event_doy)
  }
  ls_ <- rd("lfmc_season_metrics.rds")
  if (!is.null(ls_)) {
    ls_ <- ls_[ls_$usable %in% TRUE, ]
    for (fc in c("herb", "woody")) {
      a <- ls_[ls_$fuel_class == fc, ]
      add(sprintf("lfmc_%s_peak", fc), paste0("lfmc", a$series), a$lon, a$lat, a$year, a$peak_doy)
      add(sprintf("lfmc_%s_decline", fc), paste0("lfmc", a$series), a$lon, a$lat, a$year, a$rel50_down_doy)
    }
  }
  ie <- rd("inat_events.rds")
  if (!is.null(ie)) {
    a <- ie[ie$role == "flower_onset", ]
    add("inat_flower", paste0("inat", a$point_id, "|", a$taxon_group), a$grid_lon, a$grid_lat, a$year, a$doy)
  }
  E <- do.call(rbind, rows)
  if (is.null(E)) E <- data.frame(set = character(), site = character(), lon = numeric(), lat = numeric(),
                                  year = integer(), doy = numeric())
  E
}
