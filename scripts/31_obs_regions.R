#!/usr/bin/env Rscript
# Step 31 (evaluation, ground observations): do the regions separate observed phenology?
# Observations come from the GSI study exploration (config `evaluation: study_dir`,
# data/derived/*.rds; GSI_STUDY_DIR overrides). Event sets (day of year; site = plant site,
# LFMC site or iNaturalist 0.25-degree point):
#   grass_greenup   NPN graminoid % green: date green rises through 50 % (g50_up_doy)
#   grass_curing    NPN graminoid % green: date green falls through 50 % (g50_down_doy)
#   grass_cured     NPN graminoid: cured date (< 5 % green, or the first "no leaves" after green)
#   herb_onset      NPN phenometrics, herbaceous onset (first yes)
#   herb_end        NPN phenometrics, herbaceous leaves end (last yes)
#   woody_onset     NPN phenometrics, woody onset (contrast: not the target fuel)
#   lfmc_herb_peak, lfmc_herb_decline    Globe-LFMC herbaceous: peak date; 50 % decline date
#   lfmc_woody_peak, lfmc_woody_decline  Globe-LFMC woody series (contrast)
#   inat_flower     iNaturalist flowering onset (mostly shrub genera; grasses are sparse)
# A "site" is one plant type at one place: NPN site x species, LFMC series (site x species),
# iNaturalist point x genus. Each site is reduced to its median date across years (at least min_years_site years), then
# assigned to the 4 km cell it falls in and to its region in every partition (config
# `evaluation: partitions`, plus the geographic baseline) and k. Score: adjusted R^2 of site
# medians on region (one-way ANOVA; adjusted for the number of regions with sites), so a
# partition is rewarded for separating observed dates, not for having more regions.
# Outputs <run_dir>/eval/obs/: sites.csv (site medians with cell and lon/lat), obs_r2.csv,
#   sites_per_region_k<K>.csv (sampling density), figs/obs_r2.png, figs/obs_sites.png.
# Reads the .rds tables from the study folder, so run it where that folder is (the laptop).
#
# Usage: Rscript scripts/31_obs_regions.R [profile]
#   RStudio: gsi_profile <- "full"; source("scripts/31_obs_regions.R")

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

cfg <- load_cfg(profile)
log_file <- log_init("31_obs_regions", profile)
ev <- cfg$evaluation
ks <- seq(ev$k[[1]], ev$k[[2]])
os <- if (.Platform$OS.type == "windows") "windows" else "linux"
study <- Sys.getenv("GSI_STUDY_DIR", "")
if (!nzchar(study)) study <- ev$study_dir[[os]]
ddir <- file.path(path.expand(study), "data", "derived")
if (!dir.exists(ddir)) stop("no ", ddir, " (set GSI_STUDY_DIR or config evaluation: study_dir)")
out_dir <- file.path(cfg$run_dir, "eval", "obs"); fig_dir <- file.path(out_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
log_msg("=== 31_obs_regions | profile: ", profile, " | study data: ", ddir)
rd <- function(f) { p <- file.path(ddir, f); if (file.exists(p)) readRDS(p) else { log_warn("missing ", p); NULL } }

# ---- Events: one row per site-year: set, site, lon, lat, year, doy ------------------------
ev_rows <- list()
add <- function(set, site, lon, lat, year, doy) {
  d <- data.frame(set = set, site = as.character(site), lon = as.numeric(lon), lat = as.numeric(lat),
                  year = as.integer(year), doy = as.numeric(doy))
  d <- d[stats::complete.cases(d), ]
  if (nrow(d)) ev_rows[[length(ev_rows) + 1]] <<- d
  log_msg(sprintf("  %-20s %7d site-years", set, nrow(d)))
}
cm <- rd("npn_curing_metrics.rds")
if (!is.null(cm)) {
  cm <- cm[cm$usable %in% TRUE, ]
  add("grass_greenup", paste0("npn", cm$site_id, "|", cm$species_id), cm$lon, cm$lat, cm$year, cm$g50_up_doy)
  add("grass_curing", paste0("npn", cm$site_id, "|", cm$species_id), cm$lon, cm$lat, cm$year, cm$g50_down_doy)
  add("grass_cured", paste0("npn", cm$site_id, "|", cm$species_id), cm$lon, cm$lat, cm$year, cm$cured_doy)
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
E <- do.call(rbind, ev_rows)
if (is.null(E) || !nrow(E)) stop("no observations found in ", ddir)

# ---- Site medians, cells, regions -----------------------------------------------------------
S <- stats::aggregate(cbind(doy, lon, lat) ~ set + site, E, stats::median)
ny <- stats::aggregate(year ~ set + site, E, function(y) length(unique(y)))
S$n_years <- ny$year[match(paste(S$set, S$site), paste(ny$set, ny$site))]
S <- S[S$n_years >= ev$min_years_site, ]
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land_all <- which(!is.na(terra::values(grid, mat = FALSE)))
S$cell <- terra::cellFromXY(grid, as.matrix(S[, c("lon", "lat")]))
S$row <- match(S$cell, land_all)
log_msg(sum(is.na(S$row)), " of ", nrow(S), " site medians fall outside the land grid (dropped)")
S <- S[!is.na(S$row), ]
utils::write.csv(S, file.path(out_dir, "sites.csv"), row.names = FALSE)
P <- eval_partitions(cfg, grid, ks)
log_msg("partitions: ", paste(names(P), collapse = ", "), " | sites per set: ",
        paste(sprintf("%s %d", names(table(S$set)), table(S$set)), collapse = ", "))

res <- data.frame()
for (st in unique(S$set)) {
  s <- S[S$set == st, ]
  for (nm in names(P)) for (kn in colnames(P[[nm]])) {
    r <- r2_adj(s$doy, P[[nm]][s$row, kn])
    res <- rbind(res, data.frame(set = st, partition = nm, k = as.integer(sub("k", "", kn)), t(r)))
  }
}
utils::write.csv(res, file.path(out_dir, "obs_r2.csv"), row.names = FALSE)
for (k in intersect(unlist(ev$report_k), ks)) {
  log_msg(sprintf("k = %d: adjusted R^2 of site median dates on region", k))
  for (st in unique(res$set)) {
    d <- res[res$set == st & res$k == k, ]
    log_msg(sprintf("  %-20s n=%4d  %s", st, d$n[1], paste(sprintf("%s %5.2f", d$partition, d$adj_r2), collapse = "  ")))
  }
  # sampling density: sites per region and set
  for (nm in setdiff(names(P), "geo")) {
    kn <- sprintf("k%02d", k); if (!kn %in% colnames(P[[nm]])) next
    tab <- as.data.frame.matrix(table(factor(P[[nm]][S$row, kn], levels = seq_len(k)), S$set))
    tab <- cbind(region = seq_len(k), cells = tabulate(P[[nm]][, kn], k), tab)
    utils::write.csv(tab, file.path(out_dir, sprintf("sites_per_region_%s_k%02d.csv", nm, k)), row.names = FALSE)
  }
}

# ---- Figures ------------------------------------------------------------------------------
pn <- names(P); pal <- c("#2a78d6", "#1a7f37", "#eb6834", "#c0392b", "#7b4fb3")
pcol <- setNames(c(pal, "#888888")[seq_along(pn)], pn); pcol["geo"] <- "#888888"
plty <- setNames(ifelse(grepl("^B_", pn), 2, 1), pn); plty["geo"] <- 3
sets <- unique(res$set); nc <- 4; nr <- ceiling((length(sets) + 1) / nc)
png(file.path(fig_dir, "obs_r2.png"), 520 * nc, 420 * nr, res = 150)
par(mfrow = c(nr, nc), mar = c(4, 4, 2.2, 0.5), las = 1)
for (st in sets) {
  d <- res[res$set == st, ]
  yl <- range(c(0, d$adj_r2), na.rm = TRUE); yl[2] <- max(yl[2], 0.2)
  plot(NA, xlim = range(ks), ylim = yl, xlab = "k", ylab = "adjusted R^2", cex.main = 0.85,
       main = sprintf("%s (%d sites)", st, max(d$n, na.rm = TRUE)))
  abline(h = 0, col = "#cccccc")
  for (nm in pn) { e <- d[d$partition == nm, ]; lines(e$k, e$adj_r2, col = pcol[nm], lty = plty[nm], lwd = 2) }
}
plot.new(); legend("center", pn, col = pcol, lty = plty, lwd = 2, bty = "n")
dev.off()
states <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL
png(file.path(fig_dir, "obs_sites.png"), 2000, 1250, res = 150)
par(mar = c(0.5, 0.5, 2, 0.5))
terra::plot(grid, col = "#e4e3df", legend = FALSE, axes = FALSE, main = "Observation sites by set (site medians)")
if (!is.null(states)) lines(states$x, states$y, col = "#ffffff", lwd = 0.7)
sc <- setNames(grDevices::hcl.colors(length(sets), "Dark 3"), sets)
points(S$lon, S$lat, pch = 16, cex = 0.45, col = grDevices::adjustcolor(sc[S$set], 0.7))
legend("bottomleft", sets, col = sc, pch = 16, bty = "n", cex = 0.7, ncol = 2)
dev.off()
log_msg("=== done | ", nrow(S), " site medians, ", length(P), " partitions | outputs: ", out_dir, " | log: ", log_file)
