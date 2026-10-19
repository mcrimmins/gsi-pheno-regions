#!/usr/bin/env Rscript
# Step 59 (NPN grass curing study; laptop with internet, sf): EPA / CEC Level II ecoregions for
# every NPN grass site. Decided 18 Oct 2026: the NPN study is a standalone paper, so it uses a
# published regionalization (CEC 1997, Level II) instead of the fuel phenology regions (A2 k13),
# which were chosen partly on these same ground dates.
#   - downloads the Level II polygons once (config npn_study: ecoregions_url) to
#     <run_dir>/npn_study/ecoregions/
#   - assigns each site in npn_curing_metrics.rds by point in polygon; sites just outside every
#     polygon (coast, lakes) take the nearest one (flagged)
# Output data/npn_site_ecoregions.csv: site_id, eco_l2_code, eco_l2_name, eco_l1_code, eco_l1_name,
# nearest (TRUE when assigned to the nearest polygon), dist_km. Site ids and codes only (no site
# names), so it can go in git. Scripts 51, 53-55, 58 and 90 read it through npn_site_regions().
#
# Usage: Rscript npn_study/scripts/59_npn_ecoregions.R [profile]
# Needs: renv::install("sf")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("sf", quietly = TRUE)) stop("needs the sf package: renv::install('sf')")

cfg <- load_npn_cfg(profile)
log_file <- log_init("59_npn_ecoregions", profile)
ns <- cfg$npn_study %||% list()
url <- ns$ecoregions_url %||% "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip"
edir <- file.path(cfg$run_dir, "npn_study", "ecoregions"); dir.create(edir, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 59_npn_ecoregions | profile: ", profile, " | ", url)

zf <- file.path(edir, basename(url))
if (!file.exists(zf)) utils::download.file(url, zf, mode = "wb")
utils::unzip(zf, exdir = edir)
shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
shp <- shp[grepl("l2|level2|Level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
E <- sf::st_read(shp, quiet = TRUE)
gc <- attr(E, "sf_column")                                  # upper-case attribute names only, not the geometry
names(E)[names(E) != gc] <- toupper(names(E)[names(E) != gc])
E <- sf::st_make_valid(E)
log_msg(sprintf("%s: %d polygons, %d Level II regions", basename(shp), nrow(E), length(unique(E$NA_L2CODE))))

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
S <- stats::aggregate(cbind(lat, lon) ~ site_id, data = cm, FUN = stats::median)
P <- sf::st_transform(sf::st_as_sf(S, coords = c("lon", "lat"), crs = 4326), sf::st_crs(E))
hit <- sf::st_intersects(P, E)
i <- vapply(hit, function(h) if (length(h)) h[1] else NA_integer_, 0L)
near <- is.na(i)
if (any(near)) i[near] <- sf::st_nearest_feature(P[near, ], E)
dist <- rep(0, nrow(P))
if (any(near)) dist[near] <- as.numeric(sf::st_distance(P[near, ], E[i[near], ], by_element = TRUE)) / 1000
out <- data.frame(site_id = S$site_id, eco_l2_code = as.character(E$NA_L2CODE[i]), eco_l2_name = as.character(E$NA_L2NAME[i]),
                  eco_l1_code = as.character(E$NA_L1CODE[i]), eco_l1_name = as.character(E$NA_L1NAME[i]),
                  nearest = near, dist_km = round(dist, 1))
utils::write.csv(out, npn_path("data", "npn_site_ecoregions.csv"), row.names = FALSE)
tb <- sort(table(paste(out$eco_l2_code, out$eco_l2_name)), decreasing = TRUE)
log_msg(sprintf("sites: %d | %d Level II regions | %d by nearest polygon (max %.1f km)", nrow(out), length(tb), sum(near), max(c(0, dist))))
for (k in names(tb)) log_msg(sprintf("  %-50s %3d sites", k, tb[[k]]))
log_msg("=== done | data/npn_site_ecoregions.csv | log: ", log_file)
