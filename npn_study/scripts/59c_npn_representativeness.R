#!/usr/bin/env Rscript
# Step 59c (NPN grass curing study; laptop, terra + sf): how well does the NPN grass record cover
# the grass fuels of the contiguous U.S.? "Grass fuels" = cells of the LANDFIRE FBFM40 map (59a)
# in a grass (GR1-GR9) or grass-shrub (GS1-GS4) fuel model. A regular sample of the 30 m map
# (equal-area projection, so shares are area shares) is assigned EPA Level II ecoregions (59's
# polygons) and compared with the record (relaxed set; fuel model at each site from 59b).
# Outputs <run_dir>/npn_study/record/
#   representativeness_ecoregion.csv   per ecoregion: share of U.S. grass fuel area; sites,
#                                      plant-years and their shares in the record; coverage ratio
#                                      (record share / area share)
#   representativeness_fuelmodel.csv   the same per grass / grass-shrub fuel model
#   grass_fuel_sample.csv              the sampled grass-fuel cells (lon, lat, fuel model, ecoregion)
#                                      for the map
# and report/tables/npn_representativeness_{ecoregion,fuelmodel}.csv. Config npn_study:
# repr_sample (default 2e6 cells sampled before keeping grass fuels).
#
# Usage: Rscript npn_study/scripts/59c_npn_representativeness.R [profile]   (after 58, 59, 59a, 59b)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
for (p in c("terra", "sf")) if (!requireNamespace(p, quietly = TRUE)) stop("needs ", p)

cfg <- load_npn_cfg(profile)
log_file <- log_init("59c_npn_representativeness", profile)
ns <- cfg$npn_study %||% list()
base <- file.path(cfg$run_dir, "npn_study"); rdir <- file.path(base, "record")
t0 <- Sys.time()
n_samp <- as.numeric(ns$repr_sample %||% 2e6)
log_msg("=== 59c_npn_representativeness | profile: ", profile, " | sample ", format(n_samp, big.mark = ","))

# ---- LANDFIRE fuel model map ---------------------------------------------------------------------
os <- if (.Platform$OS.type == "windows") "windows" else "linux"
lf <- ns$landfire_fbfm40; if (is.list(lf)) lf <- lf[[os]]
if (is.null(lf) || !nzchar(lf)) lf <- file.path(cfg$data_root, "landfire")
lf <- path.expand(lf)
if (dir.exists(lf)) { f <- list.files(lf, "\\.tif$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[grepl("f40|fbfm40", basename(f), ignore.case = TRUE) | grepl("fbfm40", dirname(f), ignore.case = TRUE)]
  if (!length(f)) stop("no FBFM40 .tif under ", lf, " (run 59a)"); lf <- f[1] }
R <- terra::rast(lf)[[1]]; if (any(terra::is.factor(R))) levels(R) <- NULL
sb <- c(stats::setNames(paste0("GR", 1:9), 101:109), stats::setNames(paste0("GS", 1:4), 121:124))
log_msg("fuel model map: ", basename(lf), sprintf(" (%d x %d cells)", terra::nrow(R), terra::ncol(R)))
S <- terra::spatSample(R, size = n_samp, method = "regular", xy = TRUE, na.rm = FALSE)
names(S)[3] <- "code"
S <- S[!is.na(S$code), ]
n_land <- nrow(S)
S <- S[as.character(S$code) %in% names(sb), ]
S$fm <- unname(sb[as.character(S$code)]); S$group <- substr(S$fm, 1, 2)
log_msg(sprintf("sampled %s mapped cells; %s (%.0f %%) in grass or grass-shrub fuel models", format(n_land, big.mark = ","), format(nrow(S), big.mark = ","), 100 * nrow(S) / n_land))

# ---- ecoregions ------------------------------------------------------------------------------------
edir <- file.path(base, "ecoregions")
shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
shp <- shp[grepl("l2|level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
if (is.na(shp)) stop("no ecoregion shapefile in ", edir, " (run 59)")
E <- sf::st_read(shp, quiet = TRUE); gc <- attr(E, "sf_column"); names(E)[names(E) != gc] <- toupper(names(E)[names(E) != gc])
E <- sf::st_transform(sf::st_make_valid(E), terra::crs(R))
EV <- terra::vect(E)
# rasterize the ecoregion codes on a 2 km grid in the fuel map's projection and look them up
tmpl <- terra::rast(terra::ext(R), resolution = 2000, crs = terra::crs(R))
codes <- sort(unique(as.character(E$NA_L2CODE))); EV$code_i <- match(as.character(E$NA_L2CODE), codes)
ER <- terra::rasterize(EV, tmpl, field = "code_i")
ci <- terra::extract(ER, as.matrix(S[, c("x", "y")]))[, 1]
S$region <- ifelse(is.na(ci), "none", codes[ci])
ll <- terra::project(as.matrix(S[, c("x", "y")]), terra::crs(R), "EPSG:4326"); S$lon <- ll[, 1]; S$lat <- ll[, 2]
nm_tab <- unique(data.frame(region = as.character(E$NA_L2CODE), name = tools::toTitleCase(tolower(as.character(E$NA_L2NAME)))))

# ---- the record --------------------------------------------------------------------------------------
TY <- utils::read.csv(file.path(rdir, "series_typical.csv"), stringsAsFactors = FALSE, colClasses = c(region = "character"))
st <- stats::aggregate(plant_years ~ site_id + region + group, TY, sum)
st$wild <- TRUE
fm_f <- npn_path("data", "npn_site_fuelmodels.csv")
if (file.exists(fm_f)) { FM <- utils::read.csv(fm_f, stringsAsFactors = FALSE)
  FM$fm <- ifelse(!is.na(FM$fbfm_group) & FM$fbfm_group != "NB", FM$fbfm, FM$fbfm_3x3)
  st$fm <- FM$fm[match(st$site_id, FM$site_id)] } else st$fm <- NA
st$grass_site <- substr(st$fm, 1, 2) %in% c("GR", "GS")

# ---- coverage by ecoregion ------------------------------------------------------------------------
ar <- table(S$region[S$region != "none"]); ar <- ar / sum(ar)
rg <- sort(unique(c(names(ar), st$region[st$region != "none"])))
RE <- data.frame(region = rg, name = nm_tab$name[match(rg, nm_tab$region)],
                 grass_fuel_area_share = as.numeric(ar[rg]), stringsAsFactors = FALSE)
RE$grass_fuel_area_share[is.na(RE$grass_fuel_area_share)] <- 0
RE$sites <- vapply(rg, function(k) sum(st$region == k), 0); RE$sites_neon <- vapply(rg, function(k) sum(st$region == k & st$group == "NEON"), 0)
RE$plant_years <- vapply(rg, function(k) sum(st$plant_years[st$region == k]), 0)
RE$site_share <- RE$sites / sum(RE$sites); RE$plant_year_share <- RE$plant_years / sum(RE$plant_years)
RE$coverage_sites <- RE$site_share / RE$grass_fuel_area_share; RE$coverage_plant_years <- RE$plant_year_share / RE$grass_fuel_area_share
RE <- RE[order(-RE$grass_fuel_area_share), ]
# ---- coverage by fuel model -----------------------------------------------------------------------
af <- table(S$fm); af <- af / sum(af)
fms <- unname(sb)
FMt <- data.frame(fuel_model = fms, grass_fuel_area_share = as.numeric(af[fms]), stringsAsFactors = FALSE)
FMt$grass_fuel_area_share[is.na(FMt$grass_fuel_area_share)] <- 0
gs <- st[st$grass_site, ]
FMt$sites <- vapply(fms, function(k) sum(gs$fm == k, na.rm = TRUE), 0); FMt$plant_years <- vapply(fms, function(k) sum(gs$plant_years[gs$fm %in% k]), 0)
FMt$site_share <- FMt$sites / max(sum(FMt$sites), 1); FMt$plant_year_share <- FMt$plant_years / max(sum(FMt$plant_years), 1)
FMt$coverage_plant_years <- ifelse(FMt$grass_fuel_area_share > 0, FMt$plant_year_share / FMt$grass_fuel_area_share, NA)

utils::write.csv(RE, file.path(rdir, "representativeness_ecoregion.csv"), row.names = FALSE)
utils::write.csv(FMt, file.path(rdir, "representativeness_fuelmodel.csv"), row.names = FALSE)
keep <- S[, c("lon", "lat", "fm", "region")]; if (nrow(keep) > 60000) keep <- keep[sort(sample(nrow(keep), 60000)), ]
utils::write.csv(keep, file.path(rdir, "grass_fuel_sample.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(RE, file.path(tdir, "npn_representativeness_ecoregion.csv"), row.names = FALSE)
utils::write.csv(FMt, file.path(tdir, "npn_representativeness_fuelmodel.csv"), row.names = FALSE)

for (i in seq_len(min(nrow(RE), 20))) with(RE[i, ], log_msg(sprintf("  %-5s %-36s grass fuel %4.1f %% | sites %3d (NEON %2d) %4.1f %% | plant-years %5d %4.1f %% | coverage %.2f",
  region, substr(name %||% "", 1, 36), 100 * grass_fuel_area_share, sites, sites_neon, 100 * site_share, plant_years, 100 * plant_year_share, coverage_plant_years)))
for (i in seq_len(nrow(FMt))) with(FMt[i, ], log_msg(sprintf("  %-4s grass fuel %4.1f %% | grass-model sites %3d | plant-years %5d %4.1f %%", fuel_model, 100 * grass_fuel_area_share, sites, plant_years, 100 * plant_year_share)))
z <- RE[RE$grass_fuel_area_share >= 0.02, ]
log_msg(sprintf("ecoregions with >= 2 %% of U.S. grass fuel: %d, holding %.0f %% of it; %d with no site, %d with < 1 %% of plant-years",
                nrow(z), 100 * sum(z$grass_fuel_area_share), sum(z$sites == 0), sum(z$plant_year_share < 0.01)))
log_msg(sprintf("=== done | %.1f min | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), log_file))
