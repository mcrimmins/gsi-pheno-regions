#!/usr/bin/env Rscript
# Step 70 (fuel phenology by ecoregion x fuel type; P720 for CONUS, laptop for dev): LANDFIRE fuel
# models and ecoregions on the 4 km analysis grid. Scope: project doc
# claude/gsi-ecofuel-climatology-handoff.md; config `ecofuel:`. New work, separate from the NPN study
# (it reads nothing from npn_study/ except, if present, its downloaded Level II polygons).
#
# 1. LANDFIRE FBFM40 (30 m, Albers) in bands of rows, one band per worker: share of every one of the
#    45 Scott & Burgan models in ~1 km blocks (agg_factor 33 x 30 m), plus `mapped` (share of the
#    block with a valid FBFM40 code) -> <scratch>/ecofuel/tiles_<fact>/band_NNN.tif     [cached]
# 2. Mosaic the bands and area-average onto the 4 km lon/lat grid, masked to grid land cells:
#      ecofuel/fbfm_shares.tif   46 bands: GR1 ... SB4, NB1 ... NB9 (share of mapped pixels) + mapped
#      ecofuel/fuel_groups.tif   12 bands: GR GS SH TU TL SB NB_ag NB_urban NB_barren NB_snow
#                                NB_water (sums of the models above) + mapped
#    Every model is kept, so regrouping later (e.g. GR1-2 vs GR3+, agriculture) needs no rerun.
# 3. Ecoregions (config ecofuel: ecoregions; Level II by default) on the grid: each cell takes the
#    ecoregion covering most of it (sub_cells x sub_cells points per cell); land cells outside every
#    polygon (coast, lake shores) take the nearest polygon (flagged).
#      ecofuel/eco_l<L>.tif      eco (index into the lookup), eco_purity (share of the cell in that
#                                ecoregion), eco_nearest (1 = assigned to the nearest polygon)
#      ecofuel/eco_l<L>_lookup.csv  eco, code, name, l1_code, l1_name, cells, land_km2
# 4. ecofuel/shares_summary.csv: CONUS area of every model and group; log of the top units.
#
# Machine: P720 for CONUS (~15e9 LANDFIRE pixels; with 16 workers expected ~10-20 min for step 1, a
# few minutes for 2-3). The LANDFIRE GeoTIFF must be on the machine (config ecofuel: landfire_fbfm40,
# default <data_root>/landfire): copy it from the laptop, or download it with
# npn_study/scripts/59a_npn_landfire_download.R. Dev box on the laptop: ~1-2 min.
# Workers: step_workers(cfg, "ecofuel_shares") (full: workers_by_step ecofuel_shares; GSI_WORKERS).
# Temporary band tiles go to scratch (GSI_SCRATCH_ROOT; local disk on the P720), ~0.5 GB for CONUS,
# deleted at the end unless ecofuel: keep_tiles is true.
#
# Usage: Rscript scripts/70_ecofuel_shares.R [profile] [force]
#   force: redo every step (otherwise finished band tiles and outputs are reused).

suppressPackageStartupMessages({ library(future); library(furrr) })
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
force <- "force" %in% args

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "ecofuel.R"))
if (!requireNamespace("sf", quietly = TRUE)) stop("needs the sf package")

cfg <- load_cfg(profile)
log_file <- log_init("70_ecofuel_shares", profile)
notify_init(cfg)
ec <- cfg$ecofuel %||% list()
fact <- as.integer(ec$agg_factor %||% 33)
rows_per_band <- as.integer(ec$rows_per_band %||% 60)
sub_cells <- as.integer(ec$ecoregions$sub_cells %||% 10)
lev <- as.integer(ec$ecoregions$level %||% 2)
t0 <- Sys.time()

mask_file <- file.path(static_dir(cfg), "grid_mask.tif")
if (!file.exists(mask_file)) stop("static/grid_mask.tif missing; run 02_static.R first")
odir <- ecofuel_dir(cfg); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
lf_file <- ecofuel_landfire_file(cfg)
FT <- fbfm_table()
w <- step_workers(cfg, "ecofuel_shares")
log_msg("=== 70_ecofuel_shares | profile: ", profile, " | host: ", Sys.info()[["nodename"]],
        " | workers ", w$n, " (memfrac ", signif(w$memfrac, 2), ")", if (force) " | force")
log_msg("LANDFIRE: ", lf_file, " | agg_factor ", fact, " (", fact * 30, " m blocks)")
notify("70_ecofuel_shares started", sprintf("%s | %s", profile, basename(lf_file)), priority = 2)

shares_file <- file.path(odir, "fbfm_shares.tif")
groups_file <- file.path(odir, "fuel_groups.tif")
wr <- function(x, f, dt = "FLT4S") {
  tmp <- paste0(f, ".tmp.tif")
  terra::writeRaster(x, tmp, overwrite = TRUE, datatype = dt,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES", "INTERLEAVE=BAND"))
  unlink(f); invisible(file.rename(tmp, f))
}

# ---- 1. LANDFIRE bands -> ~1 km model shares (parallel) ---------------------------------------
if (force || !file.exists(shares_file)) {
  tdir <- file.path(cfg$scratch_dir, "ecofuel", sprintf("tiles_%d", fact))
  ttmp <- file.path(cfg$scratch_dir, "ecofuel", "terra_tmp")
  dir.create(tdir, showWarnings = FALSE, recursive = TRUE); dir.create(ttmp, showWarnings = FALSE, recursive = TRUE)
  bands <- ecofuel_bands(lf_file, mask_file, fact, rows_per_band)
  tiles <- file.path(tdir, sprintf("band_%03d.tif", seq_along(bands)))
  if (force) unlink(tiles)
  todo <- which(!file.exists(tiles))
  log_msg(sprintf("%d bands of %d LANDFIRE rows (%d to do) | tiles in %s", length(bands),
                  fact * rows_per_band, length(todo), tdir))
  if (length(todo)) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    t1 <- Sys.time()
    if (w$n > 1 && length(todo) > 1) plan(multisession, workers = min(w$n, length(todo))) else plan(sequential)
    memfrac <- w$memfrac; codes <- FT$code; labels <- FT$model
    res <- future_map(todo, function(b) ecofuel_band_shares(lf_file, bands[[b]], codes, labels, tiles[b],
                                                            fact, memfrac, ttmp),
                      .options = furrr_options(seed = NULL, packages = "terra"))
    plan(sequential)
    bad <- which(!vapply(res, function(r) isTRUE(r$ok), TRUE))
    for (i in bad) log_err("band ", todo[i], ": ", res[[i]]$msg)
    if (length(bad)) stop(length(bad), " band(s) failed; rerun to retry them")
    log_msg(sprintf("bands done in %.1f min", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  }

  # ---- 2. mosaic -> 4 km grid ---------------------------------------------------------------
  terra::terraOptions(memfrac = 0.5, progress = 0, tempdir = ttmp)
  mask <- terra::rast(mask_file)
  v <- terra::vrt(tiles, filename = file.path(tdir, "bands.vrt"), overwrite = TRUE)
  names(v) <- c(FT$model, "mapped")
  t1 <- Sys.time()
  s <- terra::project(v, mask, method = "average")
  s <- terra::mask(s, mask)
  names(s) <- c(FT$model, "mapped")
  # shares of mapped pixels re-normalised after averaging (blocks with partial coverage)
  tot <- sum(s[[FT$model]])
  s_mod <- terra::ifel(tot > 0, s[[FT$model]] / tot, NA)
  names(s_mod) <- FT$model
  wr(c(s_mod, s[["mapped"]]), shares_file)
  log_msg(sprintf("fbfm_shares.tif written (%.1f min)", as.numeric(difftime(Sys.time(), t1, units = "mins"))))
  if (!isTRUE(ec$keep_tiles)) unlink(tdir, recursive = TRUE)
  unlink(ttmp, recursive = TRUE)
} else log_msg("fbfm_shares.tif exists (force to redo)")

S <- terra::rast(shares_file)
grp <- ecofuel_groups()
G <- terra::rast(lapply(grp, function(g) sum(S[[FT$model[FT$group == g]]])))
names(G) <- grp
G <- c(G, S[["mapped"]])
wr(G, groups_file)
log_msg("fuel_groups.tif written: ", paste(names(G), collapse = " "))

# ---- 3. ecoregions on the grid ------------------------------------------------------------------
eco_file <- file.path(odir, sprintf("eco_l%d.tif", lev))
lookup_file <- file.path(odir, sprintf("eco_l%d_lookup.csv", lev))
mask <- terra::rast(mask_file)
if (force || !file.exists(eco_file) || !file.exists(lookup_file)) {
  ecofuel_eco_raster(cfg, ec$ecoregions, mask_file, eco_file, lookup_file,
                     water = terra::values(G[["NB_water"]], mat = FALSE))
} else log_msg(basename(eco_file), " exists (force to redo)")

# ---- 4. summary -----------------------------------------------------------------------------------
L <- utils::read.csv(lookup_file, stringsAsFactors = FALSE, colClasses = c(code = "character"))
ca <- terra::values(terra::cellSize(mask, unit = "km"), mat = FALSE)
mp <- terra::values(S[["mapped"]], mat = FALSE)
land <- which(!is.na(terra::values(mask, mat = FALSE)))
area_of <- function(r) { v <- terra::values(r, mat = FALSE); sum((v * mp * ca)[land], na.rm = TRUE) }
sm <- rbind(
  data.frame(level = "model", name = FT$model, group = FT$group,
             km2 = vapply(FT$model, function(m) area_of(S[[m]]), 0)),
  data.frame(level = "group", name = grp, group = grp, km2 = vapply(grp, function(g) area_of(G[[g]]), 0)))
sm$share <- sm$km2 / sum(sm$km2[sm$level == "group"])
sm$km2 <- round(sm$km2)
utils::write.csv(sm, file.path(odir, "shares_summary.csv"), row.names = FALSE)
gs <- sm[sm$level == "group", ]
log_msg("mapped area by group: ", paste(sprintf("%s %.1f %%", gs$name, 100 * gs$share), collapse = ", "))
log_msg(sprintf("grid land cells %d | mean mapped share %.3f | cells with mapped < 0.5: %d",
                length(land), mean(mp[land], na.rm = TRUE), sum(mp[land] < 0.5, na.rm = TRUE)))
log_msg(sprintf("ecoregions: %d (largest: %s)", nrow(L),
                paste(utils::head(sprintf("%s %s", L$code[order(-L$land_km2)], L$name[order(-L$land_km2)]), 4), collapse = "; ")))
msg <- sprintf("%d ecoregions, %d land cells in %.1f min", nrow(L), length(land),
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("70_ecofuel_shares finished", msg, priority = 3, tags = "white_check_mark")
