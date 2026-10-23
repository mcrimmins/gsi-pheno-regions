#!/usr/bin/env Rscript
# Step 55d (NPN grass curing study; P720, where the PRISM daily stacks are): daily PRISM at the
# LTER phenology sites (55c lter_sites.csv), in the same layout as 50's npn_site_daily.rds, so
# 55e computes the same site-year weather metrics (npn_site_year_metrics) as the model ladder.
# The 15 Jornada NPP sites share one point (the basin bounding-box centre; exact plot coordinates
# are not public), so they share one PRISM cell.
# Output <run_dir>/npn_study/lter/lter_site_daily.rds: list(sites, dates, ppt, tmin, tmax, vpdmax).
#
# Usage: Rscript npn_study/scripts/55d_lter_prism.R [profile]
#   needs <run_dir>/npn_study/lter/lter_sites.csv (55c; copy it to the P720 if 55c ran on the laptop)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("55d_lter_prism", profile)
ns <- cfg$npn_study %||% list()
ldir <- file.path(cfg$run_dir, "npn_study", "lter")
sf <- file.path(ldir, "lter_sites.csv")
if (!file.exists(sf)) stop("no ", sf, " (run 55c, or copy lter_sites.csv from the laptop)")
terra::terraOptions(progress = 0)
t0 <- Sys.time()
S <- utils::read.csv(sf, stringsAsFactors = FALSE)[, c("site_id", "source", "site", "lat", "lon")]
log_msg("=== 55d_lter_prism | ", nrow(S), " sites")
pfile <- function(v, y) file.path(cfg$prism_daily_dir, v, sprintf("prism_%s_%d.tif", v, y))
yr_cfg <- ns$years %||% list(start = 2001, end = 2025)
years <- seq(yr_cfg$start, yr_cfg$end)
years <- years[file.exists(vapply(years, function(y) pfile("ppt", y), ""))]
if (!length(years)) stop("no PRISM daily stacks in ", cfg$prism_daily_dir)
g <- terra::rast(pfile("ppt", years[1]))[[1]]
S$cell <- terra::cellFromXY(g, as.matrix(S[, c("lon", "lat")]))
S <- S[!is.na(S$cell) & !is.na(terra::values(g, mat = FALSE)[S$cell]), ]
vars <- unlist(ns$vars %||% list("ppt", "tmin", "tmax", "vpdmax"))
log_msg(sprintf("%d sites on the grid (%d cells); years %d-%d", nrow(S), length(unique(S$cell)), min(years), max(years)))
out <- list(sites = S)
for (v in vars) {
  out[[v]] <- do.call(cbind, lapply(years, function(y) {
    fy <- pfile(v, y); if (!file.exists(fy)) { log_warn("missing ", basename(fy)); return(NULL) }
    r <- terra::rast(fy); x <- as.matrix(terra::extract(r, S$cell)); colnames(x) <- names(r); x }))
  log_msg(sprintf("  %s: %d x %d (%.1f min)", v, nrow(out[[v]]), ncol(out[[v]]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
dd <- Reduce(intersect, lapply(vars, function(v) colnames(out[[v]])))
for (v in vars) out[[v]] <- out[[v]][, dd, drop = FALSE]
out$dates <- as.Date(dd)
saveRDS(out, file.path(ldir, "lter_site_daily.rds"))
log_msg(sprintf("=== done | %s (%s .. %s) | log: %s", file.path(ldir, "lter_site_daily.rds"), min(out$dates), max(out$dates), log_file))
