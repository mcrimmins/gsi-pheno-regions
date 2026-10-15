#!/usr/bin/env Rscript
# Step 50 (NPN grass curing study, scope: project doc gsi-npn-curing-study-scope.md): daily PRISM
# at every NPN grass site, so the climate metrics and the model ladder (51+) can run anywhere
# without the CONUS rasters.
# Sites: every site in npn_curing_metrics.rds (GSI study tables; folder from GSI_OBS_DIR or
# config evaluation: study_dir), located on the 4 km grid. Years: config npn_study: years
# (default 1990 to the last PRISM year, covering the 1991-2020 normals and the observation
# years with the year before).
# Output <run_dir>/npn_study/npn_site_daily.rds: list(sites = data.frame(site_id, lat, lon, cell,
# elev_m?), dates, ppt, tmin, tmax, vpdmax) with sites x days matrices (mm, C, C, hPa).
# Small: ~130 sites x ~13,000 days. Runs on the P720 (reads the yearly daily PRISM stacks).
#
# Usage: GSI_OBS_DIR=<folder with npn_curing_metrics.rds> Rscript scripts/50_npn_prism.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("50_npn_prism", profile)
ns <- cfg$npn_study %||% list()
out_dir <- file.path(cfg$run_dir, "npn_study"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
log_msg("=== 50_npn_prism | profile: ", profile)

ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
f <- file.path(ddir, "npn_curing_metrics.rds")
if (!file.exists(f)) stop("no ", f, " (set GSI_OBS_DIR to the folder with npn_curing_metrics.rds)")
cm <- readRDS(f)
S <- unique(cm[, c("site_id", "lat", "lon")])
S <- S[!duplicated(S$site_id), ]
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
S$cell <- terra::cellFromXY(grid, as.matrix(S[, c("lon", "lat")]))
S$on_grid <- !is.na(S$cell) & !is.na(terra::values(grid, mat = FALSE)[S$cell])
if (any(!S$on_grid)) log_warn(sum(!S$on_grid), " site(s) off the land grid: ", paste(S$site_id[!S$on_grid], collapse = ", "))
S <- S[S$on_grid, ]
pfile <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
yr_cfg <- ns$years %||% list(start = 1990, end = max(cfg$years))
years <- seq(yr_cfg$start, yr_cfg$end)
years <- years[file.exists(vapply(years, function(y) pfile("ppt", y), ""))]
vars <- unlist(ns$vars %||% list("ppt", "tmin", "tmax", "vpdmax"))
log_msg(sprintf("%d sites on the grid; years %d-%d; vars %s", nrow(S), min(years), max(years), paste(vars, collapse = ", ")))

out <- list(sites = S)
for (v in vars) {
  parts <- lapply(years, function(y) {
    fy <- pfile(v, y)
    if (!file.exists(fy)) { log_warn("missing ", basename(fy)); return(NULL) }
    r <- terra::rast(fy); x <- as.matrix(terra::extract(r, S$cell)); colnames(x) <- names(r); x
  })
  out[[v]] <- do.call(cbind, parts)
  log_msg(sprintf("  %s: %d sites x %d days (%.1f min)", v, nrow(out[[v]]), ncol(out[[v]]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
dd <- Reduce(intersect, lapply(vars, function(v) colnames(out[[v]])))
for (v in vars) out[[v]] <- out[[v]][, dd, drop = FALSE]
out$dates <- as.Date(dd)
saveRDS(out, file.path(out_dir, "npn_site_daily.rds"))
log_msg(sprintf("=== done | %.1f min | %s (%d sites, %s .. %s) | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
                file.path(out_dir, "npn_site_daily.rds"), nrow(S), min(out$dates), max(out$dates), log_file))
