#!/usr/bin/env Rscript
# Step 59b (NPN grass curing study; laptop, terra): the Scott & Burgan (2005) fire behavior fuel
# model (LANDFIRE FBFM40, 30 m) at every NPN grass site.
# Input: a local LANDFIRE FBFM40 CONUS GeoTIFF (config npn_study: landfire_fbfm40, per OS), e.g.
# "LF2023 FBFM40 CONUS" from landfire.gov (Data > Full extent downloads). Any LANDFIRE version works;
# the version is recorded from the file name.
# Each site takes the fuel model at its point, and the most common burnable fuel model within
# 90 m (3 x 3 cells) as a check against coordinate error and roads or buildings at the point.
# Output data/npn_site_fuelmodels.csv (site ids and codes only, in git): site_id, fbfm_code,
# fbfm (e.g. GR2), fbfm_group (GR grass, GS grass-shrub, SH shrub, TU timber-understory,
# TL timber litter, SB slash-blowdown, NB nonburnable), fbfm_3x3, fbfm_3x3_group, landfire_file.
#
# Usage: Rscript npn_study/scripts/59b_npn_fuelmodels.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("terra", quietly = TRUE)) stop("needs terra")

cfg <- load_npn_cfg(profile)
log_file <- log_init("59b_npn_fuelmodels", profile)
ns <- cfg$npn_study %||% list()
os <- if (.Platform$OS.type == "windows") "windows" else "linux"
lf <- ns$landfire_fbfm40; if (is.list(lf)) lf <- lf[[os]]
if (is.null(lf) || !nzchar(lf)) lf <- file.path(cfg$data_root, "landfire")
lf <- path.expand(lf)
if (dir.exists(lf)) {                                                      # a folder: find the FBFM40 GeoTIFF in it
  f <- list.files(lf, "\\.tif$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  f <- f[grepl("f40|fbfm40", basename(f), ignore.case = TRUE) | grepl("fbfm40", dirname(f), ignore.case = TRUE)]
  if (!length(f)) stop("no FBFM40 .tif under ", lf, " (set npn_study: landfire_fbfm40 in npn_study/config.yml)")
  lf <- f[1]
}
if (!file.exists(lf)) stop("LANDFIRE FBFM40 raster not found: ", lf)
log_msg("=== 59b_npn_fuelmodels | ", lf)

# Scott & Burgan (2005) codes
sb <- c(`91` = "NB1", `92` = "NB2", `93` = "NB3", `98` = "NB8", `99` = "NB9",
        stats::setNames(paste0("GR", 1:9), 101:109), stats::setNames(paste0("GS", 1:4), 121:124),
        stats::setNames(paste0("SH", 1:9), 141:149), stats::setNames(paste0("TU", 1:5), 161:165),
        stats::setNames(paste0("TL", 1:9), 181:189), stats::setNames(paste0("SB", 1:4), 201:204))
nm <- function(code) unname(ifelse(as.character(code) %in% names(sb), sb[as.character(code)], NA))
grp <- function(lab) ifelse(is.na(lab), NA, substr(lab, 1, 2))

R <- terra::rast(lf)[[1]]
if (any(terra::is.factor(R))) levels(R) <- NULL                       # LANDFIRE ships a category table: use the codes
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
S <- stats::aggregate(cbind(lat, lon) ~ site_id, data = cm, FUN = stats::median)
xy <- terra::project(as.matrix(S[, c("lon", "lat")]), "EPSG:4326", terra::crs(R))   # one row per site, kept in order
at_xy <- function(x, y) {                                                  # value at points; NA outside the raster
  cl <- terra::cellFromXY(R, cbind(x, y)); v <- rep(NA_real_, length(x)); ok <- !is.na(cl)
  if (any(ok)) v[ok] <- as.numeric(terra::extract(R, cl[ok])[, 1]); v }
pt <- at_xy(xy[, 1], xy[, 2])
# 3 x 3 window: all nine cells around each point
rs <- terra::res(R); off <- expand.grid(dx = c(-1, 0, 1) * rs[1], dy = c(-1, 0, 1) * rs[2])
W <- sapply(seq_len(nrow(off)), function(k) at_xy(xy[, 1] + off$dx[k], xy[, 2] + off$dy[k]))
if (is.null(dim(W))) W <- matrix(W, nrow = 1)
mode_burn <- apply(W, 1, function(z) { z <- z[!is.na(z) & z >= 101]; if (!length(z)) NA else as.numeric(names(sort(table(z), decreasing = TRUE))[1]) })
out <- data.frame(site_id = S$site_id, fbfm_code = pt, fbfm = nm(pt), fbfm_group = grp(nm(pt)),
                  fbfm_3x3 = nm(mode_burn), fbfm_3x3_group = grp(nm(mode_burn)), landfire_file = basename(lf))
utils::write.csv(out, npn_path("data", "npn_site_fuelmodels.csv"), row.names = FALSE)
tb <- sort(table(out$fbfm, useNA = "ifany"), decreasing = TRUE)
log_msg(sprintf("sites: %d | at the point: %s", nrow(out), paste(names(tb), tb, collapse = ", ")))
log_msg(sprintf("point nonburnable or missing: %d; of those, a burnable model within 90 m: %d",
                sum(is.na(out$fbfm_group) | out$fbfm_group == "NB"), sum((is.na(out$fbfm_group) | out$fbfm_group == "NB") & !is.na(out$fbfm_3x3))))
log_msg("=== done | data/npn_site_fuelmodels.csv | log: ", log_file)
