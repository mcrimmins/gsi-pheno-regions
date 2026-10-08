#!/usr/bin/env Rscript
# Step 12b: check the grid GSI (R/gsi.R Part 2, used by 12_gsi.R) against gsi-scout's point
# model (R/gsi.R Part 1, run_gsi()) at the check sites inside the profile's extent, for one
# year and every variant. For each site x variant: the daily smoothed GSI and phase from
# run_gsi() vs the grid code on the same series, and a few per-year values recomputed from
# the point series in plain R vs the values 12_gsi.R wrote to features/gsi_<variant>/.
# Differences should be at floating-point level (< 1e-9) and 0 phase mismatches.
# Output: logs + <run_dir>/eval/gsi_check.csv. Run after 12_gsi.R (a minute or two).
#
# Usage: Rscript scripts/12b_gsi_check.R [profile] [year]
#   RStudio: gsi_profile <- "dev"; source("scripts/12b_gsi_check.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "gsi.R"))

cfg <- load_cfg(profile)
log_file <- log_init("12b_gsi_check", profile)
g <- cfg$features$gsi; variants <- gsi_variants(g)
years <- feature_years(cfg)
Y <- if (length(args) >= 2) as.integer(args[2]) else years[1]
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
mapr <- terra::rast(file.path(static_dir(cfg), "kbdi_map.tif"))
cp <- do.call(rbind, lapply(cfg$clustering$check_points, as.data.frame))
cp$cell <- terra::cellFromXY(grid, as.matrix(cp[, c("lon", "lat")]))
cp <- cp[!is.na(cp$cell), ]
cp <- cp[!is.na(terra::extract(grid, cp$cell)[, 1]), ]
if (!nrow(cp)) stop("no check sites inside the ", profile, " extent")
log_msg("=== 12b_gsi_check | profile: ", profile, " | year ", Y, " | sites: ", paste(cp$name, collapse = "; "))

pf <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
ys <- c(Y - 1, Y, if (file.exists(pf("ppt", Y + 1)) && Y + 1 <= max(cfg$years)) Y + 1)
last_day <- as.Date(sprintf("%d-03-31", Y + 1))
series <- function(v, cells) {
  do.call(cbind, lapply(ys, function(y) {
    r <- terra::rast(pf(v, y)); d <- as.Date(names(r)); k <- which(d <= last_day)
    x <- as.matrix(terra::extract(r[[k]], cells)); colnames(x) <- format(d[k]); x
  }))
}
X <- lapply(c(ppt = "ppt", tmin = "tmin", tmax = "tmax", vpdmax = "vpdmax"), series, cells = cp$cell)
dates <- as.Date(colnames(X$ppt))
dd <- as.numeric(dates - as.Date(sprintf("%d-01-01", Y))) + 1
iy <- which(format(dates, "%Y") == as.character(Y))
lat <- terra::yFromCell(grid, cp$cell)
map_mm <- terra::extract(mapr, cp$cell)[, 1]

res <- data.frame()
for (i in seq_len(nrow(cp))) {
  met <- data.frame(date = dates, doy = as.integer(format(dates, "%j")), year = as.integer(format(dates, "%Y")),
                    tmin_c = X$tmin[i, ], tmax_c = X$tmax[i, ], vpd_pa = X$vpdmax[i, ] * 100, precip_mm = X$ppt[i, ])
  D <- list(ppt = X$ppt[i, , drop = FALSE], tmin = X$tmin[i, , drop = FALSE], tmax = X$tmax[i, , drop = FALSE],
            vpd_pa = X$vpdmax[i, , drop = FALSE] * 100)
  photo <- m_photoperiod(lat[i], met$doy, g$photo_method)
  kb <- m_kbdi(D$ppt, D$tmax, map_mm[i] / 25.4, g$kbdi$interception_in)
  for (nm in names(variants)) {
    v <- variants[[nm]]; ps <- g$param_sets[[v$param_set]]
    p <- list(tmin_lo = ps$tmin[[1]], tmin_hi = ps$tmin[[2]],
              vpd_lo = if (v$vpd) ps$vpd[[1]] else 1e12, vpd_hi = if (v$vpd) ps$vpd[[2]] else 2e12,
              photo_lo = ps$photo_h[[1]] * 3600, photo_hi = ps$photo_h[[2]] * 3600, photo_method = g$photo_method,
              use_precip = v$moisture == "precip", precip_window = ps$precip$window,
              precip_lo = ps$precip$lo, precip_hi = ps$precip$hi,
              use_soilm = FALSE, use_kbdi = v$moisture == "kbdi", kbdi_lo = g$kbdi$lo, kbdi_hi = g$kbdi$hi,
              kbdi_R_in = map_mm[i] / 25.4, combine = "product", window = ps$smooth, gsi_max = 1,
              greenup = ps$greenup, persist = g$persist, lhfm_lo = 30, lhfm_hi = 250, lwfm_lo = 60,
              lwfm_hi = 200, gate = TRUE)
    pt <- run_gsi(met, lat[i], p)
    gr <- gsi_variant_metrics(gsi_subindices(D, ps, g, photo, kb), v, ps, g, dd, iy)
    d_gsi <- max(abs(pt$gsi - gr$gsi[1, ]))
    n_phase <- sum((pt$phase == "greenup") != gr$phase[1, ])
    # per-year values in plain R from the point run, vs the file written by 12_gsi.R
    gy <- pt$gsi[iy]; ph <- pt$phase == "greenup"
    sw <- which(ph[iy] & !ph[iy - 1])
    plain <- c(gsi_mean = mean(gy), gsi_max = max(gy), peak_doy = dd[iy][which(gy >= max(gy) - 1e-6)[1]],
               green_days = sum(ph[iy]), n_pulses = length(sw), gu_doy = if (length(sw)) dd[iy][sw[1]] else NA)
    f <- feature_year_file(cfg, paste0("gsi_", nm), Y)
    filev <- if (file.exists(f)) unlist(terra::extract(terra::rast(f)[[names(plain)]], cp$cell[i])) else rep(NA, length(plain))
    d_file <- max(abs(plain - filev), na.rm = TRUE)
    na_mismatch <- sum(is.na(plain) != is.na(filev))
    res <- rbind(res, data.frame(site = cp$name[i], variant = nm, max_abs_gsi = d_gsi, phase_mismatch = n_phase,
                                 max_abs_year_values = d_file, na_mismatch = na_mismatch,
                                 gsi_max = plain[["gsi_max"]], gu_doy = plain[["gu_doy"]], green_days = plain[["green_days"]]))
  }
}
dir.create(file.path(cfg$run_dir, "eval"), showWarnings = FALSE)
utils::write.csv(res, file.path(cfg$run_dir, "eval", "gsi_check.csv"), row.names = FALSE)
for (k in seq_len(nrow(res))) with(res[k, ], log_msg(sprintf(
  "  %-38s %-20s daily GSI diff %.1e | phase mismatches %d | year values diff %.1e (NA mismatches %d) | max %.2f, green-up %s, %d green days",
  substr(site, 1, 38), variant, max_abs_gsi, phase_mismatch, max_abs_year_values, na_mismatch, gsi_max,
  ifelse(is.na(gu_doy), "none", as.character(gu_doy)), green_days)))
ok <- all(res$max_abs_gsi < 1e-9 & res$phase_mismatch == 0 & res$max_abs_year_values < 1e-4 & res$na_mismatch == 0)
log_msg("=== ", if (ok) "PASS" else "FAIL", ": grid GSI ", if (ok) "matches" else "does NOT match",
        " gsi-scout run_gsi() at ", nrow(cp), " site(s) x ", length(variants), " variants | log: ", log_file)
