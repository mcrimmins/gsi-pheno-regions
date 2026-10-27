#!/usr/bin/env Rscript
# Step 75 (fuel phenology by ecoregion x fuel type; laptop, < 1 min): fire activity as a check on the
# climatology units (decided Oct 26-27 2026, Mike). Reads the fire study's FPA-FOD table read-only.
#
# Inputs: <climatology folder>/cells.csv, climatology.csv, years.csv (72; default folder config ecofuel:
# fire: climatology = climatology_eco3_herb); fire_study/fod/fod_fires (60).
# Fires = FPA-FOD wildfires >= 10 ac, debris burning removed (the fire study's primary outcome),
# discovered in fire_years (2001-2024) on a unit's cells (4 km cell of the fire point).
#
# Seasonal: each unit's typical cured window runs from its typical curing date (eos50) to its next typical
#   green-up (sos20). concentration = share of the unit's fires discovered in that window / the window's
#   share of the year (1 = no seasonality relative to the phenology; > 1 = fires gather in the cured window).
# Year to year: fire anomaly = log(1 + fires in the unit-year's typical cured window) minus the unit's mean;
#   phenology = the unit-year's median curing (and green-up) departure from 72's years.csv. Both have the
#   all-units mean for that year removed (national reporting changes, e.g. 2020-21, and continent-wide
#   years). r per unit (units with >= min_fires_per_year in the window), pooled r over units, 90 % bootstrap
#   over years. No dryness baseline here (the fire study tests phenology over dryness in its pilot regions),
#   so this is a consistency check, not an effect estimate. Expected sign: later curing -> shorter cured
#   window -> fewer fires (negative r) where fires follow cured fuel.
# Outputs <run_dir>/ecofuel/fire/<climatology folder>/: windows.csv, seasonal.csv (unit x week), yty.csv,
#   pooled.csv.
# Usage: Rscript scripts/75_ecofuel_fire.R [profile] [climatology folder]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "ecofuel.R"))

cfg <- load_cfg(profile)
log_file <- log_init("75_ecofuel_fire", profile)
ec <- cfg$ecofuel; fc <- ec$fire %||% list()
folder <- if (length(args) >= 2 && nzchar(args[2])) args[2] else (fc$climatology %||% "climatology_eco3_herb")
cdir <- ecofuel_dir(cfg, folder)
odir <- ecofuel_dir(cfg, "fire", folder); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
yrs <- ec$fire_years %||% list(start = 2001, end = 2024)
min_fy <- as.numeric(fc$min_fires_per_year %||% 5); n_boot <- as.integer(fc$n_boot %||% 2000)
log_msg("=== 75_ecofuel_fire | profile: ", profile, " | ", folder, " | years ", yrs$start, "-", yrs$end)

need <- file.path(cdir, c("cells.csv", "climatology.csv", "years.csv"))
if (!all(file.exists(need))) stop("missing in ", cdir, ": ", paste(basename(need[!file.exists(need)]), collapse = ", "), " (run 72)")
UC <- utils::read.csv(need[1], stringsAsFactors = FALSE)
CL <- utils::read.csv(need[2], stringsAsFactors = FALSE); CL <- CL[CL$version == "p50", ]
YR <- utils::read.csv(need[3], stringsAsFactors = FALSE)
fstem <- fire_dir(cfg, "fod", "fod_fires")
if (is.na(fire_table_file(fstem))) stop("no fire_study/fod/fod_fires table (60)")
F0 <- fire_read_table(fstem)
d <- as.Date(F0$date); yr <- as.integer(format(d, "%Y"))
keep <- F0$c_plus %in% TRUE & !(F0$debris %in% TRUE) & yr >= yrs$start & yr <= yrs$end
F <- data.frame(cell = F0$cell[keep], year = yr[keep], doy = as.integer(format(d[keep], "%j")))
F$unit <- UC$unit[match(F$cell, UC$cell)]
F <- F[!is.na(F$unit), ]
n_years <- yrs$end - yrs$start + 1
log_msg(sprintf("fires on unit cells: %d of %d C+ wildfires (no debris)", nrow(F), sum(keep)))

units <- unique(UC$unit)
typ <- function(u, m) { v <- CL$typical[CL$unit == u & CL$metric == m]; if (length(v)) v[1] else NA }
in_window <- function(doy, cu, gu) {                 # cured window: curing -> next green-up, on a 365-day circle
  cu <- ((round(cu) - 1) %% 365) + 1; gu <- ((round(gu) - 1) %% 365) + 1
  if (cu <= gu) doy >= cu & doy < gu else doy >= cu | doy < gu
}
WN <- do.call(rbind, lapply(units, function(u) {
  cu <- typ(u, "eos50"); gu <- typ(u, "sos20"); f <- F[F$unit == u, ]
  if (is.na(cu) || is.na(gu)) return(NULL)
  w_days <- sum(in_window(1:365, cu, gu))
  inw <- in_window(f$doy, cu, gu)
  data.frame(unit = u, cells = sum(UC$unit == u), fires = nrow(f), fires_per_yr = nrow(f) / n_years,
             fires_per_1000km2_yr = NA, curing_doy = cu, greenup_doy = gu, window_days = w_days,
             share_fires_cured = if (nrow(f)) mean(inw) else NA, share_year_cured = w_days / 365,
             concentration = if (nrow(f)) mean(inw) / (w_days / 365) else NA,
             fires_cured_per_yr = sum(inw) / n_years, stringsAsFactors = FALSE)
}))
mask <- grid_mask(cfg); ca <- terra::values(terra::cellSize(mask, unit = "km"), mat = FALSE)
WN$fires_per_1000km2_yr <- WN$fires_per_yr / (vapply(WN$unit, function(u) sum(ca[UC$cell[UC$unit == u]]), 0) / 1000)
utils::write.csv(WN, file.path(odir, "windows.csv"), row.names = FALSE)
F$week <- pmin(52L, (F$doy - 1L) %/% 7L + 1L)
SE <- as.data.frame(table(unit = factor(F$unit, levels = units), week = factor(F$week, levels = 1:52)), stringsAsFactors = FALSE)
names(SE)[3] <- "fires"; utils::write.csv(SE, file.path(odir, "seasonal.csv"), row.names = FALSE)

# ---- year to year ----------------------------------------------------------------------------------
yy <- seq(as.integer(yrs$start), as.integer(yrs$end))
use <- WN$unit[!is.na(WN$fires_cured_per_yr) & WN$fires_cured_per_yr >= min_fy]
A <- do.call(rbind, lapply(use, function(u) {
  cu <- typ(u, "eos50"); gu <- typ(u, "sos20"); f <- F[F$unit == u, ]
  n <- vapply(yy, function(y) sum(f$year == y & in_window(f$doy, cu, gu)), 0)
  dep_c <- YR$departure[YR$unit == u & YR$metric == "eos50"][match(yy, YR$year[YR$unit == u & YR$metric == "eos50"])]
  dep_g <- YR$departure[YR$unit == u & YR$metric == "sos20"][match(yy, YR$year[YR$unit == u & YR$metric == "sos20"])]
  a <- log1p(n); data.frame(unit = u, year = yy, fires = n, anom = a - mean(a), dep_curing = dep_c, dep_greenup = dep_g)
}))
if (!is.null(A) && nrow(A)) {
  ym <- function(x) x - stats::ave(x, A$year, FUN = function(z) mean(z, na.rm = TRUE))
  A$anom_y <- ym(A$anom); A$dep_c_y <- ym(A$dep_curing); A$dep_g_y <- ym(A$dep_greenup)
  rr <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (sum(ok) >= 8) stats::cor(x[ok], y[ok]) else NA_real_ }
  YT <- do.call(rbind, lapply(use, function(u) { z <- A[A$unit == u, ]
    data.frame(unit = u, years = sum(!is.na(z$dep_c_y)), fires_cured_per_yr = mean(z$fires),
               r_curing = rr(z$dep_c_y, z$anom_y), r_greenup = rr(z$dep_g_y, z$anom_y)) }))
  utils::write.csv(YT, file.path(odir, "yty.csv"), row.names = FALSE)
  set.seed(cfg$clustering$seed %||% 1)
  boot <- function(v) {
    est <- rr(A[[v]], A$anom_y)
    bs <- replicate(n_boot, { s <- sample(yy, replace = TRUE); B <- do.call(rbind, lapply(s, function(y) A[A$year == y, ])); rr(B[[v]], B$anom_y) })
    c(r = est, lo = unname(stats::quantile(bs, 0.05, na.rm = TRUE)), hi = unname(stats::quantile(bs, 0.95, na.rm = TRUE)))
  }
  PO <- rbind(data.frame(phenology = "curing", t(boot("dep_c_y"))), data.frame(phenology = "greenup", t(boot("dep_g_y"))))
  PO$units <- length(use); PO$unit_years <- sum(!is.na(A$dep_c_y))
  PO$units_negative <- c(sum(YT$r_curing < 0, na.rm = TRUE), sum(YT$r_greenup < 0, na.rm = TRUE))
  utils::write.csv(PO, file.path(odir, "pooled.csv"), row.names = FALSE)
  for (i in seq_len(nrow(PO))) with(PO[i, ], log_msg(sprintf("year to year, %s departure vs cured-window fires: pooled r %.2f [%.2f, %.2f] over %d units (%d negative)",
                                                             phenology, r, lo, hi, units, units_negative)))
} else log_warn("no unit with >= ", min_fy, " fires per year in its cured window: year-to-year check skipped")

z <- WN[order(-WN$fires_per_yr), ]
log_msg(sprintf("seasonal: median concentration %.2f over %d units with fires (> 1 in %d); share of all unit fires in cured windows %.0f %%",
                stats::median(WN$concentration, na.rm = TRUE), sum(WN$fires > 0), sum(WN$concentration > 1, na.rm = TRUE),
                100 * sum(WN$share_fires_cured * WN$fires, na.rm = TRUE) / sum(WN$fires)))
for (i in seq_len(min(12, nrow(z)))) with(z[i, ], log_msg(sprintf("  %-18s %6.1f fires/yr | cured window %3d d (%2.0f %% of year) holds %2.0f %% of fires | concentration %.2f",
                                                                     unit, fires_per_yr, window_days, 100 * share_year_cured, 100 * share_fires_cured, concentration)))
log_msg("=== done | ", odir, " | log: ", log_file)
