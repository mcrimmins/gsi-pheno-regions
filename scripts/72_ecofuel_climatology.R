#!/usr/bin/env Rscript
# Step 72 (fuel phenology by ecoregion x fuel type; laptop, a few minutes): satellite phenology
# climatologies for every ecoregion x fuel group unit, and a check against NPN ground dates.
# Decided Oct 26 2026 (Mike; config ecofuel: climatology): units = GR and GS in every Level II
# ecoregion, plus SH in Cold Deserts (10.1); a unit's cells are those where the group holds >= 50 %
# of the land (>= 75 % as a sensitivity); a unit needs >= 200 such cells; agriculture is context only.
#
# Inputs (read only):
#   ecofuel/fuel_groups.tif, eco_l2.tif, eco_l2_lookup.csv          (70)
#   cutb/metrics/metrics_<Y>.tif, 2001-2025                            (06: per-year satellite season)
#   summaries/cutb_{ndvi,ndii}_curve_median.tif, cutb_metrics_median.tif (06: across-year medians)
#   static/herb_mask.tif                                               (03)
#   npn_study/data/npn_site_ecoregions.csv, npn_site_fuelmodels.csv    (NPN study, in git)
#   <run_dir>/npn_study/record/series_typical.csv                      (NPN study; laptop; optional)
#
# Per unit and metric (dates are days from Jan 1, so curing can run past 365):
#   typical      median over the unit's cells of each cell's median over the years (the rules
#                package's typical date, aggregated)
#   q25 / q75    spread of those cell typicals across the unit (spatial spread)
#   two_season   share of cells whose typical is > two_season_days from the unit's typical (a unit
#                mixing two seasons, e.g. spring vs monsoon green-up)
#   yty_sd / yty_p10 / yty_p90  year-to-year spread: each year's median cell departure (value - the
#                cell's own median; departures beyond max_jump days left out as season jumps, as 33)
#   sync         share of cell-year departure variance shared by the unit that year (1 - within-year
#                SS / total SS): high = the unit moves together, so a unit-level adjustment could work
#   coverage     share of cell-years with a value; jumps = share left out as season jumps
# Ground check: each NPN site in the usable record whose fuel group (59b rule: point if burnable, else
# the 90 m mode) and ecoregion match a unit; site date = median over its species series (greenup vs
# sos20, curing vs eos50, cured vs eos20). Compares the unit's typical date with the site's own 4 km
# cell's typical date: is the unit table nearly as good as the per-cell map?
#
# Outputs <run_dir>/ecofuel/climatology/
#   units.csv           unit, code, name, group, cells (main / sensitivity), land_km2, herb_mask share, eligible
#   climatology.csv     unit x metric x version (p50 main, p75 sensitivity): the numbers above
#   years.csv           unit x metric x year: median departure, cells (main version)
#   curves.csv          unit x index (ndvi, ndii) x doy: median, q25, q75 of the cells' median curves
#   ground_sites.csv    site x target: unit, ground date, unit typical, own-cell typical, errors
#   ground_summary.csv  unit x target: sites, ground median / IQR, satellite typical, bias, MAE (unit vs cell)
# Usage: Rscript scripts/72_ecofuel_climatology.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "ecofuel.R"))

cfg <- load_cfg(profile)
log_file <- log_init("72_ecofuel_climatology", profile)
ec <- cfg$ecofuel; cl <- ec$climatology
lev <- as.integer(ec$ecoregions$level %||% 2)
yrs <- seq(as.integer(cl$years$start), as.integer(cl$years$end))
mets <- unlist(cl$metrics)
date_mets <- intersect(mets, c("sos20", "sos50", "ndvi_peak_doy", "eos50", "eos20", "cure50", "ndii_peak_doy"))
max_jump <- as.numeric(cl$max_jump %||% 90); two_d <- as.numeric(cl$two_season_days %||% 60)
min_cells <- as.integer(cl$min_cells %||% 200)
odir <- ecofuel_dir(cfg, "climatology"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 72_ecofuel_climatology | profile: ", profile, " | purity ", cl$purity, " (sensitivity ", cl$purity_sensitivity,
        ") | min cells ", min_cells, " | years ", min(yrs), "-", max(yrs))

# ---- units ---------------------------------------------------------------------------------------
UM <- ecofuel_unit_cells(cfg, as.numeric(cl$purity))
US <- ecofuel_unit_cells(cfg, as.numeric(cl$purity_sensitivity))
L <- utils::read.csv(ecofuel_dir(cfg, sprintf("eco_l%d_lookup.csv", lev)), stringsAsFactors = FALSE, colClasses = c(code = "character"))
mask <- grid_mask(cfg)
ca <- terra::values(terra::cellSize(mask, unit = "km"), mat = FALSE)
hm_f <- file.path(static_dir(cfg), "herb_mask.tif")
hm <- if (file.exists(hm_f)) { h <- terra::values(terra::rast(hm_f), mat = FALSE); !is.na(h) & h > 0 } else rep(NA, terra::ncell(mask))
units <- sort(unique(UM$unit))
UT <- do.call(rbind, lapply(units, function(u) {
  i <- UM$unit == u; j <- US$unit == u
  data.frame(unit = u, code = UM$code[i][1], name = L$name[match(UM$code[i][1], L$code)], group = UM$group[i][1],
             cells = sum(i), cells_p75 = sum(j), land_km2 = round(sum(ca[UM$cell[i]])),
             herb_mask_share = mean(hm[UM$cell[i]]), stringsAsFactors = FALSE)
}))
UT$eligible <- UT$cells >= min_cells
num <- suppressWarnings(as.numeric(UT$code)); UT <- UT[order(num, UT$group), ]
utils::write.csv(UT, file.path(odir, "units.csv"), row.names = FALSE)
log_msg(sprintf("units: %d with any cell, %d eligible (>= %d cells at >= %.0f %%): %s", nrow(UT), sum(UT$eligible), min_cells,
                100 * as.numeric(cl$purity), paste(UT$unit[UT$eligible], collapse = ", ")))
keep_units <- UT$unit[UT$eligible]
UM <- UM[UM$unit %in% keep_units, ]; US <- US[US$unit %in% keep_units, ]
cells <- sort(unique(UM$cell))                         # p75 cells are a subset of p50 cells

# ---- per-year metrics at unit cells ----------------------------------------------------------------
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
have <- yrs[file.exists(vapply(yrs, mfile, ""))]
if (!length(have)) stop("no cutb/metrics/metrics_<Y>.tif for ", min(yrs), "-", max(yrs), " under ", cfg$run_dir)
if (length(have) < length(yrs)) log_warn("missing metric years: ", paste(setdiff(yrs, have), collapse = ","))
X <- setNames(lapply(mets, function(m) matrix(NA_real_, length(cells), length(have), dimnames = list(NULL, have))), mets)
for (k in seq_along(have)) {
  r <- terra::rast(mfile(have[k]))
  miss <- setdiff(mets, names(r)); if (length(miss)) stop("metrics_", have[k], " lacks ", paste(miss, collapse = ", "))
  v <- terra::values(r[[mets]])[cells, , drop = FALSE]
  for (m in mets) X[[m]][, k] <- v[, m]
}
log_msg(sprintf("read %d years x %d metrics at %d unit cells", length(have), length(mets), length(cells)))

# cell typicals and departures
TYP <- sapply(mets, function(m) apply(X[[m]], 1, stats::median, na.rm = TRUE))
if (is.null(dim(TYP))) TYP <- matrix(TYP, nrow = 1, dimnames = list(NULL, mets))
TYP[!is.finite(TYP)] <- NA
DEP <- lapply(mets, function(m) {
  d <- X[[m]] - TYP[, m]
  if (m %in% date_mets) d[abs(d) > max_jump] <- NA
  d
}); names(DEP) <- mets

# ---- climatology per unit ------------------------------------------------------------------------
sync_share <- function(D) {                            # 1 - within-year SS / total SS over valid cell-years
  ok <- !is.na(D); if (sum(ok) < 10) return(NA_real_)
  yrm <- colMeans(D, na.rm = TRUE)
  tot <- sum((D[ok] - mean(D[ok]))^2)
  within <- sum((D - matrix(yrm, nrow(D), ncol(D), byrow = TRUE))[ok]^2)
  if (tot <= 0) NA_real_ else 1 - within / tot
}
summ <- function(U, version) {
  do.call(rbind, lapply(keep_units, function(u) {
    i <- match(U$cell[U$unit == u], cells)
    do.call(rbind, lapply(mets, function(m) {
      ty <- TYP[i, m]; D <- DEP[[m]][i, , drop = FALSE]; Xm <- X[[m]][i, , drop = FALSE]
      tu <- stats::median(ty, na.rm = TRUE)
      ya <- apply(D, 2, function(z) if (sum(!is.na(z)) >= 5) stats::median(z, na.rm = TRUE) else NA_real_)
      data.frame(unit = u, version = version, metric = m, cells = length(i), cells_valid = sum(!is.na(ty)),
                 typical = tu, q25 = stats::quantile(ty, 0.25, na.rm = TRUE, names = FALSE),
                 q75 = stats::quantile(ty, 0.75, na.rm = TRUE, names = FALSE),
                 two_season = if (m %in% date_mets) mean(abs(ty - tu) > two_d, na.rm = TRUE) else NA_real_,
                 yty_sd = stats::sd(ya, na.rm = TRUE), yty_p10 = stats::quantile(ya, 0.1, na.rm = TRUE, names = FALSE),
                 yty_p90 = stats::quantile(ya, 0.9, na.rm = TRUE, names = FALSE), sync = sync_share(D),
                 coverage = mean(!is.na(Xm)), jumps = if (m %in% date_mets) mean(is.na(D) & !is.na(Xm)) else 0,
                 stringsAsFactors = FALSE)
    }))
  }))
}
CL <- rbind(summ(UM, "p50"), summ(US, "p75"))
utils::write.csv(CL, file.path(odir, "climatology.csv"), row.names = FALSE)

YR <- do.call(rbind, lapply(keep_units, function(u) {
  i <- match(UM$cell[UM$unit == u], cells)
  do.call(rbind, lapply(mets, function(m) {
    D <- DEP[[m]][i, , drop = FALSE]
    data.frame(unit = u, metric = m, year = have, departure = apply(D, 2, function(z) if (sum(!is.na(z)) >= 5) stats::median(z, na.rm = TRUE) else NA_real_),
               cells = colSums(!is.na(D)))
  }))
}))
utils::write.csv(YR, file.path(odir, "years.csv"), row.names = FALSE)

# ---- median curves -------------------------------------------------------------------------------
CV <- NULL
for (ix in c("ndvi", "ndii")) {
  f <- file.path(cfg$run_dir, "summaries", sprintf("cutb_%s_curve_median.tif", ix))
  if (!file.exists(f)) { log_warn("no ", basename(f), ": ", ix, " curves skipped"); next }
  r <- terra::rast(f)
  doy <- suppressWarnings(as.numeric(gsub("[^0-9]", "", names(r))))
  if (anyNA(doy) || length(unique(doy)) != terra::nlyr(r)) doy <- seq(1, by = 8, length.out = terra::nlyr(r))
  V <- terra::values(r)[cells, , drop = FALSE]
  for (u in keep_units) {
    i <- match(UM$cell[UM$unit == u], cells)
    q <- apply(V[i, , drop = FALSE], 2, stats::quantile, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    CV <- rbind(CV, data.frame(unit = u, index = ix, doy = doy, q25 = q[1, ], median = q[2, ], q75 = q[3, ]))
  }
}
if (!is.null(CV)) utils::write.csv(CV, file.path(odir, "curves.csv"), row.names = FALSE)

# ---- NPN ground check (read only) ----------------------------------------------------------------
rec_f <- file.path(cfg$run_dir, "npn_study", "record", "series_typical.csv")
fe <- here::here("npn_study", "data", "npn_site_ecoregions.csv"); ff <- here::here("npn_study", "data", "npn_site_fuelmodels.csv")
gmap <- unlist(cl$ground %||% list(greenup = "sos20", curing = "eos50", cured = "eos20"))
gmap <- gmap[gmap %in% mets]
if (file.exists(rec_f) && file.exists(fe) && file.exists(ff) && lev == 2) {
  ST <- utils::read.csv(rec_f, stringsAsFactors = FALSE)
  SE <- utils::read.csv(fe, stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
  SF <- utils::read.csv(ff, stringsAsFactors = FALSE)
  burn <- !is.na(SF$fbfm_group) & SF$fbfm_group != "NB"
  SF$group <- ifelse(burn, SF$fbfm_group, SF$fbfm_3x3_group)
  sites <- unique(ST[, c("site_id", "lat", "lon")]); sites <- sites[!duplicated(sites$site_id), ]
  sites$code <- SE$eco_l2_code[match(sites$site_id, SE$site_id)]
  sites$group <- SF$group[match(sites$site_id, SF$site_id)]
  sites$unit <- ifelse(paste(sites$code, sites$group) %in% keep_units, paste(sites$code, sites$group), NA)
  sites$cell <- terra::cellFromXY(mask, cbind(sites$lon, sites$lat))
  med_f <- file.path(cfg$run_dir, "summaries", "cutb_metrics_median.tif")
  CM <- if (file.exists(med_f)) terra::values(terra::rast(med_f)[[unname(gmap)]])[sites$cell, , drop = FALSE] else NULL
  if (is.null(CM)) log_warn("no summaries/cutb_metrics_median.tif: own-cell comparison left NA")
  GS <- list()
  for (tg in names(gmap)) {
    if (!tg %in% names(ST)) next
    a <- ST[!is.na(ST[[tg]]), c("site_id", "pathway", tg)]
    if (!nrow(a)) next
    s <- do.call(rbind, lapply(split(a, a$site_id), function(z)
      data.frame(site_id = z$site_id[1], series = nrow(z), ground = stats::median(z[[tg]]),
                 pathway = if (length(unique(z$pathway)) == 1) z$pathway[1] else "mixed")))
    k <- match(s$site_id, sites$site_id)
    s$unit <- sites$unit[k]; s$target <- tg; s$metric <- gmap[[tg]]
    ty <- CL[CL$version == "p50" & CL$metric == gmap[[tg]], ]
    s$unit_typical <- ty$typical[match(s$unit, ty$unit)]
    s$cell_typical <- if (is.null(CM)) NA else CM[k, gmap[[tg]]]
    s$err_unit <- cdiff(s$unit_typical, s$ground); s$err_cell <- cdiff(s$cell_typical, s$ground)
    GS[[tg]] <- s
  }
  GS <- do.call(rbind, GS)
  utils::write.csv(GS, file.path(odir, "ground_sites.csv"), row.names = FALSE)
  G1 <- GS[!is.na(GS$unit), ]
  GG <- do.call(rbind, lapply(split(G1, list(G1$unit, G1$target), drop = TRUE), function(z) {
    both <- !is.na(z$err_unit) & !is.na(z$err_cell)
    data.frame(unit = z$unit[1], target = z$target[1], metric = z$metric[1], sites = nrow(z),
               ground_median = stats::median(z$ground), ground_q25 = stats::quantile(z$ground, 0.25, names = FALSE),
               ground_q75 = stats::quantile(z$ground, 0.75, names = FALSE), sat_typical = z$unit_typical[1],
               bias = stats::median(z$err_unit, na.rm = TRUE), mae_unit = mean(abs(z$err_unit[both])),
               mae_cell = mean(abs(z$err_cell[both])), pairs = sum(both))
  }))
  allrow <- do.call(rbind, lapply(split(G1, G1$target), function(z) {
    both <- !is.na(z$err_unit) & !is.na(z$err_cell)
    data.frame(unit = "all units", target = z$target[1], metric = z$metric[1], sites = nrow(z), ground_median = NA,
               ground_q25 = NA, ground_q75 = NA, sat_typical = NA, bias = stats::median(z$err_unit, na.rm = TRUE),
               mae_unit = mean(abs(z$err_unit[both])), mae_cell = mean(abs(z$err_cell[both])), pairs = sum(both))
  }))
  GG <- rbind(GG, allrow)
  utils::write.csv(GG, file.path(odir, "ground_summary.csv"), row.names = FALSE)
  log_msg(sprintf("NPN record sites: %d | in a climatology unit: %d (%s)", nrow(sites), sum(!is.na(sites$unit)),
                  paste(names(table(sites$unit)), table(sites$unit), collapse = ", ")))
  for (i in which(GG$unit == "all units")) with(GG[i, ], log_msg(sprintf("  ground %-7s vs %-6s: %d sites | bias %+.0f d | MAE unit typical %.0f d vs own cell %.0f d (%d pairs)",
    target, metric, sites, bias, mae_unit, mae_cell, pairs)))
} else log_warn("NPN record (", rec_f, ") or site tables not here: ground check skipped (run on the laptop)")

# ---- log -----------------------------------------------------------------------------------------
fmt_d <- function(x) ifelse(is.na(x), "  -", sprintf("%3.0f", x))
for (u in keep_units) {
  z <- CL[CL$unit == u & CL$version == "p50", ]; g <- function(m, f) z[[f]][z$metric == m]
  log_msg(sprintf("  %-8s %5d cells | green-up %s [%s-%s] sd %s | peak %s | curing %s [%s-%s] sd %s sync %.2f | NDII cure %s | two-season gu %.0f%% cu %.0f%%",
    u, z$cells[1], fmt_d(g("sos20", "typical")), fmt_d(g("sos20", "q25")), fmt_d(g("sos20", "q75")), fmt_d(g("sos20", "yty_sd")),
    fmt_d(g("ndvi_peak_doy", "typical")), fmt_d(g("eos50", "typical")), fmt_d(g("eos50", "q25")), fmt_d(g("eos50", "q75")),
    fmt_d(g("eos50", "yty_sd")), g("eos50", "sync"), fmt_d(g("cure50", "typical")), 100 * g("sos20", "two_season"), 100 * g("eos50", "two_season")))
}
log_msg(sprintf("=== done | %d units | %.1f min | %s | log: %s", length(keep_units),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
