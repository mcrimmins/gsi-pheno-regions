#!/usr/bin/env Rscript
# Step 90: figures for the project page (report/index.qmd).
# Run wherever the data for a profile lives (laptop for dev, P720 for full). Each figure
# is skipped, with a log line, when its inputs don't exist yet, so this can be rerun at
# any stage of the pipeline. Outputs are small PNGs that go in git; the page embeds them,
# so rendering the page needs no data.
#
# Writes <fig_dir>/<figure>_<profile>.png and updates <fig_dir>/manifest.csv, where
# fig_dir is report/figs (dev) or <run_dir>/report_figs when config report:
# figs_in_run_dir is true (full, on the P720; 91_sync_from_p720.R copies them to the repo)
# (figure, profile, extent, years, made). The page shows the full-profile version of a
# figure when one exists, else the dev version.
#
# Figures:
#   domain          herbaceous share (core and open), strict herbaceous mask
#   block1_medians  six Block 1 features, across-year median
#   block1_iqr      year-to-year spread (IQR) of three Block 1 features
#   cutb_curves     NDVI and NDII7 composites at sample points (config report:
#                   sample_points), Terra and Aqua, last `curve_years` years
#   cutb_peak       Terra NDVI seasonal amplitude and timing of peak, latest year
#   regions         region maps (config report: region_maps; k-means rasters under clusters/)
#   eval_curves     share of the satellite-curve variance explained by each partition vs k
#                   (eval/ab/r2_curves.csv from 30_compare_ab.R)
#   eval_obs        adjusted R^2 of observed site dates on region vs k (eval/obs/obs_r2.csv
#                   from 31_obs_regions.R; sets in config report: eval_sets)
#   obs_sites       observation sites (eval/obs/sites.csv)
# eval_obs adds bootstrap ranges at report_k when eval/obs/obs_r2_boot.csv exists, and
# eval_diff_<profile>.csv (paired differences, from obs_r2_diff.csv) is written for the page.
# Also writes <fig_dir>/eval_summary_<profile>.csv (scores at evaluation: report_k), which the
# page shows as a table.
# The evaluation figures need 31's output, which is made on the laptop: run the full profile
# there (config report: figs_in_run_dir is false on Windows, so figures go to report/figs).
#
# Usage: Rscript scripts/90_report_figs.R [profile] [figure ...]
#   figure: any of the names above (default: all)
#   RStudio: gsi_profile <- "dev"; gsi_figs <- c("regions", "eval_obs"); source("scripts/90_report_figs.R")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cutb.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

all_figs <- c("domain", "block1_medians", "block1_iqr", "cutb_curves", "cutb_peak",
              "regions", "eval_curves", "eval_obs", "obs_sites", "profiles")
want <- if (length(args) >= 2) args[-1] else if (exists("gsi_figs", envir = globalenv())) {
  base::get("gsi_figs", envir = globalenv())
} else all_figs
bad <- setdiff(want, all_figs)
if (length(bad)) stop("unknown figure(s): ", paste(bad, collapse = ", "), "; use: ", paste(all_figs, collapse = ", "))
do <- function(name) name %in% want

cfg <- load_cfg(profile)
log_file <- log_init("90_report_figs", profile)
rep <- cfg$report
extent_label <- rep$extent_label %||% profile
in_run <- rep$figs_in_run_dir
if (is.list(in_run)) in_run <- in_run[[if (.Platform$OS.type == "windows") "windows" else "linux"]]
fig_dir <- if (isTRUE(in_run)) file.path(cfg$run_dir, "report_figs") else
  here::here("report", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
log_msg("=== 90_report_figs | profile: ", profile, " | extent: ", extent_label, " | ",
        paste(want, collapse = ", "), " -> ", fig_dir)

# ---- Style ------------------------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; bg <- "#fcfcfb"; landgray <- "#e4e3df"
terra_col <- "#2a78d6"; aqua_col <- "#eb6834"          # categorical slots 1 and 2
seq_pal <- function(name, n) grDevices::hcl.colors(n, name, rev = TRUE)   # light -> dark
# Share maps: blue ramp whose lightest step still reads against the page (low share is
# real land, not missing data).
blue_ramp <- function(n) grDevices::colorRampPalette(
  c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"))(n)
states <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL
land <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
asp <- 1 / cos(mean(as.vector(terra::ext(land))[3:4]) * pi / 180)  # lon/lat aspect at mid-latitude
# Figure size: panels sized to the map's aspect plus fixed title/legend margins.
# Margins: 3.2 lines top and bottom at cex 1 (~34 px per line at res 170).
dims <- function(ncol, nrow, panel_w = 1100) {
  e <- unname(as.vector(terra::ext(land)))                  # xmin, xmax, ymin, ymax
  map_h <- (panel_w - 41) * (e[4] - e[3]) / (e[2] - e[1]) * asp
  c(w = ncol * panel_w, h = unname(nrow * (map_h + 225)))
}

# Titles and colour bar are placed relative to the map extent (not the plot region, which
# can be taller than the map), in units of text-line height.
map_extent <- function() unname(as.vector(terra::ext(land)))
line_h <- function() graphics::strheight("M", cex = 1) * 1.5

# Horizontal colour bar under the lower-left of the map.
colorbar <- function(cols, labels, title) {
  e <- map_extent(); w <- e[2] - e[1]; lh <- line_h()
  x0 <- e[1] + 0.01 * w; x1 <- x0 + 0.45 * w; y1 <- e[3] - 0.5 * lh; y0 <- y1 - 0.7 * lh
  xs <- seq(x0, x1, length.out = length(cols) + 1)
  rect(xs[-length(xs)], y0, xs[-1], y1, col = cols, border = bg, lwd = 1.5, xpd = NA)
  at <- if (length(labels) == length(xs)) xs else seq(x0, x1, length.out = length(labels))
  keep <- nzchar(labels)
  # Thin labels until they fit their spacing (every 2nd, 3rd, ... non-empty label).
  if (any(keep)) {
    idx <- which(keep); gap <- (x1 - x0) / max(1, length(xs) - 1)
    wmax <- max(graphics::strwidth(labels[idx], cex = 0.8))
    k <- 1; while (wmax * 1.25 > gap * k && k < length(idx)) k <- k + 1
    keep[idx[-seq(1, length(idx), by = k)]] <- FALSE
  }
  if (any(keep)) text(at[keep], y0, labels[keep], pos = 1, cex = 0.8, col = ink2, xpd = NA)
  if (nzchar(title)) text(x1 + 0.015 * w, (y0 + y1) / 2, title, adj = 0, cex = 0.85,
                          col = ink2, xpd = NA)
}

# Text just under the map's lower-left corner (for categorical legends and notes).
under_map <- function() { e <- map_extent(); c(x = e[1] + 0.01 * (e[2] - e[1]), y = e[3] - 0.3 * line_h()) }

# One map panel: gray land, classified values, state lines, title + subtitle, colour bar.
map_panel <- function(r, title, sub, breaks, cols, labels, legend_title = "") {
  r <- terra::mask(r, land)
  terra::plot(terra::ext(land), type = "n", axes = FALSE, xlab = "", ylab = "", asp = asp,
              main = "", mar = NA)
  terra::plot(land, col = landgray, legend = FALSE, add = TRUE, axes = FALSE)
  terra::plot(r, col = cols, breaks = breaks, legend = FALSE, add = TRUE, axes = FALSE)
  if (!is.null(states)) lines(states$x, states$y, col = "#ffffff", lwd = 0.7)
  e <- map_extent(); lh <- line_h()
  text(e[1], e[4] + 1.75 * lh, title, adj = c(0, 0), cex = 1.15, font = 2, col = ink, xpd = NA)
  text(e[1], e[4] + 0.45 * lh, sub, adj = c(0, 0), cex = 0.85, col = ink2, xpd = NA)
  if (!is.null(labels)) colorbar(cols, labels, legend_title)   # NULL: caller draws a legend
}

# Breaks: linear 'pretty' over the 1-99 % range, or a 1-2-5 log sequence.
lin_breaks <- function(r, n = 8) {
  q <- terra::global(r, function(v) stats::quantile(v, c(0.01, 0.99), na.rm = TRUE))
  b <- pretty(unlist(q), n)
  list(breaks = c(-Inf, b[-c(1, length(b))], Inf), labels = c("", num_lab(b[-c(1, length(b))]), ""))
}
log_breaks <- function(r) {
  q <- unlist(terra::global(r, function(v) stats::quantile(v, c(0.01, 0.99), na.rm = TRUE)))
  q <- pmax(q, 1e-3)
  # 1-2-5 steps; if that gives fewer than 7 classes, use 1-1.5-2-3-5-7 steps
  for (m in list(c(1, 2, 5), c(1, 1.5, 2, 3, 5, 7))) {
    s <- as.vector(outer(m, 10^(-3:4)))
    b <- s[s > q[1] & s < q[2]]
    if (length(b) >= 6) break
  }
  if (length(b) < 3) return(lin_breaks(r))
  list(breaks = c(-Inf, b, Inf), labels = c("", num_lab(b), ""))
}
num_lab <- function(b) trimws(formatC(b, format = "fg", digits = 4, big.mark = ","))
fix_breaks <- function(b) list(breaks = c(-Inf, b, Inf), labels = c("", num_lab(b), ""))

# Record a figure in <fig_dir>/manifest.csv (one row per figure and profile).
record <- function(name, years) {
  mf <- file.path(fig_dir, "manifest.csv")
  m <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else
    data.frame(figure = character(), profile = character(), extent = character(),
               years = character(), made = character())
  m <- m[!(m$figure == name & m$profile == profile), , drop = FALSE]
  m <- rbind(m, data.frame(figure = name, profile = profile, extent = extent_label,
                           years = years, made = format(Sys.Date())))
  utils::write.csv(m[order(m$figure, m$profile), ], mf, row.names = FALSE)
  log_msg("wrote ", name, "_", profile, ".png")
}

# Write a map-figure PNG (ncol x nrow panels), run the plotting code, record it.
fig <- function(name, years, ncol, nrow, code, panel_w = 1100) {
  f <- file.path(fig_dir, sprintf("%s_%s.png", name, profile))
  d <- dims(ncol, nrow, panel_w)
  grDevices::png(f, width = d[["w"]], height = d[["h"]], res = 170, bg = bg)
  ok <- tryCatch({
    par(mfrow = c(nrow, ncol))                      # (mfrow shrinks cex; reset it)
    par(cex = 1, mar = c(3.2, 0.6, 3.2, 0.6), bg = bg, family = "sans")
    code(); TRUE
  }, error = function(e) { log_err(name, ": ", conditionMessage(e)); FALSE },
  finally = grDevices::dev.off())
  if (ok) record(name, years) else unlink(f)
  invisible(ok)
}
have <- function(...) { f <- c(...); ok <- all(file.exists(f))
  if (!ok) log_msg("skip: missing ", paste(basename(f[!file.exists(f)]), collapse = ", "))
  ok }
yrs_label <- function(y) if (length(y) == 1) as.character(y) else paste0(min(y), "-", max(y))

# ---- domain -------------------------------------------------------------------------------
hs_file <- file.path(static_dir(cfg), "herb_share.tif")
mask_file <- file.path(static_dir(cfg), "herb_mask.tif")
if (do("domain") && have(hs_file, mask_file)) {
  hs <- terra::rast(hs_file); hm <- terra::rast(mask_file)
  share_b <- seq(0, 1, by = 0.1)
  share_lab <- c("0", "", "0.2", "", "0.4", "", "0.6", "", "0.8", "", "1")
  fig("domain", "NLCD 2001, 2012, 2024", 3, 1, function() {
    map_panel(hs[["herb_share"]], "Herbaceous share (core)",
              "Grassland + shrub/scrub share of land", share_b,
              blue_ramp(10), share_lab)
    if ("open_herb_share" %in% names(hs)) {
      map_panel(hs[["open_herb_share"]], "Open herbaceous share",
                "Core + pasture/hay share of land", share_b,
                blue_ramp(10), share_lab)
    } else plot.new()
    map_panel(hm, "Strict herbaceous mask",
              "Core >= 0.5 and crops < 0.25", c(0.5, 1.5), terra_col, NULL)
    n <- format(as.integer(terra::global(hm, "notNA")[1, 1]), big.mark = ",")
    p <- under_map()
    text(p[["x"]], p[["y"]], paste(n, "cells in the mask"), adj = c(0, 1), cex = 0.85,
         col = ink2, xpd = NA)
  })
}

# ---- Block 1 ------------------------------------------------------------------------------
fy <- feature_years(cfg)
b1_years <- fy[file.exists(vapply(fy, function(y) feature_year_file(cfg, "block1", y), ""))]
med_file <- summary_file(cfg, "block1", "median"); iqr_file <- summary_file(cfg, "block1", "iqr")
if (do("block1_medians") && have(med_file) && length(b1_years)) {
  wy <- paste("water years", yrs_label(b1_years))
  med <- terra::rast(med_file)
  fig("block1_medians", yrs_label(b1_years), 3, 2, function() {
    b <- lin_breaks(med[["t_ann"]])
    map_panel(med[["t_ann"]], "Mean annual temperature", wy, b$breaks,
              seq_pal("YlOrRd", length(b$breaks) - 1), b$labels, "°C")
    b <- log_breaks(med[["p_ann"]])
    map_panel(med[["p_ann"]], "Annual precipitation", wy, b$breaks,
              seq_pal("YlGnBu", length(b$breaks) - 1), b$labels, "mm")
    b <- fix_breaks(seq(0.1, 0.6, by = 0.1))
    map_panel(med[["pf_jas"]], "Jul-Sep share of precipitation",
              "Monsoon / summer-rain fraction", b$breaks,
              seq_pal("BuPu", length(b$breaks) - 1), b$labels)
    # Breaks follow the UNEP aridity classes (hyper-arid < 0.05, arid < 0.2, semi-arid
    # < 0.5, dry subhumid < 0.65, humid above) and split the humid range, where the East and
    # the Plains gradient sit (most of the East is 0.9-1.5). Diverging at 0.65: six brown
    # steps below, six teal above.
    b <- fix_breaks(c(0.05, 0.1, 0.2, 0.3, 0.5, 0.65, 0.8, 1, 1.25, 1.5, 2))
    map_panel(med[["ai_ann"]], "Aridity index (P / PET)",
              "Hargreaves PET; brown = arid to dry subhumid (< 0.65)", b$breaks,
              grDevices::hcl.colors(length(b$breaks) - 1, "BrBG"), b$labels)
    b <- fix_breaks(seq(100, 350, by = 50))
    map_panel(med[["ffs_0"]], "Freeze-free season",
              "Days between last spring and first fall tmin <= 0 °C", b$breaks,
              seq_pal("Greens", length(b$breaks) - 1), b$labels, "days")
    b <- lin_breaks(med[["gdd10"]])
    map_panel(med[["gdd10"]], "Growing degree days, base 10 °C",
              "Calendar year; warm-season grass heat units", b$breaks,
              seq_pal("OrRd", length(b$breaks) - 1), b$labels, "°C-days")
  })
}
if (do("block1_iqr") && have(med_file, iqr_file) && length(b1_years) >= 2) {
  med <- terra::rast(med_file); iqr <- terra::rast(iqr_file)
  fig("block1_iqr", yrs_label(b1_years), 3, 1, function() {
    rel <- iqr[["p_ann"]] / med[["p_ann"]]
    b <- fix_breaks(seq(0.1, 0.6, by = 0.1))
    map_panel(rel, "Annual precipitation", "IQR / median across years", b$breaks,
              seq_pal("Blues", length(b$breaks) - 1), b$labels)
    b <- fix_breaks(seq(0.05, 0.35, by = 0.05))
    map_panel(iqr[["pf_jas"]], "Jul-Sep share of precipitation", "IQR across years",
              b$breaks, seq_pal("BuPu", length(b$breaks) - 1), b$labels)
    b <- fix_breaks(c(10, 20, 30, 45, 60))
    map_panel(iqr[["ffs_0"]], "Freeze-free season", "IQR across years", b$breaks,
              seq_pal("Greens", length(b$breaks) - 1), b$labels, "days")
  })
}

# ---- Cut B --------------------------------------------------------------------------------
prods <- unlist(cfg$cutb$products)
cb_years <- function(product, vars = c("ndvi", "nir", "mir")) {
  cfg$years[vapply(cfg$years, function(y)
    all(file.exists(vapply(vars, function(v) cutb_out_file(cfg, product, v, y), ""))), TRUE)]
}
# Values of one product-year at grid cells: list(dates, ndvi, ndii7) as cells x dates.
cb_cells <- function(product, year, cells) {
  rd <- function(v) as.matrix(terra::rast(cutb_out_file(cfg, product, v, year))[cells])
  r <- terra::rast(cutb_out_file(cfg, product, "ndvi", year))
  nir <- rd("nir"); mir <- rd("mir")
  d <- as.Date(names(r)); keep <- format(d, "%Y") == as.character(year)
  ndvi <- rd("ndvi")
  list(dates = d[keep], ndvi = ndvi[, keep, drop = FALSE],
       ndii7 = ((nir - mir) / (nir + mir))[, keep, drop = FALSE])
}

pts <- do.call(rbind, lapply(rep$sample_points, as.data.frame))
if (!is.null(pts) && nrow(pts)) {
  e <- as.vector(terra::ext(land))
  pts <- pts[pts$lon > e[1] & pts$lon < e[2] & pts$lat > e[3] & pts$lat < e[4], ]
  pts$cell <- terra::cellFromXY(land, as.matrix(pts[, c("lon", "lat")]))
  pts <- pts[!is.na(pts$cell) & !is.na(land[pts$cell][, 1]), ]
}
yrs_t <- if (do("cutb_curves") || do("cutb_peak")) cb_years(prods[1]) else integer()
if ((do("cutb_curves") || do("cutb_peak")) && !length(yrs_t)) log_msg("skip: no Cut B years for ", prods[1])
if (do("cutb_curves") && length(yrs_t) && !is.null(pts) && nrow(pts)) {
  cy <- utils::tail(yrs_t, rep$curve_years %||% 3)
  pts <- utils::head(pts, 4)
  series <- lapply(prods, function(p) {
    yy <- intersect(cy, cb_years(p))
    if (!length(yy)) return(NULL)
    parts <- lapply(yy, cb_cells, product = p, cells = pts$cell)
    list(product = p, dates = do.call(c, lapply(parts, `[[`, "dates")),
         ndvi = do.call(cbind, lapply(parts, `[[`, "ndvi")),
         ndii7 = do.call(cbind, lapply(parts, `[[`, "ndii7")))
  })
  series <- Filter(Negate(is.null), series)
  cols <- setNames(c(terra_col, aqua_col)[seq_along(prods)], prods)
  labs <- c(MOD13A2.061 = "Terra (MOD13A2)", MYD13A2.061 = "Aqua (MYD13A2)")
  xr <- as.Date(c(sprintf("%d-01-01", min(cy)), sprintf("%d-12-31", max(cy))))
  grDevices::png(file.path(fig_dir, sprintf("cutb_curves_%s.png", profile)),
                 width = 2400, height = 420 + 520 * nrow(pts), res = 170, bg = bg)
  ok <- tryCatch({
    par(mfrow = c(nrow(pts), 2), mar = c(2.6, 4.2, 3.4, 1), oma = c(0, 0, 2.4, 0),
        bg = bg, family = "sans", col.axis = ink2, fg = ink2)
    for (i in seq_len(nrow(pts))) for (v in c("ndvi", "ndii7")) {
      ylim <- range(unlist(lapply(series, function(s) s[[v]][i, ])), na.rm = TRUE)
      plot(xr, ylim, type = "n", xlab = "", ylab = toupper(v), las = 1, bty = "n",
           cex.axis = 0.85, xaxt = "n")
      ticks <- seq(xr[1], xr[2] + 1, by = "6 months")
      axis(1, at = ticks, labels = format(ticks, "%b %Y"), cex.axis = 0.85)
      abline(v = as.Date(sprintf("%d-01-01", cy)), col = "#e4e3df", lwd = 1)
      for (s in series) {
        lines(s$dates, s[[v]][i, ], col = cols[[s$product]], lwd = 2)
        points(s$dates, s[[v]][i, ], col = cols[[s$product]], pch = 16, cex = 0.55)
      }
      mtext(if (v == "ndvi") pts$name[i] else "", side = 3, line = 1.6, adj = 0,
            cex = 0.95, font = 2, col = ink)
      mtext(if (v == "ndvi") "NDVI: greenness" else "NDII7: canopy water / curing",
            side = 3, line = 0.4, adj = 0, cex = 0.8, col = ink2)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("top", legend = labs[names(cols)], col = cols, lwd = 2, pch = 16, horiz = TRUE,
           bty = "n", text.col = ink2, cex = 0.95)
    TRUE
  }, error = function(e) { log_err("cutb_curves: ", conditionMessage(e)); FALSE },
  finally = grDevices::dev.off())
  if (ok) record("cutb_curves", yrs_label(cy)) else
    unlink(file.path(fig_dir, sprintf("cutb_curves_%s.png", profile)))
} else if (do("cutb_curves") && length(yrs_t)) log_msg("skip cutb_curves: no sample points inside the extent")

if (do("cutb_peak") && length(yrs_t)) {
  y <- max(yrs_t)
  ndvi <- terra::rast(cutb_out_file(cfg, prods[1], "ndvi", y))
  dates <- as.Date(names(ndvi))
  # only composites starting in year y (AppEEARS adds the previous December composite)
  ndvi <- ndvi[[format(dates, "%Y") == as.character(y)]]
  dates <- as.Date(names(ndvi))
  amp <- terra::app(ndvi, function(v) {
    if (sum(!is.na(v)) < 12) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  })
  idx <- terra::which.max(ndvi)
  # Peak timing in three classes; weak seasonality (amplitude < 0.1) shown separately.
  mon <- as.integer(format(dates, "%m"))
  cls <- ifelse(mon <= 5, 1, ifelse(mon <= 7, 2, 3))
  peak <- terra::classify(idx, cbind(seq_along(cls), cls))
  peak <- terra::ifel(is.na(amp), NA, terra::ifel(amp < 0.1, 4, peak))
  fig("cutb_peak", as.character(y), 2, 1, function() {
    b <- fix_breaks(c(0.1, 0.2, 0.3, 0.4, 0.5))
    map_panel(amp, "NDVI seasonal amplitude", paste("Terra, max - min of 16-day composites,", y),
              b$breaks, seq_pal("Greens", length(b$breaks) - 1), b$labels)
    pc <- c(terra_col, aqua_col, "#1baf7a", "#b8b7b2")
    map_panel(peak, "Timing of peak NDVI", paste("Terra 16-day composites,", y),
              c(0.5, 1.5, 2.5, 3.5, 4.5), pc, NULL)
    p <- under_map()
    legend(p[["x"]], p[["y"]], xpd = NA, horiz = TRUE, bty = "n", x.intersp = 0.5,
           fill = pc, border = NA, cex = 0.8, text.col = ink2,
           legend = c("Jan-May", "Jun-Jul", "Aug-Dec", "Weak season (< 0.1)"))
  })
}

# ---- Regions and evaluation ----------------------------------------------------------------
ev <- cfg$evaluation
clus_dir <- file.path(cfg$run_dir, "clusters")
eval_dir <- file.path(cfg$run_dir, "eval")
part_file <- function(nm) file.path(clus_dir, ev$partitions[[nm]] %||% "none")
part_lab <- c(A1_groups = "Climate A1 (seasonal climate)", A2_groups = "Climate A2 (+ timing, variability)",
              A2_herbw = "Climate A2, herbaceous-weighted",
              B_raw = "Satellite curves (level + shape)", B_shape = "Satellite curves (timing only)",
              geo = "Location only (baseline)")
part_col <- c(A1_groups = terra_col, A2_groups = "#1a7f37", A2_herbw = "#7b4fb3", B_raw = aqua_col, B_shape = "#c0392b",
              geo = "#8a8985")
part_lty <- c(A1_groups = 1, A2_groups = 1, A2_herbw = 1, B_raw = 2, B_shape = 2, geo = 3)
set_lab <- c(grass_greenup = "NPN grasses: green-up (50 % green)",
             grass_curing = "NPN grasses: curing (below 50 % green)",
             grass_cured = "NPN grasses: cured",
             herb_onset = "NPN herbaceous: leaf onset", herb_end = "NPN herbaceous: leaves end",
             woody_onset = "NPN woody: leaf onset",
             lfmc_herb_peak = "Live fuel moisture, herbaceous: peak",
             lfmc_herb_decline = "Live fuel moisture, herbaceous: 50 % decline",
             lfmc_woody_peak = "Live fuel moisture, woody: peak",
             lfmc_woody_decline = "Live fuel moisture, woody: 50 % decline",
             inat_flower = "iNaturalist: flowering onset")
lab_of <- function(x, tab) ifelse(x %in% names(tab), tab[x], x)
reg_pal <- region_pal
# Line-plot figure (not a map): fixed size, the usual style, recorded in the manifest.
plot_png <- function(name, years, w, h, code) {
  f <- file.path(fig_dir, sprintf("%s_%s.png", name, profile))
  grDevices::png(f, width = w, height = h, res = 170, bg = bg)
  ok <- tryCatch({
    par(bg = bg, family = "sans", col.axis = ink2, fg = ink2, las = 1)
    code(); TRUE
  }, error = function(e) { log_err(name, ": ", conditionMessage(e)); FALSE },
  finally = grDevices::dev.off())
  if (ok) record(name, years) else unlink(f)
  invisible(ok)
}
k_panel <- function(d, ycol, main, sub, ylim, ylab, parts) {
  plot(NA, xlim = range(d$k), ylim = ylim, xlab = "number of regions", ylab = ylab, bty = "n",
       cex.axis = 0.85)
  abline(h = pretty(ylim), col = "#ecebe7", lwd = 1)
  if (ylim[1] < 0) abline(h = 0, col = "#b8b7b2")
  abline(v = unlist(ev$report_k), col = "#ecebe7", lty = 2)
  for (nm in parts) {
    e <- d[d$partition == nm, ]; e <- e[order(e$k), ]
    lines(e$k, e[[ycol]], col = part_col[[nm]] %||% ink2, lty = part_lty[[nm]] %||% 1, lwd = 2.2)
  }
  mtext(main, side = 3, line = 1.4, adj = 0, cex = 0.9, font = 2, col = ink)
  mtext(sub, side = 3, line = 0.3, adj = 0, cex = 0.75, col = ink2)
}
legend_row <- function(parts) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("top", legend = lab_of(parts, part_lab), col = part_col[parts], lty = part_lty[parts],
         lwd = 2.2, horiz = TRUE, bty = "n", text.col = ink2, cex = 0.8, seg.len = 2.5)
}
all_years <- yrs_label(cfg$years)

# regions: maps of k-means partitions, colours matched to the first panel
rm_ <- Filter(function(m) file.exists(part_file(m$partition)), rep$region_maps %||% list())
if (do("regions") && length(rep$region_maps) && !length(rm_)) log_msg("skip regions: no partition rasters")
if (do("regions") && length(rm_)) {
  land_cells <- which(!is.na(terra::values(land, mat = FALSE)))
  # shared numbering/colours (R/evaluate.R): reference = first map, regions matched to it
  RD <- region_display(cfg, land, rm_, ref = rm_[[1]])
  keep <- !vapply(RD, is.null, TRUE); rm_ <- rm_[keep]
  labs <- lapply(RD[keep], `[[`, "col")
  cp <- do.call(rbind, lapply(cfg$clustering$check_points, as.data.frame))
  nc <- min(3, length(rm_)); nr <- ceiling(length(rm_) / nc)
  fig("regions", all_years, nc, nr, function() {
    for (i in seq_along(rm_)) {
      r <- terra::rast(land); v <- rep(NA_real_, terra::ncell(land)); v[land_cells] <- labs[[i]]
      terra::values(r) <- v
      kk <- max(labs[[i]], na.rm = TRUE)
      map_panel(r, rm_[[i]]$title, "k-means; similar colours = overlapping regions",
                seq(0.5, kk + 0.5), rep_len(reg_pal, kk), NULL)
      if (!is.null(cp)) points(cp$lon, cp$lat, pch = 21, bg = "#ffffff", col = ink, cex = 0.7, lwd = 0.8)
    }
    p <- under_map()
    text(p[["x"]], p[["y"]], "White dots: grassland check sites", adj = c(0, 1), cex = 0.8, col = ink2, xpd = NA)
  })
}

# eval_curves: share of satellite curve variance explained (30_compare_ab.R)
rc_file <- file.path(eval_dir, "ab", "r2_curves.csv")
if (do("eval_curves") && have(rc_file)) {
  rc <- utils::read.csv(rc_file, stringsAsFactors = FALSE)
  if (is.null(rc$features)) rc$features <- "level_shape"
  parts <- intersect(names(part_lab), unique(rc$partition))
  fsets <- intersect(c("level_shape", "timing"), unique(rc$features))
  fs_lab <- c(level_shape = "Curves as observed (level + shape)", timing = "Timing only (curves scaled by amplitude)")
  plot_png("eval_curves", all_years, 2400, 230 + 770 * length(fsets), function() {
    par(mfrow = c(length(fsets), 2), mar = c(4.2, 4.4, 3.6, 1), oma = c(0, 0, 2.2, 0))
    for (fs in fsets) {
      d <- rc[rc$features == fs, ]
      k_panel(d[d$domain == "all", ], "r2", paste0(fs_lab[[fs]], ": all land cells"),
              "NDVI + NDII7 median seasonal curves", c(0, 1), "share of variance explained", parts)
      k_panel(d[d$domain == "herb", ], "r2", paste0(fs_lab[[fs]], ": herbaceous mask"),
              "Grass- and shrub-dominated cells", c(0, 1), "share of variance explained", parts)
    }
    legend_row(parts)
  })
}

# eval_obs: adjusted R^2 of observed site dates (31_obs_regions.R)
or_file <- file.path(eval_dir, "obs", "obs_r2.csv")
if (do("eval_obs") && have(or_file)) {
  orr <- utils::read.csv(or_file, stringsAsFactors = FALSE)
  bfile <- file.path(eval_dir, "obs", "obs_r2_boot.csv")
  bo <- if (file.exists(bfile)) utils::read.csv(bfile, stringsAsFactors = FALSE) else NULL
  sets <- intersect(unlist(rep$eval_sets %||% unique(orr$set)), unique(orr$set))
  parts <- intersect(names(part_lab), unique(orr$partition))
  yl <- range(c(0, 0.2, orr$adj_r2[orr$set %in% sets], if (!is.null(bo)) unlist(bo[bo$set %in% sets, c("lo", "hi")])), na.rm = TRUE)
  yl <- c(min(0, floor(yl[1] * 10) / 10), ceiling(yl[2] * 10) / 10)
  nc <- min(3, length(sets)); nr <- ceiling(length(sets) / nc)
  plot_png("eval_obs", "site medians across years", 800 * nc, 230 + 640 * nr, function() {
    par(mfrow = c(nr, nc), mar = c(4.2, 4.4, 3.6, 1), oma = c(0, 0, 2.2, 0))
    for (st in sets) {
      d <- orr[orr$set == st, ]
      k_panel(d, "adj_r2", lab_of(st, set_lab), sprintf("%d sites", max(d$n, na.rm = TRUE)),
              yl, "adjusted R-squared (dates ~ region)", parts)
      if (!is.null(bo)) {                      # bootstrap ranges at report_k, offset per partition
        b <- bo[bo$set == st & bo$partition %in% parts, ]
        off <- setNames(seq(-0.36, 0.36, length.out = length(parts)), parts)
        segments(b$k + off[b$partition], b$lo, b$k + off[b$partition], b$hi,
                 col = grDevices::adjustcolor(part_col[b$partition], 0.8), lwd = 2)
        points(b$k + off[b$partition], b$adj_r2, pch = 16, cex = 0.6, col = part_col[b$partition])
      }
    }
    legend_row(parts)
  })
}

# obs_sites: where the observations are
si_file <- file.path(eval_dir, "obs", "sites.csv")
if (do("obs_sites") && have(si_file)) {
  si <- utils::read.csv(si_file, stringsAsFactors = FALSE)
  grp <- c(grass = "NPN grasses (% green)", npn = "NPN other phenophases",
           lfmc = "Live fuel moisture", inat = "iNaturalist flowering")
  si$g <- ifelse(grepl("^grass_", si$set), "grass", ifelse(grepl("^lfmc_", si$set), "lfmc",
                 ifelse(si$set == "inat_flower", "inat", "npn")))
  si <- si[!duplicated(paste(si$g, round(si$lon, 2), round(si$lat, 2))), ]
  gc <- c(inat = "#b8b7b2", npn = "#7b4fb3", lfmc = aqua_col, grass = terra_col)
  fig("obs_sites", "site medians", 1, 1, function() {
    map_panel(land * NA, "Ground observation sites", "One point per place and data set", c(0, 1), landgray, NULL)
    for (g in names(gc)) { s <- si[si$g == g, ]
      points(s$lon, s$lat, pch = 16, cex = c(inat = 0.35, npn = 0.5, lfmc = 0.55, grass = 0.85)[[g]],
             col = grDevices::adjustcolor(gc[[g]], 0.8)) }
    p <- under_map()
    n <- table(factor(si$g, names(gc)))
    legend(p[["x"]], p[["y"]], xpd = NA, horiz = TRUE, bty = "n", pch = 16, col = gc,
           legend = sprintf("%s (%d)", grp[names(gc)], n), cex = 0.8, text.col = ink2, x.intersp = 0.6)
  }, panel_w = 1600)
}

# ---- Region profiles (32_region_profiles.R) ------------------------------------------------
# Per region set (config region_profiles: sets): prof_map_<id> (numbered region map),
# prof_heat_<id> (features x regions, coloured by standardized median, labelled with the
# median), prof_curves_<id> (median NDVI / NDII7 curves per region with 25-75 % bands),
# prof_obs_<id> (ground-observation dates per region) and prof_table_<id>_<profile>.csv
# (summary table for the page); crosswalk CSVs are copied as prof_crosswalk_*_<profile>.csv.
prof_root <- file.path(eval_dir, "profiles")
doy_lab <- function(d) ifelse(is.na(d), "", trimws(format(as.Date("2001-01-01") + round(d) - 1, "%b %e")))
fmt_val <- function(v, unit) {
  if (is.na(v)) return("")
  if (identical(unit, "doy")) return(doy_lab(v))
  if (identical(unit, "")) return(if (abs(v) < 10) formatC(v, format = "f", digits = 2) else num_lab(signif(v, 3)))
  num_lab(if (abs(v) >= 100) round(v) else signif(v, 3))
}
mon_axis <- function(side = 1, cex = 0.75) {
  at <- as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1
  axis(side, at = at, labels = month.abb, cex.axis = cex, tick = TRUE, las = 1)
}
pr_sets <- cfg$region_profiles$sets %||% list()
if (do("profiles") && !dir.exists(prof_root)) log_msg("skip profiles: no ", prof_root, " (run 32_region_profiles.R)")
if (do("profiles") && dir.exists(prof_root)) for (s_ in pr_sets) {
  id <- region_set_id(s_); pd <- file.path(prof_root, id)
  if (!have(file.path(pd, "regions.csv"), file.path(pd, "features.csv"))) next
  key <- utils::read.csv(file.path(pd, "regions.csv"), stringsAsFactors = FALSE)
  ft <- utils::read.csv(file.path(pd, "features.csv"), stringsAsFactors = FALSE)
  mt <- if (file.exists(file.path(pd, "metrics.csv"))) utils::read.csv(file.path(pd, "metrics.csv"), stringsAsFactors = FALSE) else NULL
  cv <- if (file.exists(file.path(pd, "curves.csv"))) utils::read.csv(file.path(pd, "curves.csv"), stringsAsFactors = FALSE) else NULL
  os <- if (file.exists(file.path(pd, "obs_sites.csv"))) utils::read.csv(file.path(pd, "obs_sites.csv"), stringsAsFactors = FALSE) else NULL
  k <- nrow(key); title <- s_$label %||% id
  min_herb <- cfg$region_profiles$min_herb_cells %||% 200
  dom_of <- function(j) if (key$herb_cells[j] >= min_herb) "herb" else "all"

  # numbered map
  RD <- region_display(cfg, land, list(s_))[[1]]
  if (!is.null(RD)) fig(paste0("prof_map_", id), all_years, 1, 1, function() {
    land_cells <- which(!is.na(terra::values(land, mat = FALSE)))
    r <- terra::rast(land); v <- rep(NA_real_, terra::ncell(land)); v[land_cells] <- RD$num
    terra::values(r) <- v
    map_panel(r, title, "Region numbers match the tables and charts below", seq(0.5, k + 0.5), key$colour, NULL)
    shadowtext <- function(x, y, lab) {
      for (dx in c(-1, 1)) for (dy in c(-1, 1)) text(x + dx * 0.06, y + dy * 0.06, lab, col = "#ffffff", font = 2, cex = 1.05)
      text(x, y, lab, col = ink, font = 2, cex = 1.05)
    }
    shadowtext(key$lon, key$lat, key$region)
  }, panel_w = 1500)

  # heatmap: features (all cells) + season metrics (all cells)
  rows <- unique(ft[, c("var", "label", "unit")])
  vals <- t(sapply(rows$var, function(v) ft$med[ft$var == v][order(ft$region[ft$var == v])]))
  if (!is.null(mt)) {
    ma <- mt[mt$domain == "all", ]
    mrows <- unique(ma[, c("var", "label", "unit")])
    mvals <- t(sapply(mrows$var, function(v) ma$med[ma$var == v][order(ma$region[ma$var == v])]))
    rows <- rbind(rows, transform(mrows, label = paste("Satellite:", label))); vals <- rbind(vals, mvals)
  }
  vals <- matrix(vals, nrow = nrow(rows))
  Z <- t(apply(vals, 1, function(v) { s <- stats::sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) v * 0 else (v - mean(v, na.rm = TRUE)) / s }))
  zc <- grDevices::colorRampPalette(c("#2a78d6", "#9ec5f4", "#f4f3ef", "#f5b386", "#c0582a"))(11)
  zb <- c(-Inf, seq(-2, 2, length.out = 10), Inf)
  nr_ <- nrow(rows)
  plot_png(paste0("prof_heat_", id), all_years, 900 + 150 * k, 260 + 62 * nr_, function() {
    par(mar = c(0.5, 22, 7, 1), xaxs = "i", yaxs = "i")
    plot(NA, xlim = c(0.5, k + 0.5), ylim = c(nr_ + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(nr_)) for (j in seq_len(k)) {
      zcol <- if (is.na(Z[i, j])) "#ffffff" else zc[findInterval(Z[i, j], zb, all.inside = TRUE)]
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = zcol, border = bg, lwd = 1.5)
      lab_ <- if (rows$var[i] == "b2med:onset_doy" && !is.na(vals[i, j]) && vals[i, j] >= 274) "none" else fmt_val(vals[i, j], rows$unit[i])
      text(j, i, lab_, cex = 0.62, col = ink)
    }
    ul <- ifelse(rows$unit %in% c("", "doy"), "", paste0(" (", rows$unit, ")"))
    axis(2, at = seq_len(nr_), labels = paste0(rows$label, ul), tick = FALSE, las = 1, cex.axis = 0.72, line = -0.6)
    rect(seq_len(k) - 0.42, 0.5 - 0.95, seq_len(k) + 0.42, 0.5 - 0.15, col = key$colour, border = NA, xpd = NA)
    text(seq_len(k), 0.5 - 0.55, key$region, col = "#ffffff", font = 2, cex = 0.8, xpd = NA)
    text(seq_len(k), 0.5 - 1.25, sprintf("%.0f%%", 100 * key$share), cex = 0.62, col = ink2, xpd = NA)
    mtext(paste0(title, ": region medians"), side = 3, line = 5.4, adj = 0, cex = 0.95, font = 2, col = ink, at = 0.5)
    mtext("Colour: below (blue) or above (orange) the average across regions; % = share of land; satellite dates for all cells",
          side = 3, line = 4.4, adj = 0, cex = 0.7, col = ink2, at = 0.5)
  })

  # season curves per region
  if (!is.null(cv) && nrow(cv)) {
    nc <- min(4, k); nr <- ceiling(k / nc)
    yl <- range(c(cv$q25, cv$q75), na.rm = TRUE)
    plot_png(paste0("prof_curves_", id), all_years, 620 * nc, 200 + 470 * nr, function() {
      par(mfrow = c(nr, nc), mar = c(2.4, 3.2, 3.2, 0.6), oma = c(0, 0, 2.6, 0))
      for (j in seq_len(k)) {
        dm <- dom_of(j)
        plot(NA, xlim = c(1, 365), ylim = yl, axes = FALSE, xlab = "", ylab = "")
        abline(h = pretty(yl), col = "#ecebe7"); mon_axis(); axis(2, cex.axis = 0.72, las = 1)
        for (ix in c("ndii", "ndvi")) {
          e <- cv[cv$region == j & cv$domain == dm & cv$index == ix, ]; e <- e[order(e$doy), ]
          if (!nrow(e)) next
          col_ <- if (ix == "ndvi") key$colour[j] else "#52514e"
          polygon(c(e$doy, rev(e$doy)), c(e$q25, rev(e$q75)), col = grDevices::adjustcolor(col_, 0.22), border = NA)
          lines(e$doy, e$med, col = col_, lwd = 2.2, lty = if (ix == "ndvi") 1 else 2)
        }
        if (dm == "herb") {                     # all-cells NDVI for comparison
          e <- cv[cv$region == j & cv$domain == "all" & cv$index == "ndvi", ]; e <- e[order(e$doy), ]
          lines(e$doy, e$med, col = ink2, lwd = 1, lty = 3)
        }
        n <- if (dm == "herb") key$herb_cells[j] else key$cells[j]
        mtext(sprintf("Region %d", j), side = 3, line = 1.3, adj = 0, cex = 0.82, font = 2, col = key$colour[j])
        mtext(sprintf("%s, %s cells", if (dm == "herb") "herbaceous cells" else "all cells (few herbaceous)",
                      format(n, big.mark = ",")), side = 3, line = 0.3, adj = 0, cex = 0.62, col = ink2)
      }
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend("top", c("NDVI (greenness), region colour", "NDII7 (canopy water / curing)", "25-75 % of cells",
                      "NDVI, all cells of the region"),
             lty = c(1, 2, NA, 3), lwd = c(2.2, 2.2, NA, 1), pch = c(NA, NA, 15, NA), pt.cex = 2, ncol = 4,
             col = c(ink, "#52514e", grDevices::adjustcolor(ink2, 0.3), ink2), bty = "n", cex = 0.85,
             text.col = ink2, x.intersp = 0.8, text.width = c(0.22, 0.2, 0.13, 0.2))
    })
  }

  # ground-observation dates per region
  if (!is.null(os) && nrow(os)) {
    sets <- intersect(unlist(cfg$region_profiles$obs_sets), unique(os$set))
    nc <- min(3, length(sets)); nr <- ceiling(length(sets) / nc)
    plot_png(paste0("prof_obs_", id), "site medians across years", 760 * nc, 200 + 560 * nr, function() {
      par(mfrow = c(nr, nc), mar = c(3.2, 4, 3.4, 0.6), oma = c(0, 0, 0, 0))
      set.seed(1)
      for (st in sets) {
        e <- os[os$set == st, ]
        yl <- range(e$doy, na.rm = TRUE); yl <- yl + c(-10, 10)
        plot(NA, xlim = c(0.5, k + 0.5), ylim = yl, axes = FALSE, xlab = "region", ylab = "")
        at <- as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1
        at <- c(at - 365, at, at + 365); lb <- rep(month.abb, 3)
        keep_ <- at >= yl[1] & at <= yl[2]
        abline(h = at[keep_], col = "#ecebe7")
        axis(2, at = at[keep_], labels = lb[keep_], cex.axis = 0.72, las = 1)
        axis(1, at = seq_len(k), cex.axis = 0.72, tick = FALSE, line = -0.8)
        for (j in seq_len(k)) {
          y <- e$doy[e$region == j]; if (!length(y)) { text(j, yl[1] + 5, "0", cex = 0.6, col = ink2); next }
          points(j + stats::runif(length(y), -0.28, 0.28), y, pch = 16, cex = 0.45,
                 col = grDevices::adjustcolor(key$colour[j], 0.55))
          segments(j - 0.35, stats::median(y), j + 0.35, stats::median(y), lwd = 2.5, col = ink)
          text(j, yl[2] - 3, length(y), cex = 0.6, col = ink2)
        }
        mtext(lab_of(st, set_lab), side = 3, line = 1.3, adj = 0, cex = 0.8, font = 2, col = ink)
        mtext("dots: sites; bar: median; number: sites", side = 3, line = 0.3, adj = 0, cex = 0.62, col = ink2)
      }
    })
  }

  # page table
  getf <- function(v) ft$med[ft$var == v][order(ft$region[ft$var == v])]
  getm <- function(v) if (is.null(mt)) rep(NA, k) else sapply(seq_len(k), function(j) {
    x <- mt$med[mt$var == v & mt$region == j & mt$domain == dom_of(j)]; if (length(x)) x else NA })
  tb <- data.frame(Region = key$region, `Share of land` = sprintf("%.0f%%", 100 * key$share),
                   `Herbaceous cells` = format(key$herb_cells, big.mark = ","), check.names = FALSE)
  add <- function(name, v, unit) tb[[name]] <<- vapply(v, fmt_val, "", unit = unit)
  add("Mean temp. (C)", getf("b1med:t_ann"), "C"); add("Precip. (mm)", getf("b1med:p_ann"), "mm")
  add("Aridity (P/PET)", getf("b1med:ai_ann"), ""); add("Jul-Sep share", getf("b1med:pf_jas"), "")
  add("Freeze-free (days)", getf("b1med:ffs_0"), "days"); add("Rain onset", getf("b2med:onset_doy"), "doy")
  add("Precip. CV", getf("b2var:cv_p_ann"), "")
  add("Green-up", getm("sos20"), "doy"); add("NDVI peak", getm("ndvi_peak_doy"), "doy")
  add("Curing", getm("cure50"), "doy"); add("NDVI amplitude", getm("ndvi_amp"), "")
  tb$`Season from` <- ifelse(key$herb_cells >= min_herb, "herbaceous cells", "all cells")
  utils::write.csv(tb, file.path(fig_dir, sprintf("prof_table_%s_%s.csv", id, profile)), row.names = FALSE)
  log_msg("wrote prof_table_", id, "_", profile, ".csv")
}
if (do("profiles") && dir.exists(prof_root)) for (f in list.files(prof_root, "^crosswalk_.*\\.csv$", full.names = TRUE)) {
  file.copy(f, file.path(fig_dir, sub("\\.csv$", sprintf("_%s.csv", profile), paste0("prof_", basename(f)))), overwrite = TRUE)
}

# summary table for the page: scores at evaluation: report_k
if (any(c("eval_curves", "eval_obs") %in% want) && (file.exists(rc_file) || file.exists(or_file))) {
  rk <- unlist(ev$report_k)
  tab <- expand.grid(partition = names(part_lab), k = rk, stringsAsFactors = FALSE)
  if (file.exists(rc_file)) {
    rc <- utils::read.csv(rc_file, stringsAsFactors = FALSE)
    if (is.null(rc$features)) rc$features <- "level_shape"
    for (fs in intersect(c("level_shape", "timing"), unique(rc$features))) for (dn in c("all", "herb")) {
      d <- rc[rc$domain == dn & rc$features == fs, ]
      cn <- if (fs == "level_shape") paste0("curves_", dn) else paste0("timing_", dn)
      tab[[cn]] <- d$r2[match(paste(tab$partition, tab$k), paste(d$partition, d$k))]
    }
  }
  if (file.exists(or_file)) {
    orr <- utils::read.csv(or_file, stringsAsFactors = FALSE)
    for (st in intersect(unlist(rep$eval_sets), unique(orr$set))) {
      d <- orr[orr$set == st, ]
      tab[[st]] <- d$adj_r2[match(paste(tab$partition, tab$k), paste(d$partition, d$k))]
    }
  }
  tab <- tab[rowSums(!is.na(tab[, -(1:2), drop = FALSE])) > 0, ]
  tab <- tab[order(tab$k, match(tab$partition, names(part_lab))), ]
  utils::write.csv(tab, file.path(fig_dir, sprintf("eval_summary_%s.csv", profile)), row.names = FALSE)
  log_msg("wrote eval_summary_", profile, ".csv")
  dfile <- file.path(eval_dir, "obs", "obs_r2_diff.csv")
  if (file.exists(dfile)) {
    dd <- utils::read.csv(dfile, stringsAsFactors = FALSE)
    dd <- dd[dd$set %in% unlist(rep$eval_sets), ]
    short <- c(A1_groups = "A1", A2_groups = "A2", A2_herbw = "A2 herb-weighted", B_raw = "Satellite (level + shape)",
               B_shape = "Satellite (timing)", geo = "location only")
    dd$comparison <- paste(lab_of(dd$first, short), "minus", lab_of(dd$second, short))
    dd$set <- lab_of(dd$set, set_lab)
    utils::write.csv(dd[, c("set", "comparison", "k", "diff", "lo", "hi", "share_above_0")],
                     file.path(fig_dir, sprintf("eval_diff_%s.csv", profile)), row.names = FALSE)
    log_msg("wrote eval_diff_", profile, ".csv")
  }
}

log_msg("=== done | figures in ", fig_dir)
