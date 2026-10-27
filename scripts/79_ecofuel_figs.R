#!/usr/bin/env Rscript
# Step 79 (fuel phenology by ecoregion x fuel type; laptop, seconds): page figures and tables for
# report/ecofuel.qmd from 70 / 71 outputs. Same conventions as 69 / 90: <fig_dir>/<figure>_<profile>.png
# + manifest.csv rows (fig_dir = report/figs on the laptop; <run_dir>/report_figs where config
# report: figs_in_run_dir). Copies the small inventory tables to report/tables/ecofuel_*.csv (counts and
# areas only; no site ids or names) so the page renders without the data root.
#   ecofuel_units      2 maps: dominant fuel group per 4 km cell (largest share of land) with ecoregion
#                      outlines; share of the land held by that group (purity)
#   ecofuel_inventory  per ecoregion: land area split by fuel group (bars), NPN grass sites in GR / GS
#                      fuel and FPA-FOD fires >= 10 ac per year (columns of numbers)
#   ecofuel_climatology  per climatology unit (72): typical green-up, NDVI peak and curing (eos50) with the
#                      spread of the unit's cells (25-75 %), NPN ground medians where >= 3 sites, and the
#                      year-to-year spread of curing; also copies 72's unit, climatology and ground-summary
#                      tables (no site ids) to report/tables
#   ecofuel_compare    73: how well each grouping (ecoregion, ecoregion x fuel, + season split, + elevation,
#                      Level III, climate regions) explains satellite typical dates (block CV), shares
#                      year-to-year variation, and predicts NPN ground dates (clustered CV)
#   ecofuel_climatology_<scheme>  the climatology figure for every other 72 scheme folder (eco3_herb[_season])
#   ecofuel_fire_<folder>         75: fire concentration in each unit's typical cured window, and the
#                      year-to-year correlation of curing with cured-window fires (pooled band)
# Usage: Rscript scripts/79_ecofuel_figs.R [profile] [figure ...]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "ecofuel.R"))

all_figs <- c("ecofuel_units", "ecofuel_inventory", "ecofuel_climatology", "ecofuel_compare", "ecofuel_fire")
want <- if (length(args) >= 2) args[-1] else all_figs
cfg <- load_cfg(profile)
log_file <- log_init("79_ecofuel_figs", profile)
rep <- cfg$report; ec <- cfg$ecofuel %||% list()
lev <- as.integer(ec$ecoregions$level %||% 2)
in_run <- rep$figs_in_run_dir
if (is.list(in_run)) in_run <- in_run[[if (.Platform$OS.type == "windows") "windows" else "linux"]]
fig_dir <- if (isTRUE(in_run)) file.path(cfg$run_dir, "report_figs") else here::here("report", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
idir <- ecofuel_dir(cfg, "inventory")
if (!file.exists(file.path(idir, "units.csv"))) stop("no ", file.path(idir, "units.csv"), " (run 71 first)")

# ---- tables for the page ---------------------------------------------------------------------------
tdir <- here::here("report", "tables"); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
for (f in c("units", "ecoregions", "models", "purity")) {
  src <- file.path(idir, paste0(f, ".csv"))
  if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_%s_%s.csv", f, profile)), overwrite = TRUE)
}
ss <- ecofuel_dir(cfg, "shares_summary.csv")
if (file.exists(ss)) file.copy(ss, file.path(tdir, sprintf("ecofuel_shares_%s.csv", profile)), overwrite = TRUE)
cdir <- ecofuel_dir(cfg, "climatology")
for (f in c("satellite", "ground", "groupings", "season_split")) {
  src <- ecofuel_dir(cfg, "compare", paste0(f, ".csv"))
  if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_cmp_%s_%s.csv", f, profile)), overwrite = TRUE)
}
for (f in c("units", "climatology", "ground_summary")) {
  src <- file.path(cdir, paste0(f, ".csv"))
  if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_clim_%s_%s.csv", f, profile)), overwrite = TRUE)
}
# other unit schemes (72 eco3_herb [+ season]): ecofuel_<folder>_<table>_<profile>.csv; fire check (75)
clim_dirs <- basename(list.dirs(ecofuel_dir(cfg), recursive = FALSE))
clim_dirs <- clim_dirs[grepl("^climatology_", clim_dirs)]
for (cd in clim_dirs) for (f in c("units", "climatology", "ground_summary")) {
  src <- ecofuel_dir(cfg, cd, paste0(f, ".csv"))
  if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_%s_%s_%s.csv", cd, f, profile)), overwrite = TRUE)
}
for (cd in c("climatology", clim_dirs)) for (f in c("windows", "yty", "pooled")) {
  src <- ecofuel_dir(cfg, "fire", cd, paste0(f, ".csv"))
  if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_fire_%s_%s_%s.csv", cd, f, profile)), overwrite = TRUE)
}
src <- ecofuel_dir(cfg, "compare", "satellite_oddeven.csv")
if (file.exists(src)) file.copy(src, file.path(tdir, sprintf("ecofuel_cmp_satellite_oddeven_%s.csv", profile)), overwrite = TRUE)
log_msg("tables -> report/tables/ecofuel_*_", profile, ".csv")

# ---- drawing helpers (as 69) ---------------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; bg <- "#fcfcfb"
# categorical slots 1-5 in fixed order (validated: CVD worst adjacent dE 9.1, normal 19.6; contrast
# relief = legend + tables), nonburnable in neutral greys, SB folded into "other burnable"
gcol <- c(GR = "#2a78d6", GS = "#eb6834", SH = "#1baf7a", TU = "#eda100", TL = "#e87ba4", SB = "#6250d6",
          NB_ag = "#d6d2c8", NB_urban = "#8e8c86", NB_barren = "#b3b0a8", NB_snow = "#b3b0a8")
glab <- ecofuel_group_labels()
record <- function(name, years) {
  mf <- file.path(fig_dir, "manifest.csv")
  m <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else
    data.frame(figure = character(), profile = character(), extent = character(), years = character(), made = character())
  m <- m[!(m$figure == name & m$profile == profile), , drop = FALSE]
  m <- rbind(m, data.frame(figure = name, profile = profile, extent = rep$extent_label %||% "CONUS",
                           years = years, made = format(Sys.Date())))
  utils::write.csv(m[order(m$figure, m$profile), ], mf, row.names = FALSE)
  log_msg("wrote ", name, "_", profile, ".png")
}
plot_png <- function(name, years, w, h, code) {
  f <- file.path(fig_dir, sprintf("%s_%s.png", name, profile))
  grDevices::png(f, width = w, height = h, res = 170, bg = bg)
  ok <- tryCatch({ par(bg = bg, family = "sans", col.axis = ink2, fg = ink2, las = 1); code(); TRUE },
                 error = function(e) { log_err(name, ": ", conditionMessage(e)); FALSE }, finally = grDevices::dev.off())
  if (ok) record(name, years) else unlink(f)
  invisible(ok)
}
lf_ver <- ec$landfire_version %||% "LANDFIRE"

# ---- maps ----------------------------------------------------------------------------------------
dfile <- file.path(idir, "dominant.tif"); efile <- ecofuel_dir(cfg, sprintf("eco_l%d.tif", lev))
if ("ecofuel_units" %in% want && file.exists(dfile) && file.exists(efile)) {
  D <- terra::rast(dfile); DG <- utils::read.csv(file.path(idir, "dominant_groups.csv"), stringsAsFactors = FALSE)
  EO <- terra::as.polygons(terra::rast(efile)[["eco"]])
  L <- utils::read.csv(ecofuel_dir(cfg, sprintf("eco_l%d_lookup.csv", lev)), stringsAsFactors = FALSE, colClasses = c(code = "character"))
  asp <- (terra::ymax(D) - terra::ymin(D)) / (terra::xmax(D) - terra::xmin(D)) / cos(mean(c(terra::ymin(D), terra::ymax(D))) * pi / 180)
  W <- 2200; H <- round(W * asp * 2 + 520)
  plot_png("ecofuel_units", paste(lf_ver, "FBFM40; EPA Level", lev, "ecoregions"), W, H, function() {
    layout(matrix(1:4, 4), heights = c(1, 0.2, 1, 0.16))
    # dominant group
    par(mar = c(0.4, 0.4, 2.6, 0.4))
    cols <- gcol[DG$group]
    terra::plot(D[["dominant"]], col = data.frame(value = DG$dominant, col = cols), type = "classes",
                legend = FALSE, axes = FALSE, box = FALSE, mar = NA, background = bg)
    terra::lines(EO, col = ink, lwd = 0.7)
    cen <- terra::centroids(EO, inside = TRUE); cxy <- terra::crds(cen)
    text(cxy[, 1], cxy[, 2], L$code[match(cen$eco, L$eco)], cex = 0.9, col = ink, font = 2)
    title(main = "Dominant fuel group (largest share of each 4 km cell's land)", adj = 0, cex.main = 1, col.main = ink, font.main = 1, line = 0.8)
    par(mar = c(0, 0.4, 0, 0.4)); plot.new()
    show <- c("GR", "GS", "SH", "TU", "TL", "NB_ag", "NB_urban", "NB_barren")
    legend("center", legend = glab[show], fill = gcol[show], border = bg, bty = "n", ncol = 4, cex = 0.85, text.col = ink2)
    # purity
    par(mar = c(0.4, 0.4, 2.6, 0.4))
    br <- c(0, 0.25, 0.5, 0.75, 1.0001); pcol <- c("#b7d3f6", "#6da7ec", "#2a78d6", "#104281")
    terra::plot(D[["dominant_share"]], breaks = br, col = pcol, legend = FALSE, axes = FALSE, box = FALSE, mar = NA, background = bg)
    terra::lines(EO, col = ink, lwd = 0.7)
    title(main = "Share of the cell's land held by its dominant group (purity)", adj = 0, cex.main = 1, col.main = ink, font.main = 1, line = 0.8)
    par(mar = c(0, 0.4, 0, 0.4)); plot.new()
    legend("center", legend = c("< 25 %", "25-50 %", "50-75 %", ">= 75 %"), fill = pcol, border = bg, bty = "n", ncol = 4, cex = 0.85, text.col = ink2)
  })
}

# ---- inventory bars ------------------------------------------------------------------------------
if ("ecofuel_inventory" %in% want) {
  U <- utils::read.csv(file.path(idir, "units.csv"), stringsAsFactors = FALSE, colClasses = c(code = "character"))
  E <- utils::read.csv(file.path(idir, "ecoregions.csv"), stringsAsFactors = FALSE, colClasses = c(code = "character"))
  E <- E[E$land_km2 > 0, ]; E <- E[order(E$land_km2), ]                     # largest on top
  segs <- c("GR", "GS", "SH", "TU", "TL", "SB", "NB_ag", "NB_urban", "NB_barren", "NB_snow")
  A <- sapply(segs, function(g) U$land_km2[match(paste(E$eco, g), paste(U$eco, U$group))]) / 1e3
  A[is.na(A)] <- 0; if (is.null(dim(A))) A <- matrix(A, nrow = 1, dimnames = list(NULL, segs))
  npn <- vapply(E$eco, function(e) sum(U$npn_sites[U$eco == e & U$group %in% c("GR", "GS")], na.rm = TRUE), 0)
  has_npn <- any(!is.na(U$npn_sites))
  fy <- vapply(E$eco, function(e) sum(U$fires_wt_per_yr[U$eco == e], na.rm = TRUE), 0)
  has_fire <- any(!is.na(U$fires_wt_per_yr))
  n <- nrow(E)
  plot_png("ecofuel_inventory", paste(lf_ver, "FBFM40; NPN grass sites; FPA-FOD 2001-2024"), 2200, 460 + 62 * n, function() {
    par(mar = c(3.2, 19, 6.5, 12), xpd = NA)
    xmax <- max(rowSums(A)) * 1.02
    plot(NA, xlim = c(0, xmax), ylim = c(0.4, n + 0.6), yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n", xaxs = "i")
    at <- pretty(c(0, xmax)); at <- at[at <= xmax]
    segments(at, 0.4, at, n + 0.5, col = grid_col, xpd = FALSE)
    axis(1, at = at, labels = format(at, big.mark = ","), tick = FALSE, line = -0.8, cex.axis = 0.8)
    mtext("Land area (thousand km2)", side = 1, line = 1.6, cex = 0.8, col = ink2)
    gap <- xmax * 0.0015
    for (i in seq_len(n)) {
      x0 <- 0
      for (g in segs) { w <- A[i, g]; if (w > 0) rect(x0, i - 0.36, x0 + max(w - gap, 0), i + 0.36, col = gcol[g], border = NA); x0 <- x0 + w }
    }
    lab <- sprintf("%s %s", E$code, E$name); lab <- ifelse(nchar(lab) > 40, paste0(substr(lab, 1, 39), "..."), lab)
    text(-xmax * 0.01, seq_len(n), lab, adj = 1, cex = 0.78, col = ink)
    xr <- xmax * 1.04
    text(xr, n + 0.95, "NPN sites\n(GR / GS)", adj = 0, cex = 0.72, col = ink2)
    text(xr + xmax * 0.13, n + 0.95, "Fires / yr\n(>= 10 ac)", adj = 0, cex = 0.72, col = ink2)
    text(xr, seq_len(n), if (has_npn) npn else "-", adj = 0, cex = 0.78, col = ink)
    text(xr + xmax * 0.13, seq_len(n), if (has_fire) format(round(fy), big.mark = ",") else "-", adj = 0, cex = 0.78, col = ink)
    show <- c("GR", "GS", "SH", "TU", "TL", "NB_ag", "NB_urban", "NB_barren")
    legend(x = 0, y = n + 1.5, legend = glab[show], fill = gcol[show], border = bg, bty = "n", ncol = 4, cex = 0.72, text.col = ink2, yjust = 0)
  })
}
# ---- climatology timeline (72) -----------------------------------------------------------------
draw_clim <- function(cdir, figname) {
  cf <- file.path(cdir, "climatology.csv"); if (!file.exists(cf)) return(invisible())
  CL <- utils::read.csv(cf, stringsAsFactors = FALSE); CL <- CL[CL$version == "p50", ]
  UT <- utils::read.csv(file.path(cdir, "units.csv"), stringsAsFactors = FALSE, colClasses = c(code = "character"))
  UT <- UT[UT$eligible, ]
  gf <- file.path(cdir, "ground_summary.csv")
  GG <- if (file.exists(gf)) utils::read.csv(gf, stringsAsFactors = FALSE) else NULL
  n <- nrow(UT)
  gv <- function(u, m, f) { v <- CL[[f]][CL$unit == u & CL$metric == m]; if (length(v)) v[1] else NA }
  gg <- function(u, t) { if (is.null(GG)) return(c(NA, NA)); z <- GG[GG$unit == u & GG$target == t, ]; if (!nrow(z)) c(NA, 0) else c(z$ground_median[1], z$sites[1]) }
  mcol <- c(gu = "#2a78d6", cu = "#eb6834")
  plot_png(figname, "MODIS 2001-2025; NPN ground 2014-2024", 2200, 360 + 64 * n, function() {
    par(mar = c(3.2, 17, 5.5, 15), xpd = NA)
    xl <- c(-92, 400)                                   # Oct of the previous year (CA green-up) to early Feb
    plot(NA, xlim = xl, ylim = c(0.4, n + 0.6), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", xaxs = "i")
    mon <- c(-91, -61, -30, as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1, 366)
    lab <- c("Oct", "Nov", "Dec", month.abb, "Jan")
    segments(mon, 0.4, mon, n + 0.5, col = grid_col, xpd = FALSE)
    axis(1, at = mon + 15, labels = lab, tick = FALSE, line = -0.8, cex.axis = 0.7)
    segments(1, 0.4, 1, n + 0.5, col = ink2, lty = 3)
    for (i in seq_len(n)) {
      y <- n + 1 - i; u <- UT$unit[i]
      for (m in c("sos20", "eos50")) {
        cc <- if (m == "sos20") mcol["gu"] else mcol["cu"]
        segments(gv(u, m, "q25"), y, gv(u, m, "q75"), y, col = cc, lwd = 2.2, lend = 1)
        points(gv(u, m, "typical"), y, pch = 21, bg = cc, col = bg, cex = 1.5, lwd = 1.5)
      }
      points(gv(u, "ndvi_peak_doy", "typical"), y, pch = 124, col = ink2, cex = 1.1)
      for (t in c("greenup", "curing")) {
        g <- gg(u, t)
        if (!is.na(g[1]) && g[2] >= 3) points(g[1], y, pch = 5, col = if (t == "greenup") mcol["gu"] else mcol["cu"], cex = 1.2, lwd = 1.6)
      }
    }
    grp <- c(GR = "grass", GS = "grass-shrub", SH = "shrub", HERB = "herbaceous")
    sea <- ifelse(grepl(" early$", UT$unit), " (early)", ifelse(grepl(" late$", UT$unit), " (late)", ""))
    labs <- sprintf("%s %s, %s%s", UT$code, substr(UT$name, 1, 24), grp[UT$group], sea)
    text(xl[1] - 6, n + 1 - seq_len(n), labs, adj = 1, cex = 0.74, col = ink)
    xr <- xl[2] + 10
    text(xr, n + 1.05, "Cells", adj = 0, cex = 0.7, col = ink2)
    text(xr + 42, n + 1.05, "NPN sites", adj = 0, cex = 0.7, col = ink2)
    text(xr + 92, n + 1.05, "Curing\nyear sd", adj = 0, cex = 0.7, col = ink2)
    ns <- vapply(UT$unit, function(u) { g <- gg(u, "greenup")[2]; c2 <- gg(u, "curing")[2]; max(c(g, c2, 0), na.rm = TRUE) }, 0)
    text(xr, n + 1 - seq_len(n), format(UT$cells, big.mark = ","), adj = 0, cex = 0.74, col = ink)
    text(xr + 42, n + 1 - seq_len(n), if (is.null(GG)) "-" else ns, adj = 0, cex = 0.74, col = ink)
    text(xr + 92, n + 1 - seq_len(n), sprintf("%.0f d", vapply(UT$unit, function(u) gv(u, "eos50", "yty_sd"), 0)), adj = 0, cex = 0.74, col = ink)
    legend(x = xl[1], y = n + 2.2, legend = c("Green-up (NDVI 20 %)", "Curing (NDVI end of season, 50 %)", "NDVI peak", "NPN ground median (>= 3 sites)"),
           pch = c(21, 21, 124, 5), pt.bg = c(mcol, NA, NA), col = c(bg, bg, ink2, ink2), pt.cex = c(1.4, 1.4, 1, 1.1), bty = "n",
           ncol = 4, cex = 0.72, text.col = ink2, yjust = 0)
    text(xl[1], n + 0.75, "Dot = typical date (median of the unit's cells); bar = 25-75 % of cells", adj = 0, cex = 0.66, col = ink2)
  })
}
for (cd in c("climatology", clim_dirs)) {
  fn <- if (cd == "climatology") "ecofuel_climatology" else paste0("ecofuel_", cd)
  if (fn %in% want || "ecofuel_climatology" %in% want) draw_clim(ecofuel_dir(cfg, cd), fn)
}

# ---- grouping comparison (73) ---------------------------------------------------------------------
sf_ <- ecofuel_dir(cfg, "compare", "satellite.csv")
if ("ecofuel_compare" %in% want && file.exists(sf_)) {
  SA <- utils::read.csv(sf_, stringsAsFactors = FALSE)
  gf2 <- ecofuel_dir(cfg, "compare", "ground.csv"); GD <- if (file.exists(gf2)) utils::read.csv(gf2, stringsAsFactors = FALSE) else NULL
  GN <- utils::read.csv(ecofuel_dir(cfg, "compare", "groupings.csv"), stringsAsFactors = FALSE)
  glab <- c(eco2 = "Level II ecoregion", eco2_fuel = "Level II x fuel (units)", eco2_fuel_season = "Level II x fuel x season",
            eco2_fuel_elev = "Level II x fuel x 500 m elevation", eco3 = "Level III ecoregion", eco3_fuel = "Level III x fuel",
            a2 = "Climate regions (A2, 13)", a2_fuel = "Climate regions x fuel")
  gs <- intersect(names(glab), unique(SA$grouping)); n <- length(gs)
  mcol <- c("#2a78d6", "#eb6834")
  panel <- function(get, title, xl, ref = 0) {
    plot(NA, xlim = xl, ylim = c(0.5, n + 0.5), yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
    at <- pretty(xl); at <- at[at >= xl[1] & at <= xl[2]]
    segments(at, 0.5, at, n + 0.5, col = grid_col, xpd = FALSE)
    if (ref >= xl[1]) segments(ref, 0.5, ref, n + 0.5, col = ink2, lty = 3)
    axis(1, at = at, tick = FALSE, line = -0.8, cex.axis = 0.75)
    for (i in seq_len(n)) {
      v <- get(gs[i]); y <- n + 1 - i
      for (j in 1:2) if (!is.na(v[j])) points(max(min(v[j], xl[2]), xl[1]), y + c(0.12, -0.12)[j], pch = if (v[j] < xl[1]) 60 else 21,
                                             bg = mcol[j], col = if (v[j] < xl[1]) mcol[j] else bg, cex = 1.4, lwd = 1.2)
    }
    mtext(title, side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink)
  }
  sv <- function(g, m, f) { v <- SA[[f]][SA$grouping == g & SA$metric == m]; if (length(v)) v[1] else NA }
  gv2 <- function(g, t) { if (is.null(GD)) return(NA); v <- GD$cv_r2[GD$grouping == g & GD$target == t]; if (length(v)) v[1] else NA }
  plot_png("ecofuel_compare", "MODIS 2001-2025; NPN ground 2014-2024", 2200, 300 + 70 * n, function() {
    layout(matrix(1:4, 1), widths = c(1.15, 1, 1, 1))
    par(mar = c(3, 0.5, 4.2, 0.5)); plot.new(); plot.window(c(0, 1), c(0.5, n + 0.5))
    lab <- sprintf("%s (%s groups)", glab[gs], GN$groups_domain[match(gs, GN$grouping)])
    text(0.98, n + 1 - seq_len(n), lab, adj = 1, cex = 0.8, col = ink, xpd = NA)
    legend("topleft", inset = c(0, -0.09), legend = c("Green-up", "Curing"), pch = 21, pt.bg = mcol, col = bg, pt.cex = 1.3,
           bty = "n", horiz = TRUE, cex = 0.78, text.col = ink2, xpd = NA)
    par(mar = c(3, 1, 4.2, 1.2))
    panel(function(g) c(sv(g, "sos20", "cv_r2"), sv(g, "eos50", "cv_r2")), "Satellite typical dates\nexplained (R2, held-out blocks)", c(0, 1))
    panel(function(g) c(sv(g, "sos20", "sync"), sv(g, "eos50", "sync")), "Year-to-year variation\nshared within groups", c(0, 1))
    panel(function(g) c(gv2(g, "greenup"), gv2(g, "curing")), "NPN ground dates predicted\n(R2, nearby sites held out)", c(-0.5, 1))
  })
}
# ---- fire check (75) ------------------------------------------------------------------------------
draw_fire <- function(cd, figname) {
  fd <- ecofuel_dir(cfg, "fire", cd); wf <- file.path(fd, "windows.csv"); if (!file.exists(wf)) return(invisible())
  WN <- utils::read.csv(wf, stringsAsFactors = FALSE); WN <- WN[WN$fires_per_yr >= 1, ]
  YT <- if (file.exists(file.path(fd, "yty.csv"))) utils::read.csv(file.path(fd, "yty.csv"), stringsAsFactors = FALSE) else NULL
  PO <- if (file.exists(file.path(fd, "pooled.csv"))) utils::read.csv(file.path(fd, "pooled.csv"), stringsAsFactors = FALSE) else NULL
  WN <- WN[order(WN$concentration), ]; n <- nrow(WN); if (!n) return(invisible())
  plot_png(figname, "FPA-FOD 2001-2024; MODIS 2001-2025", 2200, 380 + 30 * n, function() {
    layout(matrix(1:2, 1), widths = c(1.35, 1))
    par(mar = c(3.4, 15, 4.5, 1), xpd = NA)
    xl <- c(0, max(2.5, max(WN$concentration, na.rm = TRUE) * 1.05))
    plot(NA, xlim = xl, ylim = c(0.5, n + 0.5), yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
    at <- pretty(xl); segments(at, 0.5, at, n + 0.5, col = grid_col); segments(1, 0.5, 1, n + 0.5, col = ink2, lty = 3)
    axis(1, at = at, tick = FALSE, line = -0.8, cex.axis = 0.75)
    cex_p <- 0.6 + 1.6 * sqrt(WN$fires_per_yr / max(WN$fires_per_yr))
    points(WN$concentration, seq_len(n), pch = 21, bg = "#eb6834", col = bg, cex = cex_p)
    uf <- ecofuel_dir(cfg, cd, "units.csv"); lab <- WN$unit
    if (file.exists(uf)) { UU <- utils::read.csv(uf, stringsAsFactors = FALSE); k <- match(WN$unit, UU$unit)
      sea <- ifelse(grepl(" early$", WN$unit), " (early)", ifelse(grepl(" late$", WN$unit), " (late)", ""))
      lab <- ifelse(is.na(k), WN$unit, sprintf("%s %s%s", UU$code[k], substr(UU$name[k], 1, 26), sea)) }
    text(-0.03 * xl[2], seq_len(n), lab, adj = 1, cex = 0.62, col = ink)
    mtext("Fires in the typical cured window, relative to its share of the year\n(1 = no concentration; dot size = fires per year)", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink)
    par(mar = c(3.4, 1, 4.5, 1))
    plot(NA, xlim = c(-1, 1), ylim = c(0.5, n + 0.5), yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
    segments(c(-1, -0.5, 0.5, 1), 0.5, c(-1, -0.5, 0.5, 1), n + 0.5, col = grid_col); segments(0, 0.5, 0, n + 0.5, col = ink2, lty = 3)
    axis(1, at = c(-1, -0.5, 0, 0.5, 1), tick = FALSE, line = -0.8, cex.axis = 0.75)
    if (!is.null(PO)) { p <- PO[PO$phenology == "curing", ]
      rect(p$lo, 0.5, p$hi, n + 0.5, col = "#2a78d61f", border = NA); segments(p$r, 0.5, p$r, n + 0.5, col = "#2a78d6", lwd = 2) }
    if (!is.null(YT)) { r <- YT$r_curing[match(WN$unit, YT$unit)]; ok <- !is.na(r)
      points(r[ok], seq_len(n)[ok], pch = 21, bg = "#2a78d6", col = bg, cex = 1.2) }
    mtext("Year to year: later curing vs fires in the cured window (r)\nline and band = pooled over units, 90 % interval over years", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink)
  })
}
for (cd in c("climatology", clim_dirs)) {
  fn <- paste0("ecofuel_fire_", cd)
  if (fn %in% want || "ecofuel_fire" %in% want) draw_fire(cd, fn)
}
log_msg("=== done | log: ", log_file)
