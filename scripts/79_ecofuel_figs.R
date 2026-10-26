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
# Usage: Rscript scripts/79_ecofuel_figs.R [profile] [figure ...]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "ecofuel.R"))

all_figs <- c("ecofuel_units", "ecofuel_inventory")
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
log_msg("=== done | log: ", log_file)
