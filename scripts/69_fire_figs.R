#!/usr/bin/env Rscript
# Step 69 (fire activity study; laptop): figures for report/fire_study.qmd from the synced CSVs.
# Same conventions as 90_report_figs.R: <fig_dir>/<figure>_<profile>.png + manifest.csv rows
# (fig_dir = report/figs on the laptop; <run_dir>/report_figs where config report: figs_in_run_dir).
#   fire_seasonal  per pilot region: mean weekly fires (as a share of the region's peak week) and
#                  the typical share of herbaceous fuel cured, one 0-1 axis, by week (65 / 66)
#   fire_skill     held-out gain from adding phenology, per region and outcome, over the simple
#                  (60-day water balance) and the 1000-h baselines, with bootstrap 90 % intervals (66)
#   fire_checks    pooled phenology gain per year in every check and variant (66 skill.csv, shuffle.csv):
#                  fire starts, human + unknown starts, fire growth
#   fire_maps      2 x 2 maps: pilot regions; fire density (FPA-FOD fires >= 10 ac, no debris, 2001-2024,
#                  per 1,000 km2 per year, ~20 km blocks); phenology gain per year by region for fire growth
#                  and for human + unknown starts (66 skill.csv)
#   fire_maps_gsi  where the satellite and the GSI each add: gain of adding the satellite to the GSI and the
#                  GSI to the satellite, fire growth, by region (66 loyo_scores.csv, *_both_nfdrs_none)
# Usage: Rscript scripts/69_fire_figs.R [profile] [figure ...]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

all_figs <- c("fire_seasonal", "fire_skill", "fire_checks", "fire_maps", "fire_maps_gsi")
want <- if (length(args) >= 2) args[-1] else all_figs
cfg <- load_cfg(profile)
log_file <- log_init("69_fire_figs", profile)
rep <- cfg$report; fs <- cfg$fire_study
in_run <- rep$figs_in_run_dir
if (is.list(in_run)) in_run <- in_run[[if (.Platform$OS.type == "windows") "windows" else "linux"]]
fig_dir <- if (isTRUE(in_run)) file.path(cfg$run_dir, "report_figs") else here::here("report", "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e4e3df"; bg <- "#fcfcfb"
c1 <- "#2a78d6"; c2 <- "#eb6834"                 # categorical slots 1, 2 (validated pair)
reg_lab <- c(`3` = "3 Monsoon Southwest", `6` = "6 Interior West", `12` = "12 Interior West",
             `11` = "11 California annual grasslands", `10` = "10 Northern and central Great Plains",
             `13` = "13 Upper Midwest")
regs <- as.integer(unlist(fs$pilot_regions))

record <- function(name, years) {
  mf <- file.path(fig_dir, "manifest.csv")
  m <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else
    data.frame(figure = character(), profile = character(), extent = character(), years = character(), made = character())
  m <- m[!(m$figure == name & m$profile == profile), , drop = FALSE]
  m <- rbind(m, data.frame(figure = name, profile = profile, extent = rep$extent_label %||% "CONUS pilot regions",
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
mon_at <- (as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.numeric(as.Date("2001-01-01")))) / 7 + 1
mon_lab <- substr(month.abb, 1, 1)

# ---- seasonal curves ------------------------------------------------------------------------------
sf <- fire_dir(cfg, "test", "seasonal.csv")
if ("fire_seasonal" %in% want && file.exists(sf)) {
  S <- utils::read.csv(sf)
  plot_png("fire_seasonal", "2001-2024", 2000, 1500, function() {
    par(mfrow = c(3, 2), mar = c(2.6, 3.2, 2.6, 0.8), oma = c(0, 0, 2.2, 0))
    for (r in regs) {
      s <- S[S$region == r, ]; s <- s[order(s$week), ]
      f <- s$n_fires / max(s$n_fires)
      plot(NA, xlim = c(1, 52), ylim = c(0, 1.05), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
      abline(h = c(0, 0.5, 1), col = grid_col, lwd = 1)
      axis(1, at = mon_at, labels = mon_lab, tick = FALSE, cex.axis = 0.85, line = -0.6)
      axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), tick = FALSE, cex.axis = 0.85, line = -0.6)
      lines(s$week, s$cured_share_typ, col = c2, lwd = 2)
      lines(s$week, f, col = c1, lwd = 2)
      mtext(reg_lab[as.character(r)], side = 3, adj = 0, line = 0.6, col = ink, font = 2, cex = 0.85)
      mtext(sprintf("peak %.0f fires per week", max(s$n_fires)), side = 3, adj = 1, line = 0.6, col = ink2, cex = 0.7)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    legend("top", inset = 0.005, horiz = TRUE, bty = "n", lwd = 2, col = c(c1, c2), text.col = ink, cex = 0.9,
           legend = c("Fires of 10 acres or more, share of the peak week (2001-2024 mean)",
                      "Herbaceous fuel cured, typical year (share of cells)"))
  })
}

# ---- held-out skill ---------------------------------------------------------------------------------
kf <- fire_dir(cfg, "test", "skill.csv")
if ("fire_skill" %in% want && file.exists(kf)) {
  K <- utils::read.csv(kf, stringsAsFactors = FALSE)
  outs <- c(all = "Fire starts, all causes\n(2001-2024)", human_unk = "Fire starts, human + unknown\n(2001-2024)",
            lightning = "Fire starts, lightning\n(2001-2019)", px = "Fire growth, satellite pixel-days\n(2012-2024)")
  rows <- c(as.character(regs), "pooled")
  rlab <- c(reg_lab, pooled = "All pilot regions")
  plot_png("fire_skill", "2001-2024", 2200, 1150, function() {
    layout(matrix(1:5, 1), widths = c(1.55, 1, 1, 1, 1))
    par(mar = c(3.4, 0.5, 3.6, 0.5), oma = c(0, 0, 0, 0))
    plot(NA, xlim = c(0, 1), ylim = c(0.5, length(rows) + 0.5), axes = FALSE, xlab = "", ylab = "")
    text(1, rev(seq_along(rows)), rlab[rows], adj = 1, col = ink, cex = 0.85, font = ifelse(rows == "pooled", 2, 1))
    for (o in names(outs)) {
      d <- K[K$spec == o, ]
      pick <- function(cmp) { x <- d[d$comparison == cmp, ]; x[match(rows, x$region), ] }
      a <- pick("pheno_simple"); b <- pick("pheno_dfm")
      per <- function(x, v) x[[v]] / x$n_years
      lim <- range(c(per(a, "lo"), per(a, "hi"), per(b, "lo"), per(b, "hi"), 0), na.rm = TRUE)
      lim <- lim + c(-1, 1) * 0.05 * diff(lim)
      y <- rev(seq_along(rows))
      plot(NA, xlim = lim, ylim = c(0.5, length(rows) + 0.5), yaxt = "n", xlab = "", ylab = "", bty = "n", cex.axis = 0.8, tcl = -0.2, mgp = c(2, 0.4, 0))
      abline(v = 0, col = ink2, lwd = 1)
      abline(h = y, col = grid_col, lwd = 0.8)
      ok <- !is.na(a$gain)
      segments(per(a, "lo")[ok], y[ok] + 0.13, per(a, "hi")[ok], y[ok] + 0.13, col = c1, lwd = 2)
      points(per(a, "gain")[ok], y[ok] + 0.13, pch = 21, bg = c1, col = bg, cex = 1.5, lwd = 1.5)
      ok <- !is.na(b$gain)
      segments(per(b, "lo")[ok], y[ok] - 0.13, per(b, "hi")[ok], y[ok] - 0.13, col = c2, lwd = 2)
      points(per(b, "gain")[ok], y[ok] - 0.13, pch = 22, bg = c2, col = bg, cex = 1.4, lwd = 1.5)
      mtext(outs[[o]], side = 3, line = 0.4, col = ink, cex = 0.72, font = 2)
      mtext("gain per year (log score)", side = 1, line = 2, col = ink2, cex = 0.62)
    }
    par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    legend("bottomleft", inset = c(0.005, 0.005), bty = "n", pch = c(21, 22), pt.bg = c(c1, c2), col = bg, pt.cex = 1.3,
           text.col = ink, cex = 0.75, legend = c("Phenology added to the 60-day water balance baseline",
                                                 "Phenology added to the 1000-hour fuel moisture baseline"))
  })
}
# ---- checks -------------------------------------------------------------------------------------------
if ("fire_checks" %in% want && file.exists(kf)) {
  K <- utils::read.csv(kf, stringsAsFactors = FALSE)
  K <- K[K$region == "pooled" & K$comparison == "pheno_simple", ]
  sf2 <- fire_dir(cfg, "test", "shuffle.csv")
  SH <- if (file.exists(sf2)) { x <- utils::read.csv(sf2); x[x$region == "pooled", ] } else NULL
  rows <- list(c("", "Satellite cured share (main result)"), c("_nb", "Burned cells left out"), c("_hs", "This week's VPD added"),
               c("SHUFFLE", "Years shuffled (20 times)"), c("_dates", "Green-up and curing dates instead"),
               c("_gsi_nfdrs_none", "GSI (NFDRS2016) instead of satellite"), c("_gsi_nfdrs_none_hs", "GSI, with this week's VPD"),
               c("_both_nfdrs_none", "Satellite and GSI together"), c("_fc", "Curing forecast 2 weeks ahead"),
               c("_lag2", "Satellite curing 2 weeks late"), c("_rt", "Real-time, raw composites"), c("_rts", "Real-time, smoothed to date"))
  outs <- c(all = "Fire starts, all causes", human_unk = "Human + unknown starts", px = "Fire growth (2012-2024)")
  plot_png("fire_checks", "2001-2024", 2200, 1250, function() {
    layout(matrix(1:4, 1), widths = c(1.6, 1, 1, 1))
    par(mar = c(3.6, 0.5, 3.4, 0.8))
    n <- length(rows); y <- rev(seq_len(n))
    plot(NA, xlim = c(0, 1), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "")
    text(1, y, vapply(rows, `[`, "", 2), adj = 1, col = ink, cex = 0.9, font = ifelse(seq_len(n) == 1, 2, 1))
    for (o in names(outs)) {
      v <- t(vapply(rows, function(r) {
        if (r[1] == "SHUFFLE") { if (is.null(SH) || !any(SH$spec == o)) return(c(NA, NA, NA)); g <- SH$gain[SH$spec == o]
          return(c(stats::median(g), min(g), max(g)) / ifelse(o == "px", 13, 24)) }
        x <- K[K$spec == paste0(o, r[1]), ]; if (!nrow(x)) return(c(NA, NA, NA))
        c(x$gain, x$lo, x$hi) / x$n_years }, numeric(3)))
      lim <- range(c(v, 0), na.rm = TRUE); lim <- lim + c(-1, 1) * 0.05 * diff(lim)
      plot(NA, xlim = lim, ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "", ylab = "", bty = "n", cex.axis = 0.85, tcl = -0.2, mgp = c(2, 0.4, 0))
      abline(v = 0, col = ink2); abline(h = y, col = grid_col, lwd = 0.8)
      if (!is.na(v[1, 1])) abline(v = v[1, 1], col = c1, lty = 3)
      ok <- !is.na(v[, 1]); col <- ifelse(seq_len(n) == 1, c1, ifelse(rows_kind <- vapply(rows, `[`, "", 1) == "SHUFFLE", "#9a9893", ink2))
      segments(v[ok, 2], y[ok], v[ok, 3], y[ok], col = col[ok], lwd = 2)
      points(v[ok, 1], y[ok], pch = 21, bg = col[ok], col = bg, cex = 1.6, lwd = 1.5)
      if (any(!ok)) text(mean(lim), y[!ok], "not run", col = ink2, cex = 0.75)
      mtext(outs[[o]], side = 3, line = 0.6, col = ink, cex = 0.8, font = 2)
      mtext("gain per year (held-out log score)", side = 1, line = 2.2, col = ink2, cex = 0.65)
    }
  })
}
# ---- maps ------------------------------------------------------------------------------------------------
map_setup <- function() {
  grid <- terra::rast(file.path(cfg$run_dir, "static", "grid_mask.tif"))
  RG <- fire_region_values(cfg, grid, fs$partition %||% "A2_groups", as.integer(fs$k %||% 13))
  R <- terra::rast(grid); terra::values(R) <- ifelse(RG$v %in% regs, RG$v, NA)
  P <- terra::as.polygons(R, dissolve = TRUE); names(P) <- "region"
  pal <- tryCatch({ D <- region_display(cfg, grid, list(list(partition = fs$partition %||% "A2_groups", k = as.integer(fs$k %||% 13))))[[1]]
    stats::setNames(D$pal, seq_along(D$pal)) }, error = function(e) NULL)
  st <- if (requireNamespace("maps", quietly = TRUE)) maps::map("state", plot = FALSE) else NULL
  list(grid = grid, R = R, P = P, pal = pal, st = st, land = grid, asp = 1 / cos(mean(as.vector(terra::ext(grid))[3:4]) * pi / 180))
}
base_map <- function(M, title, sub) {
  e <- as.vector(terra::ext(M$grid))
  plot(NA, xlim = e[1:2], ylim = e[3:4], asp = M$asp, axes = FALSE, xlab = "", ylab = "")
  terra::plot(M$land, col = grid_col, legend = FALSE, add = TRUE, axes = FALSE)
  mtext(title, side = 3, adj = 0, line = 1.2, col = ink, font = 2, cex = 0.95)
  mtext(sub, side = 3, adj = 0, line = 0.2, col = ink2, cex = 0.72)
}
outline <- function(M, lwd = 1.4) {
  if (!is.null(M$st)) lines(M$st$x, M$st$y, col = "#ffffff", lwd = 0.6)
  terra::lines(M$P, col = ink, lwd = lwd)
}
lab_regions <- function(M, labs) {
  ct <- terra::centroids(M$P, inside = TRUE); xy <- terra::crds(ct)
  for (k in seq_len(nrow(xy))) { rg <- as.character(M$P$region[k]); if (is.na(labs[rg])) next
    text(xy[k, 1], xy[k, 2], labs[rg], col = ink, cex = 0.8, font = 2) }
}
gain_cols <- function(v) {                       # diverging: orange below zero, blue above, grey near zero
  br <- c(-Inf, -1, -0.25, 0.25, 1, 2, Inf)
  pal <- c("#c4541f", "#f3b58f", "#d9d7d2", "#9ec5f4", "#3987e5", "#184f95")
  list(col = pal[findInterval(v, br, left.open = TRUE)], br = br, pal = pal,
       lab = c("below -1", "-1 to -0.25", "about 0", "0.25 to 1", "1 to 2", "above 2"))
}
gain_map <- function(M, g, title, sub) {
  base_map(M, title, sub)
  gc <- gain_cols(g)
  rv <- terra::values(M$R, mat = FALSE); cv <- rep(NA_character_, length(rv))
  for (rg in names(g)) cv[rv %in% as.integer(rg)] <- if (is.na(g[[rg]])) "#f2f1ee" else gc$col[match(rg, names(g))]
  r2 <- M$R; terra::values(r2) <- ifelse(is.na(cv), NA, match(cv, unique(stats::na.omit(cv))))
  uc <- unique(stats::na.omit(cv))
  if (length(uc)) terra::plot(r2, col = uc, legend = FALSE, add = TRUE, axes = FALSE)
  outline(M)
  lab_regions(M, stats::setNames(ifelse(is.na(g), "not tested", sprintf("%+.1f", g)), names(g)))
  legend("bottomleft", fill = gc$pal, border = NA, legend = gc$lab, bty = "n", cex = 0.7, text.col = ink2, title = "gain per year", title.col = ink2)
}
if ("fire_maps" %in% want && file.exists(kf)) {
  M <- map_setup()
  K <- utils::read.csv(kf, stringsAsFactors = FALSE)
  per <- function(spec) { x <- K[K$spec == spec & K$comparison == "pheno_simple" & K$region != "pooled", ]
    stats::setNames(sapply(as.character(regs), function(r) { y <- x[x$region == r, ]; if (nrow(y)) y$gain / y$n_years else NA }), regs) }
  F <- tryCatch(fire_read_table(fire_dir(cfg, "fod", "fod_fires")), error = function(e) NULL)
  w <- 1300; h <- round(w / 2 * diff(as.vector(terra::ext(M$grid))[3:4]) / diff(as.vector(terra::ext(M$grid))[1:2]) * M$asp * 2 + 260)
  plot_png("fire_maps", "2001-2024", 2 * w, h, function() {
    par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2.6, 0.5))
    # (a) pilot regions
    base_map(M, "Pilot regions", "A2 regions (13-region map numbers); herbaceous phenology from each region's grass and shrub cells")
    pc <- if (!is.null(M$pal)) M$pal[as.character(regs)] else rep(c1, length(regs))
    r2 <- M$R; terra::values(r2) <- match(terra::values(M$R, mat = FALSE), regs)
    terra::plot(r2, col = grDevices::adjustcolor(pc, 0.85), legend = FALSE, add = TRUE, axes = FALSE)
    outline(M)
    lab_regions(M, stats::setNames(c("3 Monsoon SW", "6 Interior West", "12 Interior West", "11 California", "10 Great Plains", "13 Upper Midwest")[match(regs, c(3, 6, 12, 11, 10, 13))], regs))
    # (b) fire density
    base_map(M, "Wildfires of 10 acres or more", "FPA-FOD, 2001-2024, debris burning removed; per 1,000 km2 per year, ~20 km blocks")
    if (!is.null(F)) {
      F <- F[F$c_plus & !F$debris & F$year %in% seq(as.integer(fs$fod$years$start), as.integer(fs$fod$years$end)), ]
      ag <- terra::aggregate(M$grid, fact = 5, fun = "max", na.rm = TRUE)
      cnt <- terra::rasterize(terra::xyFromCell(M$grid, F$cell), ag, fun = "length", background = 0)
      ar <- terra::cellSize(ag, unit = "km")
      dens <- terra::mask(cnt / ar * 1000 / length(unique(F$year)), ag)
      br <- c(0, 0.1, 0.3, 1, 3, 10, Inf); pal <- grDevices::colorRampPalette(c("#e6eefa", "#9ec5f4", "#3987e5", "#184f95", "#0b2a57"))(6)
      terra::plot(terra::classify(dens, cbind(br[-7], br[-1], 1:6), include.lowest = TRUE), col = pal, legend = FALSE, add = TRUE, axes = FALSE)
      legend("bottomleft", fill = pal, border = NA, legend = c("< 0.1", "0.1-0.3", "0.3-1", "1-3", "3-10", "> 10"), bty = "n", cex = 0.7, text.col = ink2, title = "fires / 1,000 km2 / yr", title.col = ink2)
    }
    outline(M, 1.1)
    # (c), (d) phenology gain
    gain_map(M, per("px"), "Phenology's gain: fire growth", "Held-out gain per year from the weekly cured share, 2012-2024 (upper Midwest not tested)")
    gain_map(M, per("human_unk"), "Phenology's gain: human-caused starts", "Human + unknown-cause fires, 2001-2024 (California not split by cause)")
  })
}
if ("fire_maps_gsi" %in% want) {
  lf <- fire_dir(cfg, "test", "loyo_scores.csv")
  if (file.exists(lf)) {
    M <- map_setup()
    L <- utils::read.csv(lf)
    sc <- function(spec, model) { x <- L[L$spec == spec & L$model == model, ]; stats::aggregate(ll ~ region + year, x, sum) }
    both <- sc("px_both_nfdrs_none", "M2"); gsi <- sc("px_gsi_nfdrs_none", "M2"); sat <- sc("px", "M2")
    k <- function(d) paste(d$region, d$year)
    yrs_g <- unique(gsi$year)
    gain_of <- function(a, b) stats::setNames(sapply(as.character(regs), function(r) {
      i <- a$region == as.integer(r) & a$year %in% yrs_g; if (!any(i)) return(NA)
      sum(a$ll[i]) - sum(b$ll[match(k(a[i, ]), k(b))]) }) / length(yrs_g), regs)
    s_over_g <- gain_of(both, gsi); g_over_s <- gain_of(both, sat)
    w <- 1300; h <- round(w * diff(as.vector(terra::ext(M$grid))[3:4]) / diff(as.vector(terra::ext(M$grid))[1:2]) * M$asp + 200)
    plot_png("fire_maps_gsi", "2012-2024", 2 * w, h, function() {
      par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.6, 0.5))
      gain_map(M, s_over_g, "What the satellite adds to the GSI", "Fire growth: gain per year from adding the satellite cured share to a model with the GSI one")
      gain_map(M, g_over_s, "What the GSI adds to the satellite", "Fire growth: gain per year from adding the NFDRS2016 GSI cured share to a model with the satellite one")
    })
  }
}
log_msg("=== done | figures in ", fig_dir)
