#!/usr/bin/env Rscript
# Step 90 (NPN grass curing study): figures for the study page (npn_study/report/npn_study.qmd).
# Standalone (Oct 18 2026): moved out of gsi-pheno-regions scripts/90_report_figs.R; draws state
# outlines (maps package) instead of the regions grid. Each figure is skipped, with a log line,
# when its inputs don't exist yet. Outputs are small PNGs that go in git; the page embeds them.
#
# Writes npn_study/report/figs/<figure>_<profile>.png and updates manifest.csv there
# (figure, profile, extent, years, made).
#   npn_sites        study sites: grass type observed, NEON vs other, plant-years (58) or series-years (51)
#   npn_dates        green-up, 50 % curing and cured dates by ecoregion and grass type (58, else 51)
#   npn_ladder       skill of the model ladder against the typical date (52)
#   npn_variance     what the year-to-year variation is made of (53)
#   npn_signal       shared site-year curing signal vs summer VPD, within years (54)
#   npn_record       percent-green bins used, readings per plant-year, second green-ups by ecoregion (58)
#   npn_spatial      typical dates vs site climate: leave-one-site-out R2 by model; green-up vs
#                    Mar-Apr temperature; curing vs aridity (58)
#   npn_sensitivity  2016-17 check: NEON plants never above 50 % green by year (58); shared year
#                    signal with and without 2016-17 (53, sens/no1617)
#
# Usage: Rscript npn_study/scripts/90_npn_figs.R [profile] [figure ...]   (default: all)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

all_figs <- c("npn_sites", "npn_dates", "npn_ladder", "npn_variance", "npn_signal", "npn_record", "npn_spatial", "npn_sensitivity")
want <- if (length(args) >= 2) args[-1] else all_figs
bad <- setdiff(want, all_figs)
if (length(bad)) stop("unknown figure(s): ", paste(bad, collapse = ", "), "; use: ", paste(all_figs, collapse = ", "))
do <- function(name) name %in% want

cfg <- load_npn_cfg(profile)
log_file <- log_init("90_npn_figs", profile)
extent_label <- "CONUS"
fig_dir <- npn_path("report", "figs"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 90_npn_figs | profile: ", profile, " | ", paste(want, collapse = ", "), " -> ", fig_dir)

# ---- style (same as the gsi-pheno-regions page figures) ----------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; bg <- "#fcfcfb"; landgray <- "#e4e3df"
terra_col <- "#2a78d6"; aqua_col <- "#eb6834"

# Record a figure in <fig_dir>/manifest.csv (one row per figure and profile).
record <- function(name, years) {
  mf <- file.path(fig_dir, "manifest.csv")
  m <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else
    data.frame(figure = character(), profile = character(), extent = character(), years = character(), made = character())
  m <- m[!(m$figure == name & m$profile == profile), , drop = FALSE]
  m <- rbind(m, data.frame(figure = name, profile = profile, extent = extent_label, years = years, made = format(Sys.Date())))
  utils::write.csv(m[order(m$figure, m$profile), ], mf, row.names = FALSE)
  log_msg("wrote ", name, "_", profile, ".png")
}
# Fixed-size PNG, the usual style, recorded in the manifest.
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
have <- function(...) { f <- c(...); ok <- all(file.exists(f))
  if (!ok) log_msg("skip: missing ", paste(basename(f[!file.exists(f)]), collapse = ", "))
  ok }
yrs_label <- function(y) if (length(y) == 1) as.character(y) else paste0(min(y), "-", max(y))
mon_axis <- function(side = 1, cex = 0.75) {
  at <- as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1
  axis(side, at = at, labels = month.abb, cex.axis = cex, tick = TRUE, las = 1)
}
# CONUS base map: gray states (maps package), lon/lat at the mid-latitude aspect, title + subtitle.
if (!requireNamespace("maps", quietly = TRUE)) stop("needs the maps package: renv::install('maps')")
st_map <- function(sub, title = "NPN grass curing study sites") {
  s <- maps::map("state", fill = TRUE, plot = FALSE)
  plot(NA, xlim = c(-124.8, -66.9), ylim = c(24.4, 49.4), asp = 1 / cos(37 * pi / 180), axes = FALSE, xlab = "", ylab = "")
  graphics::polygon(s$x, s$y, col = landgray, border = "#ffffff", lwd = 0.7)
  mtext(title, side = 3, line = 1.9, adj = 0, cex = 1.15, font = 2, col = ink)
  mtext(sub, side = 3, line = 0.6, adj = 0, cex = 0.85, col = ink2)
}

npn_dir <- file.path(cfg$run_dir, "npn_study")
rec_dir <- file.path(npn_dir, "record"); st_csv <- file.path(rec_dir, "series_typical.csv")
rnm <- npn_region_names()
reg_order <- function(codes) { codes <- unique(codes[!is.na(codes) & codes != "none"]); codes[order(as.numeric(codes))] }
reg_lab <- function(codes, width = 20) vapply(codes, function(k) {
  nm <- if (k %in% names(rnm)) rnm[[k]] else ""; paste(strwrap(paste(k, nm), width), collapse = "\n") }, "")
TY <- if (file.exists(st_csv)) { x <- utils::read.csv(st_csv, stringsAsFactors = FALSE, colClasses = c(region = "character")); x[x$pathway %in% c("C3", "C4"), ] } else NULL
pw_col <- c(C3 = terra_col, C4 = aqua_col)
date_lab <- c(greenup = "Green-up (50 % green)", curing = "Curing (below 50 % green)", cured = "Cured")
sy_csv <- file.path(npn_dir, "npn_series_years.csv")
npn_neon <- function(site_id) {
  nf <- npn_path("data", "npn_site_names.csv")
  if (!file.exists(nf)) return(rep(NA, length(site_id)))
  n <- utils::read.csv(nf, stringsAsFactors = FALSE)
  grepl("\\.phenology\\.phe", n$site_name[match(site_id, n$site_id)])
}
if (any(c("npn_sites", "npn_dates") %in% want) && have(sy_csv)) {
  SY <- utils::read.csv(sy_csv, stringsAsFactors = FALSE, colClasses = c(region = "character"))
  SY <- SY[!(SY$crop %in% "yes") & SY$pathway %in% c("C3", "C4"), ]
  yl <- yrs_label(SY$year)
  if (do("npn_sites")) {
    if (!is.null(TY)) {                                   # relaxed set (58): plant-years per site
      st <- do.call(rbind, lapply(split(TY, TY$site_id), function(d) data.frame(
        site_id = d$site_id[1], lon = d$lon[1], lat = d$lat[1], n = sum(d$plant_years),
        type = if (all(d$pathway == "C3")) "C3" else if (all(d$pathway == "C4")) "C4" else "both")))
      st$neon <- st$site_id %in% TY$site_id[TY$group == "NEON"]; unit <- "plant-years"
    } else {
      st <- do.call(rbind, lapply(split(SY, SY$site_id), function(d) data.frame(
        site_id = d$site_id[1], lon = d$lon[1], lat = d$lat[1], n = nrow(d),
        type = if (all(d$pathway == "C3")) "C3" else if (all(d$pathway == "C4")) "C4" else "both")))
      st$neon <- npn_neon(st$site_id); unit <- "series-years"
    }
    tc <- c(C3 = terra_col, C4 = aqua_col, both = "#7b4fb3")
    plot_png("npn_sites", yl, 1600, 1080, function() {
      par(mar = c(3.4, 0.6, 3.6, 0.6))
      st_map(sprintf("%d sites, %s %s; point size = %s; filled = NEON site", nrow(st), format(sum(st$n), big.mark = ","), unit, unit))
      cx <- 0.5 + 1.6 * sqrt(st$n / max(st$n))
      filled <- st$neon %in% TRUE | is.na(st$neon)
      points(st$lon, st$lat, pch = ifelse(filled, 21, 1), bg = grDevices::adjustcolor(tc[st$type], 0.75),
             col = ifelse(filled, bg, tc[st$type]), cex = cx, lwd = ifelse(filled, 0.6, 1.6))
      k <- table(factor(st$type, names(tc)))
      legend("bottomleft", inset = c(0, -0.09), xpd = NA, horiz = TRUE, bty = "n", pch = 21, pt.bg = tc, col = bg, pt.cex = 1.4,
             legend = sprintf(c("cool-season (C3) only (%d)", "warm-season (C4) only (%d)", "both (%d)"), k),
             cex = 0.8, text.col = ink2, x.intersp = 0.6)
    })
  }
  if (do("npn_dates")) {
    regs <- if (!is.null(TY)) reg_order(TY$region) else reg_order(as.character(SY$region))
    cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
    plot_png("npn_dates", yl, 2600, 260 + 70 * length(regs) * 2, function() {
      par(mfrow = c(1, 3), mar = c(3.6, 9.2, 3.4, 0.8), oma = c(2.6, 0, 0, 0))
      for (tg in names(date_lab)) {
        M <- if (!is.null(TY)) data.frame(series = TY$series, region = TY$region, pathway = TY$pathway, x = TY[[tg]]) else
          stats::aggregate(SY[[tg]], list(series = SY$series, region = SY$region, pathway = SY$pathway), cmed)
        M$region <- as.character(M$region); M <- M[!is.na(M$x) & !is.na(M$region) & M$region != "none", ]
        rows <- expand.grid(pathway = c("C4", "C3"), region = rev(regs), stringsAsFactors = FALSE)
        rows$y <- seq_len(nrow(rows))
        plot(NA, xlim = c(1, 366), ylim = c(0.5, nrow(rows) + 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
        abline(v = as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1, col = "#ecebe7")
        abline(h = seq(2.5, nrow(rows), by = 2), col = "#e4e3df")
        mon_axis(1, 0.85)
        axis(2, at = seq(1.5, nrow(rows), by = 2), labels = reg_lab(rev(regs)), tick = FALSE, cex.axis = 0.78, padj = 0.5)
        for (i in seq_len(nrow(rows))) {
          d <- M[M$region == rows$region[i] & M$pathway == rows$pathway[i], ]
          if (!nrow(d)) next
          points(d$x, rows$y[i] + stats::runif(nrow(d), -0.22, 0.22), pch = 16, cex = 0.8,
                 col = grDevices::adjustcolor(pw_col[[rows$pathway[i]]], 0.55))
          if (nrow(d) >= 3) segments(cmed(d$x), rows$y[i] - 0.38, cmed(d$x), rows$y[i] + 0.38, lwd = 3, col = pw_col[[rows$pathway[i]]])
        }
        mtext(date_lab[[tg]], side = 3, line = 1.6, adj = 0, cex = 0.95, font = 2, col = ink)
        mtext("Points: one site and species; bar: region median", side = 3, line = 0.4, adj = 0, cex = 0.75, col = ink2)
      }
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend("bottom", legend = c("warm-season (C4)", "cool-season (C3)"), pch = 16, col = pw_col[c("C4", "C3")],
             horiz = TRUE, bty = "n", text.col = ink2, cex = 1.15, pt.cex = 1.5)
    })
  }
}
ld_csv <- file.path(npn_dir, "ladder", "scores.csv")
if (do("npn_ladder") && have(ld_csv)) {
  L <- utils::read.csv(ld_csv, stringsAsFactors = FALSE)
  L <- L[L$set %in% c("all", "lead") & L$level != "L0" & L$group %in% c("all", "pathway"), ]
  lv <- c(L1 = "L1 one driver", L2 = "L2 process rule", L3 = "L3 regression", L4 = "L4 random forest")
  plot_png("npn_ladder", "series-years", 2200, 900, function() {
    par(mfrow = c(1, 3), mar = c(4.2, 9.4, 3.4, 1), oma = c(2.2, 0, 0, 0))
    xl <- c(-0.5, 0.15)
    for (tg in names(date_lab)) {
      d <- L[L$target == tg, ]
      plot(NA, xlim = xl, ylim = c(0.5, 4.5), yaxt = "n", xlab = "skill against the typical date", ylab = "", bty = "n", cex.axis = 0.85)
      abline(v = pretty(xl), col = "#ecebe7"); abline(v = 0, col = ink2, lwd = 1.5)
      axis(2, at = 4:1, labels = lv, tick = FALSE, cex.axis = 0.85)
      for (k in seq_along(lv)) {
        y <- 5 - k; e <- d[d$level == names(lv)[k], ]
        for (g in c("C3", "C4")) { s <- e[e$value == g, "skill"]; if (length(s)) points(max(s, xl[1]), y + ifelse(g == "C3", 0.15, -0.15), pch = 16, cex = 1.2, col = pw_col[[g]]) }
        s <- e[e$group == "all", "skill"]
        if (length(s)) { points(max(s, xl[1]), y, pch = 21, bg = ink, col = bg, cex = 1.8)
          if (s < xl[1]) text(xl[1], y, sprintf("%.2f", s), pos = 4, cex = 0.75, col = ink2, offset = 0.8) }
      }
      n <- d$n[d$group == "all" & d$level == "L1"]
      mtext(date_lab[[tg]], side = 3, line = 1.6, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext(sprintf("%d series-years; typical-date error %.1f days", n, d$mae_L0[d$group == "all" & d$level == "L1"]), side = 3, line = 0.4, adj = 0, cex = 0.72, col = ink2)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = c("all grasses", "cool-season (C3)", "warm-season (C4)"), pch = 16, col = c(ink, pw_col),
           horiz = TRUE, bty = "n", text.col = ink2, cex = 0.95)
  })
}
vc_csv <- file.path(npn_dir, "variance", "components.csv")
if (do("npn_variance") && have(vc_csv)) {
  V <- utils::read.csv(vc_csv, stringsAsFactors = FALSE)
  parts <- list(`region-wide` = c("year", "region:year"), `site-wide` = "site:year", `one species` = "site:species:year", `single plant` = "plant")
  pc <- c("#184f95", "#6da7ec", "#b8b7b2", "#e4e3df")
  rows <- expand.grid(group = c("all", "C3", "C4"), target = names(date_lab), stringsAsFactors = FALSE)
  rows <- rows[nrow(rows):1, ]
  plot_png("npn_variance", "plant-years", 2000, 1050, function() {
    par(mar = c(4.4, 13, 3.4, 1.5))
    plot(NA, xlim = c(0, 1), ylim = c(0.5, nrow(rows) + 0.5), yaxt = "n", xaxt = "n", xlab = "share of the year-to-year variation", ylab = "", bty = "n")
    axis(1, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), " %"), cex.axis = 0.85)
    gl <- c(all = "all grasses", C3 = "cool-season (C3)", C4 = "warm-season (C4)")
    axis(2, at = seq_len(nrow(rows)), labels = paste0(c(greenup = "Green-up", curing = "Curing", cured = "Cured")[rows$target], ": ", gl[rows$group]),
         tick = FALSE, cex.axis = 0.82)
    for (i in seq_len(nrow(rows))) {
      d <- V[V$target == rows$target[i] & V$group == rows$group[i], ]; if (!nrow(d)) next
      sh <- vapply(parts, function(p) sum(d$share_year_to_year[d$component %in% p], na.rm = TRUE), 0)
      x0 <- c(0, cumsum(sh))[seq_along(sh)]
      rect(x0, i - 0.36, x0 + sh, i + 0.36, col = pc, border = bg, lwd = 1.5)
      lab <- ifelse(sh >= 0.07, paste0(round(100 * sh), "%"), "")
      text(x0 + sh / 2, i, lab, cex = 0.72, col = c("#ffffff", ink, ink, ink))
    }
    abline(h = c(3.5, 6.5), col = "#ecebe7")
    mtext("What a grass's early and late years are made of", side = 3, line = 1.9, adj = 0, cex = 1, font = 2, col = ink)
    legend("top", inset = c(0, -0.07), legend = names(parts), fill = pc, border = bg, horiz = TRUE, bty = "n", xpd = NA,
           text.col = ink2, cex = 0.82)
  })
}
sg_csv <- file.path(npn_dir, "site_signal", "signal.csv")
if (do("npn_signal") && have(sg_csv)) {
  G <- utils::read.csv(sg_csv, stringsAsFactors = FALSE)
  plot_png("npn_signal", "site-years", 2000, 900, function() {
    par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.6, 1))
    for (tg in c("curing", "cured")) {
      s <- G[G$target == tg & !is.na(G$vpd_jja), ]
      x <- s$vpd_jja - stats::ave(s$vpd_jja, s$year); y <- s$signal - stats::ave(s$signal, s$year)
      plot(x, y, pch = 16, cex = 0.8, col = grDevices::adjustcolor(terra_col, 0.5), bty = "n", cex.axis = 0.85,
           xlab = "June-August VPD, departure (hPa)", ylab = "shared signal, days later (+) or earlier (-)")
      abline(h = 0, v = 0, col = "#ecebe7")
      f <- stats::lm(y ~ x); abline(f, col = ink, lwd = 2)
      mtext(date_lab[[tg]], side = 3, line = 1.8, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext(sprintf("%d site-years; within each year; r = %.2f, slope %.1f days per hPa", nrow(s), stats::cor(x, y), stats::coef(f)[2]),
            side = 3, line = 0.5, adj = 0, cex = 0.72, col = ink2)
    }
  })
}

# ---- NPN study: record (58) ----
bn_csv <- file.path(rec_dir, "bins.csv"); rd_csv <- file.path(rec_dir, "readings.csv"); s2_csv <- file.path(rec_dir, "second.csv")
if (do("npn_record") && have(bn_csv, rd_csv, s2_csv)) {
  BN <- utils::read.csv(bn_csv, stringsAsFactors = FALSE); RD <- utils::read.csv(rd_csv, stringsAsFactors = FALSE)
  S2 <- utils::read.csv(s2_csv, stringsAsFactors = FALSE, colClasses = c(value = "character"))
  plot_png("npn_record", "plant-years", 2400, 950, function() {
    par(mfrow = c(1, 3), mar = c(4.4, 9.6, 3.8, 1.2), oma = c(0, 0, 0, 0))
    # 1. bins
    bl <- c("Less than 5%", "5-24%", "25-49%", "50-74%", "75-94%", "95% or more")
    lab <- c("< 5 %", "5-24 %", "25-49 %", "50-74 %", "75-94 %", ">= 95 %")
    gc <- c(NEON = ink, other = "#b8b7b2")
    plot(NA, xlim = c(0, 0.4), ylim = c(0.5, 6.5), yaxt = "n", xaxt = "n", xlab = "share of readings", ylab = "", bty = "n")
    axis(1, at = seq(0, 0.4, 0.1), labels = paste0(seq(0, 40, 10), " %"), cex.axis = 0.85)
    axis(2, at = 1:6, labels = lab, tick = FALSE, cex.axis = 0.85)
    abline(v = seq(0, 0.4, 0.1), col = "#ecebe7")
    for (g in names(gc)) { d <- BN[BN$group == g, ]; y <- match(d$bin, bl) + ifelse(g == "NEON", 0.17, -0.17)
      rect(0, y - 0.15, d$share, y + 0.15, col = gc[[g]], border = bg) }
    mtext("Percent-green class of each reading", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext(sprintf("%s readings; NEON (dark) and other sites (gray)", format(sum(BN$readings[BN$group == "all"]), big.mark = ",")),
          side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    # 2. readings per plant-year
    RD$lab <- sprintf("%s\n(%s plant-years)", c(NEON = "NEON", other = "other sites", managed = "gardens, yards,\ncampuses, parks")[RD$group],
                      format(RD$plant_years, big.mark = ","))
    plot(NA, xlim = c(0, 52), ylim = c(0.5, nrow(RD) + 0.5), yaxt = "n", xlab = "readings per plant-year", ylab = "", bty = "n", cex.axis = 0.85)
    abline(v = seq(0, 50, 10), col = "#ecebe7"); abline(v = c(6, 12), col = ink2, lty = 2)
    axis(2, at = rev(seq_len(nrow(RD))), labels = RD$lab, tick = FALSE, cex.axis = 0.8)
    for (i in seq_len(nrow(RD))) { y <- nrow(RD) + 1 - i
      segments(RD$n_obs_p10[i], y, RD$n_obs_p90[i], y, col = ink2, lwd = 1.5)
      rect(RD$n_obs_p25[i], y - 0.2, RD$n_obs_p75[i], y + 0.2, col = "#9ec5f4", border = bg)
      segments(RD$n_obs_median[i], y - 0.26, RD$n_obs_median[i], y + 0.26, lwd = 3, col = ink) }
    mtext("How densely plants are read", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Median, quartiles, 10-90 %; dashed: 6 and 12 readings", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    # 3. second green-ups by region
    r <- S2[S2$by == "region" & S2$value != "none" & S2$plant_years >= 30, ]
    r <- r[order(as.numeric(r$value), decreasing = TRUE), ]
    plot(NA, xlim = c(0, 0.6), ylim = c(0.5, nrow(r) + 0.5), yaxt = "n", xaxt = "n", xlab = "share of plant-years", ylab = "", bty = "n")
    axis(1, at = seq(0, 0.6, 0.2), labels = paste0(seq(0, 60, 20), " %"), cex.axis = 0.85)
    abline(v = seq(0, 0.6, 0.1), col = "#ecebe7")
    axis(2, at = seq_len(nrow(r)), labels = reg_lab(r$value, 24), tick = FALSE, cex.axis = 0.72, padj = 0.5)
    rect(0, seq_len(nrow(r)) - 0.32, pmin(r$share, 0.6), seq_len(nrow(r)) + 0.32, col = "#6da7ec", border = bg)
    text(pmin(r$share, 0.6), seq_len(nrow(r)), sprintf(" %d", r$plant_years), adj = 0, cex = 0.68, col = ink2)
    a <- S2[S2$by == "all", ]
    mtext("Second green-ups", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext(sprintf("Below 25 %%, then above 50 %%: %.0f %% overall (n = plant-years)", 100 * a$share),
          side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
  })
}

# ---- NPN study: typical dates vs site climate (58) ----
sp_csv <- file.path(rec_dir, "spatial.csv"); nm_rds <- file.path(npn_dir, "npn_metrics.rds")
if (do("npn_spatial") && have(sp_csv, st_csv, nm_rds)) {
  SPt <- utils::read.csv(sp_csv, stringsAsFactors = FALSE); SPt <- SPt[SPt$set == "relaxed", ]
  Nn <- readRDS(nm_rds)$normals
  X <- merge(TY, Nn[, c("site_id", "t_marapr", "aridity", "p_jas_share")], by = "site_id")
  cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
  unwrap <- function(x) { m <- cmed(x); d <- (x - m) %% 365.25; m + ifelse(d > 182.625, d - 365.25, d) }
  mc <- c(species = "#b8b7b2", region = "#1a7f37", climate = ink, both = "#7b4fb3")   # kept apart from the C3 / C4 colours
  ml <- c(species = "species only", region = "ecoregion", climate = "site climate", both = "region + climate")
  plot_png("npn_spatial", "site and species typical dates", 2400, 1000, function() {
    par(mfrow = c(1, 3), mar = c(4.6, 9.4, 3.8, 1.2), oma = c(2.2, 0, 0, 0))
    rows <- expand.grid(pathway = c("C4", "C3", "all"), target = rev(names(date_lab)), stringsAsFactors = FALSE)
    plot(NA, xlim = c(-0.4, 0.8), ylim = c(0.5, nrow(rows) + 0.5), yaxt = "n", xlab = "R2, leaving one site out", ylab = "", bty = "n", cex.axis = 0.85)
    abline(v = seq(-0.4, 0.8, 0.2), col = "#ecebe7"); abline(v = 0, col = ink2); abline(h = c(3.5, 6.5), col = "#e4e3df")
    pl <- c(all = "all grasses", C3 = "C3", C4 = "C4")
    axis(2, at = seq_len(nrow(rows)), labels = paste0(c(greenup = "Green-up", curing = "Curing", cured = "Cured")[rows$target], ", ", pl[rows$pathway]), tick = FALSE, cex.axis = 0.8)
    for (i in seq_len(nrow(rows))) for (k in seq_along(mc)) {
      v <- SPt$r2_loso[SPt$target == rows$target[i] & SPt$pathway == rows$pathway[i] & SPt$model == names(mc)[k]]
      if (length(v)) points(max(v, -0.4), i + (k - 2.5) * 0.17, pch = 16, cex = 1.3, col = mc[[k]])
    }
    mtext("How well typical dates are predicted at a new site", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Mixed models with a species effect; values below -0.4 drawn at -0.4", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    sc <- function(xv, tg, xlab, main, sub, logx = FALSE) {
      d <- X[!is.na(X[[tg]]) & !is.na(X[[xv]]), ]; d$y <- unwrap(d[[tg]])
      plot(d[[xv]], d$y, log = if (logx) "x" else "", pch = 16, cex = 0.5 + 0.3 * sqrt(pmin(d$years, 9)),
           col = grDevices::adjustcolor(pw_col[d$pathway], 0.6), yaxt = "n", xlab = xlab, ylab = "", bty = "n", cex.axis = 0.85)
      at <- as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1
      at <- c(at - 365, at, at + 365); lb <- rep(month.abb, 3); k <- at >= min(d$y) - 20 & at <= max(d$y) + 20
      axis(2, at = at[k], labels = lb[k], cex.axis = 0.8)
      for (p in c("C3", "C4")) { e <- d[d$pathway == p, ]; if (nrow(e) < 8) next
        xx <- if (logx) log(e[[xv]]) else e[[xv]]; f <- stats::lm(e$y ~ xx); o <- order(xx)
        lines(e[[xv]][o], stats::fitted(f)[o], col = pw_col[[p]], lwd = 2.5)
        r <- stats::cor(xx, e$y, method = "spearman")
        text(par("usr")[2], par("usr")[4] - (if (p == "C3") 0.04 else 0.10) * diff(par("usr")[3:4]),
             sprintf("%s r = %.2f (%d)", p, r, nrow(e)), adj = 1, cex = 0.78, col = pw_col[[p]]) }
      mtext(main, side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext(sub, side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    }
    sc("t_marapr", "greenup", "site March-April mean temperature (deg C)", "Green-up and spring temperature",
       "One point per site and species; size = years observed")
    sc("aridity", "curing", "site aridity, P / PET (log scale)", "Curing and site wetness",
       "Spearman r (series); lines: least squares", logx = TRUE)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = c(ml, "cool-season (C3)", "warm-season (C4)"), pch = 16, col = c(mc, pw_col[["C3"]], pw_col[["C4"]]),
           horiz = TRUE, bty = "n", text.col = ink2, cex = 0.95)
  })
}

# ---- NPN study: 2016-17 check (58; 53 main and sens/no1617) ----
pk_csv <- file.path(rec_dir, "peak_lt50.csv")
vs_main <- file.path(npn_dir, "variance", "summary.csv"); vs_sens <- file.path(npn_dir, "sens", "no1617", "variance", "summary.csv")
if (do("npn_sensitivity") && have(pk_csv, vs_main, vs_sens)) {
  PK <- utils::read.csv(pk_csv, stringsAsFactors = FALSE); PK <- PK[PK$level == "NEON_year", ]
  V0 <- utils::read.csv(vs_main, stringsAsFactors = FALSE); V1 <- utils::read.csv(vs_sens, stringsAsFactors = FALSE)
  plot_png("npn_sensitivity", "plant-years", 2200, 900, function() {
    par(mfrow = c(1, 2), mar = c(4.4, 5, 3.8, 1.2))
    odd <- PK$year %in% c(2016, 2017)
    b <- graphics::barplot(PK$share, names.arg = PK$year, col = ifelse(odd, aqua_col, "#b8b7b2"), border = bg, ylim = c(0, 0.35),
                           yaxt = "n", cex.names = 0.8, las = 2)
    axis(2, at = seq(0, 0.3, 0.1), labels = paste0(seq(0, 30, 10), " %"), cex.axis = 0.85)
    text(b, PK$share, sprintf("%d", PK$site_years), pos = 3, cex = 0.62, col = ink2)
    mtext("NEON plants that never pass 50 % green", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Share of plant-years with >= 6 readings, by year (numbers: sites)", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    tl <- c(greenup = "Green-up", curing = "Curing", cured = "Cured")
    sh <- function(V) vapply(names(tl), function(t) V$shared_share[V$target == t & V$group == "all"][1], 0)
    m <- rbind(sh(V0), sh(V1))
    b <- graphics::barplot(m, beside = TRUE, names.arg = tl, col = c(ink, "#9ec5f4"), border = bg, ylim = c(0, 0.5), yaxt = "n", cex.names = 0.9)
    axis(2, at = seq(0, 0.5, 0.1), labels = paste0(seq(0, 50, 10), " %"), cex.axis = 0.85)
    text(b, m, sprintf("%.0f %%", 100 * m), pos = 3, cex = 0.72, col = ink2)
    mtext("Year signal shared by a site's grasses", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Share of the year-to-year variation; all years (dark) and without 2016-17 (light)", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
  })
}


log_msg(sprintf("=== done | figures in %s | log: %s", fig_dir, log_file))
