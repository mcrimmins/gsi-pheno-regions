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
#   npn_species      typical season of each species with >= 3 site series: green-up, 50 % curing
#                    (with quartiles across sites) and cured (58c, report/tables/npn_species.csv)
#   npn_species_sites   the five species with most sites: typical green-up vs March-April temperature
#                       and 50 % curing vs aridity, one point per site (58, 51 normals)
#   npn_species_spread  per species: SD of the date between years vs between plants in the same
#                       site and year (58c, report/tables/npn_species_spread.csv)
#   npn_fuelmodels      Scott & Burgan fuel model (LANDFIRE FBFM40) at the sites: sites and plant-years per
#                       model, and typical 50 % curing by model and grass type (59b, 58)
#   npn_lter_c3c4       LTER cross-check: monthly share of plants green / with new leaves, C3 vs C4 at
#                       Sevilleta and C4 by vegetation zone at Jornada (55e)
#   npn_lter_agree      LTER vs NPN: same species (monthly shares), same place (Jornada black grama),
#                       and the 2016-17 check by year (55e)
#   npn_lter_weather    LTER year to year: skill of the best single weather driver, and the share of
#                       monsoon green-ups against Dec-May rain (55e with 55d)
#   npn_scale           NPN year to year at three scales (series, site, ecoregion): shared signal and
#                       leave-one-year-out skill of the best single lead driver (55f)
#   npn_power           monitoring design: reliability of a regional year mean and years needed to
#                       detect a weather driver, by sites and plants (58h)
#   npn_regions      EPA Level II ecoregions (CEC polygons cached by 59; sf) with the study sites (58)
#   npn_region_dates typical green-up and 50 % curing by ecoregion and grass type, as maps (58)
#   npn_region_variance year-to-year variance shares fitted per ecoregion, maps and bars (53b)
#   npn_region_c3c4  warm- minus cool-season typical green-up and curing by ecoregion, maps (58)
#   npn_inseason     in-season status: reversible dips by month and grass type, and weekly greenness
#                    change against VPD and rain anomalies (58i)
#   npn_lfmc         NPN study sites and herbaceous live fuel moisture sites (Globe-LFMC 2.0), distance to
#                    the nearest one, and paired readings where they meet (59d)
#   npn_repr            how the record covers U.S. grass fuels: LANDFIRE grass and grass-shrub fuel
#                       models with the sites, and shares by ecoregion (59c)
#
# Usage: Rscript npn_study/scripts/90_npn_figs.R [profile] [figure ...]   (default: all)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

all_figs <- c("npn_sites", "npn_dates", "npn_ladder", "npn_variance", "npn_signal", "npn_record", "npn_spatial", "npn_sensitivity", "npn_species",
              "npn_species_sites", "npn_species_spread", "npn_fuelmodels",
              "npn_lter_c3c4", "npn_lter_agree", "npn_lter_weather", "npn_scale", "npn_power", "npn_repr",
              "npn_regions", "npn_region_dates", "npn_region_variance", "npn_region_c3c4",
              "npn_inseason", "npn_lfmc")
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

sp_csv <- npn_path("report", "tables", "npn_species.csv")
if (do("npn_species") && have(sp_csv)) {
  P <- utils::read.csv(sp_csv, stringsAsFactors = FALSE)
  P <- P[P$pathway %in% c("C3", "C4") & P$curing_series >= 3, ]
  P <- P[order(P$pathway == "C3", -P$curing_doy), ]            # C3 block on top, each block early curing on top
  n <- nrow(P); P$y <- seq_len(n)
  plot_png("npn_species", paste0(min(P$first_year), "-", max(P$last_year)), 2200, 330 + 62 * n, function() {
    par(mar = c(4.2, 17, 3.8, 1))
    plot(NA, xlim = c(1, 366), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    abline(v = as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1, col = "#ecebe7")
    k <- sum(P$pathway == "C4"); if (k > 0 && k < n) abline(h = k + 0.5, col = ink2, lty = 3)
    mon_axis(1, 0.85)
    inv <- if ("fire_invader" %in% names(P)) P$fire_invader %in% "yes" else rep(FALSE, n)
    lab <- sprintf("%s. %s (%d sites)%s%s", substr(P$genus, 1, 1), P$species, P$sites, ifelse(P$life_cycle == "annual", " *", ""), ifelse(inv, " \u2020", ""))
    axis(2, at = P$y, labels = lab, tick = FALSE, cex.axis = 0.78, font.axis = 3)
    seg <- function(a, b, y, ...) {                                # season line, wrapping past Dec 31
      ok <- !is.na(a) & !is.na(b) & (b - a) %% 365.25 < 300     # medians out of order: no line
      a <- a[ok]; b <- b[ok]; y <- y[ok]; w <- b < a
      if (any(!w)) segments(a[!w], y[!w], b[!w], y[!w], ...)
      if (any(w)) { segments(a[w], y[w], 366, y[w], ...); segments(1, y[w], b[w], y[w], ...) }
    }
    col <- pw_col[P$pathway]
    seg(P$greenup_doy, P$curing_doy, P$y, col = grDevices::adjustcolor(col, 0.45), lwd = 5, lend = 1)
    seg(P$curing_doy, P$cured_doy, P$y, col = grDevices::adjustcolor(ink2, 0.5), lwd = 2, lend = 1)
    seg(P$curing_q25, P$curing_q75, P$y + 0.22, col = col, lwd = 1.2)
    points(P$greenup_doy, P$y, pch = 21, bg = bg, col = col, cex = 1.1, lwd = 1.5)
    points(P$curing_doy, P$y, pch = 21, bg = col, col = bg, cex = 1.4)
    points(P$cured_doy, P$y, pch = 15, col = ink2, cex = 0.9)
    mtext("Typical season by species", side = 3, line = 2.2, adj = 0, cex = 1.05, font = 2, col = ink)
    mtext("Species with 3 or more site series; median over sites. Thin line above: middle half of sites' curing dates. * annual; \u2020 invader that raises fire occurrence (Fusco et al. 2019)",
          side = 3, line = 0.8, adj = 0, cex = 0.72, col = ink2)
    legend("bottom", inset = c(0, -0.1 * 30 / max(n, 30)), xpd = NA, horiz = TRUE, bty = "n", cex = 0.78, text.col = ink2,
           legend = c("green-up", "50 % curing (C4)", "50 % curing (C3)", "cured"), pch = c(21, 21, 21, 15),
           pt.bg = c(bg, pw_col[["C4"]], pw_col[["C3"]], NA), col = c(ink2, bg, bg, ink2), pt.cex = c(1.1, 1.4, 1.4, 0.9))
  })
}

# ---- NPN study: the same species at different sites (58 series_typical, 51 normals) ----------
nm_rds2 <- file.path(npn_dir, "npn_metrics.rds")
if (do("npn_species_sites") && !is.null(TY) && have(nm_rds2)) {
  Nn <- readRDS(nm_rds2)$normals
  X <- merge(TY, Nn[, c("site_id", "t_marapr", "aridity")], by = "site_id")
  cm2 <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
  unw <- function(x) { m <- cm2(x); d <- (x - m) %% 365.25; m + ifelse(d > 182.625, d - 365.25, d) }
  ns_cur <- sort(table(paste(X$genus, X$species)[!is.na(X$curing)]), decreasing = TRUE)
  top <- names(ns_cur)[seq_len(min(5, length(ns_cur)))]
  plot_png("npn_species_sites", "site and species typical dates", 2700, 1250, function() {
    par(mfrow = c(2, length(top)), mar = c(4.2, 4.2, 3.6, 0.8), oma = c(0, 0, 2.4, 0))
    panel <- function(sp, xv, tg, xlab, logx) {
      d <- X[paste(X$genus, X$species) == sp & !is.na(X[[tg]]) & !is.na(X[[xv]]), ]
      d$y <- unw(d[[tg]]); colp <- pw_col[[d$pathway[1]]]
      at <- as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1
      at <- c(at - 365, at, at + 365); lb <- rep(month.abb, 3)
      yl <- range(d$y) + c(-15, 15); k <- at >= yl[1] & at <= yl[2]
      plot(d[[xv]], d$y, log = if (logx) "x" else "", ylim = yl, pch = 16, cex = 0.7 + 0.35 * sqrt(pmin(d$years, 9)),
           col = grDevices::adjustcolor(colp, 0.65), yaxt = "n", xlab = xlab, ylab = "", bty = "n", cex.axis = 0.85, cex.lab = 0.9)
      axis(2, at = at[k], labels = lb[k], cex.axis = 0.8)
      sub <- sprintf("%d sites", nrow(d))
      if (nrow(d) >= 5) {
        xx <- if (logx) log(d[[xv]]) else d[[xv]]; f <- stats::lm(d$y ~ xx); o <- order(xx)
        lines(d[[xv]][o], stats::fitted(f)[o], col = colp, lwd = 2)
        sub <- sprintf("%d sites; Spearman r = %.2f", nrow(d), stats::cor(xx, d$y, method = "spearman"))
      }
      mtext(sprintf("%s. %s", substr(sp, 1, 1), sub(".* ", "", sp)), side = 3, line = 1.6, adj = 0, cex = 0.9, font = 4, col = ink)
      mtext(paste0(c(greenup = "green-up", curing = "50 % curing")[[tg]], ": ", sub), side = 3, line = 0.4, adj = 0, cex = 0.7, col = ink2)
    }
    for (sp in top) panel(sp, "t_marapr", "greenup", "March-April mean temperature (deg C)", FALSE)
    for (sp in top) panel(sp, "aridity", "curing", "aridity, P / PET (log scale)", TRUE)
    mtext("The same species at different sites: typical green-up against spring temperature (top) and curing against site wetness (bottom)",
          side = 3, line = 0.6, adj = 0, outer = TRUE, cex = 0.95, font = 2, col = ink)
  })
}

# ---- NPN study: year-to-year vs plant-to-plant spread by species (58c) ------------------------
spr_csv <- npn_path("report", "tables", "npn_species_spread.csv")
if (do("npn_species_spread") && have(spr_csv)) {
  R <- utils::read.csv(spr_csv, stringsAsFactors = FALSE); R <- R[R$series_years >= 10, ]
  tgs <- intersect(c("greenup", "curing", "cured"), unique(R$target))
  ref <- R[R$target == "curing", ]; spp <- unique(c(ref$species_id[order(-ref$sd_year)], R$species_id))
  inv_ids <- if (file.exists(sp_csv)) { Q <- utils::read.csv(sp_csv, stringsAsFactors = FALSE); if ("fire_invader" %in% names(Q)) Q$species_id[Q$fire_invader %in% "yes"] else integer() } else integer()
  lab <- function(id) { r <- R[R$species_id == id, ][1, ]; sprintf("%s. %s%s", substr(r$genus, 1, 1), r$species, if (id %in% inv_ids) " \u2020" else "") }
  n <- length(spp)
  plot_png("npn_species_spread", "plant-years", 2500, 300 + 60 * n, function() {
    par(mfrow = c(1, length(tgs)), mar = c(4.4, 9.6, 3.6, 0.8), oma = c(2.2, 0, 0, 0))
    xm <- max(c(R$sd_year, R$sd_plant), na.rm = TRUE) * 1.05
    for (tg in tgs) {
      plot(NA, xlim = c(0, xm), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "standard deviation (days)", ylab = "", bty = "n", cex.axis = 0.85)
      abline(v = seq(0, xm, 10), col = "#ecebe7")
      axis(2, at = seq_len(n), labels = vapply(spp, lab, ""), tick = FALSE, cex.axis = 0.8, font.axis = 3)
      for (i in seq_len(n)) { r <- R[R$species_id == spp[i] & R$target == tg, ]; if (!nrow(r)) next
        colp <- pw_col[[r$pathway]]
        segments(min(r$sd_year, r$sd_plant, na.rm = TRUE), i, max(r$sd_year, r$sd_plant, na.rm = TRUE), i, col = "#d6d5d0", lwd = 2)
        points(r$sd_plant, i, pch = 21, bg = bg, col = ink2, cex = 1.2, lwd = 1.3)
        points(r$sd_year, i, pch = 21, bg = colp, col = bg, cex = 1.5)
        text(xm, i, r$series_years, adj = 1, cex = 0.65, col = ink2) }
      mtext(c(greenup = "Green-up", curing = "50 % curing", cured = "Cured")[[tg]], side = 3, line = 1.8, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext("Numbers at right: series-years", side = 3, line = 0.5, adj = 0, cex = 0.7, col = ink2)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", text.col = ink2, cex = 1, pch = 21, pt.cex = c(1.5, 1.5, 1.2),
           pt.bg = c(pw_col[["C4"]], pw_col[["C3"]], bg), col = c(bg, bg, ink2),
           legend = c("between years, warm-season (C4)", "between years, cool-season (C3)", "between plants in the same site and year"))
  })
}

# ---- NPN study: Scott & Burgan fuel models at the sites (59b, 58) -----------------------------
fm_csv <- npn_path("data", "npn_site_fuelmodels.csv")
if (do("npn_fuelmodels") && !is.null(TY) && have(fm_csv)) {
  FM <- utils::read.csv(fm_csv, stringsAsFactors = FALSE)
  FM$fm <- ifelse(!is.na(FM$fbfm_group) & FM$fbfm_group != "NB", FM$fbfm, FM$fbfm_3x3)    # point, else most common within 90 m
  FM$fm[is.na(FM$fm)] <- "nonburnable"
  X <- merge(TY, FM[, c("site_id", "fm")], by = "site_id")
  ord <- function(f) { g <- substr(f, 1, 2); k <- match(g, c("GR", "GS", "SH", "TU", "TL", "SB", "no")); f[order(k, f)] }
  lv <- ord(unique(X$fm)); X$fm <- factor(X$fm, lv)
  gcol <- c(GR = "#c9a227", GS = "#8a9a3b", SH = "#6b8e6b", TU = "#3f6e5a", TL = "#2f4f4f", SB = "#7a5c3e", no = "#b8b7b2")
  cmf <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
  st <- do.call(rbind, lapply(split(X, X$fm), function(d) if (nrow(d)) data.frame(fm = d$fm[1], sites = length(unique(d$site_id)),
          neon = length(unique(d$site_id[d$group == "NEON"])), py = sum(d$plant_years))))
  st$fm <- factor(as.character(st$fm), lv); st <- st[order(st$fm), ]
  n <- length(lv)
  plot_png("npn_fuelmodels", "LANDFIRE FBFM40 at the sites", 2400, 320 + 62 * n, function() {
    layout(matrix(1:2, 1), widths = c(1, 1.6)); par(mar = c(4.4, 7, 3.6, 1), oma = c(2, 0, 0, 0))
    y <- rev(seq_len(n)); names(y) <- lv
    plot(NA, xlim = c(0, max(st$sites) * 1.15), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "sites", ylab = "", bty = "n", cex.axis = 0.85)
    axis(2, at = y, labels = lv, tick = FALSE, las = 1, cex.axis = 0.85)
    cl <- gcol[substr(as.character(st$fm), 1, 2)]
    rect(0, y[as.character(st$fm)] - 0.35, st$sites, y[as.character(st$fm)] + 0.35, col = grDevices::adjustcolor(cl, 0.45), border = NA)
    rect(0, y[as.character(st$fm)] - 0.35, st$neon, y[as.character(st$fm)] + 0.35, col = cl, border = NA)
    text(st$sites, y[as.character(st$fm)], sprintf("  %s plant-yrs", format(st$py, big.mark = ",")), adj = 0, cex = 0.68, col = ink2, xpd = NA)
    mtext("Fuel model at the study sites", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Scott & Burgan (2005), LANDFIRE; dark = NEON sites", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    plot(NA, xlim = c(1, 366), ylim = c(0.5, n + 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    abline(v = as.numeric(as.Date(sprintf("2001-%02d-01", 1:12)) - as.Date("2001-01-01")) + 1, col = "#ecebe7")
    mon_axis(1, 0.85); axis(2, at = y, labels = lv, tick = FALSE, las = 1, cex.axis = 0.85)
    for (f in lv) for (p in c("C3", "C4")) {
      d <- X[X$fm == f & X$pathway == p & !is.na(X$curing), ]; if (!nrow(d)) next
      yy <- y[[f]] + if (p == "C3") 0.17 else -0.17
      points(d$curing, yy + stats::runif(nrow(d), -0.1, 0.1), pch = 16, cex = 0.8, col = grDevices::adjustcolor(pw_col[[p]], 0.55))
      if (nrow(d) >= 3) segments(cmf(d$curing), yy - 0.15, cmf(d$curing), yy + 0.15, lwd = 3, col = pw_col[[p]])
    }
    mtext("Typical 50 % curing by fuel model", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("Points: one site and species; bar: median where 3 or more", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", text.col = ink2, cex = 0.9, pch = 16, col = c(pw_col[["C3"]], pw_col[["C4"]], gcol[c("GR", "GS", "SH", "TU", "TL")]),
           legend = c("cool-season (C3)", "warm-season (C4)", "grass (GR)", "grass-shrub (GS)", "shrub (SH)", "timber-understory (TU)", "timber litter (TL)"))
  })
}

# ---- NPN study: LTER cross-checks (55e) --------------------------------------------------------
lc_dir <- file.path(npn_dir, "lter", "compare")
mon_ax <- function(cex = 0.8) axis(1, at = 1:12, labels = substr(month.abb, 1, 1), cex.axis = cex, tick = FALSE)
band <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), col = grDevices::adjustcolor(col, 0.18), border = NA)
if (do("npn_lter_c3c4") && have(file.path(lc_dir, "pathway_monthly.csv"), file.path(lc_dir, "c3c4_monthly.csv"))) {
  PM <- utils::read.csv(file.path(lc_dir, "pathway_monthly.csv"), stringsAsFactors = FALSE)
  W3 <- utils::read.csv(file.path(lc_dir, "c3c4_monthly.csv"), stringsAsFactors = FALSE)
  qs <- function(d, v) { a <- stats::aggregate(d[[v]], list(month = d$month), function(z) c(stats::quantile(z, c(0.25, 0.5, 0.75), na.rm = TRUE), n = sum(!is.na(z))))
    data.frame(month = a$month, lo = a$x[, 1], md = a$x[, 2], hi = a$x[, 3], n = a$x[, 4]) }
  yrs <- range(PM$year)
  plot_png("npn_lter_c3c4", paste(yrs, collapse = "-"), 2600, 900, function() {
    par(mfrow = c(1, 3), mar = c(3.4, 4.2, 3.8, 1), oma = c(2, 0, 0, 0))
    for (v in c("new", "green")) {
      plot(NA, xlim = c(1, 12), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "share of plants", bty = "n")
      axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), " %"), cex.axis = 0.8); mon_ax(); abline(h = 0.5, col = "#d6d5d0", lty = 2)
      for (pw in c("C4", "C3")) { d <- PM[PM$source == "SEV" & PM$site == "all" & PM$pathway == pw & !is.na(PM[[v]]), ]; if (!nrow(d)) next
        q <- qs(d, v); band(q$month, q$lo, q$hi, pw_col[[pw]]); lines(q$month, q$md, lwd = 2.5, col = pw_col[[pw]]); points(q$month, q$md, pch = 16, col = pw_col[[pw]])
        if (pw == "C3") text(q$month, 0.03, q$n, cex = 0.7, col = pw_col[["C3"]]) }
      mtext(sprintf("Sevilleta: %s", if (v == "new") "new green leaves" else "any green leaves"), side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
      mtext("median over years, band middle half; blue numbers at the bottom: years with >= 10 C3 plants", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    }
    plot(NA, xlim = c(1, 12), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "share of plants green", bty = "n")
    axis(2, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), " %"), cex.axis = 0.8); mon_ax(); abline(h = 0.5, col = "#d6d5d0", lty = 2)
    jz <- PM[PM$source == "JRN" & PM$pathway == "C4", ]; jz$zone <- substr(jz$site, 1, 1)
    zc <- c(G = "#c9a227", P = "#8a9a3b", T = "#6b8e6b", M = "#b05f3c", C = "#52514e"); zl <- c(G = "grassland", P = "playa", T = "tarbush", M = "mesquite", C = "creosote")
    for (z in names(zc)) { d <- jz[jz$zone == z, ]; if (!nrow(d)) next; q <- qs(d, "green"); lines(q$month, q$md, lwd = 2, col = zc[[z]]) }
    legend("center", legend = zl, col = zc, lwd = 2, bty = "n", cex = 0.8, text.col = ink2)
    mtext("Jornada: warm-season grasses by vegetation zone", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("share of plants not dormant; median over years and sites", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = c("cool-season (C3)", "warm-season (C4)"), col = pw_col[c("C3", "C4")], lwd = 2.5, pch = 16, horiz = TRUE, bty = "n", text.col = ink2, cex = 1)
  })
}
if (do("npn_lter_agree") && have(file.path(lc_dir, "climatology_species.csv"), file.path(lc_dir, "years_1617.csv"))) {
  CS <- utils::read.csv(file.path(lc_dir, "climatology_species.csv"), stringsAsFactors = FALSE)
  NMm <- utils::read.csv(file.path(npn_dir, "lter", "npn_monthly.csv"), stringsAsFactors = FALSE)
  STl <- utils::read.csv(file.path(npn_dir, "lter", "lter_sites.csv"), stringsAsFactors = FALSE)
  Y5 <- utils::read.csv(file.path(lc_dir, "years_1617.csv"), stringsAsFactors = FALSE)
  ov <- list.files(lc_dir, "^overlap_JRN_Bouteloua_eriopoda\\.csv$", full.names = TRUE)
  plot_png("npn_lter_agree", "LTER and NPN", 2600, 900, function() {
    par(mfrow = c(1, 3), mar = c(4.2, 4.4, 3.8, 1))
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "LTER: share of plants green", ylab = "NPN: share of readings >= 5 % green", bty = "n", cex.axis = 0.8)
    abline(0, 1, col = "#d6d5d0", lty = 2)
    xy <- stats::aggregate(cbind(lat, lon) ~ source, STl, mean); cols <- c(JRN = "#b05f3c", SEV = "#2f4f4f"); rr <- c()
    for (s in c("JRN", "SEV")) {
      c0 <- xy[xy$source == s, ]; nm <- NMm[111.2 * sqrt((NMm$lat - c0$lat)^2 + ((NMm$lon - c0$lon) * cos(c0$lat * pi / 180))^2) <= 400 & NMm$n >= 3, ]
      a <- stats::aggregate(green_mean ~ name + month, CS[CS$source == s, ], mean); b <- stats::aggregate(share_any ~ name + month, nm, mean)
      x <- merge(a, b); if (!nrow(x)) next
      points(x$green_mean, x$share_any, pch = 16, cex = 0.9, col = grDevices::adjustcolor(cols[[s]], 0.6)); rr[s] <- stats::cor(x$green_mean, x$share_any)
    }
    legend("topleft", legend = sprintf("%s (r = %.2f)", c(JRN = "Jornada", SEV = "Sevilleta")[names(rr)], rr), pch = 16, col = cols[names(rr)], bty = "n", cex = 0.8, text.col = ink2)
    mtext("Same species, by month", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("species x month means; NPN sites within 400 km", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    if (length(ov)) {
      O <- utils::read.csv(ov[1]); O$t <- O$year + (O$month - 0.5) / 12; O <- O[order(O$t), ]
      plot(NA, xlim = range(O$t), ylim = c(0, 1), xlab = "", ylab = "share green", bty = "n", cex.axis = 0.8)
      brk <- function(x, y) { g <- c(FALSE, diff(x) > 0.2); list(x = unlist(Map(function(a, b) if (b) c(NA, a) else a, x, g)), y = unlist(Map(function(a, b) if (b) c(NA, a) else a, y, g))) }
      l1 <- brk(O$t, O$green); l2 <- brk(O$t, O$share_any)                 # no lines across months without data
      lines(l1$x, l1$y, col = "#b05f3c", lwd = 2); lines(l2$x, l2$y, col = ink, lwd = 2); points(O$t, O$share_any, pch = 16, cex = 0.5, col = ink)
      legend("bottomleft", legend = c("Jornada LTER (share of plants not dormant)", "NEON Jornada (share of readings >= 5 % green)"), col = c("#b05f3c", ink), lwd = 2, bty = "n", cex = 0.72, text.col = ink2)
      mtext("Black grama at Jornada, same years", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
      mtext(sprintf("monthly; r = %.2f, anomalies r = %.2f", stats::cor(O$green, O$share_any), stats::cor(O$green - stats::ave(O$green, O$month), O$share_any - stats::ave(O$share_any, O$month))),
            side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    } else { plot.new(); text(0.5, 0.5, "no overlap table", col = ink2) }
    Y5 <- Y5[Y5$year >= 2010, ]; sc <- c(JRN = "#b05f3c", SEV = "#2f4f4f", `NPN NEON` = ink)
    plot(NA, xlim = range(Y5$year), ylim = c(0, max(0.3, max(Y5$never_50, na.rm = TRUE))), xlab = "", ylab = "share of series", bty = "n", cex.axis = 0.8)
    rect(2015.5, 0, 2017.5, 1, col = "#f1efe9", border = NA)
    for (s in names(sc)) { d <- Y5[Y5$source == s, ]; d <- d[order(d$year), ]; lines(d$year, d$never_50, col = sc[[s]], lwd = 2); points(d$year, d$never_50, pch = 16, cex = 0.7, col = sc[[s]]) }
    legend("topright", legend = c("Jornada LTER", "Sevilleta LTER", "NPN, NEON sites"), col = sc, lwd = 2, bty = "n", cex = 0.75, text.col = ink2)
    mtext("Series that never pass 50 % green", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("LTER: half the plants green; NPN: half the readings >= 50 %; 2016-17 shaded", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
  })
}
if (do("npn_lter_weather") && have(file.path(lc_dir, "weather_skill.csv"), file.path(lc_dir, "weather_monsoon_years.csv"))) {
  K <- utils::read.csv(file.path(lc_dir, "weather_skill.csv"), stringsAsFactors = FALSE)
  YR <- utils::read.csv(file.path(lc_dir, "weather_monsoon_years.csv"), stringsAsFactors = FALSE)
  plot_png("npn_lter_weather", "LTER site-years", 2300, 950, function() {
    layout(matrix(1:2, 1), widths = c(1.2, 1)); par(mar = c(4.4, 13, 3.8, 1))
    K <- K[K$n >= 100, ]; K <- K[order(K$source, K$pathway, K$var, K$date), ]; n <- nrow(K)   # rows on >= 100 site-years
    dl <- c(up1 = "green-up", down1 = "first fall", down_last = "fall"); vl <- c(green = "green", new = "new leaves")
    plot(NA, xlim = c(min(-0.4, min(K$skill)), max(0.4, max(K$skill))), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "skill vs the series mean (0 = no better)", ylab = "", bty = "n", cex.axis = 0.8)
    abline(v = 0, col = ink2); abline(v = seq(-0.4, 0.4, 0.2), col = "#ecebe7")
    axis(2, at = seq_len(n), labels = sprintf("%s %s %s, %s", K$source, K$pathway, vl[K$var], dl[K$date]), las = 1, tick = FALSE, cex.axis = 0.72)
    points(K$skill, seq_len(n), pch = 16, cex = 1.3, col = pw_col[K$pathway])
    mtext("Best single weather driver, leave one year out", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("driver chosen within the training years, as the ladder's L1; rows with >= 100 site-years", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    par(mar = c(4.4, 4.4, 3.8, 1)); Y <- YR[YR$var == "green", ]; sc <- c(JRN = "#b05f3c", SEV = "#2f4f4f")
    plot(Y$p_dm_a, Y$monsoon, pch = 16, col = sc[Y$source], xlab = "Dec-May rain, departure from the site mean (mm)", ylab = "share of series with a monsoon green-up", bty = "n", cex.axis = 0.8, ylim = c(0, 1))
    for (s in names(sc)) { d <- Y[Y$source == s, ]; if (nrow(d) > 5) { f <- stats::lm(monsoon ~ p_dm_a, d); o <- order(d$p_dm_a); lines(d$p_dm_a[o], stats::fitted(f)[o], col = sc[[s]], lwd = 2) } }
    legend("topright", legend = c("Jornada", "Sevilleta"), pch = 16, col = sc, bty = "n", cex = 0.8, text.col = ink2)
    mtext("Spring or monsoon green-up?", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("one point per year; monsoon green-up = first above half the plants green after June 1", side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
  })
}

sk_csv <- file.path(npn_dir, "scale", "skill.csv")
if (do("npn_scale") && have(sk_csv)) {
  K <- utils::read.csv(sk_csv, stringsAsFactors = FALSE); K <- K[K$drivers == "lead", ]
  tl <- c(greenup = "Green-up", curing = "50 % curing", cured = "Cured"); sl <- c(series = "one series", site = "site mean", region = "ecoregion mean")
  sc_col <- c(series = "#b8b7b2", site = "#52514e", region = ink)
  plot_png("npn_scale", "series-years", 2400, 900, function() {
    par(mfrow = c(1, 2), mar = c(4.4, 9.5, 3.8, 1), oma = c(2, 0, 0, 0))
    rows <- expand.grid(pathway = c("C4", "C3", "all"), target = rev(names(tl)), stringsAsFactors = FALSE); rows$y <- seq_len(nrow(rows))
    lab_r <- paste0(tl[rows$target], ", ", c(all = "all", C3 = "C3", C4 = "C4")[rows$pathway])
    for (pan in c("signal_ratio", "skill")) {
      xl <- if (pan == "skill") c(-0.5, 0.3) else c(0, max(5, max(K$signal_ratio, na.rm = TRUE)))
      plot(NA, xlim = xl, ylim = c(0.5, nrow(rows) + 0.5), yaxt = "n", xlab = if (pan == "skill") "skill vs no departure (0 = no better)" else "variance of the means / noise alone",
           ylab = "", bty = "n", cex.axis = 0.85)
      abline(v = if (pan == "skill") 0 else 1, col = ink2); abline(h = c(3.5, 6.5), col = "#e4e3df")
      axis(2, at = rows$y, labels = lab_r, las = 1, tick = FALSE, cex.axis = 0.8)
      for (i in seq_len(nrow(rows))) for (k in seq_along(sc_col)) {
        v <- K[[pan]][K$target == rows$target[i] & K$pathway == rows$pathway[i] & K$scale == names(sc_col)[k]]
        if (length(v) && !is.na(v)) points(max(min(v, xl[2]), xl[1]), rows$y[i] + (k - 2) * 0.22, pch = 16, cex = 1.3, col = sc_col[[k]]) }
      mtext(if (pan == "skill") "Best single lead weather driver, leave one year out" else "Is there a shared year signal?", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
      mtext(if (pan == "skill") "drivers ending >= 14 days before the typical date, chosen within the training years" else "> 1: the averaged departures vary more than plant and series noise alone would make them",
            side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = sl, col = sc_col, pch = 16, horiz = TRUE, bty = "n", text.col = ink2, cex = 1)
  })
}

pw_csv <- npn_path("report", "tables", "npn_power.csv")
if (do("npn_power") && have(pw_csv)) {
  P <- utils::read.csv(pw_csv, stringsAsFactors = FALSE)
  plc <- c(`1` = "#d6d5d0", `3` = "#9c9b96", `5` = "#52514e", `20` = ink)
  plot_png("npn_power", "variance components", 2500, 950, function() {
    par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.8, 1), oma = c(2, 0, 0, 0))
    for (tg in c("curing", "greenup")) {
      g <- P[P$target == tg & P$n_species == 2 & P$f == 0.5, ]
      plot(NA, xlim = c(1, 20), ylim = c(0, 1), log = "x", xlab = "sites in the region", ylab = "reliability of the regional year mean", bty = "n", cex.axis = 0.85, xaxt = "n")
      axis(1, at = c(1, 2, 3, 5, 10, 20), cex.axis = 0.85)
      for (k in names(plc)) { d <- g[g$n_plants == as.numeric(k), ]; d <- d[order(d$n_sites), ]; lines(d$n_sites, d$reliability, col = plc[[k]], lwd = 2.5); points(d$n_sites, d$reliability, pch = 16, col = plc[[k]]) }
      legend("bottomright", legend = paste(names(plc), ifelse(names(plc) == "1", "plant", "plants"), "per species"), col = plc, lwd = 2.5, bty = "n", cex = 0.8, text.col = ink2)
      mtext(sprintf("%s: how well is a year measured?", c(curing = "50 % curing", greenup = "Green-up")[[tg]]), side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
      mtext(sprintf("2 species per site; regional year signal SD %.0f days", g$signal_sd_days[1]), side = 3, line = 0.6, adj = 0, cex = 0.68, col = ink2)
    }
    g <- P[P$n_species == 2 & P$n_plants == 5, ]
    plot(NA, xlim = c(1, 20), ylim = c(5, 150), log = "xy", xlab = "sites in the region", ylab = "years of monitoring for 80 % power", bty = "n", cex.axis = 0.85, xaxt = "n")
    axis(1, at = c(1, 2, 3, 5, 10, 20), cex.axis = 0.85); abline(h = 11, col = "#d6d5d0", lty = 2); text(20, 11, "this record: 11 years", adj = c(1, -0.4), cex = 0.7, col = ink2)
    tc <- c(curing = aqua_col, greenup = terra_col)
    for (tg in names(tc)) for (f in c(0.5, 0.25)) { d <- g[g$target == tg & g$f == f, ]; d <- d[order(d$n_sites), ]
      lines(d$n_sites, pmin(d$years_80, 150), col = tc[[tg]], lwd = 2.5, lty = if (f == 0.5) 1 else 2); points(d$n_sites, pmin(d$years_80, 150), pch = 16, col = tc[[tg]]) }
    legend("topright", legend = c("curing, driver explains half", "curing, a quarter", "green-up, half", "green-up, a quarter"), col = tc[c("curing", "curing", "greenup", "greenup")],
           lty = c(1, 2, 1, 2), lwd = 2.5, bty = "n", cex = 0.78, text.col = ink2)
    mtext("How many years to detect a weather driver?", side = 3, line = 2, adj = 0, cex = 0.92, font = 2, col = ink)
    mtext("2 species x 5 plants per site; the driver acts on the regional year signal", side = 3, line = 0.6, adj = 0, cex = 0.68, col = ink2)
  })
}

# ---- EPA Level II ecoregion polygons for maps: clipped to the contiguous U.S., simplified, cached ----
eco_cache <- NULL
eco_l2 <- function() {
  if (!is.null(eco_cache)) return(eco_cache)
  if (!requireNamespace("sf", quietly = TRUE)) { log_msg("skip ecoregion layer: needs sf"); return(NULL) }
  edir <- file.path(npn_dir, "ecoregions"); rds <- file.path(edir, "l2_conus_plot.rds")
  if (file.exists(rds)) { eco_cache <<- readRDS(rds); return(eco_cache) }
  shp <- list.files(edir, "\\.shp$", recursive = TRUE, full.names = TRUE)
  shp <- shp[grepl("l2|level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
  if (!length(shp) || is.na(shp)) { log_msg("skip ecoregion layer: no shapefile in ", edir, " (run 59)"); return(NULL) }
  log_msg("building the ecoregion map layer from ", basename(shp), " (once; cached as l2_conus_plot.rds)")
  sf::sf_use_s2(FALSE)
  E <- sf::st_read(shp, quiet = TRUE); gc <- attr(E, "sf_column"); names(E)[names(E) != gc] <- toupper(names(E)[names(E) != gc])
  E <- E[!as.character(E$NA_L2CODE) %in% c("0.0", "0"), ]
  E <- sf::st_simplify(sf::st_make_valid(sf::st_transform(E, 5070)), dTolerance = 1500, preserveTopology = TRUE)
  us <- sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
  us <- sf::st_union(sf::st_make_valid(sf::st_transform(us, 5070)))
  E <- suppressWarnings(sf::st_intersection(sf::st_make_valid(E), us))
  E <- E[!sf::st_is_empty(E), ]
  cd <- as.character(E$NA_L2CODE); g <- sf::st_geometry(E)
  G <- do.call(c, lapply(split(seq_along(cd), cd), function(i) sf::st_union(g[i])))
  out <- sf::st_sf(code = names(split(seq_along(cd), cd)), geometry = G)
  out <- sf::st_transform(out, 4326)
  # label point: on the largest piece of each region
  pcs <- suppressWarnings(sf::st_cast(sf::st_cast(out, "MULTIPOLYGON"), "POLYGON"))
  pcs$a <- as.numeric(sf::st_area(sf::st_transform(pcs, 5070)))
  big <- do.call(rbind, lapply(split(pcs, pcs$code), function(d) d[which.max(d$a), ]))
  xy <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_transform(big, 5070))))
  xy <- sf::st_coordinates(sf::st_transform(sf::st_sfc(lapply(seq_len(nrow(xy)), function(i) sf::st_point(xy[i, 1:2])), crs = 5070), 4326))
  out$lab_x <- xy[match(out$code, big$code), 1]; out$lab_y <- xy[match(out$code, big$code), 2]
  out$area_km2 <- as.numeric(sf::st_area(sf::st_transform(out, 5070))) / 1e6
  saveRDS(out, rds); eco_cache <<- out; out
}
l1_col <- c(`5` = "#7f9cb5", `6` = "#5f8a6a", `7` = "#4f7f8f", `8` = "#8aa66a", `9` = "#c9a227", `10` = "#c98a4b",
            `11` = "#b5776b", `12` = "#a3866b", `13` = "#8f7f9f", `15` = "#6b9f9a")
eco_draw <- function(E, fill = NA, border = "#ffffff", lwd = 0.5) plot(sf::st_geometry(E), col = fill, border = border, lwd = lwd, add = TRUE)
site_tab <- function() do.call(rbind, lapply(split(TY, TY$site_id), function(d) data.frame(
  site_id = d$site_id[1], lon = d$lon[1], lat = d$lat[1], region = d$region[1], n = sum(d$plant_years), neon = any(d$group == "NEON"),
  type = if (all(d$pathway == "C3")) "C3" else if (all(d$pathway == "C4")) "C4" else "both")))
type_col <- c(C3 = terra_col, C4 = aqua_col, both = "#7b4fb3")

# ---- NPN study: Level II ecoregions with the sites (59 polygons, 58) ------------------------------
if (do("npn_regions") && !is.null(TY) && !is.null(E <- eco_l2())) {
  st <- site_tab()
  rs <- do.call(rbind, lapply(split(st, st$region), function(d) data.frame(region = d$region[1], sites = nrow(d), neon = sum(d$neon), py = sum(d$n))))
  rs <- rs[rs$region != "none", ]; rs <- rs[order(as.numeric(rs$region)), ]
  l1 <- sub("\\..*$", "", E$code); fc <- ifelse(is.na(l1_col[l1]), "#b8b7b2", l1_col[l1])
  plot_png("npn_regions", "EPA Level II ecoregions", 2900, 1300, function() {
    layout(matrix(1:2, 1), widths = c(1.75, 1)); par(mar = c(3.4, 0.6, 3.8, 0.6), oma = c(1.6, 0, 0, 0))
    st_map(sprintf("%d Level II ecoregions in the contiguous U.S.; shaded strongly where the record has sites; %d sites, sized by plant-years", nrow(E), nrow(st)),
           "EPA Level II ecoregions and the study sites")
    eco_draw(E, mapply(grDevices::adjustcolor, fc, ifelse(E$code %in% rs$region, 0.55, 0.16), USE.NAMES = FALSE), "#ffffff", 0.8)
    s <- maps::map("state", plot = FALSE); lines(s$x, s$y, col = grDevices::adjustcolor("#ffffff", 0.6), lwd = 0.4)
    lab <- E$area_km2 > 2e4 | E$code %in% rs$region
    text(E$lab_x[lab], E$lab_y[lab], E$code[lab], cex = 0.62, font = 2, col = grDevices::adjustcolor(ink, 0.75))
    filled <- st$neon
    points(st$lon, st$lat, pch = ifelse(filled, 21, 1), bg = grDevices::adjustcolor(type_col[st$type], 0.85),
           col = ifelse(filled, "#ffffff", type_col[st$type]), cex = 0.5 + 1.5 * sqrt(st$n / max(st$n)), lwd = ifelse(filled, 0.6, 1.5))
    par(mar = c(3.4, 0.5, 3.8, 1))
    plot(NA, xlim = c(0, 1), ylim = c(nrow(rs) + 0.5, 0.3), axes = FALSE, xlab = "", ylab = "")
    mtext("Ecoregions with study sites", side = 3, line = 1.9, adj = 0, cex = 1.0, font = 2, col = ink)
    mtext("sites (NEON) and plant-years", side = 3, line = 0.6, adj = 0, cex = 0.8, col = ink2)
    y <- seq_len(nrow(rs)); l1r <- sub("\\..*$", "", rs$region)
    rect(0, y - 0.4, 0.035, y + 0.4, col = grDevices::adjustcolor(ifelse(is.na(l1_col[l1r]), "#b8b7b2", l1_col[l1r]), 0.7), border = NA)
    nm <- if (length(rnm)) unname(rnm[rs$region]) else rep("", nrow(rs)); nm[is.na(nm)] <- ""
    text(0.05, y, sprintf("%s  %s", rs$region, substr(nm, 1, 50)), adj = 0, cex = 0.68, col = ink)
    text(1, y, sprintf("%d (%d)   %s", rs$sites, rs$neon, format(rs$py, big.mark = ",")), adj = 1, cex = 0.72, col = ink2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", text.col = ink2, cex = 0.85, pch = c(21, 21, 21, 1), pt.bg = c(type_col, NA), col = c(rep("#ffffff", 3), ink2), pt.cex = 1.4,
           legend = c("cool-season (C3) only", "warm-season (C4) only", "both", "open: other NPN site; filled: NEON"))
  })
}

# ---- NPN study: typical dates by ecoregion, as maps (58) --------------------------------------------
if (do("npn_region_dates") && !is.null(TY) && !is.null(E <- eco_l2())) {
  cmn <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
  pan <- list(c("C3", "greenup"), c("C3", "curing"), c("C4", "greenup"), c("C4", "curing"))
  st <- site_tab()
  dl <- function(d) format(as.Date("2001-01-01") + round(d %% 365.25) - 1, "%b %-d")
  pal <- grDevices::hcl.colors(64, "Viridis")
  plot_png("npn_region_dates", "site and species typical dates", 2900, 1900, function() {
    par(mfrow = c(2, 2), mar = c(4.2, 0.6, 3.6, 0.6), oma = c(0, 0, 0, 0))
    for (pp in pan) {
      d <- TY[TY$pathway == pp[1] & !is.na(TY[[pp[2]]]) & TY$region != "none", ]
      M <- do.call(rbind, lapply(split(d, d$region), function(z) data.frame(region = z$region[1], n = nrow(z), x = cmn(z[[pp[2]]]))))
      M <- M[M$n >= 3, ]
      ttl <- sprintf("%s grasses: %s", c(C3 = "Cool-season (C3)", C4 = "Warm-season (C4)")[[pp[1]]], c(greenup = "green-up (50 % green)", curing = "50 % curing")[[pp[2]]])
      st_map(sprintf("region mean of site-and-species typical dates, where 3 or more (%d regions); gray: fewer", nrow(M)), ttl)
      eco_draw(E, "#ecebe7", "#ffffff", 0.6)
      if (nrow(M)) {
        ctr <- cmn(M$x); off <- ((M$x - ctr + 182.625) %% 365.25) - 182.625
        lo <- floor(min(off) / 10) * 10; hi <- ceiling(max(off) / 10) * 10; if (hi - lo < 20) { lo <- lo - 10; hi <- hi + 10 }
        k <- pmax(1, pmin(64, 1 + round(63 * (off - lo) / (hi - lo))))
        i <- match(M$region, E$code); eco_draw(E[i[!is.na(i)], ], pal[k[!is.na(i)]], "#ffffff", 0.6)
        ii <- i[!is.na(i)]; text(E$lab_x[ii], E$lab_y[ii], dl(M$x[!is.na(i)]), cex = 0.6, font = 2, col = ifelse(k[!is.na(i)] > 40, ink, "#ffffff"))
        # color bar
        u <- par("usr"); x0 <- u[1] + 0.25 * diff(u[1:2]); x1 <- u[1] + 0.75 * diff(u[1:2]); yb <- u[3] - 0.05 * diff(u[3:4])
        xs <- seq(x0, x1, length.out = 65)
        rect(xs[-65], yb - 0.022 * diff(u[3:4]), xs[-1], yb, col = pal, border = NA, xpd = NA)
        tk <- pretty(c(lo, hi), 5); tk <- tk[tk >= lo & tk <= hi]
        text(x0 + (tk - lo) / (hi - lo) * (x1 - x0), yb - 0.05 * diff(u[3:4]), dl(ctr + tk), cex = 0.7, col = ink2, xpd = NA)
      }
      have_t <- st$type %in% c(pp[1], "both")
      points(st$lon[have_t], st$lat[have_t], pch = 16, cex = 0.35, col = grDevices::adjustcolor(ink, 0.6))
    }
  })
}

# ---- NPN study: year-to-year variance shares by ecoregion (53b) ------------------------------------
vr_csv <- file.path(npn_dir, "variance", "by_region.csv")
if (do("npn_region_variance") && have(vr_csv) && !is.null(E <- eco_l2())) {
  V <- utils::read.csv(vr_csv, stringsAsFactors = FALSE, colClasses = c(region = "character"))
  pal <- grDevices::hcl.colors(64, "YlGnBu", rev = TRUE)
  cc <- c(share_region_year = "#1f5f8b", share_site_year = "#6fa8c9", share_species_year = "#c9a227", share_plant = "#d6d5d0")
  cl <- c("region-wide year signal", "one site's year signal", "one species at a site", "plant to plant")
  tl <- c(greenup = "Green-up", curing = "50 % curing")
  st <- if (!is.null(TY)) site_tab() else NULL
  plot_png("npn_region_variance", "plant-years", 2900, 2000, function() {
    layout(matrix(1:4, 2, byrow = TRUE), heights = c(1.15, 1)); par(oma = c(2.4, 0, 0, 0))
    for (tg in names(tl)) {
      d <- V[V$target == tg, ]
      par(mar = c(3.4, 0.6, 3.8, 0.6))
      st_map(sprintf("share of year-to-year variance shared across the ecoregion; fitted where %d+ plant-years and 4+ years (%d regions); gray: not fitted", 60, nrow(d)),
             sprintf("%s: how much of a year's departure is region-wide?", tl[[tg]]))
      eco_draw(E, "#ecebe7", "#ffffff", 0.6)
      i <- match(d$region, E$code); ok <- !is.na(i)
      k <- pmax(1, pmin(64, 1 + round(63 * pmin(d$share_region_year, 0.6) / 0.6)))
      eco_draw(E[i[ok], ], pal[k[ok]], "#ffffff", 0.6)
      text(E$lab_x[i[ok]], E$lab_y[i[ok]], sprintf("%.0f %%\nSD %.0f d", 100 * d$share_region_year[ok], d$sd_year_to_year[ok]), cex = 0.58, font = 2,
           col = ifelse(k[ok] > 36, "#ffffff", ink))
      if (!is.null(st)) points(st$lon, st$lat, pch = 16, cex = 0.3, col = grDevices::adjustcolor(ink, 0.5))
      u <- par("usr"); x0 <- u[1] + 0.3 * diff(u[1:2]); x1 <- u[1] + 0.7 * diff(u[1:2]); yb <- u[3] - 0.04 * diff(u[3:4]); xs <- seq(x0, x1, length.out = 65)
      rect(xs[-65], yb - 0.022 * diff(u[3:4]), xs[-1], yb, col = pal, border = NA, xpd = NA)
      text(x0 + c(0, 0.25, 0.5, 0.75, 1) * (x1 - x0), yb - 0.05 * diff(u[3:4]), c("0 %", "15 %", "30 %", "45 %", "60 %+"), cex = 0.7, col = ink2, xpd = NA)
    }
    for (tg in names(tl)) {
      d <- V[V$target == tg, ]; d <- d[order(-d$share_region_year), ]; n <- nrow(d)
      par(mar = c(4.2, 13, 2.6, 1.5))
      plot(NA, xlim = c(0, 1), ylim = c(n + 0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "share of year-to-year variance", ylab = "", bty = "n")
      axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), " %"), cex.axis = 0.8)
      axis(2, at = seq_len(n), labels = sprintf("%s (%d sites, %d yrs)", reg_lab(d$region, 34), d$sites, d$years), tick = FALSE, las = 1, cex.axis = 0.62)
      left <- rep(0, n)
      for (cn in names(cc)) { rect(left, seq_len(n) - 0.38, left + d[[cn]], seq_len(n) + 0.38, col = cc[[cn]], border = "#ffffff", lwd = 0.5); left <- left + d[[cn]] }
      mtext(sprintf("%s: where the year-to-year variation sits, by ecoregion", tl[[tg]]), side = 3, line = 0.8, adj = 0, cex = 0.85, font = 2, col = ink)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = cl, fill = cc, border = NA, horiz = TRUE, bty = "n", text.col = ink2, cex = 0.95)
  })
}

# ---- NPN study: warm- minus cool-season timing by ecoregion (58) ------------------------------------
if (do("npn_region_c3c4") && !is.null(TY) && !is.null(E <- eco_l2())) {
  cmn <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
  cdf <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
  pal <- grDevices::hcl.colors(64, "Blue-Red 3")
  st <- site_tab(); stb <- st[st$type == "both", ]
  plot_png("npn_region_c3c4", "site and species typical dates", 2900, 1050, function() {
    par(mfrow = c(1, 2), mar = c(4.2, 0.6, 3.8, 0.6))
    for (tg in c("greenup", "curing")) {
      M <- do.call(rbind, lapply(split(TY, TY$region), function(z) {
        a <- z[[tg]][z$pathway == "C3"]; b <- z[[tg]][z$pathway == "C4"]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
        data.frame(region = z$region[1], n3 = length(a), n4 = length(b), d = if (length(a) >= 3 && length(b) >= 3) cdf(cmn(b), cmn(a)) else NA) }))
      M <- M[!is.na(M$d) & M$region != "none", ]
      st_map(sprintf("warm-season (C4) minus cool-season (C3) typical date, where each has 3+ site series (%d regions); gray: fewer", nrow(M)),
             sprintf("%s: warm- minus cool-season", c(greenup = "Green-up (50 % green)", curing = "50 % curing")[[tg]]))
      eco_draw(E, "#ecebe7", "#ffffff", 0.6)
      if (nrow(M)) {
        lim <- max(30, ceiling(max(abs(M$d)) / 30) * 30)
        k <- pmax(1, pmin(64, 1 + round(63 * (M$d + lim) / (2 * lim))))
        i <- match(M$region, E$code); ok <- !is.na(i)
        eco_draw(E[i[ok], ], pal[k[ok]], "#ffffff", 0.6)
        text(E$lab_x[i[ok]], E$lab_y[i[ok]], sprintf("%+.0f d", M$d[ok]), cex = 0.65, font = 2, col = ifelse(abs(k[ok] - 32.5) > 18, "#ffffff", ink))
        u <- par("usr"); x0 <- u[1] + 0.25 * diff(u[1:2]); x1 <- u[1] + 0.75 * diff(u[1:2]); yb <- u[3] - 0.04 * diff(u[3:4]); xs <- seq(x0, x1, length.out = 65)
        rect(xs[-65], yb - 0.025 * diff(u[3:4]), xs[-1], yb, col = pal, border = NA, xpd = NA)
        tk <- pretty(c(-lim, lim), 6); tk <- tk[abs(tk) <= lim]
        text(x0 + (tk + lim) / (2 * lim) * (x1 - x0), yb - 0.055 * diff(u[3:4]), sprintf("%+d d", tk), cex = 0.7, col = ink2, xpd = NA)
        text((x0 + x1) / 2, yb - 0.1 * diff(u[3:4]), "C4 earlier  <->  C4 later", cex = 0.72, col = ink2, xpd = NA)
      }
      points(stb$lon, stb$lat, pch = 21, bg = type_col[["both"]], col = "#ffffff", cex = 0.8, lwd = 0.5)
    }
  })
}

# ---- NPN study: in-season status (58i) -----------------------------------------------------------------
dp_csv <- file.path(npn_dir, "inseason", "dips.csv"); ft_csv <- file.path(npn_dir, "inseason", "pairs_fit.csv")
if (do("npn_inseason") && have(dp_csv, ft_csv)) {
  DP <- utils::read.csv(dp_csv, stringsAsFactors = FALSE); FT <- utils::read.csv(ft_csv, stringsAsFactors = FALSE)
  plot_png("npn_inseason", "plant-years", 2700, 950, function() {
    layout(matrix(1:3, 1), widths = c(1.15, 1, 1.1)); par(mar = c(4.4, 4.6, 3.8, 1), oma = c(1.6, 0, 0, 0))
    # A: dips by month of the low, C3 / C4; dark = dips that cross 50 %
    cnt <- function(pw, x) tabulate(DP$month[DP$pathway == pw & x], 12)
    a3 <- cnt("C3", TRUE); a4 <- cnt("C4", TRUE); c3 <- cnt("C3", DP$crosses50); c4 <- cnt("C4", DP$crosses50)
    ym <- max(a3, a4) * 1.12
    plot(NA, xlim = c(0.5, 12.5), ylim = c(0, ym), xaxt = "n", xlab = "", ylab = "dips (plant-years)", bty = "n", cex.axis = 0.85)
    axis(1, at = 1:12, labels = substr(month.abb, 1, 1), cex.axis = 0.85, tick = FALSE)
    rect(1:12 - 0.4, 0, 1:12, a3, col = grDevices::adjustcolor(terra_col, 0.35), border = NA); rect(1:12 - 0.4, 0, 1:12, c3, col = terra_col, border = NA)
    rect(1:12, 0, 1:12 + 0.4, a4, col = grDevices::adjustcolor(aqua_col, 0.35), border = NA); rect(1:12, 0, 1:12 + 0.4, c4, col = aqua_col, border = NA)
    mtext("Reversible dips in percent green", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("a fall of 25+ points and back up; by month of the low; dark: crossed 50 %", side = 3, line = 0.6, adj = 0, cex = 0.7, col = ink2)
    # B: share of plant-years with a dip, by ecoregion (>= 30 plant-years)
    sm_csv <- file.path(npn_dir, "inseason", "dips_summary.csv")
    SM <- if (file.exists(sm_csv)) utils::read.csv(sm_csv, stringsAsFactors = FALSE) else NULL
    R <- if (!is.null(SM)) SM[SM$by == "region" & SM$plant_years >= 30, ] else NULL
    if (!is.null(R) && nrow(R)) {
      R <- R[order(R$share_with_dip), ]; n <- nrow(R)
      par(mar = c(4.4, 11.5, 3.8, 1))
      plot(NA, xlim = c(0, 1), ylim = c(0.5, n + 0.5), yaxt = "n", xaxt = "n", xlab = "share of plant-years with a dip", ylab = "", bty = "n")
      axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), " %"), cex.axis = 0.8)
      axis(2, at = seq_len(n), labels = reg_lab(R$value, 30), las = 1, tick = FALSE, cex.axis = 0.6)
      rect(0, seq_len(n) - 0.36, R$share_with_dip, seq_len(n) + 0.36, col = grDevices::adjustcolor("#7b4fb3", 0.55), border = NA)
      text(R$share_with_dip, seq_len(n), sprintf("  %d", R$plant_years), adj = 0, cex = 0.6, col = ink2)
      mtext("Where dips happen", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext(sprintf("ecoregions with 30+ plant-years (number shown); %.0f %% of all dips rest on one low reading", 100 * mean(DP$low_readings <= 1, na.rm = TRUE)), side = 3, line = 0.6, adj = 0, cex = 0.66, col = ink2)
      par(mar = c(4.4, 4.6, 3.8, 1))
    } else plot.new()
    # C: weather slopes with 95 % intervals, series-week scale
    f <- FT[FT$scale == "series_week" & FT$stat %in% c("b_vpd", "b_rain"), ]
    f$lab <- sprintf("%s, %s", ifelse(f$pathway == "all", "all grasses", f$pathway), c(all = "whole season", early = "first 60 days", late = "later")[f$half])
    labs <- unique(f$lab); y <- rev(seq_along(labs)); names(y) <- labs
    xl <- range(c(f$lo, f$hi, -0.1, 0.1), na.rm = TRUE)
    par(mar = c(4.4, 10.5, 3.8, 1))
    plot(NA, xlim = xl, ylim = c(0.5, length(labs) + 0.5), yaxt = "n", xlab = "standardized slope (95 % interval over sites)", ylab = "", bty = "n", cex.axis = 0.85)
    axis(2, at = y, labels = labs, las = 1, tick = FALSE, cex.axis = 0.75); abline(v = 0, col = "#d6d5d0")
    sc <- c(b_vpd = "#7b4fb3", b_rain = "#2f7f6f"); off <- c(b_vpd = 0.15, b_rain = -0.15)
    for (k in names(sc)) { g <- f[f$stat == k, ]; yy <- y[g$lab] + off[[k]]
      segments(g$lo, yy, g$hi, yy, col = sc[[k]], lwd = 2); points(g$estimate, yy, pch = 16, col = sc[[k]]) }
    mtext("Weekly greenness change vs weather", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("series-week anomalies; 0 = no response", side = 3, line = 0.6, adj = 0, cex = 0.7, col = ink2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", text.col = ink2, cex = 0.9, pch = c(15, 15, 16, 16), col = c(terra_col, aqua_col, sc),
           legend = c("cool-season (C3)", "warm-season (C4)", "VPD anomaly (7 d before to the reading)", "21-day rain anomaly"))
  })
}

# ---- NPN study: NPN sites and herbaceous live fuel moisture (59d) ------------------------------------------
ls_csv <- file.path(npn_dir, "lfmc", "lfmc_sites.csv"); co_csv <- file.path(npn_dir, "lfmc", "colocation.csv")
if (do("npn_lfmc") && !is.null(TY) && have(ls_csv, co_csv)) {
  LS <- utils::read.csv(ls_csv, stringsAsFactors = FALSE); CO <- utils::read.csv(co_csv, stringsAsFactors = FALSE)
  pa_csv <- file.path(npn_dir, "lfmc", "pairs.csv"); PA <- if (file.exists(pa_csv)) utils::read.csv(pa_csv, stringsAsFactors = FALSE) else NULL
  st <- site_tab()
  np <- if (!is.null(PA) && nrow(PA) >= 10) 3 else 2
  plot_png("npn_lfmc", "NPN and Globe-LFMC 2.0", 2900, 1150, function() {
    layout(matrix(seq_len(np), 1), widths = c(1.7, 1, 1)[seq_len(np)]); par(mar = c(3.4, 0.6, 3.8, 0.6), oma = c(1.6, 0, 0, 0))
    st_map(sprintf("circles: %d NPN study sites; squares: %d herbaceous LFMC sites (Globe-LFMC 2.0), sized by samples", nrow(st), nrow(LS)),
           "Where grass greenness and herbaceous fuel moisture are measured")
    Eo <- eco_l2(); if (!is.null(Eo)) eco_draw(Eo, NA, "#cfcdc7", 0.4)
    points(LS$lon, LS$lat, pch = 22, bg = grDevices::adjustcolor("#2f7f6f", 0.6), col = "#ffffff", lwd = 0.4, cex = 0.5 + 1.4 * sqrt(LS$samples / max(LS$samples)))
    points(st$lon, st$lat, pch = 21, bg = grDevices::adjustcolor(aqua_col, 0.7), col = "#ffffff", lwd = 0.4, cex = 0.5 + 1.4 * sqrt(st$n / max(st$n)))
    par(mar = c(4.4, 4.6, 3.8, 1))
    ks <- c(1, 5, 10, 25, 50, 100, 250); cs <- vapply(ks, function(k) sum(CO$nearest_km <= k), 0); cy <- vapply(ks, function(k) sum(CO$nearest_km <= k & CO$years_in_common > 0), 0)
    plot(NA, xlim = c(0.5, length(ks) + 0.5), ylim = c(0, nrow(CO)), xaxt = "n", xlab = "distance to the nearest herbaceous LFMC site (km)", ylab = "NPN study sites", bty = "n", cex.axis = 0.85)
    axis(1, at = seq_along(ks), labels = paste0("<= ", ks), cex.axis = 0.8)
    rect(seq_along(ks) - 0.38, 0, seq_along(ks) + 0.38, cs, col = grDevices::adjustcolor("#2f7f6f", 0.35), border = NA)
    rect(seq_along(ks) - 0.38, 0, seq_along(ks) + 0.38, cy, col = "#2f7f6f", border = NA)
    text(seq_along(ks), cs, cs, pos = 3, cex = 0.75, col = ink2, xpd = NA)
    mtext("How close are they?", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("dark: with sampling years in common", side = 3, line = 0.6, adj = 0, cex = 0.7, col = ink2)
    if (np == 3) {
      plot(PA$pct_green, PA$lfmc, pch = 16, cex = 0.8, col = grDevices::adjustcolor(ifelse(PA$same_species, aqua_col, "#9c9b96"), 0.7),
           xlab = "NPN percent green (species mean, nearest reading)", ylab = "herbaceous LFMC (%)", bty = "n", cex.axis = 0.85, xlim = c(0, 100))
      abline(h = c(30, 120), col = "#d6d5d0", lty = 2)
      r <- suppressWarnings(stats::cor(PA$pct_green, PA$lfmc, method = "spearman"))
      mtext("Where they meet", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
      mtext(sprintf("%d pairs at %d site(s), +/- 7 days; Spearman r %.2f; orange: same species", nrow(PA), length(unique(PA$site_id)), r), side = 3, line = 0.6, adj = 0, cex = 0.7, col = ink2)
    }
  })
}

# ---- NPN study: representativeness against U.S. grass fuels (59c, 58) ---------------------------
re_csv <- file.path(rec_dir, "representativeness_ecoregion.csv"); gs_csv <- file.path(rec_dir, "grass_fuel_sample.csv")
if (do("npn_repr") && !is.null(TY) && have(re_csv, gs_csv)) {
  RE <- utils::read.csv(re_csv, stringsAsFactors = FALSE, colClasses = c(region = "character"))
  G <- utils::read.csv(gs_csv, stringsAsFactors = FALSE)
  RE <- RE[RE$grass_fuel_area_share >= 0.01 | RE$plant_year_share >= 0.01, ]; RE <- RE[order(RE$grass_fuel_area_share), ]
  st <- stats::aggregate(plant_years ~ site_id + lon + lat, TY, sum)
  n <- nrow(RE)
  plot_png("npn_repr", "LANDFIRE FBFM40 and the record", 2900, max(1150, 330 + 58 * n), function() {
    layout(matrix(1:2, 1), widths = c(1.45, 1)); par(mar = c(3.4, 0.6, 3.8, 0.6), oma = c(2, 0, 0, 0))
    st_map(sprintf("dots: a sample of grass (GR) and grass-shrub (GS) fuel cells; circles: %d sites sized by plant-years", nrow(st)),
           "Where U.S. grass fuels are, and where the record is")
    Eo <- eco_l2(); if (!is.null(Eo)) eco_draw(Eo, NA, "#b8b7b2", 0.5)
    gcl <- ifelse(substr(G$fm, 1, 2) == "GR", "#c9a227", "#8a9a3b")
    points(G$lon, G$lat, pch = 16, cex = 0.16, col = grDevices::adjustcolor(gcl, 0.4))
    points(st$lon, st$lat, pch = 21, bg = grDevices::adjustcolor(aqua_col, 0.7), col = bg, lwd = 0.6, cex = 0.5 + 1.6 * sqrt(st$plant_years / max(st$plant_years)))
    par(mar = c(4.4, 12, 3.8, 1.5))
    y <- seq_len(n); xm <- max(c(RE$grass_fuel_area_share, RE$plant_year_share, RE$site_share)) * 100 * 1.08
    plot(NA, xlim = c(0, xm), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = "share (%)", ylab = "", bty = "n", cex.axis = 0.85)
    axis(2, at = y, labels = reg_lab(RE$region, 30), tick = FALSE, las = 1, cex.axis = 0.62, line = -0.4)
    rect(0, y - 0.38, 100 * RE$grass_fuel_area_share, y + 0.38, col = "#dedcd6", border = NA)
    rect(0, y - 0.16, 100 * RE$plant_year_share, y + 0.16, col = aqua_col, border = NA)
    points(100 * RE$site_share, y, pch = 21, bg = terra_col, col = bg, cex = 1.1)
    mtext("Share of U.S. grass fuel area vs share of the record", side = 3, line = 2, adj = 0, cex = 0.95, font = 2, col = ink)
    mtext("EPA Level II ecoregions with 1 % or more of either", side = 3, line = 0.6, adj = 0, cex = 0.72, col = ink2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", text.col = ink2, cex = 0.9, pch = c(16, 16, 21, 15, 15, 21), pt.cex = c(1, 1, 1.3, 1.8, 1.8, 1.3),
           col = c("#c9a227", "#8a9a3b", bg, "#dedcd6", aqua_col, bg), pt.bg = c(NA, NA, aqua_col, NA, NA, terra_col),
           legend = c("grass fuel (GR)", "grass-shrub fuel (GS)", "study site", "grass fuel area", "plant-years", "sites"))
  })
}

log_msg(sprintf("=== done | figures in %s | log: %s", fig_dir, log_file))
