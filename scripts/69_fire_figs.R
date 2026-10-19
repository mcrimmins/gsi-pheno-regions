#!/usr/bin/env Rscript
# Step 69 (fire activity study; laptop): figures for report/fire_study.qmd from the synced CSVs.
# Same conventions as 90_report_figs.R: <fig_dir>/<figure>_<profile>.png + manifest.csv rows
# (fig_dir = report/figs on the laptop; <run_dir>/report_figs where config report: figs_in_run_dir).
#   fire_seasonal  per pilot region: mean weekly fires (as a share of the region's peak week) and
#                  the typical share of herbaceous fuel cured, one 0-1 axis, by week (65 / 66)
#   fire_skill     held-out gain from adding phenology, per region and outcome, over the simple
#                  (60-day water balance) and the 1000-h baselines, with bootstrap 90 % intervals (66)
# Usage: Rscript scripts/69_fire_figs.R [profile] [figure ...]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

all_figs <- c("fire_seasonal", "fire_skill")
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
log_msg("=== done | figures in ", fig_dir)
