#!/usr/bin/env Rscript
# Step 40 (product): package the green-up and curing rules decided on Oct 13 2026 after the
# ground checks (38, 39):
#   - each cell's typical green-up and curing date (median satellite NDVI green-up sos20 and
#     end of season eos50 over the record; summaries from 06) everywhere;
#   - curing shifted by a weather adjustment in config rules: adjust_regions (monsoon Southwest,
#     California annual grasslands): cues from PRISM up to the typical date minus lead_days,
#     as in 38 (R/cues.R adj_cues), refitted here on ALL years of the 4,000-cell sample.
# Outputs <run_dir>/rules/:
#   typical_dates.tif      greenup_doy, curing_doy (median, day of year; < 1 = before Jan 1,
#                          > 365 = into the next year), greenup_spread, curing_spread (IQR across
#                          years, days), herb (1 = strict herbaceous mask), region (A2 k13 map
#                          number), adjust (1 = curing adjustment applies: herbaceous cell in an
#                          adjust region with a typical curing date)
#   adjust_coefs.csv       region, cue, b (days per unit of cue departure), r2_fit, cells, lead
#   adjust_cue_means.tif   each selected cue's climatology per cell (mean over the record, at
#                          that cell's own cutoff), needed to turn a year's cues into departures
#   README.txt             what each file holds and how 41 uses them
# Runs on the P720 (PRISM daily, satellite metrics, the 34 sample cache with eos50). The per-cell
# cue climatology reads daily PRISM for every adjusted cell (forked workers by chunk).
#
# Usage: Rscript scripts/40_rules_build.R [profile]

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
source(here::here("R", "block2.R"))
source(here::here("R", "cues.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("40_rules_build", profile)
rc <- cfg$rules %||% list()
out_dir <- file.path(cfg$run_dir, "rules"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()
lead <- rc$lead_days %||% 14; max_adj <- rc$max_adj_days %||% 60; jump <- rc$max_abs_anomaly %||% 90
adj_regions <- as.integer(unlist(rc$adjust_regions %||% list(3, 11)))
log_msg("=== 40_rules_build | profile: ", profile, " | adjust regions ", paste(adj_regions, collapse = ", "), " | lead ", lead, " d")

# ---- 1. typical dates, spread, mask, regions -----------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
med <- land_values(terra::rast(summary_file(cfg, "cutb_metrics", "median")), grid, c("sos20", "eos50"))
iqr <- land_values(terra::rast(summary_file(cfg, "cutb_metrics", "iqr")), grid, c("sos20", "eos50"))
hm <- land_values(terra::rast(file.path(static_dir(cfg), "herb_mask.tif")), grid)[, 1]
herb <- as.integer(!is.na(hm) & hm > 0)
D <- region_display(cfg, grid, list(list(partition = rc$partition %||% "A2_groups", k = rc$k %||% 13)))[[1]]
region <- D$num
adjust <- as.integer(herb == 1 & region %in% adj_regions & !is.na(med[, "eos50"]))
to_rast <- function(v, name) { r <- terra::rast(grid); x <- rep(NA_real_, terra::ncell(r)); x[land] <- v; terra::values(r) <- x; names(r) <- name; r }

# ---- 2. adjustment coefficients, refitted on all years of the sample -------------------------------
cache <- file.path(cfg$run_dir, "eval", "gsi_model", sprintf("sample_%s.rds", profile))
if (!file.exists(cache)) stop("no ", cache, ": run 34 (and 35 with GSI_CURING_TARGET=eos50) first")
S <- readRDS(cache)
if (is.null(S$obs$eos50)) stop("sample cache lacks eos50")
years <- S$years; nY <- length(years)
reg_s <- region[S$cells$row]
CO <- list()
for (r in adj_regions) {
  i <- which(reg_s == r)
  if (length(i) < 25) { log_warn("region ", r, ": only ", length(i), " sampled cells, no adjustment fitted"); next }
  Y <- S$obs$eos50[i, , drop = FALSE]
  typ <- apply(Y, 1, stats::median, na.rm = TRUE)
  CU <- lapply(stats::setNames(adj_cue_names, adj_cue_names), function(x) matrix(NA_real_, length(i), nY))
  for (j in seq_len(nY)) {
    v <- adj_cues(S$drv$ppt[i, , drop = FALSE], S$drv$tmax[i, , drop = FALSE], S$dates, years[j], typ - lead, S$obs$sos[i, j])
    for (cn in adj_cue_names) CU[[cn]][, j] <- v[, cn]
  }
  X <- lapply(CU, function(x) x - rowMeans(x, na.rm = TRUE))
  dep <- Y - typ; dep[abs(dep) > jump] <- NA
  fx <- adj_fit(dep, X, X, rc$max_cues %||% 3, rc$min_gain %||% 0.01, max_adj)
  if (!length(fx$cues)) { log_warn("region ", r, ": no cue passed the selection; no adjustment"); next }
  CO[[length(CO) + 1]] <- data.frame(region = r, cue = fx$cues, b = unname(fx$b), r2_fit = fx$r2_fit, cells = length(i), lead = lead)
  log_msg(sprintf("region %d: %d sampled cells | cues %s | b %s | fitting R2 %.2f", r, length(i),
                  paste(fx$cues, collapse = ", "), paste(sprintf("%.3f", fx$b), collapse = ", "), fx$r2_fit))
}
if (!length(CO)) stop("no adjustment fitted in any region")
CO <- do.call(rbind, CO)
utils::write.csv(CO, file.path(out_dir, "adjust_coefs.csv"), row.names = FALSE)
adjust <- as.integer(adjust == 1 & region %in% CO$region)        # only regions with a fitted adjustment
# confidence: 1 = herbaceous cell in a region whose rule was checked against NPN grass sites,
# 2 = herbaceous cell in a region with few or no ground sites, 3 = mixed cover (outside the
# strict herbaceous mask: the satellite date describes forest, crops and grass together)
tested <- as.integer(unlist(rc$tested_regions %||% list(2, 3, 5, 10, 11, 13)))
confidence <- ifelse(is.na(med[, "eos50"]) & is.na(med[, "sos20"]), NA, ifelse(herb == 1, ifelse(region %in% tested, 1, 2), 3))
log_msg(sprintf("confidence: %d ground-checked herbaceous, %d other herbaceous, %d mixed-cover cells",
                sum(confidence == 1, na.rm = TRUE), sum(confidence == 2, na.rm = TRUE), sum(confidence == 3, na.rm = TRUE)))
TD <- c(to_rast(med[, "sos20"], "greenup_doy"), to_rast(med[, "eos50"], "curing_doy"),
        to_rast(iqr[, "sos20"], "greenup_spread"), to_rast(iqr[, "eos50"], "curing_spread"),
        to_rast(herb, "herb"), to_rast(region, "region"), to_rast(adjust, "adjust"), to_rast(confidence, "confidence"))
terra::writeRaster(TD, file.path(out_dir, "typical_dates.tif"), overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
log_msg(sprintf("typical dates: %d land cells, %d herbaceous, %d with the curing adjustment", length(land), sum(herb), sum(adjust)))
rm(S); invisible(gc())

# ---- 3. per-cell cue climatology for every adjusted cell ----------------------------------------------
cues_used <- unique(CO$cue)
rows <- which(adjust == 1)
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
pfile <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
py <- seq(min(years) - 1, max(years))
chunks <- split(rows, ceiling(seq_along(rows) / (rc$chunk_cells %||% 2500)))
nw <- min(as.integer(rc$workers %||% 8), length(chunks))
fork <- .Platform$OS.type == "unix" && nw > 1
pmap <- function(X, f) if (fork) parallel::mclapply(X, f, mc.cores = nw, mc.preschedule = FALSE) else lapply(X, f)
log_msg(sprintf("cue climatology: %d cells in %d chunks, %s; cues %s", length(rows), length(chunks),
                if (fork) paste(nw, "workers") else "sequential", paste(cues_used, collapse = ", ")))
clim <- pmap(seq_along(chunks), function(k) {
  rr <- chunks[[k]]; cells <- land[rr]
  rd <- function(v) do.call(cbind, lapply(py, function(y) { r <- terra::rast(pfile(v, y)); x <- as.matrix(terra::extract(r, cells)); colnames(x) <- names(r); x }))
  P <- rd("ppt"); TX <- rd("tmax"); dates <- as.Date(colnames(P))
  sos <- if ("greenup" %in% cues_used) vapply(years, function(Y) if (file.exists(mfile(Y))) terra::extract(terra::rast(mfile(Y))[["sos20"]], cells)[, 1] else rep(NA_real_, length(cells)), numeric(length(cells))) else NULL
  cut <- med[rr, "eos50"] - lead
  sm <- matrix(0, length(rr), length(cues_used), dimnames = list(NULL, cues_used)); nn <- sm
  for (j in seq_along(years)) {
    v <- adj_cues(P, TX, dates, years[j], cut, if (is.null(sos)) NULL else sos[, j])[, cues_used, drop = FALSE]
    ok <- !is.na(v); v[!ok] <- 0; sm <- sm + v; nn <- nn + ok
  }
  out <- sm / nn; out[nn == 0] <- NA
  if (k %% 5 == 0) log_msg(sprintf("  chunk %d / %d done", k, length(chunks)))
  out
})
bad <- vapply(clim, inherits, TRUE, "try-error")
if (any(bad)) stop("cue climatology failed: ", clim[bad][[1]])
CM <- do.call(rbind, clim)
CMR <- do.call(c, lapply(cues_used, function(cn) { v <- rep(NA_real_, length(land)); v[rows] <- CM[, cn]; to_rast(v, cn) }))
terra::writeRaster(CMR, file.path(out_dir, "adjust_cue_means.tif"), overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))

# ---- 4. README -------------------------------------------------------------------------------------
cue_txt <- c(p30 = "rain in the 30 days before the cutoff (mm)", p90 = "rain in the 90 days before (mm)",
             p180 = "rain in the 180 days before (mm)", tx30 = "mean daily maximum temperature, 30 days before (C)",
             tx60 = "mean daily maximum temperature, 60 days before (C)", dry60 = "days under 1 mm in the 60 days before",
             since_rain = "days since the last day with 10 mm or more (at most 180)",
             greenup = "this year's satellite green-up date, if before the cutoff (day of year)")
writeLines(c(
  "Herbaceous green-up and curing rules (fuel phenology regions project), built by scripts/40_rules_build.R",
  sprintf("Built %s from satellite years %d-%d (MODIS NDVI curves, 4 km PRISM grid).", format(Sys.Date()), min(years), max(years)),
  "",
  "Rules (decided Oct 13 2026 after checks against NPN grass observations):",
  "  green-up, everywhere: the cell's typical green-up date (greenup_doy)",
  "  curing, everywhere except adjust cells: the cell's typical curing date (curing_doy)",
  sprintf("  curing, adjust cells (regions %s): curing_doy + adjustment, from weather up to %d days before curing_doy", paste(unique(CO$region), collapse = ", "), lead),
  "",
  "typical_dates.tif layers:",
  "  greenup_doy, curing_doy  median satellite date over the record, day of year relative to Jan 1 (every land cell)",
  "                           (values below 1 fall in the previous year, above 365 in the next)",
  "  greenup_spread, curing_spread  interquartile range across years (days)",
  "  herb      1 = strict herbaceous mask (the rules were tested on grass- and shrub-dominated cells)",
  "  region    A2 climate region, 13 regions, project-page numbering",
  "  adjust    1 = curing adjustment applies",
  "  confidence  1 = herbaceous cell in a region whose rule was checked against NPN grass sites;",
  "              2 = herbaceous cell, few or no ground sites in its region;",
  "              3 = mixed cover: the date describes forest, crops and grass together, not grass fuels",
  "",
  "adjust_coefs.csv: per region, the cues and days of curing shift per unit of cue departure (b).",
  "adjust_cue_means.tif: each cue's mean over the record at each adjusted cell's own cutoff.",
  "Cues:", paste0("  ", cues_used, ": ", cue_txt[cues_used]),
  "",
  "Adjustment for year Y = sum over the region's cues of b x (cue in Y - cue mean), each cue",
  sprintf("computed from daily PRISM up to curing_doy - %d days; a missing cue counts as 0; clamped to +/- %d days.", lead, max_adj),
  "scripts/41_rules_predict.R <year> [as-of date] writes rules/predict_<year>.tif.",
  "",
  "Caveats: the typical dates are satellite dates; NDVI end of season runs about 2-3 weeks after",
  "ground curing (50 % green) on average. The adjustment in California rests on 7 ground sites."),
  file.path(out_dir, "README.txt"))
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
