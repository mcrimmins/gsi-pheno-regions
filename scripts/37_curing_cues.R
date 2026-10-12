#!/usr/bin/env Rscript
# Step 37 (evaluation, curing cues): in the regions where no model beat each cell's typical
# curing date (35), is there a weather cue that explains early and late curing years at all?
# Uses the 4,000-cell sample from 34 (daily PRISM + yearly satellite dates; R/cues.R for the
# cues). For each cell-year: the satellite curing date (NDVI end of season eos50 and NDII7
# cure50) and green-up date (sos20) as departures from the cell's own mean, and each cue's
# departure. Per working region (A2, 13 regions, page numbering) and target:
#   r      correlation of curing and cue departures, pooled over the region's cells and years
#   cv_r2  share of the curing departures that a small linear model of up to three cues
#          predicts in held-out years (odd / even years); cues picked on the fitting years
# Green-up departures are a cue too ("greenup": does curing follow a fixed lag after
# green-up, as seed set in annual grasses would). Season jumps (curing more than 90 days from
# the cell's median) are dropped, as in 33. For comparison, the best curing model's anomaly
# correlation from 35 is in the GSI / bucket region scores.
# Outputs <run_dir>/eval/curing_cues/: cue_r.csv, cue_models.csv, cues_<profile>.rds (cache).
# Runs on the P720 (needs eval/gsi_model/sample_<profile>.rds with eos50; a few minutes).
#
# Usage: Rscript scripts/37_curing_cues.R [profile]

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
log_file <- log_init("37_curing_cues", profile)
cc <- cfg$curing_cues %||% list()
out_dir <- file.path(cfg$run_dir, "eval", "curing_cues"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 37_curing_cues | profile: ", profile)

cache <- file.path(cfg$run_dir, "eval", "gsi_model", sprintf("sample_%s.rds", profile))
if (!file.exists(cache)) stop("no ", cache, ": run 34 (and 35 with GSI_CURING_TARGET=eos50) first")
S <- readRDS(cache)
if (is.null(S$obs$eos50)) stop("sample cache lacks eos50: run 35 with GSI_CURING_TARGET=eos50 once to add it")
years <- S$years; nY <- length(years); n <- nrow(S$cells)

# ---- cues per cell-year (cached) -----------------------------------------------------------
cfile <- file.path(out_dir, sprintf("cues_%s.rds", profile))
if (file.exists(cfile)) { CU <- readRDS(cfile); log_msg("cues from cache") } else {
  CU <- lapply(stats::setNames(cue_names, cue_names), function(x) matrix(NA_real_, n, nY))
  for (j in seq_along(years)) {
    v <- cue_year(S$drv$ppt, S$drv$tmin, S$drv$tmax, S$dates, years[j], cfg$features$block2)
    for (cn in cue_names) CU[[cn]][, j] <- v[, cn]
  }
  saveRDS(CU, cfile)
  log_msg(sprintf("cues: %d cells x %d years, %.1f min", n, nY, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
for (cn in cue_names) log_msg(sprintf("  %-12s median %7.1f  (share missing %.0f %%)", cn, stats::median(CU[[cn]], na.rm = TRUE),
                                      100 * mean(is.na(CU[[cn]]))))

# ---- regions (page numbering) ------------------------------------------------------------------
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
rs <- list(partition = cc$partition %||% "A2_groups", k = cc$k %||% 13)
D <- region_display(cfg, grid, list(rs))[[1]]
if (is.null(D)) stop("region set ", rs$partition, " k", rs$k, " not found")
reg <- D$num[S$cells$row]
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))
max_dev <- cc$max_abs_anomaly %||% 90

# ---- screening ---------------------------------------------------------------------------------
targets <- list(eos50 = S$obs$eos50, cure50 = S$obs$cure, sos20 = S$obs$sos)
rows_r <- list(); rows_m <- list()
for (tn in names(targets)) {
  Y <- targets[[tn]]
  Y[abs(Y - apply(Y, 1, stats::median, na.rm = TRUE)) > max_dev] <- NA     # season jumps
  YA <- cue_anom(Y)
  C <- lapply(CU, function(x) { x[is.na(Y)] <- NA; cue_anom(x) })
  if (tn != "sos20") { g <- S$obs$sos; g[is.na(Y)] <- NA; C$greenup <- cue_anom(g) }
  for (rn in c("CONUS", sort(unique(reg[!is.na(reg)])))) {
    i <- if (rn == "CONUS") seq_len(n) else which(reg == as.integer(rn))
    if (length(i) < (cc$min_cells %||% 25)) next
    sc <- cue_screen(YA[i, , drop = FALSE], lapply(C, function(x) x[i, , drop = FALSE]), folds,
                     cc$max_cues %||% 3, cc$min_gain %||% 0.01)
    rows_r[[length(rows_r) + 1]] <- data.frame(target = tn, region = rn, cells = length(i), cue = names(sc$r), n = sc$n, r = sc$r)
    rows_m[[length(rows_m) + 1]] <- data.frame(target = tn, region = rn, cells = length(i), cv_r2 = sc$cv_r2,
                                               cues_odd = sc$cues[["odd"]], cues_even = sc$cues[["even"]],
                                               anom_sd = stats::sd(YA[i, ], na.rm = TRUE))
  }
}
R <- do.call(rbind, rows_r); Mo <- do.call(rbind, rows_m)
utils::write.csv(R, file.path(out_dir, "cue_r.csv"), row.names = FALSE)
utils::write.csv(Mo, file.path(out_dir, "cue_models.csv"), row.names = FALSE)
for (i in seq_len(nrow(Mo))) with(Mo[i, ], {
  rr <- R[R$target == target & R$region == region, ]; rr <- rr[order(-abs(rr$r)), ][1:3, ]
  log_msg(sprintf("%-6s region %-5s %4d cells | departures SD %4.1f d | CV R2 %5.2f (%s / %s) | top r: %s", target, region, cells,
                  anom_sd, cv_r2, cues_odd, cues_even, paste(sprintf("%s %+.2f", rr$cue, rr$r), collapse = ", ")))
})
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
