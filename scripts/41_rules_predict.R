#!/usr/bin/env Rscript
# Step 41 (product): this year's expected green-up and curing dates from the packaged rules
# (40). Typical dates everywhere; in the adjusted cells, the typical curing date plus the
# weather adjustment, from daily PRISM up to each cell's cutoff (typical curing date minus the
# lead). With an as-of date, cells whose cutoff is still ahead keep the typical date (flag 0).
# The green-up cue uses that year's satellite green-up when cutb/metrics/metrics_<Y>.tif exists,
# otherwise it counts as typical.
# Output <run_dir>/rules/predict_<Y>.tif: greenup_typical, curing_typical, curing_adjustment
# (days, 0 where not applied), curing_expected, adjusted (1 = adjustment applied), region.
# When the year's satellite end of season exists (a past year), also predict_<Y>_check.csv:
# per adjusted region, error of the typical and the expected date against it. The coefficients
# were fitted on all years (40), so this is an in-sample check, not a test.
#
# Usage: Rscript scripts/41_rules_predict.R <year> [as-of date YYYY-MM-DD] [profile]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript scripts/41_rules_predict.R <year> [as-of YYYY-MM-DD] [profile]")
Y <- as.integer(args[1])
asof <- if (length(args) >= 2 && grepl("^\\d{4}-\\d{2}-\\d{2}$", args[2])) as.Date(args[2]) else NULL
profile <- if (length(args) >= 2 && is.null(asof)) args[2] else if (length(args) >= 3) args[3] else
  if (exists("gsi_profile", envir = globalenv())) base::get("gsi_profile", envir = globalenv()) else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "block2.R"))
source(here::here("R", "cues.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("41_rules_predict", profile)
rc <- cfg$rules %||% list()
rdir <- file.path(cfg$run_dir, "rules")
terra::terraOptions(progress = 0)
t0 <- Sys.time()
lead <- rc$lead_days %||% 14; max_adj <- rc$max_adj_days %||% 60
log_msg("=== 41_rules_predict | ", Y, if (!is.null(asof)) paste(" as of", asof) else "", " | profile: ", profile)

for (f in c("typical_dates.tif", "adjust_coefs.csv", "adjust_cue_means.tif"))
  if (!file.exists(file.path(rdir, f))) stop("no ", file.path(rdir, f), ": run 40_rules_build.R first")
TD <- terra::rast(file.path(rdir, "typical_dates.tif"))
CO <- utils::read.csv(file.path(rdir, "adjust_coefs.csv"), stringsAsFactors = FALSE)
CMR <- terra::rast(file.path(rdir, "adjust_cue_means.tif"))
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
v <- terra::values(TD, mat = TRUE)[land, , drop = FALSE]
rows <- which(v[, "adjust"] == 1 & v[, "region"] %in% CO$region)
cells <- land[rows]
log_msg(sprintf("%d adjusted cells in regions %s", length(rows), paste(unique(CO$region), collapse = ", ")))

pfile <- function(var, y) file.path(cfg$run_dir, "prism_daily", var, sprintf("prism_%s_%d.tif", var, y))
mfile <- function(y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", y))
py <- c(Y - 1, Y)
miss <- py[!file.exists(vapply(py, function(y) pfile("ppt", y), ""))]
if (length(miss)) stop("no daily PRISM for ", paste(miss, collapse = ", "), " under ", file.path(cfg$run_dir, "prism_daily"))
rd <- function(var) do.call(cbind, lapply(py, function(y) { r <- terra::rast(pfile(var, y)); x <- as.matrix(terra::extract(r, cells)); colnames(x) <- names(r); x }))
P <- rd("ppt"); TX <- rd("tmax"); dates <- as.Date(colnames(P))
sos <- if (file.exists(mfile(Y))) terra::extract(terra::rast(mfile(Y))[["sos20"]], cells)[, 1] else NULL
if (is.null(sos)) log_msg("no satellite green-up for ", Y, ": the green-up cue counts as typical")
cut <- v[rows, "curing_doy"] - lead
cues <- adj_cues(P, TX, dates, Y, cut, sos)
cut_date <- as.Date(sprintf("%d-01-01", Y)) + round(cut) - 1
last_day <- max(dates[colSums(!is.na(P)) > 0])
ready <- cut_date <= last_day & (if (is.null(asof)) TRUE else cut_date <= asof)
ready[is.na(ready)] <- FALSE
means <- terra::extract(CMR, cells)
adj <- rep(0, length(rows))
for (r in unique(CO$region)) {
  i <- which(v[rows, "region"] == r)
  adj[i] <- adj_apply(cues[i, , drop = FALSE], CO[CO$region == r, c("cue", "b")],
                      lapply(stats::setNames(names(means), names(means)), function(cn) means[[cn]][i]), max_adj)
}
adj[!ready] <- 0
to_rast <- function(x, name) { r <- terra::rast(grid); z <- rep(NA_real_, terra::ncell(r)); z[land] <- x; terra::values(r) <- z; names(r) <- name; r }
a_all <- rep(NA_real_, length(land)); a_all[v[, "herb"] == 1] <- 0; a_all[rows] <- adj
flag <- rep(NA_real_, length(land)); flag[!is.na(v[, "curing_doy"])] <- 0; flag[rows[ready]] <- 1
OUT <- c(to_rast(v[, "greenup_doy"], "greenup_typical"), to_rast(v[, "curing_doy"], "curing_typical"),
         to_rast(a_all, "curing_adjustment"), to_rast(v[, "curing_doy"] + ifelse(is.na(a_all), 0, a_all), "curing_expected"),
         to_rast(flag, "adjusted"), to_rast(v[, "region"], "region"))
of <- file.path(rdir, sprintf("predict_%d.tif", Y))
terra::writeRaster(OUT, of, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))
utils::write.csv(data.frame(year = Y, as_of = format(if (is.null(asof)) last_day else min(asof, last_day)),
                           weather_to = format(last_day), cells_adjusted = sum(ready), cells_pending = sum(!ready)),
                file.path(rdir, sprintf("predict_%d_info.csv", Y)), row.names = FALSE)
log_msg(sprintf("adjustment applied in %d of %d cells; median %+.1f d, 10-90 %%: %+.1f to %+.1f d", sum(ready), length(rows),
                stats::median(adj[ready]), stats::quantile(adj[ready], 0.1), stats::quantile(adj[ready], 0.9)))

# in-sample check against that year's satellite end of season (past years)
if (file.exists(mfile(Y))) {
  cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
  eos <- terra::extract(terra::rast(mfile(Y))[["eos50"]], cells)[, 1]
  typ <- v[rows, "curing_doy"]; dt <- cdiff(typ, eos); de <- cdiff(typ + adj, eos)
  ok <- !is.na(dt) & abs(dt) <= 90 & ready
  ck <- do.call(rbind, lapply(unique(CO$region), function(r) {
    u <- ok & v[rows, "region"] == r
    data.frame(year = Y, region = r, cells = sum(u), mae_typical = mean(abs(dt[u])), mae_expected = mean(abs(de[u])),
               r_adj = if (sum(u) > 30) stats::cor(adj[u], -dt[u]) else NA)
  }))
  utils::write.csv(ck, file.path(rdir, sprintf("predict_%d_check.csv", Y)), row.names = FALSE)
  for (i in seq_len(nrow(ck))) with(ck[i, ], log_msg(sprintf("check vs satellite %d, region %d: %d cells | typical %.1f d, expected %.1f d | r %.2f (in-sample)",
                                                           year, region, cells, mae_typical, mae_expected, r_adj)))
}
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), of, log_file))
