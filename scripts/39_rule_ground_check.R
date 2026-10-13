#!/usr/bin/env Rscript
# Step 39 (evaluation, curing rule vs the ground): does "typical date + weather adjustment"
# (38) also move curing dates toward ground observations, or only toward the satellite target?
# For each NPN grass site (graminoid % green falling through 50 %, and cured date) the rule is
# applied at the site's 4 km cell exactly as in 38: the typical date is the cell's median
# satellite end of season (eos50) over the fitting years (odd or even), the cues are computed
# from PRISM up to the typical date minus 14 days, and the adjustment uses the coefficients 38
# fitted for the site's region and fold (coefs.csv; CONUS coefficients where the region has
# none). Each held-out site-year is then compared with the ground date:
#   mae_typical / mae_adjusted   mean absolute difference from the ground date (days, circular,
#                                capped at 90), site-years without a season mismatch
#   dev_typical / dev_adjusted   the same after removing each site's average offset from the
#                                ground (satellite end of season and ground curing differ by
#                                ~20 days on average; the adjustment cannot change that): does
#                                the adjustment move each year toward the ground's early or late?
#   r_adj                        correlation of the adjustment with the ground's departure from
#                                its own site mean
#   mismatch                     share of site-years whose typical date is > 90 days from the ground
# Runs on the P720 (PRISM daily + cutb metrics there; ground table npn_curing_metrics.rds from the
# GSI study, folder given by GSI_OBS_DIR or config evaluation: study_dir). Needs 38's coefs.csv.
# Outputs <run_dir>/eval/curing_adjust/ground_site_years.csv, ground_scores.csv.
#
# Usage: GSI_OBS_DIR=<folder with npn_curing_metrics.rds> Rscript scripts/39_rule_ground_check.R [profile]

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
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("39_rule_ground_check", profile)
ca <- cfg$curing_adjust %||% list()
out_dir <- file.path(cfg$run_dir, "eval", "curing_adjust")
terra::terraOptions(progress = 0)
t0 <- Sys.time()
log_msg("=== 39_rule_ground_check | profile: ", profile)
lead <- as.integer(ca$ground_lead_days %||% 14); target <- "eos50"
cap <- ca$cap_days %||% 90; max_adj <- ca$max_adj_days %||% 60

cf_file <- file.path(out_dir, "coefs.csv")
if (!file.exists(cf_file)) stop("no ", cf_file, ": run 38 first")
CF <- utils::read.csv(cf_file, stringsAsFactors = FALSE)
CF <- CF[CF$target == target & CF$lead == lead, ]
if (!nrow(CF)) stop("coefs.csv has no rows for ", target, " at lead ", lead)

# ---- ground site-years, cells, regions -------------------------------------------------------
ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
sets <- c("grass_curing", "grass_cured")
E <- obs_events(ddir, sets)
if (!nrow(E)) stop("no NPN grass curing observations in ", ddir, " (set GSI_OBS_DIR to the folder with npn_curing_metrics.rds)")
grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
land <- which(!is.na(terra::values(grid, mat = FALSE)))
E$cell <- terra::cellFromXY(grid, as.matrix(E[, c("lon", "lat")]))
E$row <- match(E$cell, land); E <- E[!is.na(E$row), ]
D <- region_display(cfg, grid, list(list(partition = ca$partition %||% "A2_groups", k = ca$k %||% 13)))[[1]]
E$region <- D$num[E$row]

# ---- satellite dates and PRISM at the site cells -----------------------------------------------
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
pfile <- function(v, y) file.path(cfg$run_dir, "prism_daily", v, sprintf("prism_%s_%d.tif", v, y))
sat_years <- cfg$years[file.exists(vapply(cfg$years, mfile, ""))]
years <- intersect(feature_years(cfg), sat_years); nY <- length(years)
E <- E[E$year %in% years, ]
U <- sort(unique(E$row)); nU <- length(U)
log_msg(sprintf("ground: %d site-years at %d sites in %d cells; years %d-%d", nrow(E), length(unique(E$site)), nU, min(years), max(years)))
sat <- function(b) vapply(years, function(Y) terra::extract(terra::rast(mfile(Y))[[b]], land[U])[, 1], numeric(nU))
EOS <- sat(target); SOS <- sat("sos20")
py <- seq(min(years) - 1, max(years))
read_series <- function(v) do.call(cbind, lapply(py, function(y) {
  r <- terra::rast(pfile(v, y)); x <- as.matrix(terra::extract(r, land[U])); colnames(x) <- names(r); x }))
P <- read_series("ppt"); TX <- read_series("tmax"); dates <- as.Date(colnames(P))
log_msg(sprintf("drivers: %d cells x %d days (%.1f min)", nU, length(dates), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- the rule at the site cells, per fold ------------------------------------------------------
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))
reg_u <- D$num[U]
TYP <- ADJ <- matrix(NA_real_, nU, nY)
for (f in names(folds)) {
  tr <- folds[[f]]; te <- setdiff(seq_len(nY), tr)
  typ <- apply(EOS[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
  CU <- lapply(stats::setNames(adj_cue_names, adj_cue_names), function(x) matrix(NA_real_, nU, nY))
  for (j in seq_len(nY)) {
    v <- adj_cues(P, TX, dates, years[j], typ - lead, SOS[, j])
    for (cn in adj_cue_names) CU[[cn]][, j] <- v[, cn]
  }
  X <- lapply(CU, function(x) { a <- x - rowMeans(x[, tr, drop = FALSE], na.rm = TRUE); a[is.na(a)] <- 0; a })
  for (i in seq_len(nU)) {
    cf <- CF[CF$fold == f & CF$region == as.character(reg_u[i]), ]
    if (!nrow(cf) && !(as.character(reg_u[i]) %in% CF$region)) cf <- CF[CF$fold == f & CF$region == "CONUS", ]
    a <- rep(0, length(te))
    for (k in seq_len(nrow(cf))) a <- a + cf$b[k] * X[[cf$cue[k]]][i, te]
    ADJ[i, te] <- pmax(pmin(a, max_adj), -max_adj); TYP[i, te] <- typ[i]
  }
}

# ---- compare with the ground -----------------------------------------------------------------
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
iu <- match(E$row, U); jy <- match(E$year, years)
E$typical <- TYP[cbind(iu, jy)]; E$adjust <- ADJ[cbind(iu, jy)]
E <- E[!is.na(E$typical), ]
E$d_typ <- cdiff(E$typical, E$doy); E$d_adj <- cdiff(E$typical + E$adjust, E$doy)
E$mismatch <- abs(E$d_typ) > 90
# departures from each site's own mean (season mismatches left out): ground unwrapped around the
# site median; the rule's typical date is the same every test year of a fold, so its departure
# is the adjustment minus the site's mean adjustment
k <- !E$mismatch
key <- paste(E$set, E$site)
gm <- stats::ave(E$doy, key, FUN = function(v) stats::median(v, na.rm = TRUE)); gu <- gm + cdiff(E$doy, gm)
E$g_dev <- E$a_dev <- NA_real_
nk <- stats::ave(as.numeric(k), key, FUN = sum)
w <- k & nk >= 3
E$g_dev[w] <- gu[w] - stats::ave(gu[w], key[w]); E$a_dev[w] <- E$adjust[w] - stats::ave(E$adjust[w], key[w])
utils::write.csv(E, file.path(out_dir, "ground_site_years.csv"), row.names = FALSE)

score <- function(x) {
  u <- !x$mismatch; v <- !is.na(x$g_dev)
  data.frame(site_years = nrow(x), sites = length(unique(x$site)), mismatch = mean(x$mismatch),
             mae_typical = mean(pmin(abs(x$d_typ[u]), cap)), mae_adjusted = mean(pmin(abs(x$d_adj[u]), cap)),
             dev_typical = if (sum(v)) mean(abs(x$g_dev[v])) else NA,
             dev_adjusted = if (sum(v)) mean(abs(x$g_dev[v] - x$a_dev[v])) else NA,
             r_adj = if (sum(v) >= 10 && stats::sd(x$a_dev[v]) > 0) stats::cor(x$g_dev[v], x$a_dev[v]) else NA,
             n_dev = sum(v))
}
rows <- list()
for (s in sets) {
  x <- E[E$set == s, ]
  for (rn in c("CONUS", sort(unique(x$region[!is.na(x$region)])))) {
    xr <- if (rn == "CONUS") x else x[x$region %in% as.integer(rn), ]
    if (length(unique(xr$site)) < 3) next
    rows[[length(rows) + 1]] <- data.frame(set = s, region = rn, score(xr))
  }
}
Sc <- do.call(rbind, rows)
utils::write.csv(Sc, file.path(out_dir, "ground_scores.csv"), row.names = FALSE)
for (i in seq_len(nrow(Sc))) with(Sc[i, ], log_msg(sprintf(
  "%-13s region %-5s %4d site-years %3d sites | mismatch %3.0f %% | vs ground: typical %5.1f adjusted %5.1f d | offset removed: %5.1f -> %5.1f d (n %d) | r %5.2f",
  set, region, site_years, sites, 100 * mismatch, mae_typical, mae_adjusted, dev_typical, dev_adjusted, n_dev, r_adj)))
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
