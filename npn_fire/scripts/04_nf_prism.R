#!/usr/bin/env Rscript
# Step 04 (NPN x fire study; P720, minutes; terra, arrow): the dryness baseline by Level II ecoregion
# and day. ppet60 = rain minus Hargreaves PET over the last 60 days (mm), the fire study's simple
# baseline (its decision 21), from the same PRISM daily grids.
# The ecoregion mean of a 60-day sum is the 60-day sum of the ecoregion's daily mean, so only daily
# weighted means of rain and PET are computed per cell block, then rolled:
#   PET (mm / d) = 0.0023 * 0.408 * Ra * (Tmean + 17.8) * sqrt(Tmax - Tmin)   (Hargreaves; Ra FAO-56 eq 21)
# Weights per cell (03's cell_weights.csv.gz, same grid): all = every land cell with an ecoregion;
# herb / shrub / timber = that cell's LANDFIRE fuel shares (config npn_fire: prism: weights), i.e. the
# dryness where each fuel is.
# Outputs <run_dir>/npn_fire/prism/:
#   eco_daily.parquet     eco_l2_code, weight, date, p, pet, ppet60 (mm), tmean (deg C), wsum (sum of weights)
#   eco_weekly.csv        eco_l2_code, weight, year, week, ppet60, tmean (means over the week's days)
# tmean (added 2 Nov 2026 for step 07's weather-only comparison): weighted mean of (Tmin + Tmax) / 2.
#   eco_climatology.csv   eco_l2_code, weight, week, median, p10, p90, mean across years
#   years/<Y>.rds         per-year daily sums (resume; delete or `force` to recompute)
# Usage: Rscript npn_fire/scripts/04_nf_prism.R [profile] [force]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && args[1] != "force") args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
force <- "force" %in% args
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("terra", quietly = TRUE)) stop("needs the terra package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("04_nf_prism", profile)
pc <- cfg$npn_fire$prism
pdir <- pc$prism_daily_dir %||% file.path(cfg$run_dir, "prism_daily")
yrs <- seq(pc$years$start, pc$years$end)
out <- nf_out(cfg, "prism"); ydir <- nf_out(cfg, "prism", "years")
nd <- as.integer(pc$ppet_days %||% 60)
t0 <- Sys.time()
pfile <- function(v, y) file.path(pdir, v, sprintf("prism_%s_%d.tif", v, y))
miss <- unlist(lapply(yrs, function(y) { f <- vapply(c("ppt", "tmin", "tmax"), pfile, "", y = y); f[!file.exists(f)] }))
if (length(miss)) stop("PRISM files missing: ", paste(utils::head(basename(miss), 5), collapse = ", "))
log_msg("=== 04_nf_prism | profile: ", profile, " | ", min(yrs), "-", max(yrs), " | ppet ", nd, " d | ", pdir)

# ---- weights per cell (03) ---------------------------------------------------------------------------
cwf <- file.path(cfg$out_dir, "fire", "cell_weights.csv.gz")
if (!file.exists(cwf)) stop("no ", cwf, ": run 03 first")
CW <- utils::read.csv(cwf, colClasses = c(eco_l2_code = "character"))
wdef <- pc$weights %||% list(herb = c("GR", "GS"), shrub = "SH", timber = c("TU", "TL"))
WM <- cbind(all = 1, sapply(wdef, function(g) rowSums(CW[, unlist(g), drop = FALSE])))
codes <- sort(unique(CW$eco_l2_code))
g0 <- terra::rast(pfile("ppt", yrs[1]))
nr <- terra::nrow(g0); nc <- terra::ncol(g0)
if (max(CW$cell) > terra::ncell(g0)) stop("cell numbers in cell_weights exceed the PRISM grid: grids differ")
cell_row <- (CW$cell - 1) %/% nc + 1
lat_all <- terra::yFromCell(g0, CW$cell)
log_msg(sprintf("%d cells in %d ecoregions | weights: %s", nrow(CW), length(codes), paste(colnames(WM), collapse = ", ")))

ra_matrix <- function(lat, doy) {           # extraterrestrial radiation (MJ m-2 d-1), FAO-56 eq 21
  phi <- lat * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365); dec <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmin(pmax(-tan(outer(phi, rep(1, length(doy)))) * tan(outer(rep(1, length(phi)), dec)), -1), 1))
  (24 * 60 / pi) * 0.0820 * matrix(dr, length(phi), length(doy), byrow = TRUE) *
    (ws * outer(sin(phi), sin(dec)) + outer(cos(phi), cos(dec)) * sin(ws))
}

# ---- one year: weighted daily sums of p and pet per ecoregion ------------------------------------------
do_year <- function(Y) {
  f <- file.path(ydir, sprintf("%d.rds", Y))
  if (!force && file.exists(f)) return(f)
  R <- lapply(c(ppt = "ppt", tmin = "tmin", tmax = "tmax"), function(v) terra::rast(pfile(v, Y)))
  dates <- as.Date(names(R$ppt))
  if (any(is.na(dates)) || !identical(names(R$ppt), names(R$tmin)) || !identical(names(R$ppt), names(R$tmax)))
    stop(Y, ": PRISM band names (dates) missing or not the same across variables")
  doy <- as.integer(format(dates, "%j"))
  key <- factor(CW$eco_l2_code, levels = codes)
  SP <- SE <- ST <- array(0, c(length(codes), ncol(WM), length(dates)), dimnames = list(codes, colnames(WM), format(dates)))
  WS <- matrix(0, length(codes), ncol(WM), dimnames = list(codes, colnames(WM)))
  br <- as.integer(pc$block_rows %||% 40)
  for (r in R) terra::readStart(r)
  on.exit(for (r in R) terra::readStop(r), add = TRUE)
  for (row in seq(1, nr, by = br)) {
    n <- min(br, nr - row + 1)
    i <- which(cell_row >= row & cell_row < row + n)
    if (!length(i)) next
    off <- CW$cell[i] - (row - 1) * nc
    V <- lapply(R, function(r) terra::readValues(r, row = row, nrows = n, mat = TRUE)[off, , drop = FALSE])
    ok <- stats::complete.cases(V$ppt, V$tmin, V$tmax)
    if (!any(ok)) next
    P <- V$ppt[ok, , drop = FALSE]; TN <- V$tmin[ok, , drop = FALSE]; TX <- V$tmax[ok, , drop = FALSE]
    PET <- pmax(0.0023 * 0.408 * ra_matrix(lat_all[i][ok], doy) * ((TN + TX) / 2 + 17.8) * sqrt(pmax(TX - TN, 0)), 0)
    k <- key[i][ok]; w <- WM[i[ok], , drop = FALSE]
    for (j in seq_len(ncol(WM))) {
      s1 <- rowsum(P * w[, j], k, reorder = TRUE); s2 <- rowsum(PET * w[, j], k, reorder = TRUE)
      s3 <- rowsum(((TN + TX) / 2) * w[, j], k, reorder = TRUE)
      ST[rownames(s3), j, ] <- ST[rownames(s3), j, ] + s3
      SP[rownames(s1), j, ] <- SP[rownames(s1), j, ] + s1
      SE[rownames(s2), j, ] <- SE[rownames(s2), j, ] + s2
      WS[rownames(s1), j] <- WS[rownames(s1), j] + rowsum(w[, j], k, reorder = TRUE)[, 1]
    }
  }
  tmp <- paste0(f, ".tmp"); saveRDS(list(year = Y, dates = dates, SP = SP, SE = SE, ST = ST, WS = WS), tmp); file.rename(tmp, f)
  f
}
nw <- if (.Platform$OS.type == "unix") as.integer(pc$workers %||% 8) else 1L
run1 <- function(Y) tryCatch(do_year(Y), error = function(e) paste("ERROR", Y, conditionMessage(e)))
res <- if (nw > 1) parallel::mclapply(yrs, run1, mc.cores = min(nw, length(yrs)), mc.preschedule = FALSE) else lapply(yrs, run1)
err <- unlist(res[grepl("^ERROR", unlist(res))])
if (length(err)) { for (e in err) log_err(e); stop(length(err), " year(s) failed; rerun to retry") }
log_msg(sprintf("yearly sums: %d years (%.1f min)", length(yrs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- combine, roll, write ---------------------------------------------------------------------------
Yl <- lapply(yrs, function(Y) readRDS(file.path(ydir, sprintf("%d.rds", Y))))
dates <- do.call(c, lapply(Yl, `[[`, "dates"))
if (any(diff(dates) != 1)) stop("PRISM dates are not contiguous")
WS <- Yl[[1]]$WS
if (any(vapply(Yl, function(z) is.null(z$ST), TRUE))) stop("yearly sums without tmean (made before 2 Nov): rerun with force")
bind_days <- function(L) { d <- dim(L[[1]]); array(unlist(L), c(d[1], d[2], sum(vapply(L, function(z) dim(z)[3], 0L)))) }   # days are the last dim
SP <- bind_days(lapply(Yl, `[[`, "SP")); SE <- bind_days(lapply(Yl, `[[`, "SE")); ST <- bind_days(lapply(Yl, `[[`, "ST"))
roll <- function(x, n) { cs <- cumsum(c(0, x)); out <- rep(NA_real_, length(x)); i <- seq(n, length(x)); out[i] <- cs[i + 1] - cs[i + 1 - n]; out }
D <- list()
for (a in seq_along(codes)) for (j in seq_len(ncol(WM))) {
  if (WS[a, j] <= 0) next
  p <- SP[a, j, ] / WS[a, j]; e <- SE[a, j, ] / WS[a, j]; tm <- ST[a, j, ] / WS[a, j]
  D[[length(D) + 1]] <- data.frame(eco_l2_code = codes[a], weight = colnames(WM)[j], date = dates, p = round(p, 3), pet = round(e, 3),
                                   ppet60 = round(roll(p - e, nd), 2), tmean = round(tm, 2), wsum = round(WS[a, j], 2), stringsAsFactors = FALSE)
}
D <- do.call(rbind, D)
invisible(nf_write_table(D, file.path(out, "eco_daily")))
D$year <- as.integer(format(D$date, "%Y")); D$week <- nf_week(D$date)
key <- paste(D$eco_l2_code, D$weight, D$year, D$week)
WK <- D[!duplicated(key), c("eco_l2_code", "weight", "year", "week")]
WK$ppet60 <- round(as.numeric(tapply(D$ppet60, key, mean)[paste(WK$eco_l2_code, WK$weight, WK$year, WK$week)]), 2)
WK$tmean <- round(as.numeric(tapply(D$tmean, key, mean)[paste(WK$eco_l2_code, WK$weight, WK$year, WK$week)]), 2)
utils::write.csv(WK, file.path(out, "eco_weekly.csv"), row.names = FALSE)
k2 <- paste(WK$eco_l2_code, WK$weight, WK$week)
CL <- WK[!duplicated(k2), c("eco_l2_code", "weight", "week")]
q <- function(fun) as.numeric(tapply(WK$ppet60, k2, fun)[paste(CL$eco_l2_code, CL$weight, CL$week)])
CL$median <- q(function(v) stats::median(v, na.rm = TRUE)); CL$p10 <- q(function(v) stats::quantile(v, 0.1, na.rm = TRUE))
CL$p90 <- q(function(v) stats::quantile(v, 0.9, na.rm = TRUE)); CL$mean <- q(function(v) mean(v, na.rm = TRUE))
CL[c("median", "p10", "p90", "mean")] <- round(CL[c("median", "p10", "p90", "mean")], 1)
utils::write.csv(CL, file.path(out, "eco_climatology.csv"), row.names = FALSE)
s <- CL[CL$weight == "all", ]
for (cd in unique(s$eco_l2_code)) { z <- s[s$eco_l2_code == cd, ]
  log_msg(sprintf("  %-5s ppet60 weekly median: driest %6.0f mm (week %2d), wettest %6.0f mm (week %2d)", cd,
                  min(z$median, na.rm = TRUE), z$week[which.min(z$median)], max(z$median, na.rm = TRUE), z$week[which.max(z$median)])) }
log_msg(sprintf("=== done | %d ecoregion x weight series, %s to %s | %.1f min | %s | log: %s", nrow(unique(D[, c("eco_l2_code", "weight")])),
                min(dates), max(dates), as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
