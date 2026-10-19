#!/usr/bin/env Rscript
# Step 64 (fire activity study; P720): dead-fuel dryness baselines from daily PRISM (decision 6,
# Oct 18 addition), computed per grid cell for every day of the PRISM years and aggregated to
# region x day. Equations and sources in R/fire_dfm.R; parameters in config fire_study: dfm.
#   dfm1000   NFDRS 1978 1000-h dead fuel moisture (%): daily boundary from day / night EMC
#             (RH from PRISM VPDmax, Tmax, Tmin) and rain duration estimated from the daily
#             amount (NFDRS climate-class rate from each cell's aridity), 7-day mean boundary,
#             7-day lagged update
#   ppet60    rain minus Hargreaves PET over the last 60 days (mm): the simple, readable baseline
#   vpd30     mean PRISM VPDmax over the last 30 days (kPa)
#   kbdi      Keetch-Byram drought index (as 12, mean annual precipitation from the PRISM years)
#   swe, snow_share   degree-day snowpack (mm) and share of the region's weight with snow
#                     (> swe_min); plus a yearly snow-off date per cell (decision 8)
# The series starts with a spin-up pass over the first year (its Dec 31 states start Jan 1 of
# that year), so the first year needs no trimming.
# Outputs <run_dir>/fire_study/dfm/:
#   region_daily.csv (and .parquet)  region x weight (all = equal land cells, fuel = fuel share)
#                                    x date: the metrics above, weight sum
#   region_snowoff.csv     region x weight x year: mean snow-off DOY over cells with snow, share
#                          of weight with snow
#   rain_class.csv         region: share of weight with the dry-class rain rate, median P/PET
#   dfm_climatology.csv    region x weight x week (fire_week): across-year median, p10, p90
#   dfm_drivers.csv        what the 1000-h value follows: weekly anomalies (from each week's
#                          mean) of dfm1000 vs ppet60, vpd30, kbdi (r) and R2 of dfm1000 ~
#                          ppet60 + vpd30, per region and season
#   snowoff_cells.tif (one band per year), rain_rate.tif (rate, ai_ann, ai_warm)
#   chunks_<profile>/      per row-chunk results (resume; delete to recompute, or use force)
#
# Usage: Rscript scripts/64_fire_dfm.R [profile] [force]
#   Forked workers by row chunk on Linux (config workers_by_step: fire_dfm), sequential on Windows.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
force <- "force" %in% args

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "gsi.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "fire_dfm.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("64_fire_dfm", profile)
notify_init(cfg)
fs <- cfg$fire_study; dc <- fs$dfm %||% list(); sc <- dc$snow %||% list()
sw <- step_workers(cfg, "fire_dfm")
out_dir <- fire_dir(cfg, "dfm"); chunk_dir <- file.path(out_dir, paste0("chunks_", profile))
if (force && dir.exists(chunk_dir)) unlink(chunk_dir, recursive = TRUE)
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(progress = 0)
t0 <- Sys.time()

vars <- c("ppt", "tmin", "tmax", "vpdmax")
pfile <- function(v, y) prism_out_file(cfg, v, y)
yrs <- cfg$years
have <- vapply(yrs, function(y) all(file.exists(vapply(vars, pfile, "", y = y))), TRUE)
if (!have[1]) stop("PRISM files missing for the first year ", yrs[1])
yrs <- yrs[seq_len(if (all(have)) length(have) else which(!have)[1] - 1)]    # contiguous from the start
dates <- do.call(c, lapply(yrs, function(y) as.Date(names(terra::rast(pfile("ppt", y))))))
if (any(diff(dates) != 1)) stop("PRISM daily series is not contiguous")
n_spin <- sum(format(dates, "%Y") == as.character(yrs[1]))
fr <- as.numeric(dc$fr1000 %||% (1 - 0.82 * exp(-168 / 1000)))
log_msg("=== 64_fire_dfm | profile: ", profile, " | ", min(dates), " to ", max(dates), " (", length(yrs), " years) | workers ", sw$n,
        sprintf(" | 1000-h fr %.4f, rain rate %g / %g in/h, ppet %d d, vpd %d d", fr, dc$rain_rate_in_h$dry %||% 0.25,
                dc$rain_rate_in_h$wet %||% 0.05, as.integer(dc$ppet_days %||% 60), as.integer(dc$vpd_days %||% 30)))

grid_file <- file.path(static_dir(cfg), "grid_mask.tif")
grid <- terra::rast(grid_file)
RG <- fire_region_values(cfg, grid, fs$partition %||% "A2_groups", as.integer(fs$k %||% 13))
region <- as.integer(RG$v)
fuel <- fire_fuel_share(cfg); fuel[is.na(fuel)] <- 0
log_msg("regions from ", RG$src, "; fuel share from static/nlcd_shares.tif")
metrics <- c("dfm1000", "ppet60", "vpd30", "kbdi", "swe", "snow_share")
regs <- sort(unique(region[!is.na(region)]))

nr_all <- terra::nrow(grid); nc <- terra::ncol(grid)
chunk_rows <- as.integer(dc$chunk_rows %||% 3)
rows <- seq(1, nr_all, by = chunk_rows)
chunk_file <- function(row) file.path(chunk_dir, sprintf("chunk_%04d.rds", row))

do_chunk <- function(row) {
  f <- chunk_file(row)
  if (file.exists(f)) return(list(ok = TRUE, row = row, skipped = TRUE))
  tryCatch({
    terra::terraOptions(memfrac = sw$memfrac, progress = 0)
    nr <- min(chunk_rows, nr_all - row + 1)
    cells <- (row - 1) * nc + seq_len(nr * nc)
    g <- terra::rast(grid_file)
    land <- !is.na(read_rows(g, row, nr, 1)[, 1]) & !is.na(region[cells])
    res <- list(row = row, cells = integer(0))
    if (any(land)) {
      rd <- function(v) do.call(cbind, lapply(yrs, function(y) {
        r <- terra::rast(pfile(v, y)); read_rows(r, row, nr, seq_len(terra::nlyr(r)))
      }))[land, , drop = FALSE]
      P <- rd("ppt"); TN <- rd("tmin"); TX <- rd("tmax"); VP <- rd("vpdmax")
      ok <- rowSums(is.na(P) | is.na(TN) | is.na(TX) | is.na(VP)) == 0
      cl <- cells[land][ok]
      P <- P[ok, , drop = FALSE]; TN <- TN[ok, , drop = FALSE]; TX <- TX[ok, , drop = FALSE]; VP <- VP[ok, , drop = FALSE]
      if (length(cl)) {
        lat <- terra::yFromCell(g, cl)
        doy <- as.integer(format(dates, "%j"))
        PET <- hargreaves(TN, TX, ra_matrix(lat, doy))
        rr <- fdm_rain_rate(P, PET, dates, dc$dry_class$ai_ann_lt %||% -Inf, dc$dry_class$ai_warm_lt %||% 0.5,
                            dc$rain_rate_in_h$dry %||% 0.25, dc$rain_rate_in_h$wet %||% 0.05)
        map_in <- rowSums(P) / length(yrs) / 25.4
        sp <- seq_len(n_spin)                                             # spin-up: first year again
        P <- cbind(P[, sp, drop = FALSE], P); TN <- cbind(TN[, sp, drop = FALSE], TN)
        TX <- cbind(TX[, sp, drop = FALSE], TX); VP <- cbind(VP[, sp, drop = FALSE], VP)
        PET <- cbind(PET[, sp, drop = FALSE], PET); doy2 <- c(doy[sp], doy)
        photo_h <- m_photoperiod(lat, doy2, dc$photo_method %||% "geometric") / 3600
        keep <- -sp
        D <- fdm_boundary(P, TN, TX, VP, photo_h, rr$rate, as.numeric(dc$max_rain_hours %||% 8)); rm(photo_h)
        X <- list(dfm1000 = fdm_1000h(D, fr)[, keep, drop = FALSE]); rm(D)
        X$ppet60 <- m_roll_sum(P - PET, as.integer(dc$ppet_days %||% 60))[, keep, drop = FALSE]
        X$vpd30 <- m_roll_mean(VP / 10, as.integer(dc$vpd_days %||% 30))[, keep, drop = FALSE]
        X$kbdi <- m_kbdi(P, TX, map_in, as.numeric(dc$kbdi_interception_in %||% 0.2))[, keep, drop = FALSE]
        X$swe <- fdm_snow(P, TN, TX, as.numeric(sc$t_snow_c %||% 1), as.numeric(sc$t_melt_c %||% 0),
                          as.numeric(sc$ddf_mm_c %||% 3))[, keep, drop = FALSE]
        swe_min <- as.numeric(sc$swe_min_mm %||% 5)
        X$snow_share <- (X$swe > swe_min) + 0
        rm(P, TN, TX, VP, PET)
        so <- vapply(yrs, function(y) fdm_snowoff(X$swe, dates, y, swe_min, sc$snowoff_by %||% "07-31"), integer(length(cl)))
        so <- matrix(so, length(cl))
        rg <- region[cl]
        W <- list(all = rep(1, length(cl)), fuel = fuel[cl])
        sums <- lapply(W, function(w) lapply(X, function(M) rowsum(M * w, rg, reorder = TRUE)))
        wsum <- lapply(W, function(w) rowsum(w, rg, reorder = TRUE)[, 1])
        res <- list(row = row, cells = cl, sums = sums, wsum = wsum, snowoff = so,
                    rate = rr$rate, ai_ann = rr$ai_ann, ai_warm = rr$ai_warm)
      }
    }
    tmp <- paste0(f, ".tmp"); saveRDS(res, tmp); file.rename(tmp, f)
    list(ok = TRUE, row = row, n = length(res$cells))
  }, error = function(e) list(ok = FALSE, row = row, msg = conditionMessage(e)))
}

todo <- rows[!file.exists(vapply(rows, chunk_file, ""))]
log_msg(sprintf("%d row chunks of %d rows, %d to compute", length(rows), chunk_rows, length(todo)))
if (length(todo)) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  every <- max(1L, length(todo) %/% 20L)
  batches <- split(todo, ceiling(seq_along(todo) / (sw$n * every)))
  done <- 0L; fails <- list()
  for (b in batches) {
    rr <- if (.Platform$OS.type == "unix" && sw$n > 1) parallel::mclapply(b, do_chunk, mc.cores = sw$n, mc.preschedule = FALSE) else lapply(b, do_chunk)
    for (r in rr) if (!isTRUE(r$ok)) fails[[length(fails) + 1]] <- r
    done <- done + length(b)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_msg(sprintf("  chunks %d / %d (%.0f %%), %.1f min, ETA %.0f min", done, length(todo), 100 * done / length(todo),
                    el, el / done * (length(todo) - done)))
  }
  if (length(fails)) {
    for (r in fails) log_err("chunk row ", r$row, ": ", r$msg %||% "worker failed")
    stop(length(fails), " chunk(s) failed; rerun to retry them")
  }
}

# ---- combine ------------------------------------------------------------------------------------------
log_msg("combining chunks")
nd <- length(dates); nreg <- length(regs)
S <- lapply(c(all = "all", fuel = "fuel"), function(w) lapply(setNames(metrics, metrics), function(m) matrix(0, nreg, nd)))
WS <- lapply(c(all = "all", fuel = "fuel"), function(w) numeric(nreg))
SO <- matrix(NA_integer_, terra::ncell(grid), length(yrs))
RR <- matrix(NA_real_, terra::ncell(grid), 3)
for (row in rows) {
  x <- readRDS(chunk_file(row))
  if (!length(x$cells)) next
  for (w in names(S)) {
    i <- match(as.integer(names(x$wsum[[w]])), regs)
    WS[[w]][i] <- WS[[w]][i] + x$wsum[[w]]
    for (m in metrics) S[[w]][[m]][i, ] <- S[[w]][[m]][i, ] + x$sums[[w]][[m]]
  }
  SO[x$cells, ] <- x$snowoff
  RR[x$cells, ] <- cbind(x$rate, x$ai_ann, x$ai_warm)
}
DL <- do.call(rbind, lapply(names(S), function(w) {
  do.call(rbind, lapply(seq_len(nreg), function(k) {
    d <- data.frame(region = regs[k], weight = w, date = dates, wsum = WS[[w]][k])
    for (m in metrics) d[[m]] <- if (WS[[w]][k] > 0) S[[w]][[m]][k, ] / WS[[w]][k] else NA_real_
    d
  }))
}))
DL$pilot <- DL$region %in% unlist(fs$pilot_regions)
invisible(fire_write_table(DL, file.path(out_dir, "region_daily")))
DLc <- DL; DLc$date <- format(DLc$date)
for (m in metrics) DLc[[m]] <- signif(DLc[[m]], 5)
utils::write.csv(DLc, file.path(out_dir, "region_daily.csv"), row.names = FALSE); rm(DLc)

write_feature_matrix(SO, grid, paste0("snowoff_", yrs), file.path(out_dir, "snowoff_cells.tif"))
write_feature_matrix(RR, grid, c("rain_rate_in_h", "ai_ann", "ai_warm"), file.path(out_dir, "rain_rate.tif"))

# snow-off by region and year
ok <- !is.na(region)
SOR <- do.call(rbind, lapply(c("all", "fuel"), function(w) {
  wt <- if (w == "all") as.numeric(ok) else fuel * ok
  do.call(rbind, lapply(seq_along(yrs), function(j) {
    s <- SO[, j]; has <- ok & !is.na(s)
    num <- tapply(s[has] * wt[has], region[has], sum); den <- tapply(wt[has], region[has], sum)
    tot <- tapply(wt[ok], region[ok], sum)
    data.frame(region = regs, weight = w, year = yrs[j],
               snowoff_doy = as.numeric(num[as.character(regs)] / den[as.character(regs)]),
               snow_weight_share = as.numeric(den[as.character(regs)] / tot[as.character(regs)]))
  }))
}))
SOR$snow_weight_share[is.na(SOR$snow_weight_share)] <- 0
utils::write.csv(SOR, file.path(out_dir, "region_snowoff.csv"), row.names = FALSE)

RC <- do.call(rbind, lapply(regs, function(k) {
  i <- which(region == k & !is.na(RR[, 1]))
  data.frame(region = k, cells = length(i),
             dry_rate_share_all = mean(RR[i, 1] == (dc$rain_rate_in_h$dry %||% 0.25)),
             dry_rate_share_fuel = sum(fuel[i] * (RR[i, 1] == (dc$rain_rate_in_h$dry %||% 0.25))) / max(sum(fuel[i]), 1e-9),
             ai_ann_median = stats::median(RR[i, 2]), ai_warm_median = stats::median(RR[i, 3]))
}))
utils::write.csv(RC, file.path(out_dir, "rain_class.csv"), row.names = FALSE)

# ---- weekly climatology and what the 1000-h value follows -------------------------------------------
DL$year <- as.integer(format(DL$date, "%Y")); DL$week <- fire_week(DL$date)
WK <- stats::aggregate(DL[, metrics], DL[, c("region", "weight", "year", "week")], mean)
CL <- do.call(rbind, lapply(metrics, function(m) {
  a <- stats::aggregate(WK[[m]], WK[, c("region", "weight", "week")], function(v) stats::quantile(v, c(0.1, 0.5, 0.9), na.rm = TRUE))
  data.frame(a[, 1:3], metric = m, p10 = a$x[, 1], median = a$x[, 2], p90 = a$x[, 3])
}))
utils::write.csv(CL, file.path(out_dir, "dfm_climatology.csv"), row.names = FALSE)

season_of <- function(week) {                    # season of the week's middle day
  m <- as.integer(format(as.Date("2001-01-01") + (week - 1) * 7 + 3, "%m"))
  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")[m]
}
WA <- WK[WK$weight == "fuel", ]
for (m in c("dfm1000", "ppet60", "vpd30", "kbdi")) WA[[paste0(m, "_a")]] <- WA[[m]] - stats::ave(WA[[m]], WA$region, WA$week, FUN = function(v) mean(v, na.rm = TRUE))
WA$season <- season_of(WA$week)
DR <- do.call(rbind, lapply(regs, function(k) do.call(rbind, lapply(c("all", "DJF", "MAM", "JJA", "SON"), function(s) {
  d <- WA[WA$region == k & (s == "all" | WA$season == s), ]
  d <- d[stats::complete.cases(d[, c("dfm1000_a", "ppet60_a", "vpd30_a", "kbdi_a")]), ]
  if (nrow(d) < 30 || stats::sd(d$dfm1000_a) == 0) return(NULL)
  r <- function(v) suppressWarnings(stats::cor(d$dfm1000_a, d[[v]]))
  fit <- stats::lm(dfm1000_a ~ ppet60_a + vpd30_a, d)
  data.frame(region = k, season = s, n_weeks = nrow(d), sd_dfm1000 = stats::sd(d$dfm1000_a),
             r_ppet60 = r("ppet60_a"), r_vpd30 = r("vpd30_a"), r_kbdi = r("kbdi_a"),
             r2_ppet60 = r("ppet60_a")^2, r2_ppet60_vpd30 = summary(fit)$r.squared)
}))))
utils::write.csv(DR, file.path(out_dir, "dfm_drivers.csv"), row.names = FALSE)

# ---- log a readable summary for the pilot regions -------------------------------------------------------
for (k in intersect(unlist(fs$pilot_regions), regs)) {
  c1 <- CL[CL$region == k & CL$weight == "fuel" & CL$metric == "dfm1000", ]
  c2 <- CL[CL$region == k & CL$weight == "fuel" & CL$metric == "ppet60", ]
  d <- DR[DR$region == k & DR$season == "all", ]
  rc <- RC[RC$region == k, ]
  log_msg(sprintf("region %2d: 1000-h median %4.1f %% (driest week %2d) to %4.1f %% (wettest week %2d); ppet60 driest %5.0f mm | anomalies r ppet60 %.2f, vpd30 %.2f, kbdi %.2f, R2 ppet60+vpd30 %.2f | dry rain rate on %.0f %% of fuel",
                  k, min(c1$median), c1$week[which.min(c1$median)], max(c1$median), c1$week[which.max(c1$median)], min(c2$median),
                  d$r_ppet60, d$r_vpd30, d$r_kbdi, d$r2_ppet60_vpd30, 100 * rc$dry_rate_share_fuel))
}
so <- SOR[SOR$weight == "fuel" & SOR$region %in% unlist(fs$snow_regions), ]
for (k in unique(so$region)) {
  s <- so[so$region == k, ]
  log_msg(sprintf("region %2d snow-off (fuel-weighted, cells with snow): median DOY %.0f, range %.0f-%.0f; snow on %.0f %% of fuel weight (median year)",
                  k, stats::median(s$snowoff_doy, na.rm = TRUE), min(s$snowoff_doy, na.rm = TRUE), max(s$snowoff_doy, na.rm = TRUE),
                  100 * stats::median(s$snow_weight_share)))
}

msg <- sprintf("%d regions x %d days, %.1f min", nreg, nd, as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("64_fire_dfm finished", msg, priority = 3, tags = "white_check_mark")
