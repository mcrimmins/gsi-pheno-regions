#!/usr/bin/env Rscript
# Step 65 (fire activity study; P720): the region x week x year table everything is fitted on.
# Pilot regions (config fire_study: pilot_regions), weeks = fixed 7-day blocks from Jan 1 (week 52
# runs to Dec 31; fire_week()), years 2001-2024 (FPA-FOD window). One row per region-year-week.
#
# Outcomes
#   n_fires            FPA-FOD wildfires >= 10 ac (C+), debris / open burning removed, all causes
#   n_lightning, n_human, n_unknown   the same by cause (cause split used for 2001-2019 only;
#                      cause_split_ok flags those years and leaves out regions in cause_split_skip_regions)
#   acres              their final sizes (FPA-FOD, credited to the discovery week)
#   n_fires_debris     C+ fires including debris / open burning (sensitivity)
#   n_fires_top3       fires in the most herbaceous third of the region's cells (cover_sensitivity
#                      regions only; NA elsewhere)
#   px_days            VIIRS unique 375 m pixel-days, 63 `keep` detections (linked wildfire, not
#                      cropland), S-NPP plus NOAA-20 on whole CONUS days with no S-NPP data
#                      (viirs_outcome years only; NA in viirs_outcome exclude_regions)
#   px_days_nofill     S-NPP only;  fill_days  days in the week filled from NOAA-20;
#   outage_days        S-NPP outage days NOAA-20 could not fill (before 2018-04)
#   px_noncrop_all     all non-cropland detections, linked or not, with the same fill (the region 10
#                      check, viirs_unlinked_sensitivity_regions; written for every region)
# Baselines (64, fuel-share weighted region means, mean over the week's days)
#   dfm1000, ppet60, vpd30, kbdi, snow_share; snowoff_doy (that year's fuel-weighted snow-off date
#   over cells with snow, snow_regions only)
# Pixels: 375 m cells of the sinusoidal grid (fw_pixel_key).
# Phenology (satellite sos20 = green-up, eos50 = curing, cutb/metrics; typical dates from
#   rules/typical_dates.tif). Cells: strict herbaceous cells (rules `herb`), or in pheno_source
#   regions cells with fuel share >= min_fuel_share; weighted by fuel share in herb_weighted_regions.
#   cured_share        share of those cells not green at mid-week that year (green = between that
#                      year's green-up and curing, or inside the previous / next year's season)
#   cured_share_typ    the same with the typical dates (the climatological window)
#   *_nb               the same leaving out cells burned that year or the year before (MTBS any type,
#                      non-crop S-NPP detections from 2012): the burn check (Oct 19)
#   gu_week_typ, cu_week_typ  the region's typical green-up / curing week (weighted median of the
#                      typical dates; for the dates model's seasonal windows)
#   greenup_dep, curing_dep   that year's median departure (days) from the typical dates, cells with
#                      |departure| > 90 d left out as season mismatches (mismatch shares in the
#                      summary); constant within a region-year
# Flags: cause_split_ok, era2020 (era_check), viirs_ok.
# Outputs <run_dir>/fire_study/weekly/: fire_weekly.csv (+ .parquet), weekly_summary.csv
#   (per region: totals, cells, mismatch shares, fill weeks, and a first seasonal look: share of
#   fires in weeks when most fuel is typically cured vs the share of such weeks).
#
# Usage: Rscript scripts/65_fire_weekly.R [profile]
#   Needs 60 (fod_fires), 63 (link tables), 64 (dfm), cutb/metrics, rules/typical_dates.tif.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "fire.R"))
source(here::here("R", "fire_weekly.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("65_fire_weekly", profile)
notify_init(cfg)
fs <- cfg$fire_study; fc <- fs$firms
t0 <- Sys.time()
terra::terraOptions(progress = 0)
out_dir <- fire_dir(cfg, "weekly"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
regs <- as.integer(unlist(fs$pilot_regions))
yrs <- seq(as.integer(fs$fod$years$start), as.integer(fs$fod$years$end))
nw <- as.integer(fs$week$n_weeks %||% 52)
vy <- seq(as.integer(fs$viirs_outcome$years$start), as.integer(fs$viirs_outcome$years$end))
v_excl <- as.integer(unlist(fs$viirs_outcome$exclude_regions))
log_msg("=== 65_fire_weekly | profile: ", profile, " | regions ", paste(regs, collapse = ", "), " | years ", min(yrs), "-", max(yrs),
        " | VIIRS ", min(vy), "-", max(vy), " (not in ", paste(v_excl, collapse = ", "), ")")

W <- expand.grid(week = seq_len(nw), year = yrs, region = regs)[, c("region", "year", "week")]
W$week_start <- as.Date(sprintf("%d-01-01", W$year)) + (W$week - 1) * 7
key <- function(r, y, w) paste(r, y, w)
W$k <- key(W$region, W$year, W$week)
cnt <- function(k) { t <- table(k); v <- as.numeric(t[W$k]); v[is.na(v)] <- 0; v }

grid <- terra::rast(file.path(cfg$run_dir, "static", "grid_mask.tif"))
fuel <- fire_fuel_share(cfg)
TD <- terra::rast(file.path(cfg$run_dir, "rules", "typical_dates.tif"))
tdv <- terra::values(TD[[c("greenup_doy", "curing_doy", "herb", "region")]])
region <- as.integer(tdv[, "region"])

# ---- 1. FPA-FOD outcomes -------------------------------------------------------------------------------
FOD <- fire_read_table(fire_dir(cfg, "fod", "fod_fires"))
FOD$date <- as.Date(FOD$date)
FOD <- FOD[FOD$c_plus & FOD$region %in% regs & FOD$year %in% yrs, ]
FOD$k <- key(FOD$region, FOD$year, fire_week(FOD$date))
wf <- !FOD$debris
W$n_fires <- cnt(FOD$k[wf])
for (cz in c("lightning", "human", "unknown")) W[[paste0("n_", cz)]] <- cnt(FOD$k[wf & FOD$cause == cz])
ac <- tapply(FOD$size_ac[wf], FOD$k[wf], sum); W$acres <- as.numeric(ac[W$k]); W$acres[is.na(W$acres)] <- 0
W$n_fires_debris <- cnt(FOD$k)
cs <- fs$cover_sensitivity %||% list()
W$n_fires_top3 <- NA_real_
for (r in as.integer(unlist(cs$regions))) {
  thr <- stats::quantile(fuel[region %in% r], 1 - as.numeric(cs$top_fraction %||% 1 / 3), na.rm = TRUE)
  fz <- fuel[FOD$cell]
  i <- wf & FOD$region == r & !is.na(fz) & fz >= thr
  W$n_fires_top3[W$region == r] <- cnt(FOD$k[i])[W$region == r]
  log_msg(sprintf("region %d top third by fuel share (>= %.2f): %d of %d fires", r, thr, sum(i, na.rm = TRUE), sum(wf & FOD$region == r)))
}
cy <- fs$outcome$cause_split_years
W$cause_split_ok <- W$year >= as.integer(cy$start) & W$year <= as.integer(cy$end) &
  !(W$region %in% as.integer(unlist(fs$outcome$cause_split_skip_regions)))
W$era2020 <- W$year >= as.integer(fs$era_check$from_year %||% 2020)
log_msg(sprintf("FPA-FOD: %s C+ wildfires in %s region-weeks", format(sum(W$n_fires), big.mark = ","), format(sum(W$n_fires > 0), big.mark = ",")))

# ---- 2. VIIRS outcomes ------------------------------------------------------------------------------------
link_file <- function(sen, y) fire_table_file(fire_dir(cfg, "link", sen, sprintf("link_%s_%d", sen, y)))
cols <- c("lat", "lon", "date_utc", "date_solar", "region", "keep", "crop", "cell")
BURN_V <- list()                                    # cells with non-crop S-NPP detections, per year (burn check)
read_link <- function(sen, y) {
  f <- link_file(sen, y)
  if (is.na(f)) return(NULL)
  x <- if (grepl("parquet$", f)) as.data.frame(arrow::read_parquet(f, col_select = tidyselect::all_of(cols))) else readRDS(f)[, cols]
  x$date_utc <- as.Date(x$date_utc); x$date_solar <- as.Date(x$date_solar); x
}
W$px_days <- NA_real_; W$px_days_nofill <- NA_real_; W$px_noncrop_all <- NA_real_
W$fill_days <- NA_real_; W$outage_days <- NA_real_
W$viirs_ok <- W$year %in% vy & !(W$region %in% v_excl)
fill_from <- fs$viirs_fill$fill_from %||% fc$check
VS <- list()
for (y in vy) {
  P <- read_link(fc$primary, y)
  if (is.null(P)) { log_warn(fc$primary, " ", y, ": no link table, skip"); next }
  days <- seq(as.Date(sprintf("%d-01-01", y)), as.Date(sprintf("%d-12-31", y)), by = "day")
  first <- min(P$date_utc)
  out <- days[!(days %in% unique(P$date_utc)) & days >= first]
  Fd <- NULL; filled <- as.Date(character())
  if (length(out)) {
    N <- read_link(fill_from, y)
    if (!is.null(N)) { filled <- intersect(out, unique(N$date_utc)); filled <- as.Date(filled, origin = "1970-01-01"); Fd <- N[N$date_utc %in% filled, ] }
  }
  unfilled <- as.Date(setdiff(out, filled), origin = "1970-01-01")
  P <- P[P$region %in% regs, ]; if (!is.null(Fd)) Fd <- Fd[Fd$region %in% regs, ]
  BURN_V[[as.character(y)]] <- unique(P$cell[!P$crop])
  pk <- function(d) { d$k <- key(d$region, y, fire_week(d$date_solar)); d$px <- fw_pixel_key(d$lon, d$lat, d$date_solar); d }
  P <- pk(P); A <- if (!is.null(Fd) && nrow(Fd)) rbind(P, pk(Fd)) else P
  upx <- function(d) { d <- d[!duplicated(d$px), ]; cnt(d$k) }
  iy <- W$year == y
  W$px_days[iy] <- upx(A[A$keep, ])[iy]
  W$px_days_nofill[iy] <- upx(P[P$keep, ])[iy]
  W$px_noncrop_all[iy] <- upx(A[!A$crop, ])[iy]
  W$fill_days[iy] <- tabulate(fire_week(filled), nw)[W$week[iy]]
  W$outage_days[iy] <- tabulate(fire_week(unfilled), nw)[W$week[iy]]
  VS[[length(VS) + 1]] <- data.frame(year = y, outage_days = length(out), filled = length(filled), unfilled = length(unfilled))
  log_msg(sprintf("VIIRS %d: S-NPP outage days %d (filled from %s %d); kept pilot pixel-days %s (no fill %s)", y, length(out),
                  fill_from, length(filled), format(sum(W$px_days[iy]), big.mark = ","), format(sum(W$px_days_nofill[iy]), big.mark = ",")))
}
W$px_days[!W$viirs_ok] <- NA; W$px_days_nofill[!W$viirs_ok] <- NA     # px_noncrop_all kept in every region (checks)
if (length(VS)) utils::write.csv(do.call(rbind, VS), file.path(out_dir, "viirs_fill.csv"), row.names = FALSE)

# ---- 3. dryness baselines (64) -----------------------------------------------------------------------------
D <- fire_read_table(fire_dir(cfg, "dfm", "region_daily"))
D <- D[D$weight == "fuel" & D$region %in% regs, ]
D$date <- as.Date(D$date); D$year <- as.integer(format(D$date, "%Y")); D$week <- fire_week(D$date)
bm <- c("dfm1000", "ppet60", "vpd30", "kbdi", "snow_share")
A <- stats::aggregate(D[, bm], D[, c("region", "year", "week")], mean)
i <- match(W$k, key(A$region, A$year, A$week))
for (m in bm) W[[m]] <- A[[m]][i]
SO <- utils::read.csv(fire_dir(cfg, "dfm", "region_snowoff.csv"))
SO <- SO[SO$weight == "fuel", ]
W$snowoff_doy <- SO$snowoff_doy[match(paste(W$region, W$year), paste(SO$region, SO$year))]
W$snowoff_doy[!(W$region %in% as.integer(unlist(fs$snow_regions)))] <- NA
log_msg(sprintf("baselines: %d region-weeks without dfm1000", sum(is.na(W$dfm1000))))

# ---- 4. phenology -------------------------------------------------------------------------------------------
ps <- fs$pheno_source %||% list()
mfile <- function(Y) file.path(cfg$run_dir, "cutb", "metrics", sprintf("metrics_%d.tif", Y))
cells_of <- function(r) {
  if (r %in% as.integer(unlist(ps$regions))) which(region == r & !is.na(fuel) & fuel >= as.numeric(ps$min_fuel_share %||% 0.35))
  else which(region == r & tdv[, "herb"] %in% 1)
}
CELLS <- lapply(setNames(regs, regs), cells_of)
allc <- sort(unique(unlist(CELLS)))
MY <- list()
for (Y in seq(min(yrs) - 1, max(yrs) + 1)) {
  if (!file.exists(mfile(Y))) next
  r <- terra::rast(mfile(Y))
  MY[[as.character(Y)]] <- as.matrix(terra::extract(r[[c(fs$greenup_metric %||% "sos20", fs$curing_metric %||% "eos50")]], allc))
}
log_msg("satellite metrics years: ", paste(range(as.integer(names(MY))), collapse = "-"), " (", length(MY), ")")
# Burn check (decided Oct 19): cells that burned in year Y or Y-1 (any MTBS perimeter, any type; or a
# non-crop S-NPP detection, 2012 on) are left out of the *_nb phenology columns, since fire and
# post-fire regrowth change the satellite green-up and curing dates.
MC <- fire_read_table(fire_dir(cfg, "mtbs", "mtbs_cells")); MF <- fire_read_table(fire_dir(cfg, "mtbs", "mtbs_fires"))
MC$year <- as.integer(format(as.Date(MF$ig_date[match(MC$event_id, MF$event_id)]), "%Y"))
burned <- function(Y) unique(c(MC$cell[MC$year %in% c(Y, Y - 1)], unlist(BURN_V[as.character(c(Y, Y - 1))])))
mid <- fw_mid_doy(nw)
W$cured_share <- NA_real_; W$cured_share_typ <- NA_real_; W$greenup_dep <- NA_real_; W$curing_dep <- NA_real_
W$cured_share_nb <- NA_real_; W$cured_share_typ_nb <- NA_real_; W$greenup_dep_nb <- NA_real_; W$curing_dep_nb <- NA_real_
W$gu_week_typ <- NA_integer_; W$cu_week_typ <- NA_integer_
PH <- list()
for (r in regs) {
  cr <- CELLS[[as.character(r)]]; ix <- match(cr, allc)
  wt <- if (r %in% as.integer(unlist(fs$herb_weighted_regions))) fuel[cr] else rep(1, length(cr))
  tg <- tdv[cr, "greenup_doy"]; tc <- tdv[cr, "curing_doy"]
  wk <- function(doy) fire_week(as.Date("2001-01-01") + ((round(doy) - 1) %% 365))   # typical date -> week (wrapped)
  W$gu_week_typ[W$region == r] <- wk(fw_wmedian(tg, wt)); W$cu_week_typ[W$region == r] <- wk(fw_wmedian(tc, wt))
  for (Y in yrs) {
    get <- function(yy, j) { m <- MY[[as.character(yy)]]; if (is.null(m)) rep(NA_real_, length(cr)) else m[ix, j] }
    S <- cbind(get(Y - 1, 1), get(Y, 1), get(Y + 1, 1)); E <- cbind(get(Y - 1, 2), get(Y, 2), get(Y + 1, 2))
    lp <- as.numeric(as.Date(sprintf("%d-01-01", Y)) - as.Date(sprintf("%d-01-01", Y - 1)))
    lt <- as.numeric(as.Date(sprintf("%d-01-01", Y + 1)) - as.Date(sprintf("%d-01-01", Y)))
    okY <- !is.na(S[, 2]) & !is.na(E[, 2])
    G <- fw_green(mid, S, E, lp, lt)
    okT <- !is.na(tg) & !is.na(tc)
    GT <- fw_green(mid, cbind(tg, tg, tg), cbind(tc, tc, tc), lp, lt)
    iw <- W$region == r & W$year == Y
    W$cured_share[iw] <- vapply(seq_len(nw), function(j) fw_wmean(as.numeric(!G[okY, j]), wt[okY]), 0)[W$week[iw]]
    W$cured_share_typ[iw] <- vapply(seq_len(nw), function(j) fw_wmean(as.numeric(!GT[okT, j]), wt[okT]), 0)[W$week[iw]]
    dg <- S[, 2] - tg; dc <- E[, 2] - tc
    mg <- abs(dg) > 90; mc <- abs(dc) > 90
    W$greenup_dep[iw] <- fw_wmedian(ifelse(mg, NA, dg), wt)
    W$curing_dep[iw] <- fw_wmedian(ifelse(mc, NA, dc), wt)
    nb <- !(cr %in% burned(Y))
    W$cured_share_nb[iw] <- vapply(seq_len(nw), function(j) fw_wmean(as.numeric(!G[okY & nb, j]), wt[okY & nb]), 0)[W$week[iw]]
    W$cured_share_typ_nb[iw] <- vapply(seq_len(nw), function(j) fw_wmean(as.numeric(!GT[okT & nb, j]), wt[okT & nb]), 0)[W$week[iw]]
    W$greenup_dep_nb[iw] <- fw_wmedian(ifelse(mg | !nb, NA, dg), wt)
    W$curing_dep_nb[iw] <- fw_wmedian(ifelse(mc | !nb, NA, dc), wt)
    PH[[length(PH) + 1]] <- data.frame(region = r, year = Y, cells = length(cr), valid = sum(okY), burned = sum(!nb),
                                       greenup_mismatch = mean(mg, na.rm = TRUE), curing_mismatch = mean(mc, na.rm = TRUE))
  }
}
PH <- do.call(rbind, PH)
utils::write.csv(PH, file.path(out_dir, "pheno_cells.csv"), row.names = FALSE)

# ---- 5. write and summarise ----------------------------------------------------------------------------------
W$k <- NULL
invisible(fire_write_table(W, file.path(out_dir, "fire_weekly")))
Wc <- W; Wc$week_start <- format(Wc$week_start)
utils::write.csv(Wc, file.path(out_dir, "fire_weekly.csv"), row.names = FALSE)

SM <- do.call(rbind, lapply(regs, function(r) {
  w <- W[W$region == r, ]; ph <- PH[PH$region == r, ]
  cw <- w$cured_share_typ >= 0.5
  data.frame(region = r, fires = sum(w$n_fires), acres = sum(w$acres), weeks_with_fires = mean(w$n_fires > 0),
             pheno_cells = ph$cells[1], valid_share = mean(ph$valid / ph$cells),
             greenup_mismatch = mean(ph$greenup_mismatch), curing_mismatch = mean(ph$curing_mismatch),
             cured_weeks_share = mean(cw), fires_in_cured_weeks = sum(w$n_fires[cw]) / max(sum(w$n_fires), 1),
             px_days = sum(w$px_days, na.rm = TRUE), fill_weeks = sum(w$fill_days > 0, na.rm = TRUE),
             r_curing_dep_fires = {cd <- tapply(w$curing_dep, w$year, `[`, 1); nf <- tapply(w$n_fires, w$year, sum)
               ok <- !is.na(cd) & !is.na(nf); if (sum(ok) >= 3 && stats::sd(cd[ok]) > 0 && stats::sd(nf[ok]) > 0) stats::cor(cd[ok], nf[ok]) else NA_real_})
}))
utils::write.csv(SM, file.path(out_dir, "weekly_summary.csv"), row.names = FALSE)
for (i in seq_len(nrow(SM))) with(SM[i, ], log_msg(sprintf(
  "region %2d: %5d fires, %s ac; phenology cells %d (valid %.0f %%, mismatch green-up %.0f %% curing %.0f %%); typically-cured weeks %.0f %% of the year hold %.0f %% of fires",
  region, fires, format(round(acres), big.mark = ","), pheno_cells, 100 * valid_share, 100 * greenup_mismatch, 100 * curing_mismatch,
  100 * cured_weeks_share, 100 * fires_in_cured_weeks)))

msg <- sprintf("%d rows (%d regions x %d years x %d weeks), %.1f min", nrow(W), length(regs), length(yrs), nw,
               as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("65_fire_weekly finished", msg, priority = 3, tags = "white_check_mark")
