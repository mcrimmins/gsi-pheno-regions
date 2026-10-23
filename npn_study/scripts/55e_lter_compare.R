#!/usr/bin/env Rscript
# Step 55e (NPN grass curing study; laptop): independent cross-checks of the NPN results with the
# Jornada and Sevilleta LTER phenology records (55c tables; weather from 55d if present).
# The LTER records give, per site x species x month, the share of plants that are green (and at
# Sevilleta, with new green leaves). Dates below are where that share crosses 50 %, interpolated
# between monthly visits: a population date ("half the plants green"), not the NPN plant date
# ("a plant half green"); part 7 calibrates one against the other on the NPN record itself.
#   1 climatology      monthly share green / new by site x grass type and by species (mean, SD over years)
#   2 season dates     per site x species x year: up1 (first rise above 50 %), down1 (first fall
#                      after it), up2 (next rise), down_last (last fall); for green and new
#   3 C3 vs C4         Sevilleta: months where both types are below 50 % new / green, and per
#                      site-year the gap C4 up1 - C3 down1 (as 58b: positive = neither type active)
#   4 vs NPN           (a) same species: LTER monthly share green vs NPN share of readings >= 5 %
#                      green at NPN sites within 400 km; (b) same place, same years: Jornada LTER vs
#                      NEON Jornada (sites within 30 km), monthly values and anomalies
#   5 2016-17          share green in Jul-Sep and share of series never above 50 %, by year
#   6 one climate      Jornada: the same species at several of the 15 sites (one weather):
#                      spread of typical dates across sites vs from year to year
#   7 calibration      NPN record: (a) dates from one reading a month vs all readings (NEON plants
#                      with >= 20 readings); (b) population date (half the plants >= 5 % / >= 50 %
#                      green) vs the median plant 50 % date, per site-species-year
#   8 weather          (needs 55d) year-to-year anomalies of the season dates vs the ladder's
#                      site-year weather metrics: Spearman r, and leave-one-year-out skill of the
#                      best single driver (chosen within the training years) vs the series mean;
#                      and spring vs monsoon green-up (first above 50 % after June 1) against
#                      Dec-May rain
# Outputs <run_dir>/npn_study/lter/compare/*.csv; a summary in the log. Page figures: 90 npn_lter_*.
#
# Usage: Rscript npn_study/scripts/55e_lter_compare.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("55e_lter_compare", profile)
ldir <- file.path(cfg$run_dir, "npn_study", "lter"); odir <- file.path(ldir, "compare"); dir.create(odir, showWarnings = FALSE)
wcsv <- function(x, f) utils::write.csv(x, file.path(odir, f), row.names = FALSE)
log_msg("=== 55e_lter_compare | ", odir)
min_n <- 10                                                     # plants per site-species-month
lab <- function(d) ifelse(is.na(d), "-", format(as.Date("2001-01-01") + round(d) - 1, "%b %d"))

M <- utils::read.csv(file.path(ldir, "lter_grass_monthly.csv"), stringsAsFactors = FALSE)
M <- M[M$pathway %in% c("C3", "C4"), ]
Mall <- M                                                       # for the grass-type pools: plants pooled before the minimum
M <- M[M$n >= min_n, ]
NM <- utils::read.csv(file.path(ldir, "npn_monthly.csv"), stringsAsFactors = FALSE)
ST <- utils::read.csv(file.path(ldir, "lter_sites.csv"), stringsAsFactors = FALSE)
log_msg(sprintf("LTER site-species-months with >= %d plants: %d (JRN %d, SEV %d)", min_n, nrow(M), sum(M$source == "JRN"), sum(M$source == "SEV")))

# ---- helpers ----------------------------------------------------------------------------------
xing <- function(d, y, th = 0.5) {                            # all 50 % crossings of a monthly series
  k <- order(d); d <- d[k]; y <- y[k]; k <- !is.na(y) & !is.na(d); d <- d[k]; y <- y[k]
  if (length(y) < 2) return(data.frame(dir = character(), doy = numeric()))
  a <- y[-length(y)] - th; b <- y[-1] - th
  f <- function(i) d[i] + (th - y[i]) / (y[i + 1] - y[i]) * (d[i + 1] - d[i])
  up <- which(a < 0 & b >= 0); dn <- which(a >= 0 & b < 0)
  rbind(data.frame(dir = rep("up", length(up)), doy = vapply(up, f, 0)), data.frame(dir = rep("down", length(dn)), doy = vapply(dn, f, 0)))
}
season <- function(d, y, th = 0.5) {
  x <- xing(d, y, th); x <- x[order(x$doy), ]
  up1 <- x$doy[x$dir == "up"][1]
  dn1 <- x$doy[x$dir == "down" & x$doy > (if (is.na(up1)) -Inf else up1)][1]
  up2 <- if (is.na(dn1)) NA else x$doy[x$dir == "up" & x$doy > dn1][1]
  dnl <- utils::tail(x$doy[x$dir == "down"], 1)
  c(up1 = up1, down1 = dn1, up2 = up2, down_last = if (length(dnl)) dnl else NA)
}
wmean <- function(v, w) sum(v * w) / sum(w)
agg <- function(f, d, FUN = mean) if (nrow(d)) stats::aggregate(f, d, FUN) else NULL          # aggregate that allows no rows
dist_km <- function(la1, lo1, la2, lo2) 111.2 * sqrt((la1 - la2)^2 + ((lo1 - lo2) * cos(la1 * pi / 180))^2)
lter_xy <- stats::aggregate(cbind(lat, lon) ~ source, ST, mean)

# ---- 1 climatology ----------------------------------------------------------------------------
vars <- c("green", "new")
# grass-type pools per site, and for Sevilleta also over all its sites (C3 grasses are sparse:
# usually fewer than 10 plants of a species in a month); the minimum applies to the pooled plants
Mp <- rbind(Mall, transform(Mall[Mall$source == "SEV", ], site = "all"))
P <- do.call(rbind, lapply(split(Mp, paste(Mp$source, Mp$site, Mp$pathway, Mp$year, Mp$month)), function(d)
  data.frame(source = d$source[1], site = d$site[1], pathway = d$pathway[1], year = d$year[1], month = d$month[1],
             doy = stats::median(d$doy), n = sum(d$n), species = nrow(d), green = wmean(d$green, d$n),
             new = if (all(is.na(d$new))) NA else wmean(d$new, d$n))))
P <- P[P$n >= min_n, ]
clim <- function(X, by) do.call(rbind, lapply(split(X, X[, by], drop = TRUE), function(d) {
  o <- d[1, by, drop = FALSE]; o$years <- length(unique(d$year)); o$doy <- stats::median(d$doy)
  for (v in vars) { o[[paste0(v, "_mean")]] <- mean(d[[v]], na.rm = TRUE); o[[paste0(v, "_sd")]] <- stats::sd(d[[v]], na.rm = TRUE) }
  o }))
CL_p <- clim(P, c("source", "site", "pathway", "month")); CL_p$level <- "pathway"
CL_s <- clim(M, c("source", "site", "species", "name", "pathway", "month")); CL_s$level <- "species"
CL_s <- CL_s[CL_s$years >= 5, ]
wcsv(P, "pathway_monthly.csv"); wcsv(CL_p, "climatology_pathway.csv"); wcsv(CL_s, "climatology_species.csv")
for (s in c("SEV", "JRN")) for (pw in c("C3", "C4")) {
  cx <- CL_p[CL_p$source == s & CL_p$pathway == pw, ]; if (!nrow(cx)) next
  x <- stats::aggregate(green_mean ~ month, cx, mean)
  log_msg(sprintf("  %s %s share green by month: %s", s, pw, paste(sprintf("%s %.0f%%", month.abb[x$month], 100 * x$green_mean), collapse = ", ")))
}

# ---- 2 season dates ------------------------------------------------------------------------------
sd_fun <- function(X, keys, v) do.call(rbind, lapply(split(X, X[, keys], drop = TRUE), function(d) {
  if (nrow(d) < 7 || all(is.na(d[[v]]))) return(NULL)
  o <- d[1, keys, drop = FALSE]; o$var <- v; o$months <- nrow(d); s <- season(d$doy, d[[v]])
  cbind(o, as.data.frame(as.list(s))) }))
SDs <- rbind(sd_fun(M, c("source", "site", "species", "name", "pathway", "year"), "green"),
             sd_fun(M[M$source == "SEV", ], c("source", "site", "species", "name", "pathway", "year"), "new"))
SDp <- rbind(sd_fun(P, c("source", "site", "pathway", "year"), "green"), sd_fun(P[P$source == "SEV", ], c("source", "site", "pathway", "year"), "new"))
wcsv(SDs, "season_dates_species.csv"); wcsv(SDp, "season_dates_pathway.csv")
TYP <- do.call(rbind, lapply(split(SDs, SDs[, c("source", "species", "name", "pathway", "var")], drop = TRUE), function(d) {
  if (length(unique(d$year)) < 5) return(NULL)
  o <- d[1, c("source", "species", "name", "pathway", "var")]; o$site_years <- nrow(d); o$sites <- length(unique(d$site))
  for (k in c("up1", "down1", "up2", "down_last")) { o[[k]] <- stats::median(d[[k]], na.rm = TRUE); o[[paste0(k, "_sd")]] <- stats::sd(d[[k]] - stats::ave(d[[k]], d$site, FUN = function(z) mean(z, na.rm = TRUE)), na.rm = TRUE) }
  o }))
wcsv(TYP, "typical_species.csv")
log_msg(sprintf("season dates: %d species site-years; typical dates for %d species x variable", nrow(SDs), nrow(TYP)))

# ---- 3 C3 vs C4 at Sevilleta -----------------------------------------------------------------
sev <- P[P$source == "SEV", ]
W3 <- merge(sev[sev$pathway == "C3", c("site", "year", "month", "green", "new", "n")], sev[sev$pathway == "C4", c("site", "year", "month", "green", "new", "n")],
            by = c("site", "year", "month"), suffixes = c("_c3", "_c4"))
W3$both_inactive <- W3$new_c3 < 0.5 & W3$new_c4 < 0.5; W3$both_brown <- W3$green_c3 < 0.5 & W3$green_c4 < 0.5
wcsv(W3, "c3c4_monthly.csv")
G3 <- merge(SDp[SDp$source == "SEV" & SDp$pathway == "C3", c("site", "year", "var", "up1", "down1", "up2")],
            SDp[SDp$source == "SEV" & SDp$pathway == "C4", c("site", "year", "var", "up1", "down1", "down_last")], by = c("site", "year", "var"), suffixes = c("_c3", "_c4"))
G3$gap <- G3$up1_c4 - G3$down1_c3                                    # positive: C3 inactive before C4 starts
wcsv(G3, "c3c4_gap_years.csv")
if (nrow(W3)) {
  a <- stats::aggregate(cbind(new_c3, new_c4, both_inactive, both_brown) ~ month, W3, mean)
  log_msg("  SEV C3 vs C4, share with new green leaves by month (C3 / C4; share of site-years both below 50 %): ",
          paste(sprintf("%s %.0f/%.0f%% (%.0f%%)", month.abb[a$month], 100 * a$new_c3, 100 * a$new_c4, 100 * a$both_inactive), collapse = ", "))
  for (v in c("new", "green")) { g <- G3[G3$var == v & !is.na(G3$gap), ]
    if (nrow(g)) log_msg(sprintf("  SEV gap (C4 up1 - C3 down1, %s): median %+.0f d, %.0f %% of %d site-years > 0; C3 down %s, C4 up %s", v,
                                 stats::median(g$gap), 100 * mean(g$gap > 0), nrow(g), lab(stats::median(g$down1_c3)), lab(stats::median(g$up1_c4)))) }
}

# ---- 4 agreement with NPN ------------------------------------------------------------------------
shared <- intersect(unique(M$name), unique(NM$name))
A4 <- list(); O4 <- list()
for (s in unique(M$source)) {
  xy <- lter_xy[lter_xy$source == s, ]
  nm <- NM[dist_km(xy$lat, xy$lon, NM$lat, NM$lon) <= 400 & NM$n >= 3, ]
  for (sp in intersect(shared, unique(M$name[M$source == s]))) {
    a <- agg(green ~ month, M[M$source == s & M$name == sp, ])
    b <- agg(share_any ~ month, nm[nm$name == sp, ])
    if (is.null(a) || is.null(b) || nrow(b) < 6 || nrow(a) < 6) next
    x <- merge(a, b, by = "month")
    sa <- season(15 + 30.4 * (a$month - 1), a$green); sb <- season(15 + 30.4 * (b$month - 1), b$share_any)
    A4[[length(A4) + 1]] <- data.frame(source = s, name = sp, npn_sites = length(unique(nm$site_id[nm$name == sp])), months = nrow(x),
      r = if (nrow(x) >= 4 && stats::sd(x$green) > 0 && stats::sd(x$share_any) > 0) stats::cor(x$green, x$share_any) else NA, mean_abs_diff = mean(abs(x$green - x$share_any)),
      lter_up = sa[["up1"]], npn_up = sb[["up1"]], lter_down = sa[["down_last"]], npn_down = sb[["down_last"]])
  }
  # same place, same years
  near <- NM[dist_km(xy$lat, xy$lon, NM$lat, NM$lon) <= 30 & NM$n >= 3, ]
  for (sp in intersect(unique(near$name), unique(M$name[M$source == s]))) {
    a <- agg(green ~ year + month, M[M$source == s & M$name == sp, ])
    b <- agg(share_any ~ year + month, near[near$name == sp, ])
    if (is.null(a) || is.null(b)) next
    x <- merge(a, b, by = c("year", "month")); if (nrow(x) < 12) next
    x$ga <- x$green - stats::ave(x$green, x$month); x$ba <- x$share_any - stats::ave(x$share_any, x$month)
    O4[[length(O4) + 1]] <- data.frame(source = s, name = sp, years = paste(range(x$year), collapse = "-"), months = nrow(x),
      r = stats::cor(x$green, x$share_any), r_anom = stats::cor(x$ga, x$ba), mean_abs_diff = mean(abs(x$green - x$share_any)))
    wcsv(x, sprintf("overlap_%s_%s.csv", s, gsub(" ", "_", sp)))
  }
}
A4 <- do.call(rbind, A4); O4 <- do.call(rbind, O4)
if (!is.null(A4)) { wcsv(A4, "npn_same_species.csv")
  for (i in seq_len(nrow(A4))) with(A4[i, ], log_msg(sprintf("  vs NPN %s %-26s %2d NPN sites: r %.2f, |diff| %.2f; up %s vs %s, down %s vs %s",
                                                               source, name, npn_sites, r, mean_abs_diff, lab(lter_up), lab(npn_up), lab(lter_down), lab(npn_down)))) }
if (!is.null(O4)) { wcsv(O4, "npn_same_place.csv")
  for (i in seq_len(nrow(O4))) with(O4[i, ], log_msg(sprintf("  same place %s %s %s (%d months): r %.2f, anomalies r %.2f", source, name, years, months, r, r_anom))) }

# ---- 5 2016-17 --------------------------------------------------------------------------------------
Y5 <- do.call(rbind, lapply(split(M, paste(M$source, M$year)), function(d) {
  js <- d[d$month %in% 7:9 & d$pathway == "C4", ]
  mx <- agg(green ~ site + species, d, max)
  data.frame(source = d$source[1], year = d$year[1], green_jas_c4 = if (nrow(js)) wmean(js$green, js$n) else NA,
             series = nrow(mx), never_50 = mean(mx$green < 0.5))
}))
nb <- NM[NM$neon %in% TRUE, ]
Y5n <- do.call(rbind, lapply(split(nb, nb$year), function(d) data.frame(source = "NPN NEON", year = d$year[1],
  green_jas_c4 = { e <- d[d$month %in% 7:9 & d$pathway %in% "C4", ]; if (nrow(e)) wmean(e$share50, e$n) else NA },
  series = length(unique(paste(d$site_id, d$species_id))), never_50 = { mx <- agg(share50 ~ site_id + species_id, d, max); mean(mx$share50 < 0.5) })))
Y5 <- rbind(Y5, Y5n); wcsv(Y5, "years_1617.csv")
for (s in unique(Y5$source)) { y <- Y5[Y5$source == s, ]; i <- y$year %in% 2016:2017
  if (any(i) && any(!i & y$year >= 2010)) log_msg(sprintf("  2016-17 check %s: Jul-Sep C4 share green %.2f vs %.2f other years; series never above 50 %%: %.0f %% vs %.0f %%",
    s, mean(y$green_jas_c4[i], na.rm = TRUE), mean(y$green_jas_c4[!i & y$year >= 2010], na.rm = TRUE), 100 * mean(y$never_50[i]), 100 * mean(y$never_50[!i & y$year >= 2010]))) }

# ---- 6 one climate: Jornada sites ------------------------------------------------------------------
J6 <- list()
jr <- SDs[SDs$source == "JRN" & SDs$var == "green", ]
for (sp in unique(jr$name)) for (k in c("up1", "down_last")) {
  d <- jr[jr$name == sp & !is.na(jr[[k]]), ]; ns <- table(d$site); d <- d[d$site %in% names(ns)[ns >= 5], ]
  if (length(unique(d$site)) < 3) next
  sm <- tapply(d[[k]], d$site, mean)
  J6[[length(J6) + 1]] <- data.frame(name = sp, date = k, sites = length(sm), site_years = nrow(d), typical = stats::median(d[[k]]),
    sd_between_sites = stats::sd(sm), sd_within_site_years = stats::sd(d[[k]] - sm[d$site]))
}
J6 <- do.call(rbind, J6)
if (!is.null(J6)) { wcsv(J6, "jrn_sites_one_climate.csv")
  log_msg(sprintf("  Jornada, one climate: %d species x date; SD across sites median %.0f d vs year to year %.0f d", nrow(J6),
                  stats::median(J6$sd_between_sites), stats::median(J6$sd_within_site_years))) }

# ---- 7 calibration on the NPN record ----------------------------------------------------------------
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds"))
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
go <- go[!is.na(go$pct_green), ]; go$py <- paste(go$individual_id, go$year)
cmu <- cm[cm$usable %in% TRUE, ]; cmu$py <- paste(cmu$individual_id, cmu$year)
ty_f <- file.path(cfg$run_dir, "npn_study", "record", "series_typical.csv")
neon_sites <- if (file.exists(ty_f)) { t <- utils::read.csv(ty_f); unique(t$site_id[t$group == "NEON"]) } else unique(go$site_id)
# (a) monthly thinning
g7 <- go[go$py %in% cmu$py[cmu$n_obs >= 20] & go$site_id %in% neon_sites, ]
C7 <- do.call(rbind, lapply(split(g7, g7$py), function(d) {
  d <- d[order(d$doy), ]; full <- season(d$doy, d$pct_green / 100)
  m <- as.integer(substr(d$observation_date, 6, 7)); pick <- tapply(seq_len(nrow(d)), m, function(i) i[which.min(abs(as.integer(substr(d$observation_date[i], 9, 10)) - 15))])
  thin <- season(d$doy[unlist(pick)], d$pct_green[unlist(pick)] / 100)
  data.frame(up_full = full[["up1"]], up_thin = thin[["up1"]], down_full = full[["down_last"]], down_thin = thin[["down_last"]])
}))
cal <- data.frame(test = c("monthly vs all readings: green-up", "monthly vs all readings: curing"), n = c(sum(!is.na(C7$up_full) & !is.na(C7$up_thin)), sum(!is.na(C7$down_full) & !is.na(C7$down_thin))),
  bias = c(mean(C7$up_thin - C7$up_full, na.rm = TRUE), mean(C7$down_thin - C7$down_full, na.rm = TRUE)),
  mae = c(mean(abs(C7$up_thin - C7$up_full), na.rm = TRUE), mean(abs(C7$down_thin - C7$down_full), na.rm = TRUE)),
  lost = c(mean(!is.na(C7$up_full) & is.na(C7$up_thin)), mean(!is.na(C7$down_full) & is.na(C7$down_thin))))
# (b) population date vs median plant date, per site-species-year (>= 3 plants)
g7b <- go[go$py %in% cmu$py, ]; g7b$m <- as.integer(substr(g7b$observation_date, 6, 7)); g7b$ssy <- paste(g7b$site_id, g7b$species_id, g7b$year)
pm <- do.call(rbind, lapply(split(g7b, g7b$ssy), function(d) {
  if (length(unique(d$individual_id)) < 3) return(NULL)
  a <- stats::aggregate(cbind(any = pct_green >= 5, half = pct_green >= 50, doy) ~ m, d, mean)
  if (nrow(a) < 7) return(NULL)
  c_ <- cmu[cmu$site_id == d$site_id[1] & cmu$species_id == d$species_id[1] & cmu$year == d$year[1], ]
  s1 <- season(a$doy, a$any); s2 <- season(a$doy, a$half)
  data.frame(plant_up = stats::median(c_$g50_up_doy, na.rm = TRUE), plant_down = stats::median(c_$g50_down_doy, na.rm = TRUE),
             any_up = s1[["up1"]], any_down = s1[["down_last"]], half_up = s2[["up1"]], half_down = s2[["down_last"]])
}))
cal <- rbind(cal,
  data.frame(test = c("half the plants >= 5 % green vs median plant 50 %: green-up", "half the plants >= 5 % green vs median plant 50 %: curing",
                      "half the plants >= 50 % green vs median plant 50 %: green-up", "half the plants >= 50 % green vs median plant 50 %: curing"),
             n = c(sum(!is.na(pm$any_up - pm$plant_up)), sum(!is.na(pm$any_down - pm$plant_down)), sum(!is.na(pm$half_up - pm$plant_up)), sum(!is.na(pm$half_down - pm$plant_down))),
             bias = c(mean(pm$any_up - pm$plant_up, na.rm = TRUE), mean(pm$any_down - pm$plant_down, na.rm = TRUE), mean(pm$half_up - pm$plant_up, na.rm = TRUE), mean(pm$half_down - pm$plant_down, na.rm = TRUE)),
             mae = c(mean(abs(pm$any_up - pm$plant_up), na.rm = TRUE), mean(abs(pm$any_down - pm$plant_down), na.rm = TRUE), mean(abs(pm$half_up - pm$plant_up), na.rm = TRUE), mean(abs(pm$half_down - pm$plant_down), na.rm = TRUE)),
             lost = NA))
wcsv(cal, "calibration.csv")
for (i in seq_len(nrow(cal))) with(cal[i, ], log_msg(sprintf("  calibration %-62s n %5d: bias %+5.1f d, MAE %4.1f d%s", test, n, bias, mae, if (is.na(lost)) "" else sprintf(", lost %.0f %%", 100 * lost))))

# ---- 8 weather --------------------------------------------------------------------------------------
df <- file.path(ldir, "lter_site_daily.rds")
if (file.exists(df)) {
  D <- readRDS(df)
  WM <- npn_site_year_metrics(D, cfg$climate)
  WM$site_key <- WM$site_id
  drv <- setdiff(names(WM)[vapply(WM, is.numeric, TRUE)], c("year"))
  SDs$site_key <- paste(SDs$source, SDs$site, sep = "_")
  R8 <- list(); K8 <- list()
  for (s in unique(SDs$source)) for (pw in c("C3", "C4")) for (v in c("green", "new")) for (k in c("up1", "down1", "down_last")) {
    d <- SDs[SDs$source == s & SDs$pathway == pw & SDs$var == v & !is.na(SDs[[k]]), ]
    d$series <- paste(d$site, d$species); ny <- table(d$series); d <- d[d$series %in% names(ny)[ny >= 8], ]
    if (nrow(d) < 30) next
    d$y <- d[[k]] - stats::ave(d[[k]], d$series)
    X <- merge(d[, c("site_key", "year", "series", "y")], WM[, c("site_key", "year", drv)], by = c("site_key", "year"))
    if (nrow(X) < 30) next
    for (z in drv) X[[z]] <- X[[z]] - stats::ave(X[[z]], X$site_key, FUN = function(q) mean(q, na.rm = TRUE))
    rr <- vapply(drv, function(z) suppressWarnings(stats::cor(X$y, X[[z]], method = "spearman", use = "complete.obs")), 0)
    R8[[length(R8) + 1]] <- data.frame(source = s, pathway = pw, var = v, date = k, driver = drv, r = rr, n = nrow(X), series = length(unique(X$series)))
    pred <- rep(NA_real_, nrow(X))
    for (yy in unique(X$year)) { tr <- X$year != yy
      r_tr <- vapply(drv, function(z) suppressWarnings(abs(stats::cor(X$y[tr], X[[z]][tr], method = "spearman", use = "complete.obs"))), 0)
      b <- names(which.max(r_tr)); f <- stats::lm(X$y[tr] ~ X[[b]][tr])
      pred[!tr] <- stats::coef(f)[1] + stats::coef(f)[2] * X[[b]][!tr] }
    ok <- !is.na(pred)
    K8[[length(K8) + 1]] <- data.frame(source = s, pathway = pw, var = v, date = k, n = sum(ok), years = length(unique(X$year)), sd_days = stats::sd(X$y),
      best_driver = names(which.max(abs(rr))), best_r = rr[which.max(abs(rr))], skill = 1 - mean(abs(X$y[ok] - pred[ok])) / mean(abs(X$y[ok])))
  }
  # spring vs monsoon green-up: warm-season grasses that first pass 50 % green after June 1 vs earlier,
  # against winter-spring rain (Dec-May) and the 60-day moisture memory on May 1 (site-year anomalies)
  B8 <- list(); YR8 <- list()
  for (s in unique(SDs$source)) for (v in c("green", "new")) {
    d <- SDs[SDs$source == s & SDs$pathway == "C4" & SDs$var == v & !is.na(SDs$up1), ]
    if (nrow(d) < 30) next
    X <- merge(d[, c("site_key", "year", "up1")], WM[, c("site_key", "year", "p_djf", "p_mam", "ewma60_may1", "onset")], by = c("site_key", "year"))
    X$p_dm <- X$p_djf + X$p_mam
    for (z in c("p_dm", "ewma60_may1", "onset")) X[[paste0(z, "_a")]] <- X[[z]] - stats::ave(X[[z]], X$site_key, FUN = function(q) mean(q, na.rm = TRUE))
    X$monsoon <- X$up1 > 152
    if (sum(X$monsoon) < 5 || sum(!X$monsoon) < 5) next
    Yr <- stats::aggregate(cbind(monsoon, p_dm_a) ~ year, X, mean)             # one value per year (site-years are not independent)
    YR8[[length(YR8) + 1]] <- cbind(source = s, var = v, Yr)
    ry <- suppressWarnings(stats::cor.test(Yr$monsoon, Yr$p_dm_a, method = "spearman", exact = FALSE))
    B8[[length(B8) + 1]] <- data.frame(source = s, var = v, site_years = nrow(X), years = nrow(Yr), share_monsoon = mean(X$monsoon),
      r_years = unname(ry$estimate), p_years = ry$p.value,
      p_dm_spring = mean(X$p_dm[!X$monsoon]), p_dm_monsoon = mean(X$p_dm[X$monsoon]),
      p_dm_anom_monsoon = mean(X$p_dm_a[X$monsoon]), p_value = tryCatch(stats::wilcox.test(X$p_dm_a ~ X$monsoon)$p.value, error = function(e) NA),
      auc_dry_spring = tryCatch({ r <- rank(-X$p_dm_a); (sum(r[X$monsoon]) - sum(X$monsoon) * (sum(X$monsoon) + 1) / 2) / (sum(X$monsoon) * sum(!X$monsoon)) }, error = function(e) NA))
  }
  B8 <- do.call(rbind, B8); if (length(YR8)) wcsv(do.call(rbind, YR8), "weather_monsoon_years.csv")
  if (!is.null(B8)) { wcsv(B8, "weather_monsoon_greenup.csv")
    for (i in seq_len(nrow(B8))) with(B8[i, ], log_msg(sprintf("  monsoon green-up %s %-5s: %.0f %% of %d site-years green up after Jun 1; Dec-May rain %.0f vs %.0f mm in spring green-up years; across %d years, share of monsoon green-ups vs Dec-May rain r %+.2f (p %.3f)",
                                                                 source, var, 100 * share_monsoon, site_years, p_dm_monsoon, p_dm_spring, years, r_years, p_years))) }
  R8 <- do.call(rbind, R8); K8 <- do.call(rbind, K8)
  if (!is.null(R8)) wcsv(R8, "weather_r.csv")
  if (!is.null(K8)) { wcsv(K8, "weather_skill.csv")
    for (i in seq_len(nrow(K8))) with(K8[i, ], log_msg(sprintf("  weather %s %s %-5s %-9s %4d site-yrs %2d yrs, SD %3.0f d: best %-16s r %+.2f; LOYO skill %+.2f",
                                                                 source, pathway, var, date, n, years, sd_days, best_driver, best_r, skill))) }
} else log_msg("  weather: no lter_site_daily.rds yet (run 55d on the P720 and copy it to ", ldir, ")")
# small tables for the page (npn_study/report/tables, in git)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(cal, file.path(tdir, "lter_calibration.csv"), row.names = FALSE)
if (exists("K8") && !is.null(K8)) utils::write.csv(K8, file.path(tdir, "lter_weather_skill.csv"), row.names = FALSE)
if (exists("B8") && !is.null(B8)) utils::write.csv(B8, file.path(tdir, "lter_monsoon_greenup.csv"), row.names = FALSE)
if (!is.null(O4)) utils::write.csv(O4, file.path(tdir, "lter_same_place.csv"), row.names = FALSE)
if (!is.null(A4)) utils::write.csv(A4, file.path(tdir, "lter_same_species.csv"), row.names = FALSE)
if (!is.null(J6)) utils::write.csv(J6, file.path(tdir, "lter_jrn_one_climate.csv"), row.names = FALSE)
log_msg("=== done | ", odir, " | log: ", log_file)
