#!/usr/bin/env Rscript
# Step 58 (NPN grass curing study): the record (Q1), climatology (Q2) and climate across space
# (Q4), on the RELAXED set (decided 16 Oct): usable plant-years (>= 6 readings, peak >= 50 %
# green; the 02c `usable` flag), crops out, study years, series with >= 1 year. The year-to-year
# questions (Q3, Q5, Q6; scripts 52-54) keep >= 12 readings and >= 3 years; that set ("yty") is
# reported alongside as a sensitivity.
#
# Q1 the record   record.csv      plant-years, plants, sites, species, readings, dates, by set x
#                                 group (NEON / other / managed) x pathway
#                 readings.csv    readings per plant-year, median gap between readings, span
#                 bins.csv        share of readings in each percent-green bin
#                 coverage.csv    sites and plant-years per region (EPA Level II ecoregion code) x pathway
#                 peak_lt50.csv   plant-years with >= 6 readings that never pass 50 % green
#                                 (dropped by `usable`): by group, and per NEON year (share of
#                                 that year's plants) - possible drought-suppressed years
#                 second.csv      second green-ups: after falling below 25 %, back above 50 %
#                                 on two readings >= 21 days later; share by group, pathway, region; whether
#                                 the scored curing date belongs to the first or second cycle
# Q2 climatology  series_typical.csv  per site x species: typical (circular median) green-up,
#                                 50 % curing, cured, green-season length; years and plant-years
#                 climatology.csv region x pathway x date: series, sites, median, quartiles
#                 c3c4.csv        C4 minus C3 within regions with >= 3 series of each
#                                 (also annual vs perennial within C3)
# Q4 space        spatial_r.csv   Spearman r of typical dates with each site normal, by pathway
#                 spatial.csv     mixed models of typical dates, species random effect, weights
#                                 = years (capped 5): species only / region / climate / both;
#                                 marginal and conditional R2, leave-one-site-out R2, per
#                                 pathway; relaxed and yty sets
#                 spatial_coef.csv standardised climate coefficients (days per SD)
# Climate predictors (a priori, site normals from 51): green-up: Mar-Apr temperature, aridity
# (P / PET), warm-season rain share, latitude; curing / cured: Jul-Aug temperature, aridity,
# warm-season rain share, latitude.
# Names (data/npn_site_names.csv, local only) are used only for the NEON / managed flags.
# Inputs: npn_curing_metrics.rds and npn_green_obs.rds (GSI_OBS_DIR or config study_dir),
# data/npn_species_pathway.csv, <run_dir>/npn_study/npn_metrics.rds (51), regions = EPA Level II
# ecoregions (data/npn_site_ecoregions.csv, 59). Needs lme4. Runs on the laptop (a few minutes: leave-one-site-out refits).
# Output <run_dir>/npn_study/record/.
#
# Usage: Rscript npn_study/scripts/58_npn_record_climatology.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "full")

npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("58_npn_record_climatology", profile)
ns <- cfg$npn_study %||% list()
sdir <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(sdir, "record")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 58_npn_record_climatology | profile: ", profile)
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs the lme4 package")
set.seed(58)
yty_n <- ns$min_readings %||% 12; yty_years <- ns$min_years %||% 3
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)
wcsv <- function(x, f) utils::write.csv(x, file.path(odir, f), row.names = FALSE)

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cq <- function(x, p) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); m <- cmed(x); (m + stats::quantile(cdiff(x, m), p, names = FALSE)) %% 365.25 }

# ---- plant-years ------------------------------------------------------------------------------
ddir <- npn_obs_dir(cfg)
cm <- readRDS(file.path(ddir, "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id)
cm$pathway <- sp$pathway[m]; cm$life_cycle <- sp$life_cycle[m]; cm$crop <- sp$crop[m] %in% "yes"
nf <- npn_path("data", "npn_site_names.csv")
managed_re <- "garden|campus|universit|college|school|park|yard|home|house|residen|arboret|cemeter|golf|lawn|church|office|museum|zoo|nursery|farm|ranch hq|backyard|front|street|apartment|condo|library|center"
if (file.exists(nf)) {
  sn <- utils::read.csv(nf, stringsAsFactors = FALSE)
  nm <- tolower(sn$site_name[match(cm$site_id, sn$site_id)]); nm[is.na(nm)] <- ""
  cm$group <- ifelse(grepl("\\.phenology\\.phe", nm), "NEON", ifelse(grepl(managed_re, nm), "managed", "other"))
} else { log_msg("no data/npn_site_names.csv: NEON / managed flags unavailable"); cm$group <- "unknown" }

# regions: EPA / CEC Level II ecoregions (59), for every site
cm$region <- npn_site_regions(cm$site_id)
log_msg(sprintf("regions: EPA Level II ecoregions; region known for %.0f %% of plant-years, %d regions", 100 * mean(cm$region != "none"), length(unique(cm$region[cm$region != "none"]))))

u <- cm[cm$usable %in% TRUE & !cm$crop & cm$year %in% yrs, ]
u$series <- paste(u$site_id, u$species_id)
targets <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
# season mismatches (|departure from series circular median| > 90 d) set to NA per date
for (tg in names(targets)) {
  col <- targets[[tg]]; med <- c(tapply(u[[col]], u$series, cmed))
  d <- cdiff(u[[col]], med[u$series]); mis <- !is.na(d) & abs(d) > 90
  u[[tg]] <- ifelse(mis, NA, u[[col]])
  log_msg(sprintf("%-8s %d plant-years dated, %d season mismatches dropped", tg, sum(!is.na(u[[col]])), sum(mis)))
}
# yty flag: >= yty_n readings and series with >= yty_years years of any date
ny <- tapply(u$year[u$n_obs >= yty_n], u$series[u$n_obs >= yty_n], function(v) length(unique(v)))
u$yty <- u$n_obs >= yty_n & as.numeric(ny[u$series]) %in% seq(yty_years, 1000)
u$yty[is.na(u$yty)] <- FALSE
log_msg(sprintf("relaxed set: %d plant-years, %d plants, %d sites, %d species | yty subset: %d plant-years, %d sites",
                nrow(u), length(unique(u$individual_id)), length(unique(u$site_id)), length(unique(u$species_id)),
                sum(u$yty), length(unique(u$site_id[u$yty]))))

# ---- Q1 the record ------------------------------------------------------------------------------
grp_sets <- list(all = function(x) rep(TRUE, nrow(x)), NEON = function(x) x$group == "NEON",
                 other = function(x) x$group %in% c("other", "managed", "unknown"), managed = function(x) x$group == "managed")
rec <- list()
for (st in c("relaxed", "yty")) for (gn in names(grp_sets)) for (pw in c("all", "C3", "C4")) {
  x <- u[if (st == "yty") u$yty else rep(TRUE, nrow(u)), ]; x <- x[grp_sets[[gn]](x), ]
  if (pw != "all") x <- x[x$pathway %in% pw, ]
  rec[[length(rec) + 1]] <- data.frame(set = st, group = gn, pathway = pw, plant_years = nrow(x), plants = length(unique(x$individual_id)),
    sites = length(unique(x$site_id)), species = length(unique(x$species_id)), series = length(unique(x$series)), readings = sum(x$n_obs),
    first_year = if (nrow(x)) min(x$year) else NA, last_year = if (nrow(x)) max(x$year) else NA,
    greenup = sum(!is.na(x$greenup)), curing = sum(!is.na(x$curing)), cured = sum(!is.na(x$cured)),
    regions = length(unique(x$region[x$region != "none"])))
}
REC <- do.call(rbind, rec); wcsv(REC, "record.csv")
for (i in which(REC$pathway == "all" & REC$group %in% c("all", "NEON", "other")))
  with(REC[i, ], log_msg(sprintf("record %-7s %-5s %5d py %4d plants %3d sites %2d species %6d readings %d-%d | dated: green-up %d, curing %d, cured %d | %d regions",
                                 set, group, plant_years, plants, sites, species, readings, first_year, last_year, greenup, curing, cured, regions)))

# raw readings for bins, gaps and second green-ups
G <- readRDS(file.path(ddir, "npn_green_obs.rds"))
G$key <- paste(G$individual_id, G$year); u$key <- paste(u$individual_id, u$year)
G <- G[G$key %in% u$key, ]
G <- G[order(G$key, G$doy), ]
G$group <- u$group[match(G$key, u$key)]; G$pathway <- u$pathway[match(G$key, u$key)]
gap_med <- tapply(G$doy, G$key, function(d) if (length(d) > 1) stats::median(diff(unique(d))) else NA_real_)
u$gap_median <- as.numeric(gap_med[u$key])
RDN <- do.call(rbind, lapply(c("NEON", "other", "managed"), function(g) {
  x <- if (g == "other") u[u$group %in% c("other", "unknown"), ] else u[u$group == g, ]; if (!nrow(x)) return(NULL)
  q <- stats::quantile(x$n_obs, c(0.1, 0.25, 0.5, 0.75, 0.9))
  data.frame(group = g, plant_years = nrow(x), n_obs_p10 = q[[1]], n_obs_p25 = q[[2]], n_obs_median = q[[3]], n_obs_p75 = q[[4]], n_obs_p90 = q[[5]],
             gap_median_days = stats::median(x$gap_median, na.rm = TRUE), gap_p90_days = stats::quantile(x$gap_median, 0.9, na.rm = TRUE, names = FALSE),
             span_median_days = stats::median(x$span_days, na.rm = TRUE))
}))
wcsv(RDN, "readings.csv")
for (i in seq_len(nrow(RDN))) with(RDN[i, ], log_msg(sprintf("readings %-8s %5d py | per plant-year median %.0f (IQR %.0f-%.0f) | median gap %.0f d | span %.0f d",
                                                          group, plant_years, n_obs_median, n_obs_p25, n_obs_p75, gap_median_days, span_median_days)))
bin_lv <- c("Less than 5%", "5-24%", "25-49%", "50-74%", "75-94%", "95% or more")
G$grp2 <- ifelse(G$group == "NEON", "NEON", "other")
BN <- do.call(rbind, lapply(c("all", "NEON", "other"), function(g) {
  x <- if (g == "all") G else G[G$grp2 == g, ]; t <- table(factor(x$bin, bin_lv))
  data.frame(group = g, bin = bin_lv, readings = as.integer(t), share = as.numeric(t) / sum(t))
}))
wcsv(BN, "bins.csv")
log_msg(sprintf("bins (all readings): %s", paste(sprintf("%s %.0f %%", bin_lv, 100 * BN$share[BN$group == "all"]), collapse = ", ")))

CV <- stats::aggregate(list(plant_years = u$key), list(region = u$region, pathway = u$pathway), length)
CV$sites <- stats::aggregate(u$site_id, list(u$region, u$pathway), function(z) length(unique(z)))$x
CV$non_neon_sites <- stats::aggregate(u$site_id * (u$group != "NEON"), list(u$region, u$pathway), function(z) length(unique(z[z > 0])))$x
wcsv(CV[order(CV$region, CV$pathway), ], "coverage.csv")

# plant-years that never pass 50 % green (dropped by usable): count, and share by NEON year
pk <- cm[!cm$crop & cm$year %in% yrs & cm$n_obs >= 6, ]
pk$lt50 <- !(pk$peak_green >= 50)
PK <- do.call(rbind, lapply(c("NEON", "other", "managed"), function(g) {
  x <- if (g == "other") pk[pk$group %in% c("other", "unknown"), ] else pk[pk$group == g, ]
  data.frame(level = "group", group = g, year = NA, site_years = length(unique(paste(x$site_id, x$year))),
             plant_years = nrow(x), lt50 = sum(x$lt50), share = mean(x$lt50), site_years_with_lt50 = length(unique(paste(x$site_id, x$year)[x$lt50])))
}))
ne <- pk[pk$group == "NEON", ]
PY <- do.call(rbind, lapply(sort(unique(ne$year)), function(y) { x <- ne[ne$year == y, ]
  data.frame(level = "NEON_year", group = "NEON", year = y, site_years = length(unique(x$site_id)), plant_years = nrow(x), lt50 = sum(x$lt50),
             share = mean(x$lt50), site_years_with_lt50 = length(unique(x$site_id[x$lt50]))) }))
# site-years where most plants never passed 50 %
sy <- stats::aggregate(list(share = ne$lt50), list(site_id = ne$site_id, year = ne$year, pathway = ne$pathway), mean)
sy$n <- stats::aggregate(ne$lt50, list(ne$site_id, ne$year, ne$pathway), length)$x
wcsv(rbind(PK, PY), "peak_lt50.csv"); wcsv(sy[sy$share >= 0.5 & sy$n >= 3, ], "peak_lt50_site_years.csv")
log_msg(sprintf("never above 50 %% green (>= 6 readings): NEON %d of %d plant-years (%.0f %%), other %d of %d | NEON site-year x pathway with most plants below 50 %%: %d of %d",
                PK$lt50[1], PK$plant_years[1], 100 * PK$share[1], PK$lt50[2], PK$plant_years[2], sum(sy$share >= 0.5 & sy$n >= 3), sum(sy$n >= 3)))
log_msg(sprintf("  NEON share by year: %s", paste(sprintf("%d %.0f%%", PY$year, 100 * PY$share), collapse = ", ")))

# second green-ups: >= 50 %, then below 25 % (well into curing), then >= 50 % on two consecutive
# readings starting >= 21 days after the drop (16 Oct: a looser rule, back above 50 % once after
# any reading below 50 %, flagged 41 % of plant-years, mostly bin flicker around 50 %)
sec <- lapply(split(G[, c("doy", "pct_green")], G$key), function(d) {
  v <- d$pct_green; n <- length(v)
  i <- which(v >= 50)[1]; if (is.na(i)) return(c(NA, NA))
  j <- which(v < 25 & seq_len(n) > i)[1]; if (is.na(j)) return(c(0, NA))
  k <- which(v >= 50 & c(v[-1], -1) >= 50 & seq_len(n) > j & d$doy - d$doy[j] >= 21)[1]
  if (is.na(k)) c(0, NA) else c(1, d$doy[k])
})
S2 <- do.call(rbind, sec); i2 <- match(u$key, rownames(S2)); u$second <- S2[i2, 1]; u$second_doy <- S2[i2, 2]
u$curing_cycle <- ifelse(u$second %in% 1 & !is.na(u$curing), ifelse(u$curing > u$second_doy, "second", "first"), NA)
sx <- function(x, lev, val) data.frame(by = lev, value = val, plant_years = sum(!is.na(x$second)), second = sum(x$second %in% 1),
                                      share = mean(x$second[!is.na(x$second)] %in% 1), curing_on_second = sum(x$curing_cycle %in% "second"))
S2T <- rbind(sx(u, "all", "all"),
  do.call(rbind, lapply(c("NEON", "other"), function(g) sx(u[if (g == "NEON") u$group == "NEON" else u$group != "NEON", ], "group", g))),
  do.call(rbind, lapply(c("C3", "C4"), function(g) sx(u[u$pathway %in% g, ], "pathway", g))),
  do.call(rbind, lapply(sort(unique(u$region)), function(g) sx(u[u$region == g, ], "region", g))))
wcsv(S2T, "second.csv")
log_msg(sprintf("second green-ups: %d of %d plant-years (%.0f %%; C3 %.0f %%, C4 %.0f %%); curing date on the second cycle in %d | regions with share >= 10 %%: %s",
                S2T$second[1], S2T$plant_years[1], 100 * S2T$share[1], 100 * S2T$share[S2T$value == "C3"], 100 * S2T$share[S2T$value == "C4"],
                S2T$curing_on_second[1], paste(S2T$value[S2T$by == "region" & S2T$share >= 0.1 & S2T$plant_years >= 30], collapse = " ")))

# ---- Q2 climatology ------------------------------------------------------------------------------
mk_typ <- function(x) do.call(rbind, lapply(split(x, x$series), function(d) {
  o <- data.frame(series = d$series[1], site_id = d$site_id[1], species_id = d$species_id[1], genus = d$genus[1], species = d$species[1],
                  pathway = d$pathway[1], life_cycle = d$life_cycle[1], group = d$group[1], region = d$region[1],
                  lat = stats::median(d$lat), lon = stats::median(d$lon), plant_years = nrow(d), years = length(unique(d$year)))
  for (tg in names(targets)) {
    v <- d[[tg]]; o[[tg]] <- cmed(v)
    o[[paste0(tg, "_years")]] <- length(unique(d$year[!is.na(v)]))
  }
  o$season_len <- if (is.na(o$greenup) || is.na(o$curing)) NA else (o$curing - o$greenup) %% 365.25
  o
}))
TY <- mk_typ(u); TY$yty <- TY$series %in% unique(u$series[u$yty])
wcsv(TY, "series_typical.csv")
CL <- list()
for (tg in c(names(targets), "season_len")) for (rg in sort(unique(TY$region))) for (pw in c("C3", "C4")) {
  x <- TY[TY$region == rg & TY$pathway %in% pw & !is.na(TY[[tg]]), ]; if (!nrow(x)) next
  circ <- tg != "season_len"
  CL[[length(CL) + 1]] <- data.frame(target = tg, region = rg, pathway = pw, series = nrow(x), sites = length(unique(x$site_id)),
    non_neon_sites = length(unique(x$site_id[x$group != "NEON"])),
    median = if (circ) cmed(x[[tg]]) else stats::median(x[[tg]]), q25 = if (circ) cq(x[[tg]], 0.25) else stats::quantile(x[[tg]], 0.25, names = FALSE),
    q75 = if (circ) cq(x[[tg]], 0.75) else stats::quantile(x[[tg]], 0.75, names = FALSE))
}
CLT <- do.call(rbind, CL); wcsv(CLT, "climatology.csv")
for (tg in names(targets)) { z <- CLT[CLT$target == tg & CLT$series >= 3 & CLT$region != "none", ]
  log_msg(sprintf("climatology %-8s (region x pathway with >= 3 series): %s", tg,
                  paste(sprintf("%s %s %s (n %d)", z$region, z$pathway, format(as.Date("2001-01-01") + round(z$median) - 1, "%d %b"), z$series), collapse = "; "))) }
# C4 minus C3 within region; annual minus perennial within C3
cc <- list()
for (tg in c(names(targets), "season_len")) for (rg in sort(unique(TY$region[TY$region != "none"]))) {
  for (cmp in list(c("pathway", "C4", "C3"), c("life_cycle", "annual", "perennial"))) {
    x <- TY[TY$region == rg & !is.na(TY[[tg]]), ]
    if (cmp[1] == "life_cycle") x <- x[x$pathway %in% "C3", ]
    a <- x[[tg]][x[[cmp[1]]] %in% cmp[2]]; b <- x[[tg]][x[[cmp[1]]] %in% cmp[3]]
    if (length(a) < 3 || length(b) < 3) next
    dd <- if (tg == "season_len") stats::median(a) - stats::median(b) else cdiff(cmed(a), cmed(b))
    p <- tryCatch(stats::wilcox.test(a, b, exact = FALSE)$p.value, error = function(e) NA)
    cc[[length(cc) + 1]] <- data.frame(target = tg, region = rg, contrast = paste(cmp[2], "minus", cmp[3]), n_a = length(a), n_b = length(b), diff_days = dd, p = p)
  }
}
if (length(cc)) { CC <- do.call(rbind, cc); wcsv(CC, "c3c4.csv")
  for (i in seq_len(nrow(CC))) with(CC[i, ], log_msg(sprintf("  %-10s %-5s %-22s %+4.0f d (n %d vs %d, p %.2f)", target, region, contrast, diff_days, n_a, n_b, p))) }

# ---- Q4 climate across space ------------------------------------------------------------------------
M <- readRDS(file.path(sdir, "npn_metrics.rds")); N <- M$normals
preds <- list(greenup = c("t_marapr", "aridity", "p_jas_share", "lat"), curing = c("t_julaug", "aridity", "p_jas_share", "lat"),
              cured = c("t_julaug", "aridity", "p_jas_share", "lat"))
TY$wt <- pmin(TY$years, 5)
nm_cols <- setdiff(names(N), c("site_id", "lon"))
TN <- merge(TY, N[, c("site_id", setdiff(nm_cols, "lat"))], by = "site_id", all.x = TRUE)
log_msg(sprintf("normals found for %d of %d series", sum(!is.na(TN$aridity)), nrow(TN)))
# univariate Spearman r with every normal
sr <- list()
for (tg in names(targets)) for (pw in c("all", "C3", "C4")) {
  x <- TN[!is.na(TN[[tg]]) & (pw == "all" | TN$pathway %in% pw), ]
  x[[tg]] <- cmed(x[[tg]]) + cdiff(x[[tg]], cmed(x[[tg]]))           # unwrap around the circular median
  for (v in setdiff(nm_cols, c("year"))) { ok <- !is.na(x[[v]]); if (sum(ok) < 15 || stats::sd(x[[v]][ok]) == 0) next
    sr[[length(sr) + 1]] <- data.frame(target = tg, pathway = pw, normal = v, n = sum(ok), r = suppressWarnings(stats::cor(x[[tg]][ok], x[[v]][ok], method = "spearman"))) }
}
SR <- do.call(rbind, sr); SR <- SR[order(SR$target, SR$pathway, -abs(SR$r)), ]; wcsv(SR, "spatial_r.csv")
for (tg in names(targets)) for (pw in c("C3", "C4")) { z <- utils::head(SR[SR$target == tg & SR$pathway == pw, ], 4)
  log_msg(sprintf("spatial r %-8s %s: %s", tg, pw, paste(sprintf("%s %.2f", z$normal, z$r), collapse = ", "))) }

r2_mixed <- function(f) {
  vf <- stats::var(stats::predict(f, re.form = NA)); vc <- as.data.frame(lme4::VarCorr(f))
  vr <- sum(vc$vcov[vc$grp != "Residual"]); ve <- vc$vcov[vc$grp == "Residual"]
  c(marginal = vf / (vf + vr + ve), conditional = (vf + vr) / (vf + vr + ve))
}
fit_set <- function(x, tg, pw, setname) {
  x <- x[!is.na(x[[tg]]) & stats::complete.cases(x[, preds[[tg]]]) & x$region != "none" & x$pathway %in% c("C3", "C4"), ]
  if (pw != "all") x <- x[x$pathway %in% pw, ]
  if (nrow(x) < 20 || length(unique(x$site_id)) < 10) return(NULL)
  x$y <- cmed(x[[tg]]) + cdiff(x[[tg]], cmed(x[[tg]]))
  for (v in preds[[tg]]) x[[paste0("z_", v)]] <- (x[[v]] - mean(x[[v]])) / stats::sd(x[[v]])
  x$species <- factor(x$species_id); x$regionf <- factor(x$region)
  zc <- paste0("z_", preds[[tg]]); pwt <- if (pw == "all" && length(unique(x$pathway)) > 1) " + pathway" else ""
  forms <- list(species = paste0("y ~ 1", pwt), region = paste0("y ~ regionf", pwt), climate = paste0("y ~ ", paste(zc, collapse = " + "), pwt),
                both = paste0("y ~ regionf + ", paste(zc, collapse = " + "), pwt))
  out <- list(); coefs <- NULL
  for (mn in names(forms)) {
    fm <- stats::as.formula(paste(forms[[mn]], "+ (1 | species)"))
    f <- tryCatch(lme4::lmer(fm, data = x, weights = wt, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
    if (is.null(f)) next
    r2 <- r2_mixed(f)
    # leave one site out: refit without the site, predict with species effects (new species -> 0)
    pr <- rep(NA_real_, nrow(x))
    for (s in unique(x$site_id)) {
      te <- x$site_id == s; tr <- x[!te, ]
      if (mn %in% c("region", "both")) { te <- te & x$regionf %in% unique(tr$regionf) }
      if (!any(te)) next
      g <- tryCatch(lme4::lmer(fm, data = droplevels(tr), weights = wt, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE)), error = function(e) NULL)
      if (!is.null(g)) pr[te] <- tryCatch(stats::predict(g, newdata = x[te, ], allow.new.levels = TRUE), error = function(e) NA)
    }
    ok <- !is.na(pr)
    out[[mn]] <- data.frame(set = setname, target = tg, pathway = pw, model = mn, series = nrow(x), sites = length(unique(x$site_id)),
      species = length(unique(x$species)), r2_marginal = r2[["marginal"]], r2_conditional = r2[["conditional"]],
      r2_loso = 1 - sum((x$y[ok] - pr[ok])^2) / sum((x$y[ok] - mean(x$y[ok]))^2), mae_loso = mean(abs(x$y[ok] - pr[ok])), n_loso = sum(ok))
    if (mn == "climate") { cf <- lme4::fixef(f); se <- sqrt(diag(as.matrix(stats::vcov(f))))
      coefs <- data.frame(set = setname, target = tg, pathway = pw, term = sub("^z_", "", names(cf)), days_per_sd = unname(cf), se = unname(se))[names(cf) %in% zc, ] }
  }
  list(scores = do.call(rbind, out), coefs = coefs)
}
SP <- list(); SC <- list()
for (setname in c("relaxed", "yty")) for (tg in names(targets)) for (pw in c("all", "C3", "C4")) {
  x <- if (setname == "yty") TN[TN$yty, ] else TN
  r <- fit_set(x, tg, pw, setname); if (is.null(r)) next
  SP[[length(SP) + 1]] <- r$scores; SC[[length(SC) + 1]] <- r$coefs
  s <- r$scores; log_msg(sprintf("spatial %-7s %-8s %-3s %3d series %3d sites | marginal R2: region %.2f, climate %.2f, both %.2f | LOSO R2: species %.2f, region %.2f, climate %.2f, both %.2f",
    setname, tg, pw, s$series[1], s$sites[1], s$r2_marginal[s$model == "region"], s$r2_marginal[s$model == "climate"], s$r2_marginal[s$model == "both"],
    s$r2_loso[s$model == "species"], s$r2_loso[s$model == "region"], s$r2_loso[s$model == "climate"], s$r2_loso[s$model == "both"]))
}
wcsv(do.call(rbind, SP), "spatial.csv"); SCo <- do.call(rbind, SC); wcsv(SCo, "spatial_coef.csv")
z <- SCo[SCo$set == "relaxed" & SCo$pathway != "all", ]
for (i in seq_len(nrow(z))) with(z[i, ], log_msg(sprintf("  coef %-8s %s %-12s %+5.1f d per SD (se %.1f)", target, pathway, term, days_per_sd, se)))

log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
