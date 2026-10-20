#!/usr/bin/env Rscript
# Step 58b (NPN grass curing study): ground timing of the gap between cool-season (C3) and
# warm-season (C4) grass growth. Request from the fire study (project doc gsi-npn-c3c4-request.md,
# 19 Oct 2026): does fire activity dip between the C3 and C4 pulses in mixed prairie?
# Two gaps, both C4 green-up (50 % green going up) minus a C3 date, in days (positive = a gap
# with neither type green; negative = overlap, C4 is green before C3 has dried). Plain within-calendar-
# year differences, not circular: a circular difference flips sign past +/- 182 d (fixed 20 Oct):
#   gap_cured   C3 cured date  -> C4 green-up   (the request's definition)
#   gap_curing  C3 50 % curing -> C4 green-up   (both types below 50 % green in between)
# Two views per region:
#   same_site   sites that observe both pathways: each site's typical dates (circular median over
#               its C3 / C4 plant-years); year-to-year spread = SD of the site-year gap around the
#               site's mean gap (site-years with both a C3 and a C4 date), pooled in the region
#   regional    all sites: region medians of the site typical dates; year-to-year spread = SD over
#               region-years of (median C4 green-up) - (median C3 date) across that year's sites
# Plant-years: the relaxed set of 58 (02c usable: >= 6 readings, peak > 50 % green), crops out,
# study years, season mismatches (> 90 d from the series' circular median) dropped.
# Regions: EPA Level II ecoregions (default, data/npn_site_ecoregions.csv) or any CSV with
# site_id + region, or with site = "npn<site_id>|..." + region (e.g. the fuel phenology regions'
# eval/target_check/site_years.csv, A2 13-region map numbers, for the fire page).
# Output <run_dir>/npn_study/record/c3c4_gap_<tag>/: sites.csv (per site), site_years.csv,
# regions.csv (the table). tag = "eco" or the second argument's tag.
#
# Usage: Rscript npn_study/scripts/58b_npn_c3c4_gap.R [profile] [region_csv] [tag]
#   e.g. Rscript npn_study/scripts/58b_npn_c3c4_gap.R full C:/Users/Crimmins/data/gsi-pheno-regions/full/eval/target_check/site_years.csv a2

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
reg_csv <- if (length(args) >= 2) args[2] else ""
tag <- if (length(args) >= 3) args[3] else if (nzchar(reg_csv)) "custom" else "eco"
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("58b_npn_c3c4_gap", profile)
ns <- cfg$npn_study %||% list()
odir <- file.path(cfg$run_dir, "npn_study", "record", paste0("c3c4_gap_", tag)); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 58b_npn_c3c4_gap | profile: ", profile, " | regions: ", if (nzchar(reg_csv)) reg_csv else "EPA Level II ecoregions", " | tag ", tag)
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
doy_lab <- function(d) ifelse(is.na(d), "", format(as.Date("2001-01-01") + round(d) - 1, "%d %b"))

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id); cm$pathway <- sp$pathway[m]; cm$crop <- sp$crop[m] %in% "yes"
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$year %in% yrs & cm$pathway %in% c("C3", "C4"), ]
u$series <- paste(u$site_id, u$species_id)
cols <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
for (tg in names(cols)) {
  med <- c(tapply(u[[cols[[tg]]]], u$series, cmed)); d <- cdiff(u[[cols[[tg]]]], med[u$series])
  u[[tg]] <- ifelse(!is.na(d) & abs(d) > 90, NA, u[[cols[[tg]]]])
}

# regions
if (nzchar(reg_csv)) {
  r <- utils::read.csv(reg_csv, stringsAsFactors = FALSE)
  if (!"site_id" %in% names(r)) r$site_id <- as.integer(sub("^npn(\\d+)\\|.*$", "\\1", r$site))
  r <- r[!is.na(r$site_id) & !is.na(r$region), ]
  rr <- tapply(as.character(r$region), r$site_id, function(z) names(sort(table(z), decreasing = TRUE))[1])
  u$region <- as.character(rr[as.character(u$site_id)]); u$region[is.na(u$region)] <- "none"
  rname <- stats::setNames(paste("Region", sort(unique(u$region))), sort(unique(u$region)))
} else {
  u$region <- npn_site_regions(u$site_id); rname <- npn_region_names()
}
log_msg(sprintf("plant-years: %d (C3 %d, C4 %d) at %d sites; region known for %.0f %% of sites", nrow(u), sum(u$pathway == "C3"), sum(u$pathway == "C4"),
                length(unique(u$site_id)), 100 * mean(tapply(u$region != "none", u$site_id, any))))

# per site typical dates by pathway
S <- do.call(rbind, lapply(split(u, u$site_id), function(d) {
  f <- function(pw, tg) cmed(d[[tg]][d$pathway == pw])
  data.frame(site_id = d$site_id[1], region = d$region[1], lat = stats::median(d$lat), lon = stats::median(d$lon),
             c3_series = length(unique(d$series[d$pathway == "C3"])), c4_series = length(unique(d$series[d$pathway == "C4"])),
             c3_greenup = f("C3", "greenup"), c3_curing = f("C3", "curing"), c3_cured = f("C3", "cured"),
             c4_greenup = f("C4", "greenup"), c4_curing = f("C4", "curing"), c4_cured = f("C4", "cured"))
}))
S$both <- S$c3_series > 0 & S$c4_series > 0
S$gap_cured <- S$c4_greenup - S$c3_cured; S$gap_curing <- S$c4_greenup - S$c3_curing
utils::write.csv(S, file.path(odir, "sites.csv"), row.names = FALSE)

# per site-year (medians over that year's plants of each pathway)
SY <- do.call(rbind, lapply(split(u, paste(u$site_id, u$year)), function(d) {
  f <- function(pw, tg) cmed(d[[tg]][d$pathway == pw])
  data.frame(site_id = d$site_id[1], year = d$year[1], region = d$region[1],
             c3_curing = f("C3", "curing"), c3_cured = f("C3", "cured"), c4_greenup = f("C4", "greenup"))
}))
SY$gap_cured <- SY$c4_greenup - SY$c3_cured; SY$gap_curing <- SY$c4_greenup - SY$c3_curing
utils::write.csv(SY, file.path(odir, "site_years.csv"), row.names = FALSE)

# region table
R <- list()
for (rg in setdiff(sort(unique(S$region)), "none")) {
  s <- S[S$region == rg, ]; b <- s[s$both, ]; y <- SY[SY$region == rg, ]
  # same-site year-to-year: departures of site-year gaps from the site's mean gap
  dep <- function(col) { z <- y[!is.na(y[[col]]) & y$site_id %in% b$site_id, ]
    if (nrow(z) < 3) return(c(n = nrow(z), sd = NA)); z$d <- z[[col]] - stats::ave(z[[col]], z$site_id); c(n = nrow(z), sd = stats::sd(z$d) * sqrt(nrow(z) / max(nrow(z) - length(unique(z$site_id)), 1))) }
  dc <- dep("gap_cured"); dg <- dep("gap_curing")
  # regional: per year, median C4 green-up minus median C3 date across that year's sites
  ry <- do.call(rbind, lapply(split(y, y$year), function(z) data.frame(year = z$year[1],
          g = cmed(z$c4_greenup) - cmed(z$c3_cured), h = cmed(z$c4_greenup) - cmed(z$c3_curing))))
  R[[length(R) + 1]] <- data.frame(region = rg, region_name = if (rg %in% names(rname)) rname[[rg]] else "",
    sites = nrow(s), sites_c3 = sum(s$c3_series > 0), sites_c4 = sum(s$c4_series > 0), sites_both = nrow(b),
    # same-site view
    ss_c3_cured = cmed(b$c3_cured), ss_c3_curing = cmed(b$c3_curing), ss_c4_greenup = cmed(b$c4_greenup),
    ss_gap_cured = stats::median(b$gap_cured, na.rm = TRUE), ss_gap_curing = stats::median(b$gap_curing, na.rm = TRUE),
    ss_share_gap_cured_pos = mean(b$gap_cured > 0, na.rm = TRUE), ss_site_years = dc[["n"]], ss_sd_gap_cured = dc[["sd"]], ss_sd_gap_curing = dg[["sd"]],
    # regional view
    rg_c3_cured = cmed(s$c3_cured), rg_c3_curing = cmed(s$c3_curing), rg_c4_greenup = cmed(s$c4_greenup),
    rg_gap_cured = cmed(s$c4_greenup) - cmed(s$c3_cured), rg_gap_curing = cmed(s$c4_greenup) - cmed(s$c3_curing),
    rg_years = sum(!is.na(ry$g)), rg_sd_gap_cured = if (sum(!is.na(ry$g)) >= 3) stats::sd(ry$g, na.rm = TRUE) else NA,
    rg_sd_gap_curing = if (sum(!is.na(ry$h)) >= 3) stats::sd(ry$h, na.rm = TRUE) else NA)
}
RT <- do.call(rbind, R)
utils::write.csv(RT, file.path(odir, "regions.csv"), row.names = FALSE)
for (i in seq_len(nrow(RT))) with(RT[i, ], log_msg(sprintf(
  "%-5s %-28s sites %2d (C3 %2d, C4 %2d, both %2d) | same-site: C3 cured %s, C3 50%% %s, C4 green-up %s, gap cured %+4.0f d (%.0f %% of sites > 0), gap 50%% %+4.0f d, yr-to-yr SD %s / %s d (%d site-years) | regional: gap cured %+4.0f, gap 50%% %+4.0f, SD %s / %s d (%d years)",
  region, substr(region_name, 1, 28), sites, sites_c3, sites_c4, sites_both, doy_lab(ss_c3_cured), doy_lab(ss_c3_curing), doy_lab(ss_c4_greenup),
  ss_gap_cured, 100 * ss_share_gap_cured_pos, ss_gap_curing, format(round(ss_sd_gap_cured)), format(round(ss_sd_gap_curing)), as.integer(ss_site_years),
  rg_gap_cured, rg_gap_curing, format(round(rg_sd_gap_cured)), format(round(rg_sd_gap_curing)), as.integer(rg_years))))
log_msg(sprintf("=== done | %s | log: %s", odir, log_file))
