#!/usr/bin/env Rscript
# Step 58c (NPN grass curing study): the species list for the page, with the metadata needed to
# read the results: grass type (C3 / C4), life cycle, origin, how sure the C3/C4 label is, where
# each species was observed (EPA Level II ecoregions, states, NEON vs other sites), how much record
# it has (plant-years, plants, readings, years, series in the year-to-year set) and its typical
# green-up, 50 % curing and cured dates (circular median over its site series, with quartiles
# across sites).
# Same relaxed set as 58: 02c usable plant-years, crops out, study years. Site groups (NEON /
# other / managed) and series typical dates come from 58's record/series_typical.csv, so no site
# names are read here.
# Outputs
#   <run_dir>/npn_study/record/species.csv     one row per species, all columns
#   npn_study/report/tables/npn_species.csv     the same, for the page (no site ids or names)
#   npn_study/report/tables/npn_species_spread.csv  year-to-year vs plant-to-plant SD per species and
#                                                date (stricter year-to-year set; see below)
#
# Usage: Rscript npn_study/scripts/58c_npn_species.R [profile]   (after 58)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))

cfg <- load_npn_cfg(profile)
log_file <- log_init("58c_npn_species", profile)
ns <- cfg$npn_study %||% list()
rdir <- file.path(cfg$run_dir, "npn_study", "record")
ty_csv <- file.path(rdir, "series_typical.csv")
if (!file.exists(ty_csv)) stop("run 58 first: missing ", ty_csv)
log_msg("=== 58c_npn_species | profile: ", profile)
yrs <- seq(ns$years$start %||% 2001, ns$years$end %||% 2025)

cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
# quartile across sites, unwrapped around the circular median (returned in day of year)
cq <- function(x, p) { x <- x[!is.na(x)]; if (length(x) < 3) return(NA_real_); m <- cmed(x); d <- (x - m) %% 365.25
  d <- ifelse(d > 182.625, d - 365.25, d); (m + stats::quantile(d, p, names = FALSE)) %% 365.25 }

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
m <- match(cm$species_id, sp$species_id)
cm$crop <- sp$crop[m] %in% "yes"
u <- cm[cm$usable %in% TRUE & !cm$crop & cm$year %in% yrs, ]
ty <- utils::read.csv(ty_csv, stringsAsFactors = FALSE)
grp <- ty$group[match(u$site_id, ty$site_id)]; grp[is.na(grp)] <- "other"
u$neon <- grp == "NEON"
u$region <- npn_site_regions(u$site_id)
log_msg(sprintf("relaxed set: %d plant-years, %d species, %d sites", nrow(u), length(unique(u$species_id)), length(unique(u$site_id))))

collapse <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x) & x != "none"])), collapse = ", ")
eco_sort <- function(x) { k <- unique(x[!is.na(x) & x != "none"]); paste(k[order(as.numeric(k))], collapse = ", ") }
S <- do.call(rbind, lapply(split(u, u$species_id), function(d) {
  id <- d$species_id[1]; s <- sp[sp$species_id == id, ][1, ]; t <- ty[ty$species_id == id, ]
  out <- data.frame(species_id = id, genus = s$genus, species = s$species, common_name = s$common_name,
    pathway = s$pathway, life_cycle = s$life_cycle, origin = s$origin, pathway_confidence = s$confidence,
    sites = length(unique(d$site_id)), sites_neon = length(unique(d$site_id[d$neon])),
    plant_years = nrow(d), plant_years_neon = sum(d$neon), plants = length(unique(d$individual_id)),
    readings = sum(d$n_obs), first_year = min(d$year), last_year = max(d$year), years = length(unique(d$year)),
    ecoregions = eco_sort(d$region), states = collapse(d$state),
    series = nrow(t), series_yty = sum(t$yty %in% TRUE), stringsAsFactors = FALSE)
  for (tg in c("greenup", "curing", "cured")) {
    out[[paste0(tg, "_doy")]] <- cmed(t[[tg]]); out[[paste0(tg, "_q25")]] <- cq(t[[tg]], 0.25)
    out[[paste0(tg, "_q75")]] <- cq(t[[tg]], 0.75); out[[paste0(tg, "_series")]] <- sum(!is.na(t[[tg]]))
  }
  out
}))
S <- S[order(S$pathway, -S$plant_years), ]
utils::write.csv(S, file.path(rdir, "species.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(S, file.path(tdir, "npn_species.csv"), row.names = FALSE)

# ---- year-to-year vs plant-to-plant spread per species (year-to-year set) ----------------------
# Plant-years with >= min_readings readings; a series (site x species) needs the date in >= min_years
# years. Departures from the series' circular median, unwrapped; |departure| > 90 d dropped as a
# season mismatch (as in 52-54).
#   sd_year   SD of series-year departures (median over that year's plants), pooled over series
#   sd_plant  SD of plants around their series-year median, series-years with >= 2 plants, pooled
# Output report/tables/npn_species_spread.csv: one row per species x date with >= 3 series.
min_rd <- ns$min_readings %||% 12; min_yr <- ns$min_years %||% 3
cdev <- function(x, m) { d <- (x - m) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
v <- u[u$n_obs >= min_rd, ]; v$series <- paste(v$site_id, v$species_id)
cols <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
SPR <- list()
for (tg in names(cols)) {
  w <- v[!is.na(v[[cols[[tg]]]]), c("series", "species_id", "year", "individual_id", cols[[tg]])]; names(w)[5] <- "x"
  ny <- tapply(w$year, w$series, function(z) length(unique(z))); w <- w[w$series %in% names(ny)[ny >= min_yr], ]
  if (!nrow(w)) next
  med <- c(tapply(w$x, w$series, cmed)); w$d <- cdev(w$x, med[w$series]); w <- w[abs(w$d) <= 90, ]
  w$sy <- paste(w$series, w$year)
  w$sy_med <- stats::ave(w$d, w$sy, FUN = stats::median)
  for (id in unique(w$species_id)) {
    z <- w[w$species_id == id, ]; ns_ <- length(unique(z$series)); if (ns_ < 3) next
    yd <- stats::aggregate(d ~ series + year, z, stats::median)          # series-year departures
    yd$d <- yd$d - stats::ave(yd$d, yd$series)                            # around each series' mean
    sd_year <- sqrt(sum(yd$d^2) / max(nrow(yd) - ns_, 1))
    k <- table(z$sy); zp <- z[z$sy %in% names(k)[k >= 2], ]
    sd_plant <- if (nrow(zp)) sqrt(sum((zp$d - zp$sy_med)^2) / max(nrow(zp) - length(unique(zp$sy)), 1)) else NA
    SPR[[length(SPR) + 1]] <- data.frame(species_id = id, target = tg, series = ns_, series_years = nrow(yd),
      plant_years = nrow(z), sd_year = sd_year, sd_plant = sd_plant)
  }
}
SPR <- do.call(rbind, SPR)
if (!is.null(SPR)) {
  SPR <- merge(S[, c("species_id", "genus", "species", "common_name", "pathway", "life_cycle")], SPR, by = "species_id")
  utils::write.csv(SPR, file.path(tdir, "npn_species_spread.csv"), row.names = FALSE)
  for (tg in unique(SPR$target)) { s <- SPR[SPR$target == tg, ]
    log_msg(sprintf("spread %-8s %2d species: SD between years median %.0f d (range %.0f-%.0f), between plants median %.0f d",
                    tg, nrow(s), stats::median(s$sd_year), min(s$sd_year), max(s$sd_year), stats::median(s$sd_plant, na.rm = TRUE))) }
}

tb <- table(S$pathway, S$life_cycle)
log_msg(sprintf("species: %d (%s)", nrow(S), paste(sprintf("%s %s %d", rep(rownames(tb), ncol(tb)), rep(colnames(tb), each = nrow(tb)), as.integer(tb))[as.integer(tb) > 0], collapse = ", ")))
log_msg(sprintf("species with >= 3 sites: %d; >= 100 plant-years: %d; top 10 carry %.0f %% of plant-years",
                sum(S$sites >= 3), sum(S$plant_years >= 100), 100 * sum(sort(S$plant_years, decreasing = TRUE)[1:min(10, nrow(S))]) / sum(S$plant_years)))
log_msg(sprintf("=== done | %s, %s | log: %s", file.path(rdir, "species.csv"), file.path(tdir, "npn_species.csv"), log_file))
