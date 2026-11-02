#!/usr/bin/env Rscript
# Step 05 (NPN x fire study; laptop, < 2 min; arrow): Tier 1 inventory. Where do usable NPN dates and
# fire activity in the matching fuel coincide, by Level II ecoregion x plant group x phenophase event?
# Rules (config npn_fire: inventory) were fixed on 1 Nov 2026 before any test result was seen.
#
# Inputs (synced from the P720): clean/plant_years (02b), fire/fire_weekly (03); groups re-joined from
# data/nf_species_groups.csv (02a).
#   1. Usable dates: plant-year date of the event with its bracketing gap <= gap_days. Site-year date =
#      median over the site's plants; site typical = median of its site-years (circular-safe: dates are
#      within one year). Ecoregion typical = median of site typicals (sites weighted equally).
#   2. Fire in the matched fuel (config match): weekly sum over the plant group's fuel groups, per
#      outcome; window = typical week -window_weeks .. +window_weeks (clipped to weeks 1-52).
#   3. Per outcome: overlap years (NPN usable sites that year and the outcome's years); covered years =
#      overlap years with >= min_events_window in the window. Lead outcome = most covered years, ties by
#      the config order. Tier 2 possible: >= tier2 min_sites sites over the record and lead covered years
#      >= min_covered_years. Tier 3 possible: >= tier3 min_years years with >= min_sites_year sites that
#      are also covered by the lead outcome.
# Outputs <run_dir>/npn_fire/inventory/:
#   site_dates.parquet  eco_l2_code, group, event, site_id, year, doy, n_plants (site-year dates; 06 / 07 use them)
#   events.csv          eco x group x event: sites, site_years, years_ge5 / ge10 sites, typical_doy, typical_week, iqr_doy
#   coverage.csv        eco x group x event x outcome: overlap / covered years, window and yearly means
#   inventory.csv       eco x group x event: the above + lead outcome, tier2_ok, tier3_ok, tier3_years
#   fire_nonburnable.csv eco x outcome: share of each record in mostly nonburnable / unmapped cells (not matched)
# Usage: Rscript npn_fire/scripts/05_nf_inventory.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))

cfg <- load_nf_cfg(profile)
log_file <- log_init("05_nf_inventory", profile)
ic <- cfg$npn_fire$inventory
gap_ok <- ic$gap_days %||% 14; ww <- ic$window_weeks %||% 4; mev <- ic$min_events_window %||% 10
out <- nf_out(cfg, "inventory")
log_msg("=== 05_nf_inventory | profile: ", profile, " | gap ", gap_ok, " d | window +-", ww, " weeks | covered: >= ", mev, " events")

# ---- NPN usable dates -> site-year dates -------------------------------------------------------------------
pyf <- nf_table_path(file.path(cfg$out_dir, "clean", "plant_years"))
if (is.na(pyf)) stop("no clean/plant_years: run 02b and sync")
P <- nf_read_table(pyf); P$eco_l2_code <- as.character(P$eco_l2_code)
SG <- utils::read.csv(nf_path("data", "nf_species_groups.csv"), stringsAsFactors = FALSE)
P$group <- SG$group[match(P$species_id, SG$species_id)]
P <- P[P$group %in% names(ic$match) & !is.na(P$eco_l2_code), ]
gap_col <- c(onset = "onset_gap", end = "end_gap", up50 = "up50_gap", down50 = "down50_gap")
SD <- list()
for (ev in names(ic$events)) {
  e <- ic$events[[ev]]; dc <- paste0(e$date, "_doy"); gc <- gap_col[[e$date]]
  q <- P[P$pheno_class_id == e$pheno_class & !is.na(P[[dc]]) & !is.na(P[[gc]]) & P[[gc]] <= gap_ok, ]
  if (!nrow(q)) next
  key <- paste(q$eco_l2_code, q$group, q$site_id, q$year, sep = "|")
  med <- tapply(q[[dc]], key, stats::median); np <- tapply(q$individual_id, key, function(z) length(unique(z)))
  k <- do.call(rbind, strsplit(names(med), "|", fixed = TRUE))
  SD[[ev]] <- data.frame(eco_l2_code = k[, 1], group = k[, 2], event = ev, site_id = as.integer(k[, 3]), year = as.integer(k[, 4]),
                         doy = as.numeric(med), n_plants = as.integer(np[names(med)]), stringsAsFactors = FALSE)
}
SD <- do.call(rbind, SD); rownames(SD) <- NULL
invisible(nf_write_table(SD, file.path(out, "site_dates")))
log_msg(sprintf("site-year dates: %d (%s)", nrow(SD), paste(names(table(SD$event)), table(SD$event), collapse = ", ")))

# events: per eco x group x event
k3 <- paste(SD$eco_l2_code, SD$group, SD$event, sep = "|")
site_typ <- stats::aggregate(doy ~ eco_l2_code + group + event + site_id, SD, stats::median)
kt <- paste(site_typ$eco_l2_code, site_typ$group, site_typ$event, sep = "|")
u <- unique(k3); kk <- do.call(rbind, strsplit(u, "|", fixed = TRUE))
EV <- data.frame(eco_l2_code = kk[, 1], group = kk[, 2], event = kk[, 3], stringsAsFactors = FALSE)
EV$sites <- as.integer(tapply(SD$site_id, k3, function(z) length(unique(z)))[u])
EV$site_years <- as.integer(table(k3)[u])
spy <- tapply(seq_len(nrow(SD)), k3, function(i) table(SD$year[i]))
EV$years_ge5 <- vapply(u, function(x) sum(spy[[x]] >= 5), 0L)
EV$years_ge10 <- vapply(u, function(x) sum(spy[[x]] >= 10), 0L)
EV$typical_doy <- round(as.numeric(tapply(site_typ$doy, kt, stats::median)[u]))
EV$iqr_doy <- round(as.numeric(tapply(site_typ$doy, kt, stats::IQR)[u]))
EV$typical_week <- nf_week(as.Date(EV$typical_doy - 1, origin = "2001-01-01"))
utils::write.csv(EV, file.path(out, "events.csv"), row.names = FALSE)

# ---- fire in the matched fuel ---------------------------------------------------------------------------
ff <- nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))
if (is.na(ff)) stop("no fire/fire_weekly: run 03 and sync")
FW <- nf_read_table(ff); FW$eco_l2_code <- as.character(FW$eco_l2_code)
FW <- FW[FW$outcome %in% unlist(ic$outcomes), ]
NBt <- stats::aggregate(value ~ outcome + eco_l2_code, FW, sum)
NBn <- stats::aggregate(value ~ outcome + eco_l2_code, FW[FW$fuel %in% c("NB", "unmapped"), ], sum)
NBt$nonburnable_share <- round(NBn$value[match(paste(NBt$outcome, NBt$eco_l2_code), paste(NBn$outcome, NBn$eco_l2_code))] / NBt$value, 3)
NBt$nonburnable_share[is.na(NBt$nonburnable_share)] <- 0
utils::write.csv(NBt[, c("outcome", "eco_l2_code", "nonburnable_share")], file.path(out, "fire_nonburnable.csv"), row.names = FALSE)
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))           # years each record covers

# matched-fuel weekly series per eco x group x outcome x year: array over weeks
fire_week_mat <- function(eco, grp, o) {
  d <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% unlist(ic$match[[grp]]), ]
  yy <- oy[[o]]; M <- matrix(0, length(yy), 52, dimnames = list(yy, NULL))
  if (nrow(d)) { s <- tapply(d$value, list(factor(d$year, levels = yy), factor(d$week, levels = 1:52)), sum); s[is.na(s)] <- 0; M[] <- s }
  M
}
CV <- list()
for (r in seq_len(nrow(EV))) {
  e <- EV[r, ]
  win <- max(1, e$typical_week - ww):min(52, e$typical_week + ww)
  npn_years <- as.integer(names(spy[[paste(e$eco_l2_code, e$group, e$event, sep = "|")]]))
  for (o in unlist(ic$outcomes)) {
    M <- fire_week_mat(e$eco_l2_code, e$group, o)
    wsum <- rowSums(M[, win, drop = FALSE]); ov <- intersect(npn_years, as.integer(rownames(M)))
    CV[[length(CV) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o,
      overlap_years = length(ov), covered_years = sum(wsum[as.character(ov)] >= mev),
      window_per_year = round(mean(wsum), 2), year_per_year = round(mean(rowSums(M)), 2),
      window_share = round(sum(wsum) / max(sum(M), 1e-9), 3), stringsAsFactors = FALSE)
  }
}
CV <- do.call(rbind, CV)
utils::write.csv(CV, file.path(out, "coverage.csv"), row.names = FALSE)

# ---- lead outcome and tiers ---------------------------------------------------------------------------
ord <- unlist(ic$outcomes)
INV <- EV
kc <- paste(CV$eco_l2_code, CV$group, CV$event, sep = "|"); ke <- paste(INV$eco_l2_code, INV$group, INV$event, sep = "|")
lead <- tapply(seq_len(nrow(CV)), kc, function(i) { d <- CV[i, ]; d <- d[order(-d$covered_years, match(d$outcome, ord)), ]; d$outcome[1] })
INV$lead_outcome <- as.character(lead[ke])
lc <- CV[match(paste(ke, INV$lead_outcome), paste(kc, CV$outcome)), ]
INV$lead_covered_years <- lc$covered_years; INV$lead_window_per_year <- lc$window_per_year; INV$lead_window_share <- lc$window_share
t2 <- ic$tier2; t3 <- ic$tier3
INV$tier2_ok <- INV$sites >= t2$min_sites & INV$lead_covered_years >= t2$min_covered_years
# Tier 3: years with >= min_sites_year sites AND covered by the lead outcome
INV$tier3_years <- vapply(seq_len(nrow(INV)), function(r) {
  e <- INV[r, ]; tb <- spy[[ke[r]]]; yrs10 <- as.integer(names(tb)[tb >= t3$min_sites_year])
  if (!length(yrs10)) return(0L)
  M <- fire_week_mat(e$eco_l2_code, e$group, e$lead_outcome)
  win <- max(1, e$typical_week - ww):min(52, e$typical_week + ww)
  ws <- rowSums(M[, win, drop = FALSE]); sum(ws[as.character(intersect(yrs10, as.integer(rownames(M))))] >= mev)
}, 0L)
INV$tier3_ok <- INV$tier2_ok & INV$tier3_years >= t3$min_years
INV$fuel_match <- vapply(INV$group, function(g) paste(unlist(ic$match[[g]]), collapse = "+"), "")
INV <- INV[order(INV$eco_l2_code, INV$group, match(INV$event, names(ic$events))), ]
utils::write.csv(INV, file.path(out, "inventory.csv"), row.names = FALSE)

log_msg(sprintf("inventory: %d eco x group x event rows | Tier 2 possible: %d | Tier 3 possible: %d", nrow(INV), sum(INV$tier2_ok), sum(INV$tier3_ok)))
tl <- table(INV$lead_outcome[INV$tier2_ok]); log_msg("lead outcome where Tier 2 is possible: ", paste(names(tl), tl, collapse = ", "))
for (g in names(ic$match)) { s <- INV[INV$group == g, ]
  log_msg(sprintf("  %-16s Tier 2 in %2d ecoregion x event rows (%2d ecoregions); Tier 3 in %2d (%s)", g, sum(s$tier2_ok),
                  length(unique(s$eco_l2_code[s$tier2_ok])), sum(s$tier3_ok),
                  paste(unique(paste(s$eco_l2_code, s$event)[s$tier3_ok]), collapse = ", "))) }
log_msg("=== done | ", out, " | log: ", log_file)
