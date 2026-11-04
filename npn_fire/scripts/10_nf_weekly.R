#!/usr/bin/env Rscript
# Step 10 (NPN x fire study; P720, ~1 h; arrow): weekly status and intensity tables from the raw NPN records,
# for the intensity tests (12) and machine learning (11). Roadmap phase 1 (Mike, 3 Nov 2026; project doc
# claude/gsi-npn-fire-roadmap.md).
# Inputs: _shared/raw/npn/status/<Y>/*.parquet (01), lookups (phenophases), data/nf_species_groups.csv (02a),
#   clean/sites.csv (02b: Level II ecoregion per site; bad-coordinate sites have none).
# Per year (resumable: weekly/years/<Y>_*.parquet skipped if present; `force` redoes):
#   - plant pheno classes 1-13, status 0 / 1, same-day repeat visits merged (yes if any; highest intensity), as 02b;
#     percent classes: "no" = 0 %
#   - plant x pheno class x week: status = yes if any of the class's phenophases was yes in any visit that
#     week; pct = highest percent reading that week
#   - site x group x class x week (`site_weekly`, P720 / laptop only, never published): plants, share of
#     plants yes, plants with a percent reading, their mean percent
#   - ecoregion x group x class x year x week (`eco_weekly`): sites, mean over sites of the share yes (sites
#     weighted equally), sites with a percent reading, mean over those sites of the mean percent
# Anomalies (value minus that week's mean over exploration years) are computed by 11 / 12, not here.
# Outputs <run_dir>/npn_fire/weekly/: eco_weekly.parquet, site_weekly.parquet, weekly_summary.csv.
# Usage: Rscript npn_fire/scripts/10_nf_weekly.R [profile] [years] [force]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && !grepl("^[0-9]|^force$", args[1])) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
rest <- if (length(args) && identical(args[1], profile)) args[-1] else args
force <- "force" %in% rest; rest <- setdiff(rest, "force")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))

cfg <- load_nf_cfg(profile)
log_file <- log_init("10_nf_weekly", profile)
nf <- cfg$npn_fire; cc <- nf$clean
years <- nf_parse_years(rest[grepl("^[0-9]", rest)][1]) %||% seq(nf$years$start, nf$years$end)
out <- nf_out(cfg, "weekly"); ydir <- nf_out(cfg, "weekly", "years")
t0 <- Sys.time()
log_msg("=== 10_nf_weekly | profile: ", profile, " | years ", paste(range(years), collapse = "-"), " | force: ", force)

lk <- file.path(cfg$raw_dir, "npn", "lookups")
PH <- utils::read.csv(file.path(lk, "phenophases.csv"), stringsAsFactors = FALSE)
SG <- utils::read.csv(nf_path("data", "nf_species_groups.csv"), stringsAsFactors = FALSE)
SX <- utils::read.csv(file.path(cfg$out_dir, "clean", "sites.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
plant_classes <- unlist(cc$pheno_classes %||% 1:13)

read_year <- function(Y) {
  fs <- list.files(file.path(cfg$raw_dir, "npn", "status", Y), "^npn_status_.*\\.(parquet|csv\\.gz)$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  keep <- c("site_id", "species_id", "individual_id", "phenophase_id", "observation_date", "phenophase_status", "intensity_value")
  do.call(rbind, lapply(fs, function(f) {
    d <- nf_read_table(f)
    for (k in setdiff(keep, names(d))) d[[k]] <- NA
    d <- d[, keep]; d$intensity_value <- as.character(d$intensity_value); d$observation_date <- as.character(d$observation_date)
    d
  }))
}

do_year <- function(Y) {
  f_s <- file.path(ydir, sprintf("site_weekly_%d", Y)); f_e <- file.path(ydir, sprintf("eco_weekly_%d", Y))
  if (!force && !is.na(nf_table_path(f_s)) && !is.na(nf_table_path(f_e))) return(paste(Y, "on disk"))
  t1 <- Sys.time()
  d <- read_year(Y)
  if (is.null(d) || !nrow(d)) return(paste(Y, "no files"))
  d$pheno_class_id <- PH$pheno_class_id[match(d$phenophase_id, PH$phenophase_id)]
  d <- d[d$pheno_class_id %in% plant_classes & d$species_id %in% SG$species_id & d$phenophase_status %in% c(0, 1), ]
  d$group <- SG$group[match(d$species_id, SG$species_id)]
  d$eco_l2_code <- SX$eco_l2_code[match(d$site_id, SX$site_id)]
  d <- d[!is.na(d$eco_l2_code) & !is.na(d$group), ]
  d$date <- as.Date(d$observation_date); d$week <- nf_week(d$date)
  iv <- unique(d$intensity_value[!is.na(d$intensity_value)])
  pm <- nf_pct_mid(iv); pt <- nf_intensity_type(iv) == "percent"
  d$pct <- ifelse(d$intensity_value %in% iv[pt], pm[match(d$intensity_value, iv)], NA_real_)
  pct_ph <- unique(d$phenophase_id[!is.na(d$pct)])
  d$pct[d$phenophase_id %in% pct_ph & d$phenophase_status == 0] <- 0          # "no" = 0 % for percent phenophases
  # plant x class x week: yes if any yes; highest percent
  k <- paste(d$individual_id, d$pheno_class_id, d$week)
  st <- tapply(d$phenophase_status, k, max)
  pc <- suppressWarnings(tapply(d$pct, k, max, na.rm = TRUE)); pc[!is.finite(pc)] <- NA
  first <- !duplicated(k)
  PW <- data.frame(site_id = d$site_id[first], eco_l2_code = d$eco_l2_code[first], group = d$group[first],
                   pheno_class_id = d$pheno_class_id[first], week = d$week[first], key = k[first], stringsAsFactors = FALSE)
  PW$yes <- as.numeric(st[PW$key]); PW$pct <- as.numeric(pc[PW$key]); PW$key <- NULL
  # site x group x class x week
  ks <- paste(PW$site_id, PW$group, PW$pheno_class_id, PW$week)
  f1 <- !duplicated(ks)
  SW <- data.frame(year = Y, site_id = PW$site_id[f1], eco_l2_code = PW$eco_l2_code[f1], group = PW$group[f1],
                   pheno_class_id = PW$pheno_class_id[f1], week = PW$week[f1], stringsAsFactors = FALSE)
  kk <- ks[f1]
  SW$plants <- as.integer(table(ks)[kk]); SW$share_yes <- round(as.numeric(tapply(PW$yes, ks, mean)[kk]), 4)
  SW$plants_pct <- as.integer(tapply(!is.na(PW$pct), ks, sum)[kk])
  mp <- suppressWarnings(tapply(PW$pct, ks, mean, na.rm = TRUE)[kk]); mp[!is.finite(mp)] <- NA
  SW$mean_pct <- round(as.numeric(mp), 2)
  # ecoregion x group x class x week (sites weighted equally)
  ke <- paste(SW$eco_l2_code, SW$group, SW$pheno_class_id, SW$week)
  f2 <- !duplicated(ke); ke1 <- ke[f2]
  EW <- data.frame(year = Y, eco_l2_code = SW$eco_l2_code[f2], group = SW$group[f2], pheno_class_id = SW$pheno_class_id[f2],
                   week = SW$week[f2], stringsAsFactors = FALSE)
  EW$sites <- as.integer(table(ke)[ke1]); EW$plants <- as.integer(tapply(SW$plants, ke, sum)[ke1])
  EW$share_yes <- round(as.numeric(tapply(SW$share_yes, ke, mean)[ke1]), 4)
  hp <- !is.na(SW$mean_pct)
  EW$sites_pct <- as.integer(tapply(hp, ke, sum)[ke1])
  mp2 <- suppressWarnings(tapply(SW$mean_pct, ke, mean, na.rm = TRUE)[ke1]); mp2[!is.finite(mp2)] <- NA
  EW$mean_pct <- round(as.numeric(mp2), 2)
  invisible(nf_write_table(SW, f_s)); invisible(nf_write_table(EW, f_e))
  sprintf("%d: %d records -> %d plant-class-weeks | %d site rows | %d ecoregion rows | %.0f s", Y, nrow(d), nrow(PW), nrow(SW), nrow(EW),
          as.numeric(difftime(Sys.time(), t1, units = "secs")))
}
nw <- if (.Platform$OS.type == "unix") as.integer(nf$weekly$workers %||% 8) else 1L
run1 <- function(Y) tryCatch(do_year(Y), error = function(e) paste("ERROR", Y, conditionMessage(e)))
res <- if (nw > 1) parallel::mclapply(years, run1, mc.cores = min(nw, length(years)), mc.preschedule = FALSE) else lapply(years, run1)
for (m in unlist(res)) if (grepl("^ERROR", m)) log_err(m) else log_msg("  ", m)
if (any(grepl("^ERROR", unlist(res)))) stop("year(s) failed; rerun to retry (finished years are skipped)")

all_years <- seq(nf$years$start, nf$years$end)
rd <- function(stem) { f <- stats::na.omit(vapply(all_years, function(Y) nf_table_path(file.path(ydir, sprintf("%s_%d", stem, Y))), ""))
  do.call(rbind, lapply(f, nf_read_table)) }
EW <- rd("eco_weekly"); SW <- rd("site_weekly")
invisible(nf_write_table(EW, file.path(out, "eco_weekly"))); invisible(nf_write_table(SW, file.path(out, "site_weekly")))
SM <- stats::aggregate(cbind(sites, sites_pct) ~ group + pheno_class_id + year, EW, sum)
utils::write.csv(SM, file.path(out, "weekly_summary.csv"), row.names = FALSE)
log_msg(sprintf("eco_weekly: %d rows | site_weekly: %d rows | ecoregion-weeks with >= 5 sites: %d (%d with >= 5 percent sites)",
                nrow(EW), nrow(SW), sum(EW$sites >= 5), sum(EW$sites_pct >= 5)))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
