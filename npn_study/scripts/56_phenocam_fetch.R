#!/usr/bin/env Rscript
# Step 56 (NPN grass curing study, exploratory; laptop with internet, phenocamr): PhenoCam
# grassland greenness at the NEON sites in the NPN grass record, as an observer-free check of the
# site-year signal (54). PhenoCam is blocked from the cloud workspace, so this runs on the laptop.
#   1. NEON NPN sites: data/npn_site_names.csv (54a; NEON = "*.phenology.phe*") and the site
#      coordinates in npn_curing_metrics.rds. NEON site code = the 4-letter upper-case token in
#      the NPN site name.
#   2. PhenoCam ROIs (phenocamr::list_rois): NEON cameras ("NEON.Dxx.CODE.DP1.000nn"), veg type
#      GR (grassland). Matched on site code; sites without a code match take the nearest NEON
#      camera with a GR ROI within max_km. Both tower (DP1.00033) and understory (DP1.00042)
#      cameras are kept; 57 picks one per site.
#   3. download_phenocam(frequency 3, outlier detection, smoothing, phenophase = TRUE): 3-day
#      GCC and transition dates (rising = green-up, falling = green-down; 10 / 25 / 50 %).
# Nothing personal is written: NEON site codes and camera names only.
# Output <run_dir>/npn_study/phenocam/: matches.csv, rois_all_neon.csv, raw/ (phenocamr files).
#
# Usage: Rscript npn_study/scripts/56_phenocam_fetch.R [profile] [max_km, default 5]
# Needs: renv::install("phenocamr")

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
max_km <- if (length(args) >= 2) as.numeric(args[2]) else 5

npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("phenocamr", quietly = TRUE)) stop("needs the phenocamr package: renv::install('phenocamr')")

cfg <- load_npn_cfg(profile)
log_file <- log_init("56_phenocam_fetch", profile)
odir <- file.path(cfg$run_dir, "npn_study", "phenocam"); rdir <- file.path(odir, "raw")
dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 56_phenocam_fetch | profile: ", profile, " | max_km: ", max_km)

# 1. NEON sites in the NPN grass record
nf <- npn_path("data", "npn_site_names.csv")
if (!file.exists(nf)) stop("needs data/npn_site_names.csv: run scripts/54a_npn_site_names.R first")
sn <- utils::read.csv(nf, stringsAsFactors = FALSE)
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
u <- cm[cm$usable %in% TRUE, ]
S <- stats::aggregate(cbind(lat, lon) ~ site_id, data = u, FUN = stats::median)
S$plant_years <- as.integer(table(u$site_id)[as.character(S$site_id)])
S$site_name <- sn$site_name[match(S$site_id, sn$site_id)]
S <- S[grepl("\\.phenology\\.phe", tolower(S$site_name)), ]
tok <- regmatches(S$site_name, regexpr("(?<![A-Za-z])[A-Z]{4}(?![A-Za-z])", S$site_name, perl = TRUE))
S$neon_code <- NA_character_
S$neon_code[grepl("(?<![A-Za-z])[A-Z]{4}(?![A-Za-z])", S$site_name, perl = TRUE)] <- tok
S$site_name <- NULL                                   # nothing beyond the NEON code is kept
log_msg(sprintf("NEON NPN sites: %d (%d plant-years); code found for %d: %s", nrow(S), sum(S$plant_years), sum(!is.na(S$neon_code)),
                paste(sort(unique(stats::na.omit(S$neon_code))), collapse = " ")))

# 2. NEON PhenoCam ROIs
rois <- phenocamr::list_rois()
names(rois) <- tolower(names(rois))
sc <- grep("^site$", names(rois), value = TRUE)[1]
rois <- rois[grepl("^NEON\\.", rois[[sc]]), ]
rois$neon_code <- sub("^NEON\\.D\\d+\\.([A-Z0-9]{4})\\..*$", "\\1", rois[[sc]])
rois$camera <- ifelse(grepl("00042", rois[[sc]]), "understory", ifelse(grepl("00033", rois[[sc]]), "tower", "other"))
utils::write.csv(rois, file.path(odir, "rois_all_neon.csv"), row.names = FALSE)
gr <- rois[toupper(rois$veg_type) == "GR", ]
log_msg(sprintf("NEON PhenoCam ROIs: %d (%d sites); GR ROIs: %d at %d NEON codes", nrow(rois), length(unique(rois[[sc]])),
                nrow(gr), length(unique(gr$neon_code))))

hav <- function(la1, lo1, la2, lo2) { r <- pi / 180; a <- sin((la2 - la1) * r / 2)^2 + cos(la1 * r) * cos(la2 * r) * sin((lo2 - lo1) * r / 2)^2; 12742 * asin(sqrt(a)) }
M <- list()
for (i in seq_len(nrow(S))) {
  g <- gr[gr$neon_code %in% S$neon_code[i], ]; how <- "code"
  if (!nrow(g)) {
    d <- hav(S$lat[i], S$lon[i], gr$lat, gr$lon); g <- gr[d <= max_km, ]; how <- "distance"
  }
  if (!nrow(g)) next
  g$km <- hav(S$lat[i], S$lon[i], g$lat, g$lon)
  M[[length(M) + 1]] <- data.frame(site_id = S$site_id[i], neon_code = S$neon_code[i], plant_years = S$plant_years[i], match = how,
    phenocam_site = g[[sc]], camera = g$camera, roi_id = g$roi_id_number, km = round(g$km, 2),
    first_date = as.character(g$first_date), last_date = as.character(g$last_date), site_years = g$site_years,
    missing_pct = if ("missing_data_pct" %in% names(g)) g$missing_data_pct else NA)
}
MT <- if (length(M)) do.call(rbind, M) else stop("no NEON GR ROIs matched any NEON NPN site")
utils::write.csv(MT, file.path(odir, "matches.csv"), row.names = FALSE)
nomatch <- setdiff(S$site_id, MT$site_id)
log_msg(sprintf("matched %d of %d NEON NPN sites (%d plant-years of %d) to %d GR ROIs on %d cameras; unmatched codes: %s",
                length(unique(MT$site_id)), nrow(S), sum(S$plant_years[S$site_id %in% MT$site_id]), sum(S$plant_years), nrow(MT),
                length(unique(MT$phenocam_site)), paste(sort(stats::na.omit(S$neon_code[S$site_id %in% nomatch])), collapse = " ")))
far <- MT[MT$km > max_km, ]
if (nrow(far)) log_msg(sprintf("  note: %d code matches more than %g km from the NPN site (max %.1f km)", nrow(far), max_km, max(far$km)))

# 3. download 3-day GCC and transition dates
todo <- unique(MT[, c("phenocam_site", "roi_id")])
for (j in seq_len(nrow(todo))) {
  s <- todo$phenocam_site[j]; id <- todo$roi_id[j]
  ok <- tryCatch({
    phenocamr::download_phenocam(site = paste0(s, "$"), veg_type = "GR", roi_id = id, frequency = "3",
                                 outlier_detection = TRUE, smooth = TRUE, phenophase = TRUE, out_dir = rdir)
    TRUE
  }, error = function(e) { log_msg(sprintf("  %s GR %s: %s", s, id, conditionMessage(e))); FALSE })
  if (ok) log_msg(sprintf("  %d/%d %s GR %s", j, nrow(todo), s, id))
}
tf <- list.files(rdir, "transition_dates\\.csv$")
log_msg(sprintf("transition files: %d | GCC files: %d", length(tf), length(list.files(rdir, "_3day\\.csv$"))))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
