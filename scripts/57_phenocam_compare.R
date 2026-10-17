#!/usr/bin/env Rscript
# Step 57 (NPN grass curing study, exploratory): is the NPN site-year signal a real vegetation
# signal? Compares PhenoCam GR transition-date departures (56) with the shared site-year signal of
# the NPN dates (54) and with PRISM weather departures, on the same site-years.
#   camera events   per camera ROI and year, the rising (green-up) / falling (green-down) 50 %
#                   transition of gcc_90 nearest the site's typical NPN green-up / curing date
#                   (within +/- window_days); departure = date minus the camera's circular median
#                   over its years (>= 3 years needed)
#   one ROI per site  code match before distance match, then the most years overlapping the NPN
#                   signal; where a site has two cameras, their agreement is reported (camera noise)
#   comparisons     camera green-up vs NPN green-up signal; camera green-down vs NPN curing and
#                   cured signals. For each: r(camera, NPN), r(camera, PRISM), r(NPN, PRISM) on the
#                   same site-years, PRISM = a priori drivers (green-up: t_marapr, gdd5_apr30;
#                   green-down: vpd_jja, p_amj, events_amjjas; departures from 54's signal.csv).
#                   95 % intervals from a site-block bootstrap (sites resampled), n is small.
# Reading the result: camera ~ NPN with PRISM missing both -> real local signal PRISM misses;
# camera ~ PRISM but not NPN -> NPN site-year signal mostly crew / schedule; neither -> inconclusive.
# Inputs: <run_dir>/npn_study/phenocam/ (56), site_signal/signal.csv (54), npn_curing_metrics.rds.
# Output <run_dir>/npn_study/phenocam/: camera_dates.csv, pairs.csv, compare.csv, camera_agreement.csv.
#
# Usage: Rscript scripts/57_phenocam_compare.R [profile] [window_days, default 60]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "dev")
win_d <- if (length(args) >= 2) as.numeric(args[2]) else 60

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "obs.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("57_phenocam_compare", profile)
ns <- cfg$npn_study %||% list()
sdir <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(sdir, "phenocam"); rdir <- file.path(odir, "raw")
t0 <- Sys.time()
log_msg("=== 57_phenocam_compare | profile: ", profile, " | window: ", win_d, " d")
set.seed(57)

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }

MT <- utils::read.csv(file.path(odir, "matches.csv"), stringsAsFactors = FALSE)
SG <- utils::read.csv(file.path(sdir, "site_signal", "signal.csv"), stringsAsFactors = FALSE)
SG <- SG[SG$site_id %in% MT$site_id, ]

# typical NPN dates per site (same plant-year filter as 53 / 54)
cm <- readRDS(file.path(Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg)), "npn_curing_metrics.rds"))
u <- cm[cm$usable %in% TRUE & cm$n_obs >= (ns$min_readings %||% 12) & cm$site_id %in% MT$site_id, ]
typ <- lapply(c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy"), function(col) c(tapply(u[[col]], u$site_id, cmed)))

# transition dates
tf <- list.files(rdir, "transition_dates\\.csv$", full.names = TRUE)
if (!length(tf)) stop("no transition files in ", rdir, ": run 56_phenocam_fetch.R")
TR <- do.call(rbind, lapply(tf, function(f) {
  x <- utils::read.csv(f, comment.char = "#", stringsAsFactors = FALSE)
  names(x) <- tolower(names(x))
  if ("gcc_value" %in% names(x)) x <- x[x$gcc_value == "gcc_90", ]
  data.frame(phenocam_site = x$site, roi_id = as.integer(x$roi_id), direction = x$direction, date = as.Date(x$transition_50))
}))
TR <- TR[!is.na(TR$date), ]
TR$year <- as.integer(format(TR$date, "%Y")); TR$doy <- as.integer(format(TR$date, "%j"))
log_msg(sprintf("transitions: %d (rising %d, falling %d) on %d ROIs", nrow(TR), sum(TR$direction == "rising"), sum(TR$direction == "falling"),
                length(unique(paste(TR$phenocam_site, TR$roi_id)))))

# camera event per site x ROI x year x target, nearest the typical NPN date
pairs_def <- list(greenup = "rising", curing = "falling", cured = "falling")
ev <- list()
for (k in seq_len(nrow(MT))) for (tg in names(pairs_def)) {
  td <- typ[[tg]][as.character(MT$site_id[k])]; if (is.na(td)) next
  x <- TR[TR$phenocam_site == MT$phenocam_site[k] & TR$roi_id == MT$roi_id[k] & TR$direction == pairs_def[[tg]], ]
  if (!nrow(x)) next
  x$off <- cdiff(x$doy, td); x <- x[abs(x$off) <= win_d, ]
  if (!nrow(x)) next
  x <- x[order(x$year, abs(x$off)), ]; x <- x[!duplicated(x$year), ]
  if (nrow(x) < 3) next
  med <- cmed(x$doy)
  ev[[length(ev) + 1]] <- data.frame(site_id = MT$site_id[k], neon_code = MT$neon_code[k], match = MT$match[k], camera = MT$camera[k],
    phenocam_site = MT$phenocam_site[k], roi_id = MT$roi_id[k], target = tg, year = x$year, cam_doy = x$doy, cam_dep = cdiff(x$doy, med))
}
CE <- do.call(rbind, ev)
utils::write.csv(CE, file.path(odir, "camera_dates.csv"), row.names = FALSE)

# camera-to-camera agreement where a site has two ROIs
ag <- list()
for (tg in names(pairs_def)) for (s in unique(CE$site_id)) {
  x <- CE[CE$site_id == s & CE$target == tg, ]; rk <- unique(paste(x$phenocam_site, x$roi_id))
  if (length(rk) < 2) next
  w <- stats::reshape(data.frame(year = x$year, k = paste(x$phenocam_site, x$roi_id), d = x$cam_dep), idvar = "year", timevar = "k", direction = "wide")
  w <- w[stats::complete.cases(w[, 2:3]), ]
  if (nrow(w) >= 3) ag[[length(ag) + 1]] <- data.frame(target = tg, site_id = s, rois = paste(rk[1:2], collapse = " | "), n_years = nrow(w),
                                                         r = stats::cor(w[[2]], w[[3]]), mad_diff = stats::median(abs(w[[2]] - w[[3]])))
}
if (length(ag)) {
  AG <- do.call(rbind, ag); utils::write.csv(AG, file.path(odir, "camera_agreement.csv"), row.names = FALSE)
  for (tg in unique(AG$target)) { a <- AG[AG$target == tg, ]
    log_msg(sprintf("camera agreement %-8s: %d sites, median r %.2f, median |difference| %.0f d", tg, nrow(a), stats::median(a$r), stats::median(a$mad_diff))) }
}

# one ROI per site and target: code match first, then most years overlapping the NPN signal
sel <- list()
for (tg in names(pairs_def)) for (s in unique(CE$site_id)) {
  x <- CE[CE$site_id == s & CE$target == tg, ]; if (!nrow(x)) next
  sy <- SG$year[SG$site_id == s & SG$target == tg]
  sc <- stats::aggregate(list(n = x$year %in% sy), list(phenocam_site = x$phenocam_site, roi_id = x$roi_id, match = x$match), sum)
  sc <- sc[order(sc$match != "code", -sc$n), ][1, ]
  sel[[length(sel) + 1]] <- x[x$phenocam_site == sc$phenocam_site & x$roi_id == sc$roi_id, ]
}
CS <- do.call(rbind, sel)
P <- merge(CS, SG, by = c("site_id", "year", "target"))
drivers <- list(greenup = c("t_marapr", "gdd5_apr30"), curing = c("vpd_jja", "p_amj", "events_amjjas"), cured = c("vpd_jja", "p_amj", "events_amjjas"))
keepc <- unique(c("site_id", "neon_code", "camera", "phenocam_site", "roi_id", "target", "year", "cam_doy", "cam_dep", "signal", unlist(drivers)))
P <- P[, intersect(keepc, names(P))]
utils::write.csv(P, file.path(odir, "pairs.csv"), row.names = FALSE)

# correlations with a site-block bootstrap
boot_r <- function(a, b, site, nb = 2000) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]; site <- site[ok]
  if (length(a) < 8 || length(unique(site)) < 3) return(c(n = length(a), sites = length(unique(site)), r = NA, lo = NA, hi = NA))
  us <- unique(site); idx <- split(seq_along(site), site)
  rb <- replicate(nb, { s <- sample(us, replace = TRUE); i <- unlist(idx[as.character(s)]); suppressWarnings(stats::cor(a[i], b[i])) })
  q <- stats::quantile(rb, c(0.025, 0.975), na.rm = TRUE)
  c(n = length(a), sites = length(us), r = stats::cor(a, b), lo = unname(q[1]), hi = unname(q[2]))
}
CR <- list()
for (tg in names(pairs_def)) {
  x <- P[P$target == tg, ]; if (nrow(x) < 8) { log_msg(sprintf("%s: only %d paired site-years, skipped", tg, nrow(x))); next }
  add <- function(what, a, b) { z <- boot_r(a, b, x$site_id)
    CR[[length(CR) + 1]] <<- data.frame(target = tg, comparison = what, n = z[["n"]], sites = z[["sites"]], r = z[["r"]], lo95 = z[["lo"]], hi95 = z[["hi"]]) }
  add("camera ~ NPN signal", x$cam_dep, x$signal)
  for (v in intersect(drivers[[tg]], names(x))) { add(paste("camera ~", v), x$cam_dep, x[[v]]); add(paste("NPN signal ~", v), x$signal, x[[v]]) }
  ss <- c(stats::sd(x$cam_dep), stats::sd(x$signal))
  b <- stats::coef(stats::lm(signal ~ cam_dep, data = x))[2]
  log_msg(sprintf("%-8s %d site-years at %d sites | SD camera %.0f d, NPN signal %.0f d | slope NPN on camera %.2f", tg, nrow(x), length(unique(x$site_id)), ss[1], ss[2], b))
}
CRd <- do.call(rbind, CR)
utils::write.csv(CRd, file.path(odir, "compare.csv"), row.names = FALSE)
for (i in seq_len(nrow(CRd))) log_msg(sprintf("  %-8s %-28s n %3d (%2d sites) r %5.2f [%5.2f, %5.2f]",
  CRd$target[i], CRd$comparison[i], CRd$n[i], CRd$sites[i], CRd$r[i], CRd$lo95[i], CRd$hi95[i]))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
