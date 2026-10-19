#!/usr/bin/env Rscript
# Step 66b (fire activity study; laptop, seconds): ground check of the satellite phenology at the
# scale the fire test uses (decided Oct 19). The regions work's target check (36) holds, per NPN grass
# site-year, the ground date (grass green-up, curing = 50 % green going down, cured) and the satellite
# date at the site's cell, each as a departure from its own site median (g_anom, s_anom; season
# mismatches out). Site-year agreement is limited by plant-to-plant noise; the fire test uses regional
# shares, so here departures are also averaged over the sites of a region in a year and compared again.
# Input: <run_dir>/eval/target_check/site_years.csv (36). Output fire_study/test/ground_check.csv:
# set x metric x scale (site-year / region-year with >= 2 or >= 3 sites) x region (all + each with
# >= 5 region-years): n, r.
#
# Usage: Rscript scripts/66b_fire_ground_check.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "fire.R"))
cfg <- load_cfg(profile)
log_file <- log_init("66b_fire_ground_check", profile)
f <- file.path(cfg$run_dir, "eval", "target_check", "site_years.csv")
if (!file.exists(f)) stop("no ", f, " (run 36_curing_target_check.R)")
D <- utils::read.csv(f, stringsAsFactors = FALSE)
D <- D[!is.na(D$g_anom) & !is.na(D$s_anom), ]
pairs <- list(c("grass_curing", "eos50"), c("grass_curing", "cure50"), c("grass_cured", "eos20"), c("grass_greenup", "sos20"))
r_of <- function(a, b) if (length(a) >= 5 && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_
OUT <- list()
for (p in pairs) {
  x <- D[D$set == p[1] & D$metric == p[2], ]
  if (!nrow(x)) next
  OUT[[length(OUT) + 1]] <- data.frame(set = p[1], metric = p[2], scale = "site-year", region = "all", n = nrow(x), r = r_of(x$g_anom, x$s_anom))
  ry <- stats::aggregate(cbind(g = x$g_anom, s = x$s_anom) ~ region + year, x, mean)
  ry$sites <- stats::aggregate(site ~ region + year, x, function(v) length(unique(v)))$site
  for (k in c(2, 3)) {
    z <- ry[ry$sites >= k, ]
    OUT[[length(OUT) + 1]] <- data.frame(set = p[1], metric = p[2], scale = sprintf("region-year, >= %d sites", k), region = "all", n = nrow(z), r = r_of(z$g, z$s))
    if (k == 2) for (rg in sort(unique(z$region))) { zz <- z[z$region == rg, ]
      if (nrow(zz) >= 5) OUT[[length(OUT) + 1]] <- data.frame(set = p[1], metric = p[2], scale = "region-year, >= 2 sites", region = as.character(rg), n = nrow(zz), r = r_of(zz$g, zz$s)) }
  }
}
O <- do.call(rbind, OUT)
dir.create(fire_dir(cfg, "test"), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(O, fire_dir(cfg, "test", "ground_check.csv"), row.names = FALSE)
for (i in which(O$region == "all")) with(O[i, ], log_msg(sprintf("%-13s %-6s %-26s n %4d  r %.2f", set, metric, scale, n, r)))
log_msg("=== done | log: ", log_file)
