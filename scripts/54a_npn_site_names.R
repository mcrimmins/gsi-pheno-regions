#!/usr/bin/env Rscript
# Step 54a (NPN grass curing study, laptop with internet): NPN site names for the managed-site
# check in 54. Looks up every site in npn_curing_metrics.rds with rnpn::npn_stations() and writes
# data/npn_site_names.csv (site_id, site_name, state). Names only; nothing personal is kept.
#
# Usage: Rscript scripts/54a_npn_site_names.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R"))
source(here::here("R", "obs.R"))
if (!requireNamespace("rnpn", quietly = TRUE)) stop("needs the rnpn package: renv::install('rnpn')")
cfg <- load_cfg(profile)
cm <- readRDS(file.path(Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg)), "npn_curing_metrics.rds"))
S <- unique(cm[cm$usable %in% TRUE, c("site_id", "state")])
st <- do.call(rbind, lapply(sort(unique(S$state)), function(s) {
  x <- tryCatch(rnpn::npn_stations(state_code = s), error = function(e) { message("  ", s, ": ", conditionMessage(e)); NULL })
  if (is.null(x) || !nrow(x)) return(NULL)
  idc <- grep("station_id|site_id", names(x), value = TRUE)[1]; nmc <- grep("station_name|site_name", names(x), value = TRUE)[1]
  data.frame(site_id = as.integer(x[[idc]]), site_name = as.character(x[[nmc]]), state = s)
}))
out <- merge(S["site_id"], st[!duplicated(st$site_id), ], by = "site_id", all.x = TRUE)
utils::write.csv(out, here::here("data", "npn_site_names.csv"), row.names = FALSE)
message(sprintf("names found for %d of %d sites -> data/npn_site_names.csv", sum(!is.na(out$site_name)), nrow(out)))
