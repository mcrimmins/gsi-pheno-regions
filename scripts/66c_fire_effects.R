#!/usr/bin/env Rscript
# Step 66c (fire study, Oct 26 2026; laptop, seconds): phenology effects in usable units. 66's held-out log-score
# gains say whether the cured share helps; this says how much. From the full-data fits (66 coefs.csv, model M2 =
# season + ppet60 + cured share; M2d with dfm1000), per region and outcome:
#   pct_per10   % change in the outcome for a week 10 percentage points more cured than usual (exp(b) - 1;
#               counts: NB, so % more fires; growth: Gaussian on log(1 + pixel-days), so % more 1 + pixel-days)
#   dep_typ     typical size of the departure: median |cured share - typical| in the region's transition weeks
#               (typical share 0.1-0.9), in percentage points, from the same cured-share series the spec used
#   pct_typ     % change for a departure of that typical size
#   lo / hi     90 % Wald intervals (b +/- 1.645 se); ppet60_pct = % change per 100 mm more 60-day rain minus PET,
#               for scale
# Reads <run_dir>/fire_study/test/coefs.csv and test/<sub>/coefs.csv (66 runs with GSI_FIRE_TEST_SUB), weekly
# fire_weekly.csv, realtime/realtime_weekly*.csv -> fire_study/test/effects.csv.
# Usage: Rscript scripts/66c_fire_effects.R [profile]
args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")
source(here::here("R", "config.R")); source(here::here("R", "log.R")); source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
cfg <- load_cfg(profile); log_file <- log_init("66c_fire_effects", profile)
specs <- c("all", "human_unk", "px", "px_rt", "px_rt4", "px_rt0", "all_vd", "human_unk_vd", "px_vd", "px_vdx")
td <- fire_dir(cfg, "test")
cf <- c(main = file.path(td, "coefs.csv"), stats::setNames(file.path(list.dirs(td, recursive = FALSE), "coefs.csv"), basename(list.dirs(td, recursive = FALSE))))
cf <- cf[file.exists(cf)]
W <- utils::read.csv(fire_dir(cfg, "weekly", "fire_weekly.csv"))
dep_series <- function(spec) {          # cured share and typical share behind a spec's phenology term
  suf <- sub("^(all|human_unk|px)_?", "", spec)
  if (suf == "") return(data.frame(region = W$region, s = W$cured_share, t = W$cured_share_typ))
  f <- if (grepl("^vd", suf)) "realtime_weekly_vd.csv" else if (grepl("^rt[0-9]+$", suf)) sprintf("realtime_weekly_lag%s.csv", sub("rt", "", suf)) else "realtime_weekly.csv"
  f <- fire_dir(cfg, "realtime", f); if (!file.exists(f)) return(NULL)
  R <- utils::read.csv(f); k <- sub("[0-9]+$", "", suf); if (k == "rt" && grepl("[0-9]", suf)) k <- "rt"
  data.frame(region = R$region, s = R[[paste0("cured_share_", k)]], t = R[[paste0("cured_share_", k, "_typ")]])
}
OUT <- list()
for (src in names(cf)) {
  co <- utils::read.csv(cf[[src]])
  for (sp in intersect(specs, unique(co$spec))) for (m in c("M2", "M2d")) {
    x <- co[co$spec == sp & co$model == m, ]; ph <- x[grepl("cured_dep", x$term), ]; if (!nrow(ph)) next
    pp <- x[grepl(":ppet$", x$term), ]
    ds <- dep_series(sp)
    for (k in seq_len(nrow(ph))) {
      r <- as.integer(sub("^region_f([0-9]+):.*", "\\1", ph$term[k]))
      d <- if (is.null(ds)) numeric() else { z <- ds[ds$region == r & ds$t > 0.1 & ds$t < 0.9, ]; abs(z$s - z$t) }
      dt <- if (length(d)) 100 * stats::median(d, na.rm = TRUE) else NA_real_
      b <- ph$est[k]; se <- ph$se[k]; pr <- pp[pp$term == sprintf("region_f%d:ppet", r), ]
      OUT[[length(OUT) + 1]] <- data.frame(source = src, spec = sp, model = m, region = r, b = b, se = se,
        pct_per10 = 100 * (exp(b) - 1), lo = 100 * (exp(b - 1.645 * se) - 1), hi = 100 * (exp(b + 1.645 * se) - 1),
        dep_typ = dt, pct_typ = 100 * (exp(b * dt / 10) - 1),
        ppet60_pct = if (nrow(pr)) 100 * (exp(pr$est[1]) - 1) else NA_real_)
    }
  }
}
E <- do.call(rbind, OUT); out <- fire_dir(cfg, "test", "effects.csv"); utils::write.csv(E, out, row.names = FALSE)
for (k in which(E$model == "M2")) with(E[k, ], log_msg(sprintf("%-5s %-13s region %2d: %+5.0f %% per 10 pts [%+.0f, %+.0f]; typical departure %4.1f pts -> %+4.0f %%; ppet60 %+4.0f %% per 100 mm",
  source, spec, region, pct_per10, lo, hi, dep_typ, pct_typ, ppet60_pct)))
log_msg("=== done | ", out, " | log: ", log_file)
