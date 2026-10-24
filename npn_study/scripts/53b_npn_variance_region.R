#!/usr/bin/env Rscript
# Step 53b (NPN grass curing study; laptop): the year-to-year variance partition of 53, fitted
# separately in each EPA Level II ecoregion, for the maps on the page. Same plant-years, unwrapping
# and season-mismatch rule as 53 (all grasses, crops out). Per ecoregion:
#   date ~ 1 + (1|site) + (1|species) + (1|site:species)      place and species
#            + (1|year)                                         year signal shared across the region
#            + (1|site:year)                                    year signal of one site
#            + (1|site:species:year)                            year signal of one species at a site
#            + residual                                         plant within series-year
# With one site in a region, site terms drop out and the year term is that site's year signal.
# Shares are of the year-to-year part (all but place). Fitted where a region has at least
# npn_study: region_min_plant_years plant-years (default 60) and 4 years.
# Output <run_dir>/npn_study/variance/by_region.csv and report/tables/npn_variance_region.csv.
#
# Usage: Rscript npn_study/scripts/53b_npn_variance_region.R [profile]   (after 59)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs the lme4 package")

cfg <- load_npn_cfg(profile)
log_file <- log_init("53b_npn_variance_region", profile)
ns <- cfg$npn_study %||% list()
odir <- file.path(cfg$run_dir, "npn_study", "variance"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)
min_py <- as.numeric(ns$region_min_plant_years %||% 60)
log_msg("=== 53b_npn_variance_region | profile: ", profile, " | regions with >= ", min_py, " plant-years and >= 4 years")

cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
u <- cm[cm$usable %in% TRUE & cm$n_obs >= (ns$min_readings %||% 12), ]
m <- match(u$species_id, sp$species_id)
u$pathway <- sp$pathway[m]; u <- u[!(sp$crop[m] %in% "yes"), ]
u$region <- npn_site_regions(u$site_id)

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
targets <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
ctl <- lme4::lmerControl(calc.derivs = FALSE, check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
rows <- list()
for (tg in names(targets)) {
  x <- u[!is.na(u[[targets[[tg]]]]), ]
  x$series <- paste(x$site_id, x$species_id)
  med <- c(tapply(x[[targets[[tg]]]], x$series, cmed))
  d <- cdiff(x[[targets[[tg]]]], med[x$series]); x$y <- as.numeric(med[x$series]) + d
  x <- x[abs(d) <= 90 & x$region != "none", ]
  for (rg in sort(unique(x$region))) {
    g <- x[x$region == rg, ]
    nsit <- length(unique(g$site_id)); nyr <- length(unique(g$year)); nsp <- length(unique(g$species_id))
    if (nrow(g) < min_py || nyr < 4) next
    g$site <- factor(g$site_id); g$species <- factor(g$species_id); g$year <- factor(g$year)
    terms <- c(if (nsit > 1) "(1 | site)", if (nsp > 1) "(1 | species)", if (nsit > 1 && nsp > 1) "(1 | site:species)",
               "(1 | year)", if (nsit > 1) "(1 | site:year)", "(1 | site:species:year)")
    if (nsit == 1 && nsp == 1) terms <- setdiff(terms, "(1 | site:species:year)")   # = year then
    fit <- tryCatch(suppressMessages(lme4::lmer(stats::as.formula(paste("y ~ 1 +", paste(terms, collapse = " + "))), data = g, REML = TRUE, control = ctl)),
                    error = function(e) { log_msg("  ", tg, " ", rg, ": fit failed (", conditionMessage(e), ")"); NULL })
    if (is.null(fit)) next
    vc <- as.data.frame(lme4::VarCorr(fit)); v <- stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "plant", vc$grp))
    gv <- function(k) if (k %in% names(v)) v[[k]] else 0
    yy <- gv("year") + gv("site:year") + gv("site:species:year") + gv("plant")
    rows[[length(rows) + 1]] <- data.frame(target = tg, region = rg, plant_years = nrow(g), sites = nsit, species = nsp, years = nyr,
      sd_year_to_year = sqrt(yy), sd_region_year = sqrt(gv("year")),
      share_region_year = gv("year") / yy, share_site_year = gv("site:year") / yy,
      share_species_year = gv("site:species:year") / yy, share_plant = gv("plant") / yy, one_site = nsit == 1)
    with(rows[[length(rows)]], log_msg(sprintf("  %-7s %-5s %4d plant-years, %2d sites, %2d years | year-to-year SD %4.1f d | region-wide %3.0f %%, site %3.0f %%, species %3.0f %%, plant %3.0f %%%s",
      target, region, plant_years, sites, years, sd_year_to_year, 100 * share_region_year, 100 * share_site_year, 100 * share_species_year, 100 * share_plant,
      if (one_site) " (one site)" else "")))
  }
}
R <- do.call(rbind, rows)
utils::write.csv(R, file.path(odir, "by_region.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE)
utils::write.csv(R, file.path(tdir, "npn_variance_region.csv"), row.names = FALSE)
log_msg("=== done | ", file.path(odir, "by_region.csv"), " | log: ", log_file)
