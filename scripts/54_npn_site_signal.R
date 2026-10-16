#!/usr/bin/env Rscript
# Step 54 (NPN grass curing study): does the shared year signal follow local weather, or the site?
# 53 found a year signal shared by the grasses at a site (~18 d SD) that the ladder (52) missed.
# Here that signal is estimated per site-year from the same mixed model (year + region:year +
# site:year random effects, the part every grass at the site shares) and tested directly:
#   1. correlation with each site-year weather departure (51 metrics, departure from the site's
#      own mean), plus moisture memory, rain and heat at the site's typical date minus 14 d (lead)
#      and at the typical date itself (hindcast)
#   2. leave-one-year-out forward selection (<= 3 drivers) on those departures: cross-validated
#      R^2 of the shared signal
#   3. if data/npn_site_names.csv exists (54a_npn_site_names.R, laptop with rnpn): likely managed
#      sites (garden, campus, park, yard, ...) flagged; model refitted with and without them
# The estimates are shrunken (BLUPs), so correlations are conservative but not biased in sign.
# Inputs: npn_curing_metrics.rds (GSI_OBS_DIR or study_dir), <run_dir>/npn_study/npn_metrics.rds
# (51), npn_site_daily.rds (50), region from eval/target_check. Needs lme4.
# Output <run_dir>/npn_study/site_signal/: signal.csv (site-year signal + drivers),
# driver_r.csv, cv.csv, managed.csv (when names exist).
#
# Usage: Rscript scripts/54_npn_site_signal.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "obs.R"))
source(here::here("R", "features.R"))
source(here::here("R", "block2.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
source(here::here("R", "npn.R"))

cfg <- load_cfg(profile)
log_file <- log_init("54_npn_site_signal", profile)
ns <- cfg$npn_study %||% list()
sdir <- file.path(cfg$run_dir, "npn_study"); odir <- file.path(sdir, "site_signal")
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 54_npn_site_signal | profile: ", profile)
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs the lme4 package")
lead <- ns$lead_days %||% 14

ddir <- Sys.getenv("GSI_OBS_DIR", obs_study_dir(cfg))
cm <- readRDS(file.path(ddir, "npn_curing_metrics.rds"))
sp <- utils::read.csv(here::here("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
M <- readRDS(file.path(sdir, "npn_metrics.rds")); YM <- M$site_years
D <- readRDS(file.path(sdir, "npn_site_daily.rds"))
u <- cm[cm$usable %in% TRUE & cm$n_obs >= (ns$min_readings %||% 12), ]
m <- match(u$species_id, sp$species_id); u$pathway <- sp$pathway[m]; u <- u[!(sp$crop[m] %in% "yes"), ]
tc <- file.path(cfg$run_dir, "eval", "target_check", "site_years.csv"); u$region <- 0
if (file.exists(tc)) {
  r <- utils::read.csv(tc, stringsAsFactors = FALSE)
  r$site_id <- as.integer(sub("^npn(\\d+)\\|.*$", "\\1", r$site)); r <- r[!is.na(r$site_id) & !is.na(r$region), ]
  u$region <- r$region[match(u$site_id, r$site_id)]; u$region[is.na(u$region)] <- 0
}

# likely managed sites from names (optional)
nf <- here::here("data", "npn_site_names.csv")
managed_re <- "garden|campus|universit|college|school|park|yard|home|house|residen|arboret|cemeter|golf|lawn|church|office|museum|zoo|nursery|farm|ranch hq|backyard|front|street|apartment|condo|library|center"
u$managed <- NA
if (file.exists(nf)) {
  sn <- utils::read.csv(nf, stringsAsFactors = FALSE)
  u$site_name <- sn$site_name[match(u$site_id, sn$site_id)]
  u$managed <- grepl(managed_re, tolower(u$site_name))
  log_msg(sprintf("site names for %d of %d sites; %d flagged as likely managed", sum(unique(u$site_id) %in% sn$site_id),
                  length(unique(u$site_id)), length(unique(u$site_id[u$managed %in% TRUE]))))
} else log_msg("no data/npn_site_names.csv: managed-site check skipped (run 54a_npn_site_names.R on the laptop)")

cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
cmed <- function(x) { x <- x[!is.na(x)]; a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
fml <- y ~ 1 + (1 | site) + (1 | species) + (1 | site:species) + (1 | year) + (1 | region:year) + (1 | site:year) + (1 | site:species:year)
shared <- c("year", "region:year", "site:year"); place <- c("site", "species", "site:species")
prep <- function(col) {
  x <- u[!is.na(u[[col]]), ]; x$series <- paste(x$site_id, x$species_id)
  med <- c(tapply(x[[col]], x$series, cmed)); d <- cdiff(x[[col]], med[x$series])
  x$y <- as.numeric(med[x$series]) + d; x <- x[abs(d) <= 90, ]
  x$site <- factor(x$site_id); x$species <- factor(x$species_id); x$year <- factor(x$year); x$region <- factor(x$region); x
}
fit_vc <- function(g) {
  f <- lme4::lmer(fml, data = g, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))
  vc <- as.data.frame(lme4::VarCorr(f)); v <- stats::setNames(vc$vcov, ifelse(vc$grp == "Residual", "plant", vc$grp))
  list(fit = f, v = v, within = sum(v[setdiff(names(v), place)]))
}

# daily fields for the cutoff drivers
dates <- D$dates; doy <- as.integer(format(dates, "%j"))
PET <- hargreaves(D$tmin, D$tmax, ra_matrix(D$sites$lat, doy))
EW30 <- ewma_rows(D$ppt - PET, 30); EW60 <- ewma_rows(D$ppt - PET, 60)
cut_drivers <- function(site_id, year, day, tag) {
  srow <- match(site_id, D$sites$site_id); k <- match(as.Date(sprintf("%d-01-01", year)) + round(day) - 1, dates)
  w <- function(M, n, f) vapply(seq_along(k), function(i) if (is.na(k[i]) || is.na(srow[i]) || k[i] <= n) NA_real_ else f(M[srow[i], (k[i] - n + 1):k[i]]), 0)
  o <- data.frame(EW30[cbind(srow, k)], EW60[cbind(srow, k)], w(D$ppt, 30, sum), w(D$ppt, 90, sum), w(D$tmax, 30, mean), w(D$vpdmax, 30, mean))
  names(o) <- paste0(c("ew30", "ew60", "p30", "p90", "tx30", "vpd30"), "_", tag); o
}

targets <- c(greenup = "g50_up_doy", curing = "g50_down_doy", cured = "cured_doy")
fixed <- setdiff(names(YM), c("site_id", "year"))
sig_rows <- list(); r_rows <- list(); cv_rows <- list(); man_rows <- list()
for (tg in names(targets)) {
  x <- prep(targets[[tg]]); F <- fit_vc(x)
  re <- lme4::ranef(F$fit)
  sy <- unique(x[, c("site_id", "year", "region")])
  sy$sig_year <- re$year[as.character(sy$year), 1]
  ry <- re[["region:year"]]; sy$sig_region <- ry[paste(sy$region, sy$year, sep = ":"), 1]
  sty <- re[["site:year"]]; sy$sig_site <- sty[paste(sy$site_id, sy$year, sep = ":"), 1]
  sy[is.na(sy)] <- 0
  sy$signal <- sy$sig_year + sy$sig_region + sy$sig_site
  sy$n_series <- as.integer(table(paste(x$site_id, x$year))[paste(sy$site_id, sy$year)])
  # drivers: fixed-date departures from the site mean, and cutoff drivers at the site typical date
  mm <- match(paste(sy$site_id, sy$year), paste(YM$site_id, YM$year))
  for (v in fixed) { clim <- tapply(YM[[v]], YM$site_id, mean, na.rm = TRUE); sy[[v]] <- YM[[v]][mm] - clim[as.character(sy$site_id)] }
  typ <- c(tapply(x$y, x$site_id, cmed)); td <- as.numeric(typ[as.character(sy$site_id)])
  cd <- cbind(cut_drivers(sy$site_id, as.integer(as.character(sy$year)), td - lead, "lead"),
              cut_drivers(sy$site_id, as.integer(as.character(sy$year)), td, "atdate"))
  for (v in names(cd)) cd[[v]] <- cd[[v]] - stats::ave(cd[[v]], sy$site_id, FUN = function(z) mean(z, na.rm = TRUE))
  sy <- cbind(sy, cd); drivers <- c(fixed, names(cd))
  sy$target <- tg; sig_rows[[tg]] <- sy
  log_msg(sprintf("%s: %d site-years at %d sites | shared signal SD %.1f d (model %.1f d)", tg, nrow(sy), length(unique(sy$site_id)),
                  stats::sd(sy$signal), sqrt(sum(F$v[shared]))))
  # 1. univariate correlations
  rr <- do.call(rbind, lapply(drivers, function(v) { ok <- !is.na(sy[[v]]); if (sum(ok) < 30) return(NULL)
    ct <- stats::cor.test(sy[[v]][ok], sy$signal[ok]); data.frame(target = tg, driver = v, n = sum(ok), r = unname(ct$estimate), p = ct$p.value) }))
  rr$p_bh <- stats::p.adjust(rr$p, "BH"); rr <- rr[order(-abs(rr$r)), ]; r_rows[[tg]] <- rr
  log_msg(sprintf("  strongest drivers: %s | %d of %d with BH p < 0.05", paste(sprintf("%s %.2f", utils::head(rr$driver, 4), utils::head(rr$r, 4)), collapse = ", "),
                  sum(rr$p_bh < 0.05), nrow(rr)))
  # 2. leave-one-year-out forward selection, two sets (no at-date drivers in "lead")
  for (st in c("lead", "hindcast")) {
    dv <- if (st == "lead") setdiff(drivers, grep("_atdate$", drivers, value = TRUE)) else drivers
    X <- as.matrix(sy[, dv]); X[is.na(X)] <- 0; y <- sy$signal; pred <- rep(NA_real_, nrow(sy)); picks <- character()
    for (Y in unique(sy$year)) {
      te <- which(sy$year == Y); tr <- which(sy$year != Y); sel <- character(); best <- 0; sst <- sum(y[tr]^2)
      for (s in 1:3) {
        g <- vapply(setdiff(dv, sel), function(cn) 1 - sum(stats::lm.fit(X[tr, c(sel, cn), drop = FALSE], y[tr])$residuals^2) / sst, 0)
        if (max(g) - best < 0.01) break
        sel <- c(sel, names(which.max(g))); best <- max(g)
      }
      pred[te] <- if (length(sel)) X[te, sel, drop = FALSE] %*% stats::lm.fit(X[tr, sel, drop = FALSE], y[tr])$coefficients else 0
      picks <- c(picks, sel)
    }
    r2 <- 1 - sum((y - pred)^2) / sum(y^2)
    top <- names(sort(table(picks), decreasing = TRUE))[seq_len(min(3, length(unique(picks))))]
    cv_rows[[length(cv_rows) + 1]] <- data.frame(target = tg, set = st, site_years = nrow(sy), r2_cv = r2, r_cv = stats::cor(pred, y),
                                                 top_drivers = paste(top, collapse = "; "))
    log_msg(sprintf("  %s: cross-validated R2 of the shared signal %.2f (r %.2f) | most picked: %s", st, r2, stats::cor(pred, y), paste(top, collapse = ", ")))
  }
  # 3. managed sites
  if (any(!is.na(x$managed))) for (mg in c(TRUE, FALSE)) {
    g <- x[x$managed %in% mg, ]; if (length(unique(g$site)) < 8) next
    Fg <- fit_vc(droplevels(g))
    man_rows[[length(man_rows) + 1]] <- data.frame(target = tg, managed = mg, sites = length(unique(g$site)), plant_years = nrow(g),
      sd_year_to_year = sqrt(Fg$within), shared_share = sum(Fg$v[shared]) / Fg$within, site_year_share = Fg$v[["site:year"]] / Fg$within,
      plant_share = Fg$v[["plant"]] / Fg$within)
    with(man_rows[[length(man_rows)]], log_msg(sprintf("  managed=%s: %d sites | year-to-year SD %.0f d | shared %.0f %% (site-wide %.0f %%), plant %.0f %%",
                                                       managed, sites, sd_year_to_year, 100 * shared_share, 100 * site_year_share, 100 * plant_share)))
  }
}
keep <- unique(unlist(lapply(sig_rows, names)))
S <- do.call(rbind, lapply(sig_rows, function(d) { d[setdiff(keep, names(d))] <- NA; d[keep] }))
utils::write.csv(S, file.path(odir, "signal.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, r_rows), file.path(odir, "driver_r.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, cv_rows), file.path(odir, "cv.csv"), row.names = FALSE)
if (length(man_rows)) utils::write.csv(do.call(rbind, man_rows), file.path(odir, "managed.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), odir, log_file))
