#!/usr/bin/env Rscript
# Step 07 (NPN x fire study; laptop, minutes; MASS): Tier 3 year-to-year test. In years when a phenophase
# comes early or late, does the fire change around it move with it, beyond dryness and beyond what weather
# alone would predict? Design approved by Mike on 2 Nov 2026 and fixed in config (npn_fire: yeartoyear)
# before any result was seen.
#
# Cases: 05 tier3_ok rows (ecoregion x plant group x event), lead fire record, matched fuel (as 06).
# Cases need >= min_years (tier3 min_years, 8) kept years.
# Year departure (NPN): each site's date minus that site's own typical date (median over its years; sites
#   with >= min_site_years years), averaged over sites (sites weighted, not plants). Years kept: >= min_sites
#   sites with a departure, inside the fire record. dep_wk = round(departure / 7), clipped to +-max_shift.
# Windows as 06: typical week -8..+7 in each kept year. Weekly fires ~ dryness + step, where the step
#   ('after') starts at:
#     M0 calendar   the typical week in every year (= step 2)
#     M1 NPN        the typical week + that year's NPN departure
#     M2 weather    the typical week + a weather-predicted departure: NPN departure regressed (training years)
#                   on the year's mean temperature anomaly over the weeks before the typical date (wx_weeks),
#                   matched-fuel weight; predicted for the held-out year
# Leave-one-year-out (LOYO): each model is fitted (negative binomial, year effects + dryness + after) on the
#   other years; the held-out year is scored by the multinomial log-likelihood of how its window fires fall
#   across weeks (its own total given, so no year level needs predicting).
#   gain_npn = sum over years of LL(M1) - LL(M0); gain_npn_vs_wx = LL(M1) - LL(M2); gain_wx = LL(M2) - LL(M0).
#   Bootstrap: held-out years resampled B times; 90 % percentile interval of the mean per-year gain.
# Pass (NPN adds information): gain_npn > 0 with its interval above 0, AND M1's step in the expected direction
#   (events with one; config inventory: events). 'Beats weather' additionally: gain_npn_vs_wx interval above 0.
# Feasibility: years, sites per year, SD of the year departure (days), share of years the NPN step week
#   differs from the typical week; with few years or small departures the test cannot see much.
# Outputs <run_dir>/npn_fire/yeartoyear/: results.csv (case), years.csv (case x year: departure, weather
#   anomaly, LOYO log-likelihoods), summary.csv (group x event).
# Usage: Rscript npn_fire/scripts/07_nf_yeartoyear.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("07_nf_yeartoyear", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; yc <- cfg$npn_fire$yeartoyear
W <- as.integer(sc$window_weeks %||% 8); B <- as.integer(yc$bootstrap$B %||% 1000); lev <- yc$bootstrap$level %||% 0.9
min_sites <- as.integer(yc$min_sites %||% ic$tier3$min_sites_year %||% 10); min_sy <- as.integer(yc$min_site_years %||% 2)
max_shift <- as.integer(yc$max_shift_weeks %||% 6); wxw <- as.integer(yc$wx_weeks %||% 8)
min_years <- as.integer(yc$min_years %||% ic$tier3$min_years %||% 8); n_short <- 0
out <- nf_out(cfg, "yeartoyear")
set.seed(yc$seed %||% 20261102)
t0 <- Sys.time()
log_msg("=== 07_nf_yeartoyear | profile: ", profile, " | window +-", W, " weeks | >= ", min_sites, " sites / year | LOYO, bootstrap B ", B)

INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
SD <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "inventory", "site_dates"))); SD$eco_l2_code <- as.character(SD$eco_l2_code)
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
has_t <- "tmean" %in% names(PW)
if (!has_t) log_err("eco_weekly.csv has no tmean (rerun 04 with force on the P720 and sync): the weather comparison (M2) is skipped")
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))
cases <- INV[INV$tier3_ok %in% TRUE, ]
mx <- suppressWarnings(as.integer(Sys.getenv("NF_MAX_CASES", ""))); if (is.finite(mx)) cases <- cases[sample(nrow(cases), min(mx, nrow(cases))), ]  # testing only
log_msg(sprintf("cases: %d ecoregion x group x event (Tier 3 possible)", nrow(cases)))

wrap <- function(d) ifelse(d > 182, d - 365, ifelse(d < -182, d + 365, d))
# yearly NPN departure for one case: mean over sites of (site date - site typical date)
npn_dep <- function(eco, g, ev) {
  d <- SD[SD$eco_l2_code == eco & SD$group == g & SD$event == ev & !is.na(SD$doy), c("site_id", "year", "doy")]
  if (!nrow(d)) return(NULL)
  ny <- tapply(d$year, d$site_id, function(v) length(unique(v)))
  d <- d[d$site_id %in% names(ny)[ny >= min_sy], ]
  if (!nrow(d)) return(NULL)
  d$dep <- wrap(d$doy - stats::ave(d$doy, d$site_id, FUN = stats::median))
  a <- stats::aggregate(dep ~ year, d, mean); n <- stats::aggregate(site_id ~ year, d, function(v) length(unique(v)))
  a$sites <- n$site_id[match(a$year, n$year)]
  a[a$sites >= min_sites, ]
}
wk_val <- function(eco, wt, yy, var) {
  d <- PW[PW$eco_l2_code == eco & PW$weight == wt & PW$year %in% yy, ]
  v <- rep(NA_real_, length(yy) * 52); v[(d$year - min(yy)) * 52 + d$week] <- d[[var]]
  v
}
series <- function(eco, fuels, o, yy) {
  d <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% fuels & FW$year %in% yy, ]
  v <- numeric(length(yy) * 52)
  if (nrow(d)) { s <- tapply(d$value, (d$year - min(yy)) * 52 + d$week, sum); v[as.integer(names(s))] <- s }
  v
}
# multinomial log-likelihood of a year's weekly fires given linear predictors eta (its total given)
mn_ll <- function(y, eta) { if (sum(y) <= 0) return(0); p <- exp(eta - max(eta)); p <- p / sum(p); sum(y * log(pmax(p, 1e-12))) }
fit_step <- function(D) {
  m <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ wy + dry + after, data = D)), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- stats::coef(m); if (any(is.na(b[c("dry", "after")]))) return(NULL)
  c(dry = unname(b[["dry"]]), after = unname(b[["after"]]))
}
bci <- function(g) { if (length(g) < 3) return(c(NA, NA))
  bs <- vapply(seq_len(B), function(i) mean(sample(g, replace = TRUE)), 0)
  unname(stats::quantile(bs, c((1 - lev) / 2, 1 - (1 - lev) / 2))) }

RES <- list(); YRS <- list()
for (r in seq_len(nrow(cases))) {
  e <- cases[r, ]; o <- e$lead_outcome
  fuels <- unlist(ic$match[[e$group]]); wt <- sc$dryness_weight[[e$group]]
  exp_dir <- ic$events[[e$event]]$expected %||% "none"
  A <- npn_dep(e$eco_l2_code, e$group, e$event)
  yy_all <- oy[[o]]
  if (is.null(A)) next
  A <- A[A$year %in% yy_all, ]
  if (nrow(A) < min_years) { n_short <- n_short + 1; next }
  f <- series(e$eco_l2_code, fuels, o, yy_all); dr <- wk_val(e$eco_l2_code, wt, yy_all, "ppet60") / 100
  tm <- if (has_t) wk_val(e$eco_l2_code, wt, yy_all, "tmean") else NULL
  rel <- seq(-W, W - 1); tw <- e$typical_week
  Y <- list()
  for (k in seq_len(nrow(A))) {
    i <- match(A$year[k], yy_all); t <- (i - 1) * 52 + tw + rel
    if (min(t) < 1 || max(t) > length(f) || any(is.na(dr[t]))) next
    wx <- if (has_t) { tb <- (i - 1) * 52 + tw - seq_len(wxw); if (min(tb) < 1) NA_real_ else mean(tm[tb], na.rm = TRUE) } else NA_real_
    Y[[length(Y) + 1]] <- list(year = A$year[k], dep = A$dep[k], sites = A$sites[k], y = f[t], dry = dr[t], wx = wx,
                               wk = tw)
  }
  if (length(Y) < min_years) { n_short <- n_short + 1; next }
  yrs <- vapply(Y, `[[`, 0, "year"); dep <- vapply(Y, `[[`, 0, "dep"); wx <- vapply(Y, `[[`, 0, "wx")
  wxa <- wx - mean(wx, na.rm = TRUE)                                   # temperature anomaly (deg C)
  shift <- function(d) pmax(-max_shift, pmin(max_shift, round(d / 7)))
  long <- function(idx, sh) do.call(rbind, lapply(seq_along(idx), function(j) { z <- Y[[idx[j]]]
    data.frame(wy = z$year, y = z$y, dry = z$dry, after = as.integer(rel >= sh[j])) }))
  ll <- matrix(NA_real_, length(Y), 3, dimnames = list(NULL, c("cal", "npn", "wx")))
  for (h in seq_along(Y)) {
    tr <- setdiff(seq_along(Y), h)
    sh_tr <- list(cal = rep(0, length(tr)), npn = shift(dep[tr]), wx = NULL)
    sh_h <- c(cal = 0, npn = shift(dep[h]), wx = NA)
    if (has_t && sum(is.finite(wxa[tr])) >= 3 && is.finite(wxa[h])) {
      ok <- is.finite(wxa[tr]); lm1 <- stats::lm(dep[tr][ok] ~ wxa[tr][ok])
      cf <- stats::coef(lm1); pr <- function(x) cf[[1]] + cf[[2]] * x
      sh_tr$wx <- shift(ifelse(is.finite(wxa[tr]), pr(wxa[tr]), 0)); sh_h[["wx"]] <- shift(pr(wxa[h]))
    }
    for (mname in c("cal", "npn", "wx")) {
      if (is.null(sh_tr[[mname]]) || is.na(sh_h[[mname]])) next
      D <- long(tr, sh_tr[[mname]]); D$wy <- factor(D$wy)
      if (sum(D$y) <= 0) next
      b <- fit_step(D); if (is.null(b)) next
      eta <- b[["dry"]] * Y[[h]]$dry + b[["after"]] * as.integer(rel >= sh_h[[mname]])
      ll[h, mname] <- mn_ll(Y[[h]]$y, eta)
    }
  }
  # full-data M1 step (direction check) and feasibility numbers
  Df <- long(seq_along(Y), shift(dep)); Df$wy <- factor(Df$wy); bf <- fit_step(Df)
  ok1 <- is.finite(ll[, "cal"]) & is.finite(ll[, "npn"]); g1 <- (ll[, "npn"] - ll[, "cal"])[ok1]
  ok2 <- ok1 & is.finite(ll[, "wx"]); g2 <- (ll[, "npn"] - ll[, "wx"])[ok2]; g3 <- (ll[, "wx"] - ll[, "cal"])[ok2]
  c1 <- bci(g1); c2 <- bci(g2); c3 <- bci(g3)
  step_b <- if (is.null(bf)) NA_real_ else bf[["after"]]
  dir_ok <- switch(exp_dir, down = step_b < 0, up = step_b > 0, TRUE)
  pass <- length(g1) >= 4 && sum(g1) > 0 && is.finite(c1[1]) && c1[1] > 0 && isTRUE(dir_ok)
  RES[[length(RES) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o, expected = exp_dir,
    typical_week = tw, years = length(Y), years_scored = length(g1), sites_median = stats::median(vapply(Y, `[[`, 0, "sites")),
    dep_sd_days = round(stats::sd(dep), 1), share_shifted = round(mean(shift(dep) != 0), 2),
    wx_dep_r = if (sum(is.finite(wxa)) >= 4) round(stats::cor(dep, wxa, use = "complete.obs"), 2) else NA_real_,
    step_ratio_npn = round(exp(step_b), 3), expected_direction = if (exp_dir == "none") NA else dir_ok,
    gain_npn = round(sum(g1), 2), gain_npn_lo = round(c1[1] * length(g1), 2), gain_npn_hi = round(c1[2] * length(g1), 2),
    gain_npn_vs_wx = if (length(g2)) round(sum(g2), 2) else NA_real_,
    gain_npn_vs_wx_lo = round(c2[1] * length(g2), 2), gain_npn_vs_wx_hi = round(c2[2] * length(g2), 2),
    gain_wx = if (length(g3)) round(sum(g3), 2) else NA_real_,
    pass = pass, beats_weather = pass && length(g2) >= 4 && is.finite(c2[1]) && c2[1] > 0, stringsAsFactors = FALSE)
  YRS[[length(YRS) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, year = yrs,
    sites = vapply(Y, `[[`, 0, "sites"), dep_days = round(dep, 1), wx_anom_c = round(wxa, 2),
    ll_cal = round(ll[, "cal"], 3), ll_npn = round(ll[, "npn"], 3), ll_wx = round(ll[, "wx"], 3))
  if (r %% 25 == 0) log_msg(sprintf("  %d / %d cases (%.1f min)", r, nrow(cases), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!length(RES)) stop("no case could be fitted: check inventory.csv (tier3_ok), site_dates, fire_weekly and eco_weekly")
R <- do.call(rbind, RES); utils::write.csv(R, file.path(out, "results.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, YRS), file.path(out, "years.csv"), row.names = FALSE)

key <- paste(R$group, R$event)
S <- do.call(rbind, lapply(split(R, key), function(d) data.frame(group = d$group[1], event = d$event[1], expected = d$expected[1],
  cases = nrow(d), years_median = stats::median(d$years), dep_sd_median = stats::median(d$dep_sd_days),
  gain_positive = sum(d$gain_npn > 0), pass = sum(d$pass), beats_weather = sum(d$beats_weather, na.rm = TRUE),
  pass_ecoregions = paste(d$eco_l2_code[d$pass], collapse = ", "), stringsAsFactors = FALSE)))
S <- S[order(S$group, match(S$event, names(ic$events))), ]
utils::write.csv(S, file.path(out, "summary.csv"), row.names = FALSE)

# chance level: with no information in the departure, gain_npn > 0 about half the time or less
log_msg(sprintf("cases left out (< %d years with >= %d sites whose own typical date is known): %d", min_years, min_sites, n_short))
log_msg(sprintf("cases fitted: %d | median years %d | median departure SD %.1f days | NPN step moved in %.0f %% of case-years",
  nrow(R), as.integer(stats::median(R$years)), stats::median(R$dep_sd_days), 100 * mean(R$share_shifted)))
log_msg(sprintf("gain_npn > 0 (NPN date beats the calendar date in held-out years): %d of %d (%.0f %%) | pass %d | beats weather %d",
  sum(R$gain_npn > 0), nrow(R), 100 * mean(R$gain_npn > 0), sum(R$pass), sum(R$beats_weather, na.rm = TRUE)))
if (has_t) log_msg(sprintf("gain_wx > 0 (weather date beats the calendar): %d of %d | NPN beats weather (gain_npn_vs_wx > 0): %d of %d | median r(departure, temperature anomaly) %.2f",
  sum(R$gain_wx > 0, na.rm = TRUE), sum(!is.na(R$gain_wx)), sum(R$gain_npn_vs_wx > 0, na.rm = TRUE), sum(!is.na(R$gain_npn_vs_wx)),
  stats::median(R$wx_dep_r, na.rm = TRUE)))
for (ev in names(ic$events)) { s <- R[R$event == ev, ]; if (!nrow(s)) next
  log_msg(sprintf("  %-22s %3d cases | gain > 0 %3.0f %% | pass %2d | beats weather %2d | median dep SD %4.1f d", ev, nrow(s),
                  100 * mean(s$gain_npn > 0), sum(s$pass), sum(s$beats_weather, na.rm = TRUE), stats::median(s$dep_sd_days))) }
log_msg(sprintf("=== done | %d cases | %.1f min | %s | log: %s", nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
