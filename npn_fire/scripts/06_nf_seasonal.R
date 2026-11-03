#!/usr/bin/env Rscript
# Step 06 (NPN x fire study; laptop, minutes; MASS): Tier 2 seasonal test. Does the weekly fire rate in
# the matched fuel change at an event's typical date, beyond what the seasonal course of dryness
# explains? Design approved by Mike on 1 Nov 2026 and fixed in config (npn_fire: seasonal) before any
# result was seen.
#
# Cases: every ecoregion x plant group x event where 05 found the seasonal test possible (tier2_ok); per
# case the lead record and every other record with >= tier2 min_covered_years covered years.
# For each window-year Y (each year the record covers): weeks typical_week - W .. typical_week + W - 1 on
# a continuous weekly axis (so windows can cross Dec / Jan); rel = week - typical week; after = rel >= 0.
#   fires ~ factor(window-year) + ppet60 + after          negative binomial, log link
# fires = weekly fires (or pixel-days) in the plant group's matched fuel (03); ppet60 = the ecoregion's
# weekly 60-day rain minus PET weighted by the matched fuel type (04; config dryness_weight), in 100 mm.
# Ratio = exp(after): the fire rate from the typical date on, relative to before, dryness held fixed.
#   Placebo: the same model (dispersion fixed at the real fit's) at every typical-week shift k with
#     |k| >= placebo_min_shift (cyclic over 52 weeks). placebo_pct = share of placebo ratios below the real one.
#   Pass: expected down -> placebo_pct <= pass_tail; up -> >= 1 - pass_tail; none -> either tail.
#   Expected directions: config inventory: events (Mike, 1 Nov: all 13 pheno classes; flowering / fruit 'up').
#   Bootstrap (lead record): window-years resampled B times, 90 % percentile interval of the ratio.
# Composite curves: per case, weekly fires / the window-year's mean, averaged over window-years, by rel.
# Outputs <run_dir>/npn_fire/seasonal/: results.csv, curves.csv, summary.csv (group x event x outcome:
#   cases, share in the expected direction, share passing).
# Usage: Rscript npn_fire/scripts/06_nf_seasonal.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && !grepl("=", args[1])) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
# Exploratory override (Mike, 3 Nov 2026): outcome=<fire outcome> runs the same test with that outcome in place
# of the lead record, for the same cases, writing to <out>_<outcome>/ (cases with < 5 years holding >= 3 fires
# of the outcome, or of its count companion for area outcomes, in the window are skipped)
ov <- sub("^outcome=", "", grep("^outcome=", args, value = TRUE)); ov <- if (length(ov)) ov[1] else NA_character_
count_of <- c(fod_cplus_acres = "fod_cplus", mtbs_km2 = "mtbs_fires")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("06_nf_seasonal", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal
W <- as.integer(sc$window_weeks %||% 8); tail_p <- sc$pass_tail %||% 0.05; B <- as.integer(sc$bootstrap$B %||% 200)
lev <- sc$bootstrap$level %||% 0.9
out <- nf_out(cfg, if (is.na(ov)) "seasonal" else paste0("seasonal_", ov))
if (!is.na(ov)) B <- 0L   # exploratory runs: no bootstrap
set.seed(sc$seed %||% 20261101)
t0 <- Sys.time()
log_msg("=== 06_nf_seasonal | profile: ", profile, " | window +-", W, " weeks | placebo |shift| >= ", sc$placebo_min_shift,
        " | pass tail ", tail_p, " | bootstrap B ", B)

if (!is.na(ov)) log_msg("EXPLORATORY run: outcome ", ov, " in place of the lead record -> ", out)
INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
CV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "coverage.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))
cases <- INV[INV$tier2_ok %in% TRUE, ]
log_msg(sprintf("cases: %d ecoregion x group x event (Tier 2 possible)", nrow(cases)))

# Continuous weekly series (52 weeks a year) for one ecoregion, matched fuel and outcome over years yy
series <- function(eco, fuels, o, yy) {
  d <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% fuels & FW$year %in% yy, ]
  v <- numeric(length(yy) * 52)
  if (nrow(d)) { s <- tapply(d$value, (d$year - min(yy)) * 52 + d$week, sum); v[as.integer(names(s))] <- s }
  v
}
dry_series <- function(eco, wt, yy) {
  d <- PW[PW$eco_l2_code == eco & PW$weight == wt & PW$year %in% yy, ]
  v <- rep(NA_real_, length(yy) * 52); v[(d$year - min(yy)) * 52 + d$week] <- d$ppet60 / 100
  v
}
# Long data for windows anchored at week w0 of every year in yy (windows must lie inside the series)
win_data <- function(f, dr, yy, w0) {
  rel <- seq(-W, W - 1)
  L <- lapply(seq_along(yy), function(i) {
    t <- (i - 1) * 52 + w0 + rel
    if (min(t) < 1 || max(t) > length(f)) return(NULL)
    data.frame(wy = yy[i], rel = rel, y = f[t], dry = dr[t])
  })
  D <- do.call(rbind, L)
  if (is.null(D)) return(NULL)
  D <- D[!is.na(D$dry), ]; D$after <- as.integer(D$rel >= 0); D$wy <- factor(D$wy)
  D
}
# Fixed-dispersion NB fit on a prebuilt model matrix (fast path for placebo and bootstrap refits)
fit_fixed <- function(D, fam) {
  X <- stats::model.matrix(~ wy + dry + after, D)
  m <- tryCatch(suppressWarnings(stats::glm.fit(X, D$y, family = fam, control = stats::glm.control(maxit = 50))), error = function(e) NULL)
  if (is.null(m) || !m$converged) return(NA_real_)
  b <- m$coefficients[["after"]]; if (is.na(b)) NA_real_ else unname(b)
}

RES <- list(); CUR <- list()
for (r in seq_len(nrow(cases))) {
  e <- cases[r, ]
  fuels <- unlist(ic$match[[e$group]]); wt <- sc$dryness_weight[[e$group]]
  cv <- CV[CV$eco_l2_code == e$eco_l2_code & CV$group == e$group & CV$event == e$event, ]
  outs <- if (is.na(ov)) unique(c(e$lead_outcome, cv$outcome[cv$covered_years >= ic$tier2$min_covered_years])) else ov
  if (!is.na(ov)) e$lead_outcome <- ov
  exp_dir <- ic$events[[e$event]]$expected %||% "none"
  for (o in outs) {
    yy <- oy[[o]]; f <- series(e$eco_l2_code, fuels, o, yy); dr <- dry_series(e$eco_l2_code, wt, yy)
    D <- win_data(f, dr, yy, e$typical_week)
    if (is.null(D) || sum(D$y) <= 0 || length(unique(D$wy)) < 3) next
    if (!is.na(ov)) { oc <- if (ov %in% names(count_of)) count_of[[ov]] else ov
      Dc <- if (oc == ov) D else win_data(series(e$eco_l2_code, fuels, oc, yy), dr, yy, e$typical_week)
      if (is.null(Dc) || sum(tapply(Dc$y, Dc$wy, sum) >= 3) < 5) next
      D$y <- round(D$y) }
    m <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ wy + dry + after, data = D)), error = function(err) NULL)
    if (is.null(m)) next
    b <- stats::coef(m); theta <- m$theta; fam <- MASS::negative.binomial(theta)
    est <- unname(b[["after"]]); dry_b <- unname(b[["dry"]])
    # placebo dates
    shifts <- setdiff(-26:25, -(sc$placebo_min_shift - 1):(sc$placebo_min_shift - 1))
    pl <- vapply(shifts, function(k) { w <- ((e$typical_week - 1 + k) %% 52) + 1; Dp <- win_data(f, dr, yy, w)
      if (is.null(Dp) || sum(Dp$y) <= 0) NA_real_ else fit_fixed(Dp, fam) }, 0)
    pl <- pl[is.finite(pl)]
    pct <- if (length(pl)) mean(pl < est) else NA_real_
    pass <- if (is.na(pct)) NA else switch(exp_dir, down = pct <= tail_p, up = pct >= 1 - tail_p, pct <= tail_p | pct >= 1 - tail_p)
    dir_ok <- switch(exp_dir, down = est < 0, up = est > 0, NA)
    lo <- hi <- NA_real_
    if (o == e$lead_outcome && B > 0) {
      ys <- levels(D$wy)
      bs <- vapply(seq_len(B), function(i) { s <- sample(ys, replace = TRUE)
        Db <- do.call(rbind, lapply(seq_along(s), function(j) { d <- D[D$wy == s[j], ]; d$wy <- j; d })); Db$wy <- factor(Db$wy)
        fit_fixed(Db, fam) }, 0)
      q <- stats::quantile(bs, c((1 - lev) / 2, 1 - (1 - lev) / 2), na.rm = TRUE); lo <- exp(q[[1]]); hi <- exp(q[[2]])
    }
    RES[[length(RES) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o,
      lead = o == e$lead_outcome, typical_week = e$typical_week, window_years = length(unique(D$wy)), fires_in_windows = round(sum(D$y), 1),
      ratio = round(exp(est), 3), ratio_lo = round(lo, 3), ratio_hi = round(hi, 3), dryness_coef_per_100mm = round(dry_b, 3),
      theta = round(theta, 3), placebo_n = length(pl), placebo_pct = round(pct, 3), expected = exp_dir,
      expected_direction = dir_ok, pass = pass, stringsAsFactors = FALSE)
    # composite curve: weekly fires / window-year mean, averaged over window-years
    D$norm <- D$y / stats::ave(D$y, D$wy, FUN = function(v) if (mean(v) > 0) mean(v) else NA)
    cu <- tapply(D$norm, D$rel, mean, na.rm = TRUE)
    CUR[[length(CUR) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o,
                                         rel = as.integer(names(cu)), rel_rate = round(as.numeric(cu), 4))
  }
  if (r %% 25 == 0) log_msg(sprintf("  %d / %d cases (%.1f min)", r, nrow(cases), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!length(RES)) stop("no case could be fitted: check inventory.csv (tier2_ok), fire_weekly and eco_weekly")
R <- do.call(rbind, RES); CU <- do.call(rbind, CUR)
utils::write.csv(R, file.path(out, "results.csv"), row.names = FALSE)
utils::write.csv(CU, file.path(out, "curves.csv"), row.names = FALSE)

key <- paste(R$group, R$event, R$outcome)
S <- R[!duplicated(key), c("group", "event", "outcome", "expected")]
k <- paste(S$group, S$event, S$outcome)
S$cases <- as.integer(table(key)[k])
S$share_expected_direction <- round(as.numeric(tapply(R$expected_direction, key, function(v) mean(v, na.rm = TRUE))[k]), 3)
S$share_pass <- round(as.numeric(tapply(R$pass, key, function(v) mean(v, na.rm = TRUE))[k]), 3)
S$median_ratio <- round(as.numeric(tapply(R$ratio, key, stats::median)[k]), 3)
S <- S[order(S$outcome != "fod_all", S$group, match(S$event, names(ic$events))), ]
utils::write.csv(S, file.path(out, "summary.csv"), row.names = FALSE)

L <- R[R$lead, ]
log_msg(sprintf("lead-record cases: %d | expected direction %d of %d with a direction | pass %d of %d", nrow(L),
                sum(L$expected_direction %in% TRUE), sum(!is.na(L$expected_direction)), sum(L$pass %in% TRUE), sum(!is.na(L$pass))))
for (ev in names(ic$events)) { s <- L[L$event == ev, ]; if (!nrow(s)) next
  ed <- if (all(is.na(s$expected_direction))) "  - " else sprintf("%3.0f %%", 100 * mean(s$expected_direction, na.rm = TRUE))
  log_msg(sprintf("  %-22s (%-4s) %3d cases | median ratio %.2f | expected direction %s | pass %3.0f %%", ev, ic$events[[ev]]$expected %||% "none",
                  nrow(s), stats::median(s$ratio), ed, 100 * mean(s$pass, na.rm = TRUE))) }
log_msg(sprintf("=== done | %d results | %.1f min | %s | log: %s", nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
