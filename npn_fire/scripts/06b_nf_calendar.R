#!/usr/bin/env Rscript
# Step 06b (NPN x fire study; laptop, ~10-20 min; MASS): fixed-calendar check on the Tier 2 seasonal test.
# Asked by Mike on 2 Nov 2026 (option A, before step 3). Rules fixed in config (npn_fire: calendar) before
# any result was seen.
#
# Question: is the fire change at an ecoregion's own phenology date sharper than at the dates the same
# phenophase has elsewhere? If the ecoregion's own date does no better than borrowed dates, the step 2
# result is mostly time of year, not local phenology.
#
# Cases: as 06 (05 tier2_ok rows), lead record only. Same windows and model as 06 (fires ~ factor(window-year)
# + ppet60 + after, NB); the dispersion is estimated once at the own date and reused at every other date.
# Strength of a fire change b = log ratio: b for events expected 'up', -b for 'down', |b| for 'none'.
# Alternative dates, all for the same plant group and event:
#   national  the median typical week over all ecoregions where 05 found a typical date ("one calendar")
#   donors    each other ecoregion's typical week (distinct weeks, the own week left out)
# Per case: own_beats_national (strength own > national; NA when the national week is the own week);
#   donor_pct = share of donor dates with a weaker change than the own date (0.5 = no better than chance),
#   also within distance bands (|week difference| 1-2, 3-4, 5+).
# Units: many events share an ecoregion x group x typical week; summaries are also given with each
#   such unit counted once.
# Outputs <run_dir>/npn_fire/calendar/: results.csv (case), donors.csv (case x donor date), summary.csv.
# Usage: Rscript npn_fire/scripts/06b_nf_calendar.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("06b_nf_calendar", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; kc <- cfg$npn_fire$calendar
W <- as.integer(sc$window_weeks %||% 8)
bands <- kc$bands %||% list(near = c(1, 2), mid = c(3, 4), far = c(5, 26))
out <- nf_out(cfg, "calendar")
t0 <- Sys.time()
log_msg("=== 06b_nf_calendar | profile: ", profile, " | window +-", W, " weeks | lead record | alternative dates: national median + other ecoregions")

INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))
cases <- INV[INV$tier2_ok %in% TRUE, ]
mx <- suppressWarnings(as.integer(Sys.getenv("NF_MAX_CASES", ""))); if (is.finite(mx)) cases <- cases[sample(nrow(cases), min(mx, nrow(cases))), ]  # testing only
dates <- INV[!is.na(INV$typical_week), c("eco_l2_code", "group", "event", "typical_week")]
log_msg(sprintf("cases: %d | ecoregion typical dates available: %d", nrow(cases), nrow(dates)))

# --- same series, windows and fixed-dispersion fit as 06 ---
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
fit_fixed <- function(D, fam) {
  if (is.null(D) || sum(D$y) <= 0) return(NA_real_)
  X <- stats::model.matrix(~ wy + dry + after, D)
  m <- tryCatch(suppressWarnings(stats::glm.fit(X, D$y, family = fam, control = stats::glm.control(maxit = 50))), error = function(e) NULL)
  if (is.null(m) || !m$converged) return(NA_real_)
  b <- m$coefficients[["after"]]; if (is.na(b)) NA_real_ else unname(b)
}
strength <- function(b, dir) switch(dir, up = b, down = -b, abs(b))
wk_dist <- function(a, b) { d <- abs(a - b) %% 52; pmin(d, 52 - d) }

RES <- list(); DON <- list()
for (r in seq_len(nrow(cases))) {
  e <- cases[r, ]; o <- e$lead_outcome
  fuels <- unlist(ic$match[[e$group]]); wt <- sc$dryness_weight[[e$group]]
  exp_dir <- ic$events[[e$event]]$expected %||% "none"
  yy <- oy[[o]]; f <- series(e$eco_l2_code, fuels, o, yy); dr <- dry_series(e$eco_l2_code, wt, yy)
  D <- win_data(f, dr, yy, e$typical_week)
  if (is.null(D) || sum(D$y) <= 0 || length(unique(D$wy)) < 3) next
  m <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ wy + dry + after, data = D)), error = function(err) NULL)
  if (is.null(m)) next
  fam <- MASS::negative.binomial(m$theta)
  b_own <- unname(stats::coef(m)[["after"]]); s_own <- strength(b_own, exp_dir)
  dd <- dates[dates$group == e$group & dates$event == e$event, ]
  nat_wk <- round(stats::median(dd$typical_week))
  b_nat <- if (nat_wk == e$typical_week) NA_real_ else fit_fixed(win_data(f, dr, yy, nat_wk), fam)
  don <- sort(unique(dd$typical_week[dd$eco_l2_code != e$eco_l2_code])); don <- don[don != e$typical_week]
  b_don <- vapply(don, function(w) fit_fixed(win_data(f, dr, yy, w), fam), 0)
  s_don <- vapply(b_don, strength, 0, dir = exp_dir); dist <- wk_dist(don, e$typical_week)
  ok <- is.finite(s_don)
  bp <- vapply(bands, function(bd) { k <- ok & dist >= bd[1] & dist <= bd[2]; if (any(k)) mean(s_don[k] < s_own) else NA_real_ }, 0)
  RES[[length(RES) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o,
    expected = exp_dir, typical_week = e$typical_week, ratio_own = round(exp(b_own), 3),
    national_week = nat_wk, week_from_national = e$typical_week - nat_wk,
    ratio_national = round(exp(b_nat), 3),
    own_beats_national = if (is.finite(b_nat)) s_own > strength(b_nat, exp_dir) else NA,
    donors = sum(ok), donor_pct = if (any(ok)) round(mean(s_don[ok] < s_own), 3) else NA_real_,
    t(setNames(round(bp, 3), paste0("donor_pct_", names(bands)))),
    donors_near = sum(ok & dist <= bands[[1]][2]), stringsAsFactors = FALSE, check.names = FALSE)
  if (length(don)) DON[[length(DON) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event,
    donor_week = don, week_diff = dist, ratio = round(exp(b_don), 3), stronger_than_own = s_don >= s_own)
  if (r %% 50 == 0) log_msg(sprintf("  %d / %d cases (%.1f min)", r, nrow(cases), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!length(RES)) stop("no case could be fitted: check inventory.csv (tier2_ok), fire_weekly and eco_weekly")
R <- do.call(rbind, RES); utils::write.csv(R, file.path(out, "results.csv"), row.names = FALSE)
if (length(DON)) utils::write.csv(do.call(rbind, DON), file.path(out, "donors.csv"), row.names = FALSE)

# --- summaries: all cases, and each ecoregion x group x typical week once (first event in config order) ---
evo <- match(R$event, names(ic$events))
U <- R[order(evo), ]; U <- U[!duplicated(paste(U$eco_l2_code, U$group, U$typical_week)), ]
summ <- function(d, by, lab) {
  k <- if (is.null(by)) rep("all", nrow(d)) else d[[by]]
  do.call(rbind, lapply(split(d, as.character(k), drop = TRUE), function(s) data.frame(set = lab, level = if (is.null(by)) "all" else s[[by]][1],
    cases = nrow(s), beats_national = round(mean(s$own_beats_national, na.rm = TRUE), 3),
    national_tested = sum(!is.na(s$own_beats_national)), median_donor_pct = round(stats::median(s$donor_pct, na.rm = TRUE), 3),
    mean_donor_pct = round(mean(s$donor_pct, na.rm = TRUE), 3),
    mean_donor_pct_near = round(mean(s$donor_pct_near, na.rm = TRUE), 3),
    mean_donor_pct_mid = round(mean(s$donor_pct_mid, na.rm = TRUE), 3),
    mean_donor_pct_far = round(mean(s$donor_pct_far, na.rm = TRUE), 3), stringsAsFactors = FALSE)))
}
R$season <- cut(R$typical_week, c(0, 13, 26, 39, 53), labels = c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec"))
U$season <- cut(U$typical_week, c(0, 13, 26, 39, 53), labels = c("Jan-Mar", "Apr-Jun", "Jul-Sep", "Oct-Dec"))
S <- rbind(summ(R, NULL, "cases"), summ(U, NULL, "units"), summ(R, "season", "cases"), summ(R, "group", "cases"),
           summ(R, "event", "cases"), summ(R, "eco_l2_code", "cases"))
utils::write.csv(S, file.path(out, "summary.csv"), row.names = FALSE)

# sign tests on units (each ecoregion x group x week once): own date vs the national date; donor_pct vs 0.5
bt <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(c(n = 0, k = 0, p = NA)); c(n = length(x), k = sum(x), p = stats::binom.test(sum(x), length(x))$p.value) }
a <- bt(U$own_beats_national); b <- bt(U$donor_pct[!is.na(U$donor_pct) & U$donor_pct != 0.5] > 0.5)
log_msg(sprintf("cases fitted: %d | units (ecoregion x group x week): %d", nrow(R), nrow(U)))
log_msg(sprintf("own date beats the national date: %.0f %% of cases | units %d of %d (sign test p = %.3g)",
                100 * mean(R$own_beats_national, na.rm = TRUE), a[["k"]], a[["n"]], a[["p"]]))
log_msg(sprintf("donor_pct (share of other ecoregions' dates with a weaker change): mean %.2f cases | %.2f units | units above 0.5: %d of %d (p = %.3g)",
                mean(R$donor_pct, na.rm = TRUE), mean(U$donor_pct, na.rm = TRUE), b[["k"]], b[["n"]], b[["p"]]))
for (nm in names(bands)) log_msg(sprintf("  donors %d-%d weeks away: mean donor_pct %.2f (units %.2f)", bands[[nm]][1], bands[[nm]][2],
  mean(R[[paste0("donor_pct_", nm)]], na.rm = TRUE), mean(U[[paste0("donor_pct_", nm)]], na.rm = TRUE)))
for (s in levels(R$season)) { x <- R[R$season == s, ]; if (!nrow(x)) next
  log_msg(sprintf("  %-8s %4d cases | beats national %3.0f %% | mean donor_pct %.2f", s, nrow(x),
                  100 * mean(x$own_beats_national, na.rm = TRUE), mean(x$donor_pct, na.rm = TRUE))) }
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
