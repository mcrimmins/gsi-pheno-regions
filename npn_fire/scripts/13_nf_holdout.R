#!/usr/bin/env Rscript
# Step 13 (NPN x fire study; laptop or P720, minutes; MASS): roadmap phase 4, the ONE-TIME check on the held-out
# years (config roadmap: holdout_years, 2021-2024) of the canopy-fullness results of 12. Tests and decision
# rules fixed in config (npn_fire: holdout) on 4 Nov 2026 by Mike, before any hold-out result was seen.
# The machine-learning checks (tests 2-3) run in 11 with `final`.
#
# Test 1  deciduous-tree canopy fullness (leaves, class 3, mean percent) vs timber-fuel FPA-FOD fires in
#         8.1, 8.3, 8.4. Test 1b  the same canopy measure vs fire size given an ignition (share of fires
#         reaching 10 acres) in 8.3, 8.4.
# Same construction as 12: weeks = the weeks 12 used for the series (from the exploration years); canopy
#   anomaly = hold-out value minus the exploration-year mean for that week (per 10 points; missing = 0);
#   model = factor(year) + factor(week) + dryness + 8-week temperature anomaly + canopy anomaly (NB for counts,
#   binomial with dispersion for size). One-sided p (fewer / smaller fires): LR test, halved when the
#   coefficient is negative (1 - p/2 otherwise). With 4 years a permutation over years has only 24 orders; it is
#   reported as a check, not used for the decision.
# Decision: Fisher-combined one-sided p over the listed ecoregions < alpha.
# Outputs <run_dir>/npn_fire/holdout/: canopy.csv, decision.csv.
# Usage: Rscript npn_fire/scripts/13_nf_holdout.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("13_nf_holdout", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; rc <- cfg$npn_fire$roadmap; xc <- cfg$npn_fire$intensity; hc <- cfg$npn_fire$holdout
ey <- unlist(rc$explore_years); ho <- unlist(rc$holdout_years); ex_y <- seq(ey[1], ey[2]); ho_y <- seq(ho[1], ho[2])
min_sites <- as.integer(xc$min_sites %||% 5); share_years <- xc$share_years %||% 0.6; tw_w <- as.integer(xc$temp_weeks %||% 8)
alpha <- hc$alpha %||% 0.05; g <- hc$group %||% "deciduous_tree"
tests <- hc$tests %||% list(canopy_fires = list(outcome = "fod_all", ecoregions = list("8.1", "8.3", "8.4")),
                             canopy_size = list(outcome = "size", ecoregions = list("8.3", "8.4")))
out <- nf_out(cfg, "holdout")
t0 <- Sys.time()
log_msg("=== 13_nf_holdout | profile: ", profile, " | HOLD-OUT ", min(ho_y), "-", max(ho_y), " (run once) | canopy fullness, ", g)

EW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "weekly", "eco_weekly"))); EW$eco_l2_code <- as.character(EW$eco_l2_code)
EW <- EW[EW$group == g & EW$pheno_class_id == 3 & EW$sites_pct >= min_sites & !is.na(EW$mean_pct), ]
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
FW <- FW[FW$year %in% ho_y, ]
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
fuels <- unlist(ic$match[[g]]); wt <- sc$dryness_weight[[g]]
wkv <- function(eco, var) { d <- PW[PW$eco_l2_code == eco & PW$weight == wt, ]; setNames(d[[var]], paste(d$year, d$week)) }
fire_wk <- function(eco, o) { d <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% fuels, ]
  if (!nrow(d)) return(setNames(numeric(0), character(0))); s <- tapply(d$value, paste(d$year, d$week), sum); setNames(as.numeric(s), names(s)) }

one <- function(eco, o) {
  mx <- EW[EW$eco_l2_code == eco & EW$year %in% ex_y, ]
  ny <- length(unique(mx$year)); wk_share <- tapply(mx$year, mx$week, function(v) length(unique(v))) / ny
  wks <- as.integer(names(wk_share)[wk_share >= share_years])                  # the weeks 12 used (exploration years)
  clim <- tapply(mx$mean_pct[mx$week %in% wks], mx$week[mx$week %in% wks], mean)
  mh <- EW[EW$eco_l2_code == eco & EW$year %in% ho_y & EW$week %in% wks, ]
  G <- expand.grid(week = wks, year = ho_y); key <- paste(G$year, G$week)
  v <- setNames(mh$mean_pct, paste(mh$year, mh$week))[key]
  G$x <- (v - as.numeric(clim[as.character(G$week)])) / 10; G$has <- !is.na(G$x); G$x[is.na(G$x)] <- 0
  dr <- wkv(eco, "ppet60"); tmv <- wkv(eco, "tmean")
  G$dry <- as.numeric(dr[key]) / 100
  trail <- vapply(seq_len(nrow(G)), function(i) { w <- G$week[i] - 0:(tw_w - 1); y <- G$year[i] - (w < 1); w[w < 1] <- w[w < 1] + 52
    mean(tmv[paste(y, w)], na.rm = TRUE) }, 0)
  G$ta <- trail - stats::ave(trail, G$week, FUN = function(z) mean(z, na.rm = TRUE))
  G <- G[is.finite(G$dry) & is.finite(G$ta), ]; key <- paste(G$year, G$week)
  G$fy <- factor(G$year); G$fw <- factor(G$week)
  binom <- o == "size"
  if (binom) { a <- fire_wk(eco, "fod_all"); cp <- fire_wk(eco, "fod_cplus")
    G$n <- round(ifelse(is.na(a[key]), 0, a[key])); G$k <- pmin(round(ifelse(is.na(cp[key]), 0, cp[key])), G$n); G <- G[G$n > 0, ]
  } else { f <- fire_wk(eco, o); G$y <- round(ifelse(is.na(f[key]), 0, f[key])) }
  f1 <- if (binom) cbind(k, n - k) ~ fy + fw + dry + ta + x else y ~ fy + fw + dry + ta + x
  f0 <- stats::update(f1, . ~ . - x)
  if (binom) { m1 <- stats::glm(f1, data = G, family = stats::quasibinomial()); m0 <- stats::glm(f0, data = G, family = stats::quasibinomial())
    disp <- max(1, summary(m1)$dispersion); lr <- max(m0$deviance - m1$deviance, 0) / disp
  } else { m1 <- suppressWarnings(MASS::glm.nb(f1, data = G)); m0 <- suppressWarnings(stats::glm(f0, data = G, family = MASS::negative.binomial(m1$theta)))
    m1b <- suppressWarnings(stats::glm(f1, data = G, family = MASS::negative.binomial(m1$theta))); lr <- max(m0$deviance - m1b$deviance, 0) }
  b <- unname(stats::coef(m1)[["x"]]); p2 <- stats::pchisq(lr, 1, lower.tail = FALSE); p1 <- if (b < 0) p2 / 2 else 1 - p2 / 2
  # permutation over the 4! year orders (check only): year a gets year b's canopy anomaly, week by week
  xf <- setNames(G$x, key); perms <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
  bs <- vapply(perms(ho_y), function(pm) { map <- setNames(pm, ho_y); Gp <- G; xv <- xf[paste(map[as.character(G$year)], G$week)]; xv[is.na(xv)] <- 0; Gp$x <- xv
    mm <- tryCatch(suppressWarnings(stats::glm(f1, data = Gp, family = if (binom) stats::quasibinomial() else MASS::negative.binomial(m1$theta))), error = function(e) NULL)
    if (is.null(mm)) NA_real_ else stats::coef(mm)[["x"]] }, 0)
  data.frame(eco_l2_code = eco, outcome = o, years = length(unique(G$year)), weeks = length(unique(G$week)), share_weeks_with_canopy = round(mean(G$has), 2),
             fires = if (binom) sum(G$n) else sum(G$y), ratio_per_10 = round(exp(b), 3), p_one_sided = signif(p1, 3),
             perm_rank = sum(bs <= b, na.rm = TRUE), perm_n = sum(is.finite(bs)), stringsAsFactors = FALSE)
}

RES <- list(); DEC <- list()
for (tn in names(tests)) {
  tt <- tests[[tn]]
  rr <- do.call(rbind, lapply(unlist(tt$ecoregions), function(eco) tryCatch(one(eco, tt$outcome), error = function(e) { log_err(tn, " ", eco, ": ", conditionMessage(e)); NULL })))
  if (is.null(rr)) next
  rr$test <- tn; RES[[tn]] <- rr
  pf <- stats::pchisq(-2 * sum(log(rr$p_one_sided)), 2 * nrow(rr), lower.tail = FALSE)
  DEC[[tn]] <- data.frame(test = tn, outcome = tt$outcome, ecoregions = nrow(rr), fewer = sum(rr$ratio_per_10 < 1), combined_p = signif(pf, 3),
                          decision = if (pf < alpha) "confirmed on the hold-out years" else "not confirmed")
  for (i in seq_len(nrow(rr))) log_msg(sprintf("  %-13s %-4s %d years | ratio per 10 points %.3f | one-sided p %.3g | permutation rank %d of %d | canopy data in %.0f %% of weeks",
    tn, rr$eco_l2_code[i], rr$years[i], rr$ratio_per_10[i], rr$p_one_sided[i], rr$perm_rank[i], rr$perm_n[i], 100 * rr$share_weeks_with_canopy[i]))
  log_msg(sprintf("DECISION %-13s %d ecoregions | fewer in %d | Fisher one-sided p %.3g | %s", tn, nrow(rr), DEC[[tn]]$fewer, pf, DEC[[tn]]$decision))
}
if (!length(RES)) stop("no hold-out test could be run")
utils::write.csv(do.call(rbind, RES), file.path(out, "canopy.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, DEC), file.path(out, "decision.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
