#!/usr/bin/env Rscript
# Step 09 (NPN x fire study; laptop, minutes; MASS): CONFIRMATORY test of the one pattern the exploratory search
# (08) found: in the eastern deciduous forests, when deciduous-tree leaf-out runs ahead of schedule, timber-fuel
# fires that week are lower. Does NPN leaf-out status add anything beyond spring temperature?
# Asked by Mike on 3 Nov 2026; the test is fixed in config (npn_fire: confirm) before its result was seen.
#
# Cases: config confirm: ecoregions (8.1, 8.2, 8.3, 8.4) x deciduous trees x primary event (leaves_onset);
#   secondary events reported but not part of the decision. Lead record (FPA-FOD all), timber fuels.
# Years: as 07 / 08 (>= 10 sites with a date; >= 8 years). Weeks: typical week -8..+7.
# Progress anomaly (as 08): share of the year's sites past the event at that week minus its mean over years at
#   that week, per 10 points. Temperature anomaly: mean tmean over the 8 weeks up to and including the week
#   (timber weight), minus its mean over years at that week, per 1 deg C.
# Models (NB, dispersion from the full model, then fixed):
#   T    fires ~ factor(year) + factor(week) + dryness + temperature
#   TP   T + progress                       <- the test: progress coefficient b (expected < 0)
#   P    fires ~ factor(year) + factor(week) + dryness + progress   (08's model, for comparison)
# p (one-sided, b < 0): progress anomaly curves permuted across years (whole years) n_perm times, refit TP.
# Decision (primary event): Fisher combination of the ecoregions' one-sided p < 0.05 = NPN leaf-out status adds
#   information beyond spring temperature; each ecoregion's p and ratio reported.
# Outputs <run_dir>/npn_fire/confirm/: results.csv (ecoregion x event), decision.csv.
# Usage: Rscript npn_fire/scripts/09_nf_confirm.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("09_nf_confirm", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; yc <- cfg$npn_fire$yeartoyear; cc <- cfg$npn_fire$confirm
W <- as.integer(sc$window_weeks %||% 8); min_sites <- as.integer(yc$min_sites %||% 10); min_years <- as.integer(yc$min_years %||% 8)
ecos <- unlist(cc$ecoregions %||% list("8.1", "8.2", "8.3", "8.4")); grp <- cc$group %||% "deciduous_tree"
prim <- cc$primary_event %||% "leaves_onset"; sec <- unlist(cc$secondary_events %||% list("initial_growth_onset", "young_leaves_onset"))
tw_w <- as.integer(cc$temp_weeks %||% 8); n_perm <- as.integer(cc$n_perm %||% 999); alpha <- cc$alpha %||% 0.05
out <- nf_out(cfg, "confirm")
set.seed(cc$seed %||% 20261104)
t0 <- Sys.time()
log_msg("=== 09_nf_confirm | profile: ", profile, " | CONFIRMATORY | ", grp, " ", prim, " in ", paste(ecos, collapse = ", "),
        " | temperature: ", tw_w, "-week trailing tmean anomaly | permutations ", n_perm)

INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
SD <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "inventory", "site_dates"))); SD$eco_l2_code <- as.character(SD$eco_l2_code)
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
if (!"tmean" %in% names(PW)) stop("eco_weekly.csv has no tmean: rerun 04 with force on the P720 and sync")
fuels <- unlist(ic$match[[grp]]); wt <- sc$dryness_weight[[grp]]
rel <- seq(-W, W - 1)

wk_val <- function(eco, yy, var) {
  d <- PW[PW$eco_l2_code == eco & PW$weight == wt & PW$year %in% yy, ]
  v <- rep(NA_real_, length(yy) * 52); v[(d$year - min(yy)) * 52 + d$week] <- d[[var]]
  v
}
RES <- list()
for (eco in ecos) for (ev in c(prim, sec)) {
  e <- INV[INV$eco_l2_code == eco & INV$group == grp & INV$event == ev, ]
  if (!nrow(e) || is.na(e$typical_week)) { log_msg("  ", eco, " ", ev, ": no typical date"); next }
  o <- e$lead_outcome; tw <- e$typical_week
  yy <- sort(unique(FW$year[FW$outcome == o]))
  d <- SD[SD$eco_l2_code == eco & SD$group == grp & SD$event == ev & !is.na(SD$doy), c("site_id", "year", "doy")]
  d <- d[!duplicated(d[, c("site_id", "year")]), ]
  ns <- table(d$year); yk <- intersect(as.integer(names(ns)[ns >= min_sites]), yy)
  fr <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% fuels & FW$year %in% yy, ]
  f <- numeric(length(yy) * 52); if (nrow(fr)) { s <- tapply(fr$value, (fr$year - min(yy)) * 52 + fr$week, sum); f[as.integer(names(s))] <- s }
  dr <- wk_val(eco, yy, "ppet60") / 100; tm <- wk_val(eco, yy, "tmean")
  sw <- pmin(52L, pmax(1L, as.integer(ceiling(d$doy / 7))))
  L <- do.call(rbind, lapply(yk, function(y) { i <- match(y, yy); t <- (i - 1) * 52 + tw + rel
    if (min(t) - tw_w + 1 < 1 || max(t) > length(f) || anyNA(dr[t])) return(NULL)
    s <- sw[d$year == y]
    data.frame(year = y, rel = rel, y = round(f[t]), dry = dr[t], prog = vapply(tw + rel, function(w) mean(s <= w), 0),
               temp = vapply(t, function(j) mean(tm[(j - tw_w + 1):j], na.rm = TRUE), 0)) }))
  if (is.null(L) || length(unique(L$year)) < min_years) { log_msg("  ", eco, " ", ev, ": fewer than ", min_years, " years"); next }
  L <- L[order(L$year, L$rel), ]
  L$pa <- (L$prog - stats::ave(L$prog, L$rel)) * 10; L$ta <- L$temp - stats::ave(L$temp, L$rel)
  L$fy <- factor(L$year); L$fr <- factor(L$rel)
  mTP <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ fy + fr + dry + ta + pa, data = L)), error = function(err) NULL)
  mP <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ fy + fr + dry + pa, data = L)), error = function(err) NULL)
  mT <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ fy + fr + dry + ta, data = L)), error = function(err) NULL)
  if (is.null(mTP)) { log_msg("  ", eco, " ", ev, ": model failed"); next }
  fam <- MASS::negative.binomial(mTP$theta); X <- stats::model.matrix(~ fy + fr + dry + ta + pa, L); jp <- which(colnames(X) == "pa")
  b_obs <- unname(stats::coef(mTP)[["pa"]])
  ub <- unique(L$year); blocks <- split(seq_len(nrow(L)), L$year)
  bp <- vapply(seq_len(n_perm), function(i) { sh <- sample(ub); Xp <- X
    for (k in seq_along(ub)) Xp[blocks[[as.character(ub[k])]], jp] <- L$pa[blocks[[as.character(sh[k])]]]
    m <- tryCatch(suppressWarnings(stats::glm.fit(Xp, L$y, family = fam, control = stats::glm.control(maxit = 50))), error = function(err) NULL)
    if (is.null(m) || !m$converged) NA_real_ else m$coefficients[jp] }, 0)
  bp <- bp[is.finite(bp)]; p1 <- (1 + sum(bp <= b_obs)) / (1 + length(bp))
  RES[[length(RES) + 1]] <- data.frame(eco_l2_code = eco, group = grp, event = ev, role = if (ev == prim) "primary" else "secondary",
    typical_week = tw, years = length(unique(L$year)), fires = sum(L$y), r_progress_temp = round(stats::cor(L$pa, L$ta), 2),
    ratio_progress_given_temp = round(exp(b_obs), 3), p_one_sided = signif(p1, 3),
    ratio_progress_alone = if (is.null(mP)) NA else round(exp(stats::coef(mP)[["pa"]]), 3),
    ratio_temp_alone = if (is.null(mT)) NA else round(exp(stats::coef(mT)[["ta"]]), 3),
    ratio_temp_given_progress = round(exp(stats::coef(mTP)[["ta"]]), 3), stringsAsFactors = FALSE)
  log_msg(sprintf("  %-4s %-22s %2d years | progress ratio per 10 points ahead %.3f given temperature (alone %.3f) | one-sided p %.3g | r(progress, temperature) %.2f | temperature ratio per 1 C %.3f alone, %.3f given progress",
    eco, ev, length(unique(L$year)), exp(b_obs), RES[[length(RES)]]$ratio_progress_alone, p1, RES[[length(RES)]]$r_progress_temp,
    RES[[length(RES)]]$ratio_temp_alone, RES[[length(RES)]]$ratio_temp_given_progress))
}
if (!length(RES)) stop("no case could be fitted")
R <- do.call(rbind, RES); utils::write.csv(R, file.path(out, "results.csv"), row.names = FALSE)
fisher <- function(p) { p <- p[is.finite(p)]; if (!length(p)) return(NA_real_); stats::pchisq(-2 * sum(log(p)), 2 * length(p), lower.tail = FALSE) }
D <- do.call(rbind, lapply(split(R, R$event), function(d) data.frame(event = d$event[1], role = d$role[1], ecoregions = nrow(d),
  fewer_fires_when_ahead = sum(d$ratio_progress_given_temp < 1), combined_p = signif(fisher(d$p_one_sided), 3),
  decision = if (d$role[1] == "primary") (if (fisher(d$p_one_sided) < alpha) "NPN adds beyond temperature" else "not shown") else "secondary (not decided)")))
utils::write.csv(D, file.path(out, "decision.csv"), row.names = FALSE)
for (i in seq_len(nrow(D))) log_msg(sprintf("DECISION %-9s %-22s %d ecoregions | fewer fires when ahead in %d | Fisher combined one-sided p %.3g | %s",
  D$role[i], D$event[i], D$ecoregions[i], D$fewer_fires_when_ahead[i], D$combined_p[i], D$decision[i]))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
