#!/usr/bin/env Rscript
# Step 08 (NPN x fire study; laptop, ~30-60 min; MASS): EXPLORATORY search for NPN-fire relationships beyond
# steps 2-3. Asked by Mike on 3 Nov 2026 ("keep digging"); angles and rules chosen by Mike and fixed in config
# (npn_fire: explore) before any result was seen. Everything here is hypothesis-generating: results are
# flagged by a Benjamini-Hochberg false-discovery rate within each family of tests.
#
# Cases: 05 tier3_ok rows (as 07). Outcomes (config explore: outcomes): all-cause, lightning and human FPA-FOD
# fires, FPA-FOD 10 acres and up (fires and acres), MTBS fires and area, VIIRS pixel-days; matched fuel.
# Years: as 07 (>= 10 sites with >= 2 years each; >= 8 years); a case x outcome needs >= 5 years with >= 3
# fires of the outcome (or its count companion for area) in the +-8 week window.
# Family 'size' (does an early or late year change how much fire, not when):
#   per year: fires after the typical date (weeks 0..+7) and in the whole window (-8..+7);
#   NB: total ~ year (linear trend) + window dryness (ppet60 mean) + departure (per 7 days); LR test for the
#   departure; effect = ratio per week later. Expected sign: events expected 'down' (fewer fires after, e.g.
#   leaf-out) -> later year, more fire (ratio > 1); 'up' (curing, leaf fall, flowering) -> later, less fire.
# Family 'progress' (in-season status): weekly share of the year's sites already past the event
#   (site date's week <= week), minus its mean over years at that week = progress anomaly (per 10 points);
#   NB: weekly fires ~ factor(year) + factor(week in window) + dryness + progress anomaly; LR test.
#   Expected sign: 'down' events -> ahead of schedule, fewer fires (ratio < 1); 'up' -> more fires.
# Family 'reruns': if 06 / 07 were run with outcome=<x> (seasonal_<x>/, yeartoyear_<x>/), their pass counts
#   against the chance rate (binomial test).
# Outputs <run_dir>/npn_fire/explore/: size.csv, progress.csv (case x outcome: effect, p, q, expected,
#   direction), reruns.csv, summary.csv (family x outcome).
# Usage: Rscript npn_fire/scripts/08_nf_explore.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("08_nf_explore", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; yc <- cfg$npn_fire$yeartoyear; xc <- cfg$npn_fire$explore
W <- as.integer(sc$window_weeks %||% 8); min_sites <- as.integer(yc$min_sites %||% 10); min_sy <- as.integer(yc$min_site_years %||% 2)
min_years <- as.integer(yc$min_years %||% 8); fdr <- xc$fdr_q %||% 0.1
screen_p <- xc$screen_p %||% 0.05; n_perm <- as.integer(xc$n_perm %||% 999)
set.seed(xc$seed %||% 20261103)
outs <- unlist(xc$outcomes %||% list("fod_all", "fod_all_lightning", "fod_all_human", "fod_cplus", "fod_cplus_acres", "mtbs_fires", "mtbs_km2", "viirs_px"))
count_of <- c(fod_cplus_acres = "fod_cplus", mtbs_km2 = "mtbs_fires")
out <- nf_out(cfg, "explore")
t0 <- Sys.time()
log_msg("=== 08_nf_explore | profile: ", profile, " | EXPLORATORY | outcomes: ", paste(outs, collapse = ", "), " | FDR q ", fdr)

INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
SD <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "inventory", "site_dates"))); SD$eco_l2_code <- as.character(SD$eco_l2_code)
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))
outs <- intersect(outs, names(oy))
cases <- INV[INV$tier3_ok %in% TRUE, ]
mx <- suppressWarnings(as.integer(Sys.getenv("NF_MAX_CASES", ""))); if (is.finite(mx)) cases <- cases[sample(nrow(cases), min(mx, nrow(cases))), ]  # testing only
log_msg(sprintf("cases: %d (Tier 3 possible) x %d outcomes", nrow(cases), length(outs)))

wrap <- function(d) ifelse(d > 182, d - 365, ifelse(d < -182, d + 365, d))
rel <- seq(-W, W - 1)
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
# LR test of the extra term x (NB; dispersion from the full model, then fixed). Fires within a year are not
# independent (one large fire, many pixel-days), so tests with parametric p < screen_p are re-tested by
# permutation: x is permuted across years (whole yearly blocks), n_perm times; p = permutation p.
nb_lr <- function(f1, f0, D, perm_block) {
  m1 <- tryCatch(suppressWarnings(MASS::glm.nb(f1, data = D)), error = function(e) NULL)
  if (is.null(m1) || !is.finite(m1$theta)) return(c(b = NA, p = NA, p_param = NA))
  fam <- MASS::negative.binomial(m1$theta)
  X1 <- stats::model.matrix(f1, D); X0 <- stats::model.matrix(f0, D); y <- stats::model.response(stats::model.frame(f1, D))
  dev <- function(X) { m <- tryCatch(suppressWarnings(stats::glm.fit(X, y, family = fam, control = stats::glm.control(maxit = 50))), error = function(e) NULL)
    if (is.null(m) || !m$converged) NA_real_ else m$deviance }
  d0 <- dev(X0); d1 <- dev(X1); if (!is.finite(d0) || !is.finite(d1)) return(c(b = NA, p = NA, p_param = NA))
  b <- stats::coef(m1)[["x"]]; lr <- max(d0 - d1, 0); pp <- stats::pchisq(lr, 1, lower.tail = FALSE); p <- pp
  if (is.finite(pp) && pp < screen_p) {
    jx <- which(colnames(X1) == "x"); blk <- perm_block; ub <- unique(blk)
    lrp <- vapply(seq_len(n_perm), function(i) { map <- setNames(sample(ub), ub)
      Xp <- X1; Xp[, jx] <- D$x[order_within(blk, map[as.character(blk)])]; dp <- dev(Xp); if (is.finite(dp)) max(d0 - dp, 0) else NA_real_ }, 0)
    lrp <- lrp[is.finite(lrp)]; p <- (1 + sum(lrp >= lr)) / (1 + length(lrp))
  }
  c(b = unname(b), p = p, p_param = pp)
}
# row index that moves block a's x values onto block map[a]'s rows (blocks of equal length, same within-block order)
order_within <- function(blk, target) {
  idx <- seq_along(blk); src <- split(idx, blk)
  out <- integer(length(blk)); for (a in names(src)) out[src[[as.character(target[src[[a]][1]])]]] <- src[[a]]
  out
}
dir_ok <- function(b, expd, fam) { if (is.na(b) || expd == "none") return(NA)
  want_pos <- if (fam == "size") expd == "down" else expd == "up"; if (want_pos) b > 0 else b < 0 }

SZ <- list(); PG <- list()
for (r in seq_len(nrow(cases))) {
  e <- cases[r, ]; tw <- e$typical_week
  fuels <- unlist(ic$match[[e$group]]); wt <- sc$dryness_weight[[e$group]]
  expd <- ic$events[[e$event]]$expected %||% "none"
  d <- SD[SD$eco_l2_code == e$eco_l2_code & SD$group == e$group & SD$event == e$event & !is.na(SD$doy), c("site_id", "year", "doy")]
  d <- d[!duplicated(d[, c("site_id", "year")]), ]
  ny <- tapply(d$year, d$site_id, function(v) length(unique(v)))
  dd <- d[d$site_id %in% names(ny)[ny >= min_sy], ]
  dd$dep <- wrap(dd$doy - stats::ave(dd$doy, dd$site_id, FUN = stats::median))
  A <- stats::aggregate(dep ~ year, dd, mean); A$sites <- as.integer(table(dd$year)[as.character(A$year)])
  A <- A[A$sites >= min_sites, ]
  for (o in outs) {
    yy <- oy[[o]]; yk <- A$year[A$year %in% yy]
    if (length(yk) < min_years) next
    f <- series(e$eco_l2_code, fuels, o, yy); dr <- dry_series(e$eco_l2_code, wt, yy)
    oc <- if (o %in% names(count_of)) count_of[[o]] else o; fc <- if (oc == o) f else series(e$eco_l2_code, fuels, oc, yy)
    idx <- lapply(yk, function(y) (match(y, yy) - 1) * 52 + tw + rel)
    okw <- vapply(idx, function(t) min(t) >= 1 && max(t) <= length(f) && !anyNA(dr[t]), TRUE)
    yk <- yk[okw]; idx <- idx[okw]
    if (length(yk) < min_years) next
    if (sum(vapply(idx, function(t) sum(fc[t]) >= 3, TRUE)) < 5) next
    # --- size family ---
    Y <- data.frame(year = yk, dep = A$dep[match(yk, A$year)] / 7,
                    after = vapply(idx, function(t) round(sum(f[t[rel >= 0]])), 0), total = vapply(idx, function(t) round(sum(f[t])), 0),
                    dry = vapply(idx, function(t) mean(dr[t]), 0))
    Y$x <- Y$dep; Y$yr <- Y$year - mean(Y$year)
    for (tot in c("after", "total")) {
      Y$y <- Y[[tot]]; if (sum(Y$y) <= 0) next
      t1 <- nb_lr(y ~ yr + dry + x, y ~ yr + dry, Y, Y$year)
      SZ[[length(SZ) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o, measure = tot,
        expected = expd, years = nrow(Y), fires = sum(Y$y), ratio_per_week_later = round(exp(t1[["b"]]), 3), p = signif(t1[["p"]], 3), p_param = signif(t1[["p_param"]], 3),
        expected_direction = dir_ok(t1[["b"]], expd, "size"), stringsAsFactors = FALSE)
    }
    # --- progress family ---
    sw <- pmin(52L, pmax(1L, as.integer(ceiling(d$doy / 7))))
    L <- do.call(rbind, lapply(seq_along(yk), function(j) { s <- sw[d$year == yk[j]]; n <- length(s); if (n < min_sites) return(NULL)
      wk <- tw + rel   # same-year weeks; before week 1 no site is past, after week 52 all are
      data.frame(year = yk[j], rel = rel, y = round(f[idx[[j]]]), dry = dr[idx[[j]]], prog = vapply(wk, function(w) mean(s <= w), 0)) }))
    if (is.null(L) || length(unique(L$year)) < min_years || sum(L$y) <= 0) next
    L$x <- (L$prog - stats::ave(L$prog, L$rel)) * 10
    if (stats::sd(L$x) == 0) next
    L$fy <- factor(L$year); L$fr <- factor(L$rel)
    L <- L[order(L$year, L$rel), ]
    t2 <- nb_lr(y ~ fy + fr + dry + x, y ~ fy + fr + dry, L, L$year)
    PG[[length(PG) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, outcome = o,
      expected = expd, years = length(unique(L$year)), fires = sum(L$y), ratio_per_10pt_ahead = round(exp(t2[["b"]]), 3),
      p = signif(t2[["p"]], 3), p_param = signif(t2[["p_param"]], 3), expected_direction = dir_ok(t2[["b"]], expd, "progress"), stringsAsFactors = FALSE)
  }
  if (r %% 25 == 0) log_msg(sprintf("  %d / %d cases (%.1f min)", r, nrow(cases), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
bh <- function(p) stats::p.adjust(p, "BH")
S1 <- do.call(rbind, SZ); S1$q <- signif(bh(S1$p), 3); utils::write.csv(S1, file.path(out, "size.csv"), row.names = FALSE)
S2 <- do.call(rbind, PG); S2$q <- signif(bh(S2$p), 3); utils::write.csv(S2, file.path(out, "progress.csv"), row.names = FALSE)

# --- reruns of 06 / 07 with other outcomes (if present) ---
RR <- list()
for (o in outs) {
  f6 <- file.path(cfg$out_dir, paste0("seasonal_", o), "results.csv"); f7 <- file.path(cfg$out_dir, paste0("yeartoyear_", o), "results.csv")
  if (file.exists(f6)) { x <- utils::read.csv(f6, stringsAsFactors = FALSE); x <- x[!is.na(x$pass), ]
    ch <- ifelse(x$expected == "none", 2 * sc$pass_tail, sc$pass_tail)
    RR[[length(RR) + 1]] <- data.frame(family = "rerun_seasonal", outcome = o, tests = nrow(x), pass = sum(x$pass), chance = round(sum(ch), 1),
      p = signif(stats::poisson.test(sum(x$pass), sum(ch), alternative = "greater")$p.value, 3)) }
  if (file.exists(f7)) { x <- utils::read.csv(f7, stringsAsFactors = FALSE)
    RR[[length(RR) + 1]] <- data.frame(family = "rerun_yeartoyear", outcome = o, tests = nrow(x), pass = sum(x$pass), chance = round(0.05 * nrow(x), 1),
      p = signif(stats::poisson.test(sum(x$pass), 0.05 * nrow(x), alternative = "greater")$p.value, 3)) }
}
RR <- if (length(RR)) do.call(rbind, RR) else NULL
if (!is.null(RR)) utils::write.csv(RR, file.path(out, "reruns.csv"), row.names = FALSE)

fam_sum <- function(X, fam) do.call(rbind, lapply(split(X, X$outcome), function(x) data.frame(family = fam, outcome = x$outcome[1],
  tests = nrow(x), p05 = sum(x$p < 0.05, na.rm = TRUE), expected_p05 = round(0.05 * sum(!is.na(x$p)), 1), q_hits = sum(x$q < fdr, na.rm = TRUE),
  q_hits_expected_dir = sum(x$q < fdr & x$expected_direction %in% TRUE, na.rm = TRUE),
  q_hits_wrong_dir = sum(x$q < fdr & x$expected_direction %in% FALSE, na.rm = TRUE))))
SM <- rbind(fam_sum(transform(S1, outcome = paste(outcome, measure)), "size"), fam_sum(S2, "progress"))
utils::write.csv(SM, file.path(out, "summary.csv"), row.names = FALSE)

log_msg(sprintf("size tests: %d | progress tests: %d", nrow(S1), nrow(S2)))
for (i in seq_len(nrow(SM))) log_msg(sprintf("  %-9s %-26s %4d tests | p < 0.05: %3d (chance %5.1f) | FDR q < %.2f: %3d (expected dir %d, wrong dir %d)",
  SM$family[i], SM$outcome[i], SM$tests[i], SM$p05[i], SM$expected_p05[i], fdr, SM$q_hits[i], SM$q_hits_expected_dir[i], SM$q_hits_wrong_dir[i]))
if (!is.null(RR)) for (i in seq_len(nrow(RR))) log_msg(sprintf("  %-17s %-18s %4d tests | pass %3d (chance %5.1f) | p = %.3g", RR$family[i], RR$outcome[i],
  RR$tests[i], RR$pass[i], RR$chance[i], RR$p[i]))
top <- function(X, lab, eff) { h <- X[X$q < fdr & !is.na(X$q), ]; h <- h[order(h$q), ]
  for (i in seq_len(min(nrow(h), 25))) log_msg(sprintf("  %s hit: %s %s %s %s%s | %s %.2f | q %.3g | expected %s", lab, h$eco_l2_code[i], h$group[i], h$event[i],
    h$outcome[i], if ("measure" %in% names(h)) paste0(" ", h$measure[i]) else "", eff, h[[eff]][i], h$q[i], h$expected[i])) }
top(S1, "size", "ratio_per_week_later"); top(S2, "progress", "ratio_per_10pt_ahead")
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
