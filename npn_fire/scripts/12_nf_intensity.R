#!/usr/bin/env Rscript
# Step 12 (NPN x fire study; P720 or laptop, ~20-40 min; MASS): roadmap phase 2, intensity tests. Does the
# weekly NPN status / intensity anomaly go with weekly fire, beyond temperature and dryness? EXPLORATORY,
# on the exploration years only (config roadmap: explore_years; 2021-2024 are never read here).
# Fixed before results (Mike, 3 Nov 2026; roadmap doc claude/gsi-npn-fire-roadmap.md; config npn_fire: intensity).
#
# Measures (from 10's eco_weekly), per ecoregion x plant group x week, as anomalies (value minus that week's
# mean over the years used):
#   canopy   leaves (class 3) mean percent: deciduous / evergreen trees, conifers, deciduous / evergreen shrubs
#   green    grasses: leaves (class 3) mean percent (% green); forbs: share with leaves (no percent recorded)
#   index    share with leaves minus share with colored or falling leaves (classes 3 - 4 - 5), every group
# Rows: years x weeks where the measure has >= min_sites sites (grasses / forbs: min_sites_herb) in >= 60 % of
#   the years; missing anomalies set to 0 (no information); >= 6 years and >= 8 weeks.
# Outcomes (matched fuel, as 03 / 05): config roadmap: outcomes (fire counts / area, NB), plus 'size' = share of
#   FPA-FOD fires reaching 10 acres (binomial, fod_cplus of fod_all), if roadmap: fire_size_given_ignition.
# Model: outcome ~ factor(year) + factor(week) + dryness + 8-week temperature anomaly + NPN anomaly (per 10
#   points); LR test; tests with p < 0.05 re-tested by permuting the anomaly across whole years (999);
#   BH q within each family (canopy, green, index) x outcome type. Expected: greener than usual -> fewer fires
#   and smaller fires (ratio < 1).
# Outputs <run_dir>/npn_fire/intensity/: results.csv, summary.csv.
# Usage: Rscript npn_fire/scripts/12_nf_intensity.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("12_nf_intensity", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; rc <- cfg$npn_fire$roadmap; xc <- cfg$npn_fire$intensity
ey <- unlist(rc$explore_years %||% list(2009, 2020)); yrs_use <- seq(ey[1], ey[2])
outs <- unlist(rc$outcomes %||% list("fod_all")); do_size <- isTRUE(rc$fire_size_given_ignition %||% TRUE)
min_sites <- as.integer(xc$min_sites %||% 5); min_herb <- as.integer(xc$min_sites_herb %||% 3)
share_years <- xc$share_years %||% 0.6; min_years <- as.integer(xc$min_years %||% 6); min_weeks <- as.integer(xc$min_weeks %||% 8)
tw_w <- as.integer(xc$temp_weeks %||% 8); screen_p <- xc$screen_p %||% 0.05; n_perm <- as.integer(xc$n_perm %||% 999); fdr <- xc$fdr_q %||% 0.1
out <- nf_out(cfg, "intensity")
set.seed(xc$seed %||% 20261104)
t0 <- Sys.time()
log_msg("=== 12_nf_intensity | profile: ", profile, " | EXPLORATORY | years ", ey[1], "-", ey[2], " (hold-out never read) | outcomes: ",
        paste(c(outs, if (do_size) "size"), collapse = ", "))

EW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "weekly", "eco_weekly"))); EW$eco_l2_code <- as.character(EW$eco_l2_code)
EW <- EW[EW$year %in% yrs_use, ]
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
FW <- FW[FW$year %in% yrs_use, ]
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
if (!"tmean" %in% names(PW)) stop("prism/eco_weekly.csv has no tmean: rerun 04 with force and sync")
PW <- PW[PW$year %in% (min(yrs_use) - 1):max(yrs_use), ]

# --- NPN measures: eco x group x year x week value + sites ---
pick <- function(g, cls, var, nmin) { d <- EW[EW$group == g & EW$pheno_class_id == cls, ]
  d$n <- if (var == "mean_pct") d$sites_pct else d$sites; d$v <- d[[var]]; d[d$n >= nmin & !is.na(d$v), c("eco_l2_code", "year", "week", "v")] }
woody <- c("deciduous_tree", "evergreen_tree", "conifer", "deciduous_shrub", "evergreen_shrub")
M <- list()
for (g in woody) M[[length(M) + 1]] <- cbind(family = "canopy", group = g, pick(g, 3, "mean_pct", min_sites))
M[[length(M) + 1]] <- cbind(family = "green", group = "grass", pick("grass", 3, "mean_pct", min_herb))
M[[length(M) + 1]] <- cbind(family = "green", group = "forb", { x <- pick("forb", 3, "share_yes", min_herb); x$v <- x$v * 100; x })
for (g in c(woody, "grass", "forb")) {
  nm <- if (g %in% c("grass", "forb")) min_herb else min_sites
  a <- pick(g, 3, "share_yes", nm); b4 <- pick(g, 4, "share_yes", nm); b5 <- pick(g, 5, "share_yes", nm)
  k <- function(x) paste(x$eco_l2_code, x$year, x$week)
  a$v <- 100 * (a$v - ifelse(is.na(match(k(a), k(b4))), 0, b4$v[match(k(a), k(b4))]) - ifelse(is.na(match(k(a), k(b5))), 0, b5$v[match(k(a), k(b5))]))
  M[[length(M) + 1]] <- cbind(family = "index", group = g, a)
}
M <- do.call(rbind, M); M$family <- as.character(M$family); M$group <- as.character(M$group)
log_msg(sprintf("NPN measure rows (ecoregion x group x year x week with enough sites): %d", nrow(M)))

fuel_of <- function(g) unlist(ic$match[[g]] %||% ic$match[["deciduous_shrub"]])
wt_of <- function(g) sc$dryness_weight[[g]] %||% "all"
wkv <- function(eco, wt, var) { d <- PW[PW$eco_l2_code == eco & PW$weight == wt, ]; setNames(d[[var]], paste(d$year, d$week)) }
fire_wk <- function(eco, fuels, o) { d <- FW[FW$eco_l2_code == eco & FW$outcome == o & FW$fuel %in% fuels, ]
  if (!nrow(d)) return(setNames(numeric(0), character(0))); s <- tapply(d$value, paste(d$year, d$week), sum); setNames(as.numeric(s), names(s)) }

lr_test <- function(D, binom, xfull) {   # xfull: anomaly on the complete year x week grid, named 'year week'
  f1 <- if (binom) cbind(k, nmk) ~ fy + fw + dry + ta + x else y ~ fy + fw + dry + ta + x
  f0 <- if (binom) cbind(k, nmk) ~ fy + fw + dry + ta else y ~ fy + fw + dry + ta
  wv <- NULL
  if (binom) { fam <- stats::binomial(); yv <- D$k / D$n; wv <- D$n } else {
    m1 <- tryCatch(suppressWarnings(MASS::glm.nb(f1, data = D)), error = function(e) NULL)
    if (is.null(m1) || !is.finite(m1$theta)) return(c(b = NA, p = NA, p_param = NA)); fam <- MASS::negative.binomial(m1$theta); yv <- D$y }
  X1 <- stats::model.matrix(f1, D); X0 <- stats::model.matrix(f0, D); jx <- which(colnames(X1) == "x")
  fit <- function(X) { m <- tryCatch(suppressWarnings(stats::glm.fit(X, yv, weights = wv, family = fam, control = stats::glm.control(maxit = 50))), error = function(e) NULL)
    if (is.null(m) || !m$converged) NULL else m }
  m0 <- fit(X0); m1f <- fit(X1); if (is.null(m0) || is.null(m1f)) return(c(b = NA, p = NA, p_param = NA))
  disp <- if (binom) { mu <- m1f$fitted.values; max(1, sum(m1f$prior.weights * (m1f$y - mu)^2 / pmax(mu * (1 - mu), 1e-9)) / max(m1f$df.residual, 1)) } else 1
  lr <- max(m0$deviance - m1f$deviance, 0) / disp; pp <- stats::pchisq(lr, 1, lower.tail = FALSE); p <- pp; b <- m1f$coefficients[jx]
  if (is.finite(pp) && pp < screen_p) {
    yrs <- sort(unique(as.integer(sub(" .*", "", names(xfull)))))
    lrp <- vapply(seq_len(n_perm), function(i) { map <- setNames(sample(yrs), yrs); Xp <- X1   # year a gets year map[a]'s anomaly, week by week
      xv <- xfull[paste(map[as.character(D$year)], D$week)]; xv[is.na(xv)] <- 0; Xp[, jx] <- xv
      m <- fit(Xp); if (is.null(m)) NA_real_ else max(m0$deviance - m$deviance, 0) / disp }, 0)
    lrp <- lrp[is.finite(lrp)]; p <- (1 + sum(lrp >= lr)) / (1 + length(lrp))
  }
  c(b = unname(b), p = p, p_param = pp)
}

RES <- list()
keys <- unique(M[, c("family", "group", "eco_l2_code")])
log_msg(sprintf("measure x ecoregion series: %d", nrow(keys)))
for (r in seq_len(nrow(keys))) {
  kf <- keys$family[r]; g <- keys$group[r]; eco <- keys$eco_l2_code[r]
  m <- M[M$family == kf & M$group == g & M$eco_l2_code == eco, ]
  ny <- length(unique(m$year)); if (ny < min_years) next
  wk_share <- tapply(m$year, m$week, function(v) length(unique(v))) / ny
  wks <- as.integer(names(wk_share)[wk_share >= share_years]); if (length(wks) < min_weeks) next
  yk <- sort(unique(m$year[m$week %in% wks])); yk <- yk[vapply(yk, function(y) mean(wks %in% m$week[m$year == y]) >= share_years, TRUE)]
  if (length(yk) < min_years) next
  G <- expand.grid(week = wks, year = yk); key <- paste(G$year, G$week)
  v <- setNames(m$v, paste(m$year, m$week))[key]
  clim <- tapply(v, G$week, mean, na.rm = TRUE); G$x <- (v - clim[as.character(G$week)]) / 10; G$x[is.na(G$x)] <- 0
  if (stats::sd(G$x) == 0) next
  fuels <- fuel_of(g); wt <- wt_of(g)
  dr <- wkv(eco, wt, "ppet60"); tmv <- wkv(eco, wt, "tmean")
  G$dry <- as.numeric(dr[key]) / 100
  # 8-week trailing temperature anomaly (weeks cross the year boundary through year - 1)
  trail <- vapply(seq_len(nrow(G)), function(i) { w <- G$week[i] - 0:(tw_w - 1); y <- G$year[i] - (w < 1); w[w < 1] <- w[w < 1] + 52
    mean(tmv[paste(y, w)], na.rm = TRUE) }, 0)
  G$ta <- trail - stats::ave(trail, G$week, FUN = function(z) mean(z, na.rm = TRUE))
  G$fy <- factor(G$year); G$fw <- factor(G$week)
  xfull <- setNames(G$x, paste(G$year, G$week))
  ok <- is.finite(G$dry) & is.finite(G$ta); G <- G[ok, ]; if (nrow(G) < min_years * min_weeks / 2) next
  key <- paste(G$year, G$week)
  for (o in c(outs, if (do_size) "size")) {
    D <- G
    if (o == "size") { a <- fire_wk(eco, fuels, "fod_all"); cp <- fire_wk(eco, fuels, "fod_cplus")
      D$n <- round(ifelse(is.na(a[key]), 0, a[key])); D$k <- pmin(round(ifelse(is.na(cp[key]), 0, cp[key])), D$n); D$nmk <- D$n - D$k
      D <- D[D$n > 0, ]; if (nrow(D) < 20 || sum(D$k) < 5 || length(unique(D$year)) < min_years) next
    } else { f <- fire_wk(eco, fuels, o); D$y <- round(ifelse(is.na(f[key]), 0, f[key]))
      if (sum(D$y) < 20 || sum(tapply(D$y, D$year, sum) >= 3) < 5) next }
    D <- D[order(D$year, D$week), ]
    t1 <- lr_test(D, o == "size", xfull)
    RES[[length(RES) + 1]] <- data.frame(family = kf, group = g, eco_l2_code = eco, outcome = o, years = length(unique(D$year)),
      weeks = length(unique(D$week)), fires = if (o == "size") sum(D$n) else sum(D$y), ratio_per_10 = round(exp(t1[["b"]]), 3),
      p = signif(t1[["p"]], 3), p_param = signif(t1[["p_param"]], 3), expected_direction = t1[["b"]] < 0, stringsAsFactors = FALSE)
  }
  if (r %% 20 == 0) log_msg(sprintf("  %d / %d series (%.1f min)", r, nrow(keys), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (!length(RES)) stop("no series could be tested")
R <- do.call(rbind, RES)
R$otype <- ifelse(R$outcome == "size", "size", "counts")
R$q <- NA_real_
for (k in unique(paste(R$family, R$otype))) { i <- paste(R$family, R$otype) == k; R$q[i] <- signif(stats::p.adjust(R$p[i], "BH"), 3) }
utils::write.csv(R, file.path(out, "results.csv"), row.names = FALSE)
S <- do.call(rbind, lapply(split(R, paste(R$family, R$outcome)), function(d) data.frame(family = d$family[1], outcome = d$outcome[1], tests = nrow(d),
  p05 = sum(d$p < 0.05, na.rm = TRUE), p05_fewer = sum(d$p < 0.05 & d$expected_direction %in% TRUE, na.rm = TRUE),
  p05_more = sum(d$p < 0.05 & d$expected_direction %in% FALSE, na.rm = TRUE), chance_each_side = round(0.025 * sum(!is.na(d$p)), 1),
  q_hits = sum(d$q < fdr, na.rm = TRUE), median_ratio = round(stats::median(d$ratio_per_10, na.rm = TRUE), 3))))
utils::write.csv(S, file.path(out, "summary.csv"), row.names = FALSE)
for (i in seq_len(nrow(S))) log_msg(sprintf("  %-7s %-18s %3d tests | p < 0.05: fewer %2d, more %2d (chance %.1f each) | FDR hits %d | median ratio %.3f",
  S$family[i], S$outcome[i], S$tests[i], S$p05_fewer[i], S$p05_more[i], S$chance_each_side[i], S$q_hits[i], S$median_ratio[i]))
h <- R[!is.na(R$q) & R$q < fdr, ]; h <- h[order(h$q), ]
for (i in seq_len(min(30, nrow(h)))) log_msg(sprintf("  hit: %-7s %-16s %-5s %-18s ratio per 10 points %.3f | q %.3g | %d years", h$family[i], h$group[i],
  h$eco_l2_code[i], h$outcome[i], h$ratio_per_10[i], h$q[i], h$years[i]))
log_msg(sprintf("=== done | %d tests | %.1f min | %s | log: %s", nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
