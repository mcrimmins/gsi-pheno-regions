#!/usr/bin/env Rscript
# Step 07b (NPN x fire study; laptop, ~1 h; MASS): how many NPN sites a year does the step 3 signal need?
# Asked by Mike on 3 Nov 2026; rules fixed in config (npn_fire: thinning) before any result was seen.
#
# Cases: the 07 cases that passed, plus every Mixed Wood Plains (8.1) deciduous-tree event with a typical
#   date in the first half of the year (where step 3 found the signal), lead record.
# For each site count k (config sites_per_year; 'all' = the full record) and draw d = 1..draws: in each year
#   that has >= k sites with a departure, k sites are drawn at random (without replacement) and the year's
#   NPN departure is their mean (each site's departure is still from its own full-record typical date).
#   Years kept: those with >= max(k, yeartoyear min_sites) sites, i.e. 07's years for k <= 10 and fewer
#   years for larger k; k below 10 thins 07's years. Then 07's test, unchanged: step at the calendar date (M0)
#   vs at the year's NPN date (M1), leave-one-year-out multinomial log-likelihood, bootstrap over held-out
#   years (B here = config thinning: B), pass = gain > 0, 90 % interval above 0, M1 step in the expected
#   direction. A draw with fewer than min_years years is 'too few years'.
# Outputs <run_dir>/npn_fire/thinning/: draws.csv (case x k x draw), summary.csv (case x k: share of
#   draws with gain > 0, share passing, median gain, years), curve.csv (k: across cases).
# Usage: Rscript npn_fire/scripts/07b_nf_thinning.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("MASS", quietly = TRUE)) stop("needs the MASS package")

cfg <- load_nf_cfg(profile)
log_file <- log_init("07b_nf_thinning", profile)
ic <- cfg$npn_fire$inventory; sc <- cfg$npn_fire$seasonal; yc <- cfg$npn_fire$yeartoyear; tc <- cfg$npn_fire$thinning
W <- as.integer(sc$window_weeks %||% 8); B <- as.integer(tc$B %||% 200); lev <- yc$bootstrap$level %||% 0.9
min_sy <- as.integer(yc$min_site_years %||% 2); max_shift <- as.integer(yc$max_shift_weeks %||% 6)
min_years <- as.integer(yc$min_years %||% 8); draws <- as.integer(tc$draws %||% 20)
min_sites <- as.integer(yc$min_sites %||% 10)
nd_env <- suppressWarnings(as.integer(Sys.getenv("NF_DRAWS", ""))); if (is.finite(nd_env)) draws <- nd_env   # testing only
ks <- tc$sites_per_year %||% list(5, 10, 15, 20, 30, 40, "all")
out <- nf_out(cfg, "thinning")
set.seed(tc$seed %||% 20261103)
t0 <- Sys.time()
log_msg("=== 07b_nf_thinning | profile: ", profile, " | sites per year: ", paste(unlist(ks), collapse = ", "), " | draws ", draws, " | B ", B)

R7 <- utils::read.csv(file.path(cfg$out_dir, "yeartoyear", "results.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
INV <- utils::read.csv(file.path(cfg$out_dir, "inventory", "inventory.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
SD <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "inventory", "site_dates"))); SD$eco_l2_code <- as.character(SD$eco_l2_code)
FW <- nf_read_table(nf_table_path(file.path(cfg$out_dir, "fire", "fire_weekly"))); FW$eco_l2_code <- as.character(FW$eco_l2_code)
PW <- utils::read.csv(file.path(cfg$out_dir, "prism", "eco_weekly.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
oy <- tapply(FW$year, FW$outcome, function(z) sort(unique(z)))
fe <- tc$focus %||% list(eco_l2_code = "8.1", group = "deciduous_tree", max_week = 26)
sel <- R7$pass %in% TRUE | (R7$eco_l2_code == fe$eco_l2_code & R7$group == fe$group & R7$typical_week <= fe$max_week)
cases <- R7[sel, c("eco_l2_code", "group", "event", "outcome", "expected", "typical_week")]
log_msg(sprintf("cases: %d (07 passes %d, plus %s %s events up to week %d)", nrow(cases), sum(R7$pass %in% TRUE), fe$eco_l2_code, fe$group, fe$max_week))

wrap <- function(d) ifelse(d > 182, d - 365, ifelse(d < -182, d + 365, d))
site_deps <- function(eco, g, ev) {
  d <- SD[SD$eco_l2_code == eco & SD$group == g & SD$event == ev & !is.na(SD$doy), c("site_id", "year", "doy")]
  ny <- tapply(d$year, d$site_id, function(v) length(unique(v)))
  d <- d[d$site_id %in% names(ny)[ny >= min_sy], ]
  d$dep <- wrap(d$doy - stats::ave(d$doy, d$site_id, FUN = stats::median))
  d
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
mn_ll <- function(y, eta) { if (sum(y) <= 0) return(0); p <- exp(eta - max(eta)); p <- p / sum(p); sum(y * log(pmax(p, 1e-12))) }
fit_step <- function(D) {
  m <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ wy + dry + after, data = D)), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- stats::coef(m); if (any(is.na(b[c("dry", "after")]))) return(NULL)
  c(dry = unname(b[["dry"]]), after = unname(b[["after"]]))
}
rel <- seq(-W, W - 1)
shift <- function(d) pmax(-max_shift, pmin(max_shift, round(d / 7)))
# LOYO log-likelihoods for steps placed at shifts sh (one per year) over years Y (list of y, dry)
loyo <- function(Y, sh) {
  vapply(seq_along(Y), function(h) {
    tr <- setdiff(seq_along(Y), h)
    D <- do.call(rbind, lapply(tr, function(j) data.frame(wy = j, y = Y[[j]]$y, dry = Y[[j]]$dry, after = as.integer(rel >= sh[j]))))
    D$wy <- factor(D$wy); if (sum(D$y) <= 0) return(NA_real_)
    b <- fit_step(D); if (is.null(b)) return(NA_real_)
    mn_ll(Y[[h]]$y, b[["dry"]] * Y[[h]]$dry + b[["after"]] * as.integer(rel >= sh[h]))
  }, 0)
}

DR <- list(); n_run <- 0
for (r in seq_len(nrow(cases))) {
  e <- cases[r, ]; o <- e$outcome; tw <- e$typical_week
  fuels <- unlist(ic$match[[e$group]]); wt <- sc$dryness_weight[[e$group]]
  S <- site_deps(e$eco_l2_code, e$group, e$event)
  yy_all <- oy[[o]]; S <- S[S$year %in% yy_all, ]
  f <- series(e$eco_l2_code, fuels, o, yy_all); dr <- wk_val(e$eco_l2_code, wt, yy_all, "ppet60") / 100
  nsite <- tapply(S$site_id, S$year, function(v) length(unique(v)))
  for (k in ks) {
    kk <- if (identical(k, "all")) 0L else as.integer(k)
    yrs <- as.integer(names(nsite)[nsite >= max(kk, min_sites)])   # the years 07 used, or fewer when k > min_sites
    Y <- list(); keep <- integer(0)
    for (yv in yrs) { i <- match(yv, yy_all); t <- (i - 1) * 52 + tw + rel
      if (min(t) < 1 || max(t) > length(f) || any(is.na(dr[t]))) next
      Y[[length(Y) + 1]] <- list(y = f[t], dry = dr[t]); keep <- c(keep, yv) }
    nd <- if (kk == 0) 1L else draws
    if (length(Y) < min_years) {
      DR[[length(DR) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, k = as.character(k), draw = NA_integer_,
        years = length(Y), gain = NA_real_, gain_lo = NA_real_, pass = NA, too_few_years = TRUE)
      next
    }
    ll0 <- loyo(Y, rep(0, length(Y)))                      # calendar: the same for every draw
    for (d in seq_len(nd)) {
      dep <- vapply(keep, function(yv) { s <- S[S$year == yv, ]; s <- s[!duplicated(s$site_id), ]
        if (kk > 0) s <- s[sample(nrow(s), kk), ]; mean(s$dep) }, 0)
      sh <- shift(dep); ll1 <- loyo(Y, sh); g <- (ll1 - ll0)[is.finite(ll1) & is.finite(ll0)]
      lo <- if (length(g) >= 3) stats::quantile(vapply(seq_len(B), function(i) mean(sample(g, replace = TRUE)), 0), (1 - lev) / 2) * length(g) else NA_real_
      Df <- do.call(rbind, lapply(seq_along(Y), function(j) data.frame(wy = j, y = Y[[j]]$y, dry = Y[[j]]$dry, after = as.integer(rel >= sh[j]))))
      Df$wy <- factor(Df$wy); bf <- fit_step(Df); sb <- if (is.null(bf)) NA_real_ else bf[["after"]]
      dir_ok <- switch(e$expected, down = sb < 0, up = sb > 0, TRUE)
      DR[[length(DR) + 1]] <- data.frame(eco_l2_code = e$eco_l2_code, group = e$group, event = e$event, k = as.character(k), draw = d,
        years = length(Y), gain = round(sum(g), 2), gain_lo = round(unname(lo), 2),
        pass = length(g) >= 4 && sum(g) > 0 && is.finite(lo) && lo > 0 && isTRUE(dir_ok), too_few_years = FALSE)
      n_run <- n_run + 1
    }
  }
  log_msg(sprintf("  %d / %d cases (%s %s %s) | %.1f min", r, nrow(cases), e$eco_l2_code, e$group, e$event, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
D <- do.call(rbind, DR); utils::write.csv(D, file.path(out, "draws.csv"), row.names = FALSE)
korder <- vapply(ks, as.character, "")
key <- paste(D$eco_l2_code, D$group, D$event, D$k)
S <- do.call(rbind, lapply(split(D, key), function(d) data.frame(eco_l2_code = d$eco_l2_code[1], group = d$group[1], event = d$event[1], k = d$k[1],
  years = stats::median(d$years), too_few_years = all(d$too_few_years), draws = sum(!d$too_few_years),
  share_gain_pos = if (all(d$too_few_years)) NA_real_ else round(mean(d$gain > 0, na.rm = TRUE), 2),
  share_pass = if (all(d$too_few_years)) NA_real_ else round(mean(d$pass, na.rm = TRUE), 2),
  median_gain = if (all(d$too_few_years)) NA_real_ else round(stats::median(d$gain, na.rm = TRUE), 1), stringsAsFactors = FALSE)))
S <- S[order(S$eco_l2_code, S$group, match(S$event, names(ic$events)), match(S$k, korder)), ]
utils::write.csv(S, file.path(out, "summary.csv"), row.names = FALSE)
foc <- S$eco_l2_code == fe$eco_l2_code & S$group == fe$group
C <- do.call(rbind, lapply(korder, function(k) { s <- S[S$k == k, ]; f <- s[foc[S$k == k], ]
  data.frame(k = k, cases = nrow(s), cases_tested = sum(!s$too_few_years), share_gain_pos = round(mean(s$share_gain_pos, na.rm = TRUE), 2),
             share_pass = round(mean(s$share_pass, na.rm = TRUE), 2), focus_cases_tested = sum(!f$too_few_years),
             focus_share_gain_pos = round(mean(f$share_gain_pos, na.rm = TRUE), 2), focus_share_pass = round(mean(f$share_pass, na.rm = TRUE), 2)) }))
utils::write.csv(C, file.path(out, "curve.csv"), row.names = FALSE)
log_msg(sprintf("draws run: %d", n_run))
for (i in seq_len(nrow(C))) log_msg(sprintf("  k = %-4s | cases tested %2d / %2d | gain > 0 %3.0f %% | pass %3.0f %% || %s %s: tested %2d | gain > 0 %3.0f %% | pass %3.0f %%",
  C$k[i], C$cases_tested[i], C$cases[i], 100 * C$share_gain_pos[i], 100 * C$share_pass[i], fe$eco_l2_code, fe$group,
  C$focus_cases_tested[i], 100 * C$focus_share_gain_pos[i], 100 * C$focus_share_pass[i]))
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out, log_file))
