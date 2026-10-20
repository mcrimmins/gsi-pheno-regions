#!/usr/bin/env Rscript
# Step 66 (fire activity study; laptop or P720, minutes): the tests (scope decisions 13, 21).
# Input: fire_study/weekly/fire_weekly.csv (65).
#
# Year-to-year test: one pooled negative binomial model of weekly counts with region effects, scored
# leaving one year out (fit on the other years, predict the held-out year's 52 weeks per region,
# score = negative binomial log predictive density). Models (every term region-specific):
#   M0 season       region + region x 3 annual harmonics of week (the typical fire curve)
#   M1 simple       M0 + 60-day rain minus PET (ppet60, per 100 mm)          [decision 21: fitted first]
#   M2 simple+pheno M1 + cured-share departure (that year's cured share - typical, per 0.1)
#   M1d / M2d       the same with the 1000-h fuel moisture (dfm1000, per 1 %) instead of ppet60
# Snow regions add their weekly snow share to every model. Skill of phenology = held-out log score of
# M2 minus M1 (and M2d - M1d), summed over years, per region and pooled; uncertainty by bootstrap
# over years (resampling held-out years). Also the correlation of held-out yearly totals with
# observed, per region (what managers would see), and the full-data phenology coefficients.
# Outcome specs (the primary and the sensitivities in the scope doc):
#   all          n_fires, all regions, 2001-2024 (primary)
#   era          all + region x (year >= 2020) in every model
#   lightning / human   cause split, 2001-2019, regions without cause_split_skip_regions
#   human_unk    human + unknown, all years, same regions
#   debris       C+ fires including debris / open burning
#   top3         region 13 fires in the most herbaceous third of its cells
#   px / px_nofill   VIIRS pixel-days 2012-2024 (no region 13); px_nofill drops weeks with filled days
#   px10_all     region 10, all non-cropland detections (linked or not)
#   all_nb, human_unk_nb, px_nb   burn check: phenology from cells not burned that year or the year before
#   all_fc, human_unk_fc, px_fc   forecast test: cured share from curing dates predicted 14 days ahead
#                (67; skipped when forecast_weekly.csv is missing)
#   all_lag2, human_unk_lag2, px_lag2   nowcast: cured-share departure lagged 2 weeks (satellite latency)
#   all_rt, human_unk_rt, px_rt   real-time cured share (65b): only composites available that week
#   all_rts, human_unk_rts, px_rts   the same from the smoothed nowcast (65b *_rts)
#   all_hs, human_unk_hs, px_hs   hot-spell check: 7-day VPD (64 vpd7) in every model
#   all_gsi_<v>, human_unk_gsi_<v>, px_gsi_<v>   cured share from GSI dates (config gsi_compare)
#   all_dates, human_unk_dates, px_dates   dates version: green-up / curing departures (per 10 days) in
#                seasonal windows around the typical green-up (-4..+3 weeks) and curing (0..+7 weeks)
# Year-shuffle check (all, px): each region's cured-share departures assigned to random other years,
# GSI_FIRE_SHUFFLES times (default 20) -> shuffle.csv; the real gain should beat the shuffled ones.
#   (pixel-day specs: Gaussian model of log(1 + pixel-days); counts too heavy-tailed for the NB)
# Seasonal test: per region and week, the mean fire count / pixel-days, typical cured share, mean
# dryness -> seasonal.csv; and the weeks where the typical fire curve and the typical cured share
# cross their halfway points (window open / close) -> seasonal_windows.csv.
# Outputs <run_dir>/fire_study/test/: loyo_scores.csv (spec x model x region x year), skill.csv
# (spec x comparison x region: total gain, bootstrap 90 % interval, share of years improved),
# yearly_r.csv, coefs.csv, seasonal.csv, seasonal_windows.csv, fit_notes.csv.
#
# Usage: Rscript scripts/66_fire_test.R [profile] [spec,...]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "notify.R"))
source(here::here("R", "fire.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
suppressPackageStartupMessages(library(MASS))

cfg <- load_cfg(profile)
log_file <- log_init("66_fire_test", profile)
notify_init(cfg)
fs <- cfg$fire_study
t0 <- Sys.time()
out_dir <- fire_dir(cfg, "test"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
W <- utils::read.csv(fire_dir(cfg, "weekly", "fire_weekly.csv"))
regs <- as.integer(unlist(fs$pilot_regions))
snow_regs <- as.integer(unlist(fs$snow_regions))
skip_cause <- as.integer(unlist(fs$outcome$cause_split_skip_regions))
set.seed(66)
n_boot <- 2000

# predictors
W$cured_dep <- (W$cured_share - W$cured_share_typ) * 10            # per 0.1 of the fuel
W$ppet <- W$ppet60 / 100                                           # per 100 mm
for (k in 1:3) { W[[paste0("s", k)]] <- sin(2 * pi * k * W$week / 52); W[[paste0("c", k)]] <- cos(2 * pi * k * W$week / 52) }
for (r in snow_regs) W[[paste0("snow", r)]] <- ifelse(W$region == r, W$snow_share, 0)
# burn check (Oct 19): cured-share departure from cells not burned that year or the year before
W$cured_dep_nb <- (W$cured_share_nb - W$cured_share_typ_nb) * 10
# dates version (Oct 19): departures in days (per 10 d) acting in seasonal windows, each with its own
# coefficient per region: green-up departure in the 8 weeks around the typical green-up week (-4..+3),
# curing departure in the 8 weeks from the typical curing week (0..+7); weeks wrap around the year.
in_win <- function(w, w0, lo, hi) { d <- (w - w0) %% 52; d[d > 26] <- d[d > 26] - 52; !is.na(w0) & d >= lo & d <= hi }
gd <- W$greenup_dep; gd[is.na(gd)] <- 0; cd <- W$curing_dep; cd[is.na(cd)] <- 0
W$gdep_win <- ifelse(in_win(W$week, W$gu_week_typ, -4, 3), gd / 10, 0)
W$cdep_win <- ifelse(in_win(W$week, W$cu_week_typ, 0, 7), cd / 10, 0)
# nowcast (Oct 19): the cured-share departure as the satellite showed it 2 weeks earlier (weeks 1-2 take
# the previous year's weeks 51-52)
W <- W[order(W$region, W$year, W$week), ]
lag_key <- paste(W$region, ifelse(W$week <= 2, W$year - 1, W$year), ifelse(W$week <= 2, W$week + 50, W$week - 2))
W$cured_dep_lag2 <- W$cured_dep[match(lag_key, paste(W$region, W$year, W$week))]

specs <- list(
  all = list(y = "n_fires", rows = TRUE),
  era = list(y = "n_fires", rows = TRUE, era = TRUE),
  lightning = list(y = "n_lightning", rows = quote(cause_split_ok)),
  human = list(y = "n_human", rows = quote(cause_split_ok)),
  human_unk = list(y = "n_human_unk", rows = quote(!(region %in% skip_cause))),
  debris = list(y = "n_fires_debris", rows = TRUE),
  top3 = list(y = "n_fires_top3", rows = quote(region %in% as.integer(unlist(fs$cover_sensitivity$regions)))),
  px = list(y = "px_days", rows = quote(viirs_ok), log = TRUE),
  px_nofill = list(y = "px_days_nofill", rows = quote(viirs_ok & fill_days %in% 0), log = TRUE),
  all_nb = list(y = "n_fires", rows = TRUE, pheno = "cured_dep_nb"),
  human_unk_nb = list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = "cured_dep_nb"),
  px_nb = list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = "cured_dep_nb"),
  all_lag2 = list(y = "n_fires", rows = TRUE, pheno = "cured_dep_lag2"),
  human_unk_lag2 = list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = "cured_dep_lag2"),
  px_lag2 = list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = "cured_dep_lag2"),
  all_dates = list(y = "n_fires", rows = TRUE, pheno = c("gdep_win", "cdep_win")),
  human_unk_dates = list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = c("gdep_win", "cdep_win")),
  px_dates = list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = c("gdep_win", "cdep_win")),
  all_fc = list(y = "n_fires", rows = TRUE, pheno = "cured_dep_fc"),
  human_unk_fc = list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = "cured_dep_fc"),
  px_fc = list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = "cured_dep_fc"),
  px10_all = list(y = "px_noncrop_all", log = TRUE, rows = quote(region %in% as.integer(unlist(fs$viirs_unlinked_sensitivity_regions)) &
                                                      year >= fs$viirs_outcome$years$start & year <= fs$viirs_outcome$years$end))
)
W$n_human_unk <- W$n_human + W$n_unknown
# forecast test (67, Oct 19): cured share from predicted curing dates (lead-time weather adjustment refit
# leaving each year out), when 67 has run
ff <- fire_dir(cfg, "forecast", "forecast_weekly.csv")
if (file.exists(ff)) {
  FC <- utils::read.csv(ff)
  i <- match(paste(W$region, W$year, W$week), paste(FC$region, FC$year, FC$week))
  W$cured_dep_fc <- (FC$cured_share_fc[i] - FC$cured_share_typ_fc[i]) * 10
} else specs <- specs[!grepl("_fc$", names(specs))]
# real-time cured share (65b, Oct 19): from the composites available each week, typical = other years
rf <- fire_dir(cfg, "realtime", "realtime_weekly.csv")
if (file.exists(rf)) {
  RT <- utils::read.csv(rf)
  i <- match(paste(W$region, W$year, W$week), paste(RT$region, RT$year, RT$week))
  W$cured_dep_rt <- (RT$cured_share_rt[i] - RT$cured_share_rt_typ[i]) * 10
  specs$all_rt <- list(y = "n_fires", rows = TRUE, pheno = "cured_dep_rt")
  specs$human_unk_rt <- list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = "cured_dep_rt")
  specs$px_rt <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = "cured_dep_rt")
  if ("cured_share_rts" %in% names(RT)) {
    W$cured_dep_rts <- (RT$cured_share_rts[i] - RT$cured_share_rts_typ[i]) * 10
    specs$all_rts <- list(y = "n_fires", rows = TRUE, pheno = "cured_dep_rts")
    specs$human_unk_rts <- list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = "cured_dep_rts")
    specs$px_rts <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = "cured_dep_rts")
  }
}
for (v in c("viirs_ok", "cause_split_ok", "era2020")) W[[v]] <- as.logical(W[[v]])     # robust to "True" / "TRUE"
# hot-spell check (Oct 19): 7-day VPD added to every model, so phenology must add beyond this week's dryness
if ("vpd7" %in% names(W)) {
  specs$all_hs <- list(y = "n_fires", rows = TRUE, extra = "vpd7")
  specs$human_unk_hs <- list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), extra = "vpd7")
  specs$px_hs <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, extra = "vpd7")
}
# GSI comparison (Oct 19): cured share from GSI season dates (65 columns cured_share_gsi_<variant>)
for (v in unlist(fs$gsi_compare)) {
  cs <- paste0("cured_share_gsi_", v); ct <- paste0("cured_share_typ_gsi_", v)
  if (!all(c(cs, ct) %in% names(W))) next
  pv <- paste0("cured_dep_gsi_", v); W[[pv]] <- (W[[cs]] - W[[ct]]) * 10
  specs[[paste0("all_gsi_", v)]] <- list(y = "n_fires", rows = TRUE, pheno = pv)
  specs[[paste0("human_unk_gsi_", v)]] <- list(y = "n_human_unk", rows = quote(!(region %in% skip_cause)), pheno = pv)
  specs[[paste0("px_gsi_", v)]] <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = pv)
  # does the GSI's gain come from its VPD ramp? GSI with this week's VPD (vpd7) in every model; and
  # satellite + GSI together (the gain over each alone comes from loyo_scores: M2 of *_both vs M2 of each)
  if ("vpd7" %in% names(W)) {
    specs[[paste0("all_gsi_", v, "_hs")]] <- list(y = "n_fires", rows = TRUE, pheno = pv, extra = "vpd7")
    specs[[paste0("px_gsi_", v, "_hs")]] <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = pv, extra = "vpd7")
  }
  specs[[paste0("all_both_", v)]] <- list(y = "n_fires", rows = TRUE, pheno = c("cured_dep", pv))
  specs[[paste0("px_both_", v)]] <- list(y = "px_days", rows = quote(viirs_ok), log = TRUE, pheno = c("cured_dep", pv))
}
pick <- args[-1]
if (length(pick)) specs <- specs[intersect(names(specs), unlist(strsplit(pick, ",")))]

terms_for <- function(d, model, era, pheno = "cured_dep", extra = NULL) {
  multi <- length(unique(d$region)) > 1
  rz <- function(x) if (multi) paste0("region_f:(", x, ")") else x
  h <- "s1 + c1 + s2 + c2 + s3 + c3"
  f <- c(if (multi) "region_f", rz(h))
  sn <- intersect(paste0("snow", snow_regs), paste0("snow", unique(d$region)))
  sn <- sn[vapply(sn, function(v) stats::sd(d[[v]], na.rm = TRUE) > 0, TRUE)]
  f <- c(f, sn)
  if (era) f <- c(f, rz("era2020"))
  if (length(extra)) f <- c(f, rz(paste(extra, collapse = " + ")))
  dry <- switch(model, M1 = , M2 = "ppet", M1d = , M2d = "dfm1000", NULL)
  if (!is.null(dry)) f <- c(f, rz(dry))
  if (model %in% c("M2", "M2d")) f <- c(f, rz(paste(pheno, collapse = " + ")))
  stats::as.formula(paste("y ~", paste(f, collapse = " + ")))
}
fit_nb <- function(fm, d) {
  m <- tryCatch(suppressWarnings(MASS::glm.nb(fm, data = d, control = glm.control(maxit = 100))), error = function(e) NULL)
  if (is.null(m)) {                  # heavy-tailed counts (pixel-days): start from the Poisson fit, theta near 1
    p0 <- stats::glm(fm, data = d, family = poisson())
    m <- tryCatch(suppressWarnings(MASS::glm.nb(fm, data = d, start = stats::coef(p0), init.theta = 1,
                                                control = glm.control(maxit = 200))), error = function(e) NULL)
  }
  if (is.null(m)) { m <- stats::glm(fm, data = d, family = poisson()); m$theta <- Inf; m$note <- "poisson fallback" }
  m
}
# Pixel-days are too heavy-tailed for a count model (negative binomial fits fail in some folds):
# those specs model log(1 + pixel-days) with a Gaussian linear model (score = Gaussian log density).
fit_lm <- function(fm, d) { m <- stats::lm(fm, data = d); m$sigma <- summary(m)$sigma; m }
ll_nb <- function(y, mu, theta) if (is.finite(theta)) stats::dnbinom(y, size = theta, mu = mu, log = TRUE) else stats::dpois(y, mu, log = TRUE)
# Leave-one-year-out: per held-out year, fit on the others and score the year's weeks.
loyo <- function(d, fm, fitf, logy) {
  notes <- c(); out <- list()
  for (yo in sort(unique(d$year))) {
    tr <- d[d$year != yo, ]; te <- d[d$year == yo, ]
    m <- fitf(fm, tr)
    if (!is.null(m$note)) notes <- c(notes, yo)
    te$region_f <- factor(te$region, levels = levels(tr$region_f))
    mu <- stats::predict(m, newdata = te, type = "response")
    ll <- if (logy) stats::dnorm(te$y, mu, m$sigma, log = TRUE) else ll_nb(te$y, mu, m$theta)
    y <- te$y
    if (logy) { y <- expm1(y); mu <- expm1(mu) }                   # yearly totals on the count scale
    out[[length(out) + 1]] <- data.frame(region = te$region, year = yo, ll = ll, y = y, mu = mu, week = te$week)
  }
  list(scores = do.call(rbind, out), notes = notes)
}
models <- c("M0", "M1", "M2", "M1d", "M2d")
comps <- list(pheno_simple = c("M2", "M1"), pheno_dfm = c("M2d", "M1d"), dry_simple = c("M1", "M0"), dry_dfm = c("M1d", "M0"),
              dfm_vs_simple = c("M1d", "M1"))

SC <- list(); CO <- list(); YR <- list(); NOTE <- list()
for (sn in names(specs)) {
  sp <- specs[[sn]]
  d <- W[eval(if (is.logical(sp$rows)) quote(rep(TRUE, nrow(W))) else sp$rows, W, environment()), ]
  d$y <- if (isTRUE(sp$log)) log1p(d[[sp$y]]) else d[[sp$y]]
  ph <- sp$pheno %||% "cured_dep"
  need <- c("y", "ppet", "dfm1000", ph, sp$extra)
  d <- d[stats::complete.cases(d[, need]), ]
  d$region_f <- factor(d$region)
  yrs <- sort(unique(d$year))
  log_msg(sprintf("spec %-10s %s: %d rows, regions %s, years %d-%d, total %s", sn, sp$y, nrow(d), paste(sort(unique(d$region)), collapse = ","),
                  min(yrs), max(yrs), format(sum(d[[sp$y]]), big.mark = ",")))
  for (mo in models) {
    fm <- terms_for(d, mo, isTRUE(sp$era), ph, sp$extra)
    fitf <- if (isTRUE(sp$log)) fit_lm else fit_nb
    full <- fitf(fm, d)
    if (!is.null(full$note)) NOTE[[length(NOTE) + 1]] <- data.frame(spec = sn, model = mo, fold = "full", note = full$note)
    if (mo %in% c("M2", "M2d", "M1", "M1d")) {
      cf <- summary(full)$coefficients
      keep <- grepl("cured_dep|dep_win|ppet|dfm1000", rownames(cf))
      if (any(keep)) CO[[length(CO) + 1]] <- data.frame(spec = sn, model = mo, term = rownames(cf)[keep], est = cf[keep, 1], se = cf[keep, 2],
                                                         theta = full$theta %||% NA_real_, row.names = NULL)
    }
    lo <- loyo(d, fm, fitf, isTRUE(sp$log))
    for (nt in lo$notes) NOTE[[length(NOTE) + 1]] <- data.frame(spec = sn, model = mo, fold = nt, note = "poisson fallback")
    SC[[length(SC) + 1]] <- data.frame(spec = sn, model = mo, lo$scores)
  }
  log_msg(sprintf("  %s done, %.1f min", sn, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
S <- do.call(rbind, SC)
LY <- stats::aggregate(cbind(ll, y, mu) ~ spec + model + region + year, S, sum)
utils::write.csv(LY, file.path(out_dir, "loyo_scores.csv"), row.names = FALSE)

# skill: held-out log-score gain, per region and pooled, bootstrap over years
boot_gain <- function(g) {           # g = per-year gains
  b <- replicate(n_boot, sum(sample(g, replace = TRUE)))
  c(gain = sum(g), lo = unname(stats::quantile(b, 0.05)), hi = unname(stats::quantile(b, 0.95)), years_better = mean(g > 0), n_years = length(g))
}
SK <- list()
for (sn in unique(LY$spec)) for (cn in names(comps)) {
  a <- LY[LY$spec == sn & LY$model == comps[[cn]][1], ]; b <- LY[LY$spec == sn & LY$model == comps[[cn]][2], ]
  m <- merge(a, b, by = c("region", "year"), suffixes = c("_a", "_b"))
  m$g <- m$ll_a - m$ll_b
  for (r in c(sort(unique(m$region)), NA)) {
    x <- if (is.na(r)) stats::aggregate(g ~ year, m, sum) else m[m$region == r, c("year", "g")]
    SK[[length(SK) + 1]] <- data.frame(spec = sn, comparison = cn, region = if (is.na(r)) "pooled" else as.character(r), t(boot_gain(x$g)))
  }
}
SK <- do.call(rbind, SK)
utils::write.csv(SK, file.path(out_dir, "skill.csv"), row.names = FALSE)

# yearly totals: held-out predicted vs observed (per region)
YR <- do.call(rbind, lapply(split(LY, list(LY$spec, LY$model, LY$region), drop = TRUE), function(x)
  data.frame(spec = x$spec[1], model = x$model[1], region = x$region[1], r_yearly = suppressWarnings(stats::cor(log1p(x$y), log1p(x$mu))))))
utils::write.csv(YR, file.path(out_dir, "yearly_r.csv"), row.names = FALSE)
CO <- do.call(rbind, CO)
CO$z <- CO$est / CO$se
utils::write.csv(CO, file.path(out_dir, "coefs.csv"), row.names = FALSE)
utils::write.csv(if (length(NOTE)) do.call(rbind, NOTE) else data.frame(note = "none"), file.path(out_dir, "fit_notes.csv"), row.names = FALSE)

# ---- year-shuffle check (Oct 19) ------------------------------------------------------------------------------
# Each region's cured-share departures are given to other years at random (weeks kept); if phenology's gain
# comes from year-specific timing, the shuffled gain should fall to about zero. n_shuffle permutations,
# M2 vs M1 (60-day water balance baseline), specs all and px.
n_shuffle <- as.integer(Sys.getenv("GSI_FIRE_SHUFFLES", "20"))
SH <- list()
for (sn in intersect(c("all", "px"), names(specs))) {
  if (n_shuffle < 1) break
  sp <- specs[[sn]]
  d <- W[eval(if (is.logical(sp$rows)) quote(rep(TRUE, nrow(W))) else sp$rows, W, environment()), ]
  d$y <- if (isTRUE(sp$log)) log1p(d[[sp$y]]) else d[[sp$y]]
  d <- d[stats::complete.cases(d[, c("y", "ppet", "dfm1000", "cured_dep")]), ]
  d$region_f <- factor(d$region)
  fitf <- if (isTRUE(sp$log)) fit_lm else fit_nb
  base <- stats::aggregate(ll ~ region + year, LY[LY$spec == sn & LY$model == "M1", ], sum)
  obs <- sum(LY$ll[LY$spec == sn & LY$model == "M2"]) - sum(base$ll)
  for (k in seq_len(n_shuffle)) {
    d$cured_dep_shuf <- NA_real_
    for (r in unique(d$region)) {
      ys <- sort(unique(d$year[d$region == r])); src <- stats::setNames(sample(ys), ys)
      i <- d$region == r
      key_from <- paste(r, src[as.character(d$year[i])], d$week[i])
      d$cured_dep_shuf[i] <- d$cured_dep[match(key_from, paste(d$region, d$year, d$week))]
    }
    fm <- terms_for(d, "M2", FALSE, "cured_dep_shuf")
    lo <- loyo(d, fm, fitf, isTRUE(sp$log))$scores
    g <- stats::aggregate(ll ~ region, lo, sum)
    b <- stats::aggregate(ll ~ region, base, sum)
    SH[[length(SH) + 1]] <- data.frame(spec = sn, shuffle = k, region = c(g$region, "pooled"),
                                       gain = c(g$ll - b$ll[match(g$region, b$region)], sum(lo$ll) - sum(base$ll)), observed_pooled = obs)
  }
  x <- do.call(rbind, SH); x <- x[x$spec == sn & x$region == "pooled", ]
  log_msg(sprintf("shuffle %-4s observed pooled gain %+.1f; shuffled median %+.1f, max %+.1f (%d of %d shuffles >= observed)",
                  sn, obs, stats::median(x$gain), max(x$gain), sum(x$gain >= obs), nrow(x)))
}
if (length(SH)) utils::write.csv(do.call(rbind, SH), file.path(out_dir, "shuffle.csv"), row.names = FALSE)

# ---- seasonal test ------------------------------------------------------------------------------------------
SE <- stats::aggregate(W[, c("n_fires", "px_days", "cured_share_typ", "cured_share", "dfm1000", "ppet60", "snow_share")],
                       W[, c("region", "week")], mean, na.rm = TRUE, na.action = stats::na.pass)
utils::write.csv(SE, file.path(out_dir, "seasonal.csv"), row.names = FALSE)
half_cross <- function(x, up = TRUE) {    # first week the series crosses half its range going up / down
  h <- min(x) + (max(x) - min(x)) / 2
  w <- if (up) which(x[-length(x)] < h & x[-1] >= h) + 1 else which(x[-length(x)] >= h & x[-1] < h) + 1
  if (length(w)) paste(w, collapse = ";") else NA_character_
}
SW <- do.call(rbind, lapply(regs, function(r) {
  s <- SE[SE$region == r, ]; s <- s[order(s$week), ]
  fz <- stats::filter(s$n_fires, rep(1 / 3, 3), circular = TRUE)
  data.frame(region = r, fire_peak_week = s$week[which.max(fz)], fire_up = half_cross(as.numeric(fz), TRUE), fire_down = half_cross(as.numeric(fz), FALSE),
             cured_up = half_cross(s$cured_share_typ, TRUE), cured_down = half_cross(s$cured_share_typ, FALSE),
             r_fire_cured = stats::cor(s$n_fires, s$cured_share_typ))
}))
utils::write.csv(SW, file.path(out_dir, "seasonal_windows.csv"), row.names = FALSE)

# ---- log -----------------------------------------------------------------------------------------------------------
for (sn in unique(SK$spec)) for (cn in c("dry_simple", "pheno_simple", "dfm_vs_simple", "pheno_dfm")) {
  x <- SK[SK$spec == sn & SK$comparison == cn, ]
  log_msg(sprintf("%-10s %-13s %s", sn, cn, paste(sprintf("%s %+.1f [%+.1f, %+.1f] %.0f%%", x$region, x$gain, x$lo, x$hi, 100 * x$years_better), collapse = " | ")))
}
msg <- sprintf("%d specs, %.1f min", length(specs), as.numeric(difftime(Sys.time(), t0, units = "mins")))
log_msg("=== done | ", msg, " | log: ", log_file)
notify("66_fire_test finished", msg, priority = 3, tags = "white_check_mark")
