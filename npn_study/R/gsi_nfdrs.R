## =============================================================================
## gsi_nfdrs.R (NPN grass curing study) -- the project's existing GSI code, vendored so the
## study stays self-contained (it becomes its own repo at publication).
##   1. GrowingSeasonIndex/R/gsi.R, copied unchanged below (19,631 bytes, 2026-09).
##   2. gsi_dates() from GrowingSeasonIndex/study/R/study_utils.R (unchanged).
##   3. npn_gsi_params_v4(): the NFDRS v4 parameter set exactly as in
##      GrowingSeasonIndex/study/05_build_tables.R (Jolly et al. 2024; daily maximum VPD).
## If the originals change, re-copy them here; nothing else in npn_study edits these functions.
## =============================================================================

## -----------------------------------------------------------------------------
## gsi.R -- Growing Season Index model, as pure functions.
##
## Baseline follows Jolly, Nemani & Running (2005) Glob Change Biol 11:619-632
## and the NFDRS2016 implementation (NWCG). Optional extensions follow
## Daham et al. (2018) J Water Clim Change (precipitation control) plus
## structural variants that are the point of this tool.
##
## Nothing here knows about Shiny. Everything is vectorised over a day series.
## -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- generic linear ramp -----------------------------------------------------

#' Linear ramp bounded to [0, 1]
#' @param x numeric vector
#' @param lo value at which the index is 0 (or 1 if decreasing)
#' @param hi value at which the index is 1 (or 0 if decreasing)
#' @param decreasing if TRUE the index falls from 1 at `lo` to 0 at `hi`
gsi_ramp <- function(x, lo, hi, decreasing = FALSE) {
  if (isTRUE(all.equal(lo, hi))) {
    y <- as.numeric(x >= hi)
  } else {
    y <- pmin(1, pmax(0, (x - lo) / (hi - lo)))
  }
  if (decreasing) 1 - y else y
}

# ---- individual indicator functions -----------------------------------------

#' Minimum temperature indicator. NFDRS2016 defaults: -2 C to 5 C.
idx_tmin <- function(tmin_c, lo = -2, hi = 5) {
  gsi_ramp(tmin_c, lo, hi)
}

#' Vapour pressure deficit indicator. NFDRS2016 defaults: 900 Pa to 4100 Pa.
#' Decreasing: unconstrained below `lo`, fully limiting above `hi`.
idx_vpd <- function(vpd_pa, lo = 900, hi = 4100) {
  gsi_ramp(vpd_pa, lo, hi, decreasing = TRUE)
}

#' Photoperiod indicator. NFDRS2016 defaults: 36000 s (10 h) to 39600 s (11 h).
idx_photo <- function(photo_s, lo = 36000, hi = 39600) {
  gsi_ramp(photo_s, lo, hi)
}

#' Precipitation indicator (Daham et al. 2018), generalised.
#' Daham used a running total with PrecipMax = 10 mm; the accumulation window is
#' left open here because a 10 mm / short-window rule is far too short-memory for
#' monsoon-driven systems.
idx_precip <- function(precip_mm, window = 21, lo = 0, hi = 25) {
  acc <- roll_sum_trailing(precip_mm, window)
  gsi_ramp(acc, lo, hi)
}

# ---- rolling helpers ---------------------------------------------------------

roll_sum_trailing <- function(x, n) {
  x[is.na(x)] <- 0
  cs <- c(0, cumsum(x))
  idx <- seq_along(x)
  lo <- pmax(0L, idx - n)
  cs[idx + 1L] - cs[lo + 1L]
}

roll_mean_trailing <- function(x, n) {
  x0 <- x; x0[is.na(x0)] <- 0
  cs <- c(0, cumsum(x0))
  idx <- seq_along(x)
  lo <- pmax(0L, idx - n)
  (cs[idx + 1L] - cs[lo + 1L]) / (idx - lo)
}

roll_mean_centered <- function(x, n) {
  half <- floor(n / 2)
  n_ <- length(x)
  vapply(seq_len(n_), function(i) {
    mean(x[max(1L, i - half):min(n_, i + half)], na.rm = TRUE)
  }, numeric(1))
}

# ---- photoperiod -------------------------------------------------------------

#' Sun-elevation offsets defining the end of "day", in degrees below horizon.
PHOTO_CONVENTIONS <- c(
  geometric      =  0.000,  # sun centre on the horizon; FAO-56 / Allen et al.
  sunrise_sunset = -0.833,  # sun's upper limb + refraction: published daylength
  civil_twilight = -6.000   # usable light; arguably the biological cue
)

#' Daylength in seconds from day-of-year and latitude.
#'
#' `method = "geometric"` reproduces the FAO-56 relation used by Daham et al.
#' (2018) and equivalent to the Monteith & Unsworth formulation cited by Jolly
#' et al. (2005). It is the default so the app reproduces the published index.
#'
#' The convention is NOT a cosmetic choice. At Tucson (32.2 N) the December
#' solstice daylength is 9.89 h geometric but 10.04 h by the ordinary
#' sunrise-to-sunset definition. With the NFDRS2016 photoperiod thresholds of
#' 10-11 h that single decision determines whether iPhoto reaches zero -- and
#' therefore whether GSI is forced to zero -- across the desert Southwest every
#' December. Under civil twilight it never approaches zero anywhere south of
#' about 40 N. Same model, same station, three different answers.
photoperiod_seconds <- function(doy, lat_deg,
                                method = c("geometric", "sunrise_sunset",
                                           "civil_twilight")) {
  method <- match.arg(method)
  h0    <- PHOTO_CONVENTIONS[[method]] * pi / 180
  phi   <- lat_deg * pi / 180
  delta <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  arg   <- (sin(h0) - sin(phi) * sin(delta)) / (cos(phi) * cos(delta))
  ws    <- acos(pmin(1, pmax(-1, arg)))
  (24 / pi) * ws * 3600
}

# ---- vapour pressure deficit -------------------------------------------------

#' Saturation vapour pressure (kPa) from temperature (C), Tetens/FAO-56.
svp_kpa <- function(t_c) 0.6108 * exp(17.27 * t_c / (t_c + 237.3))

#' Daily VPD (Pa) from Tmax, Tmin and dewpoint, FAO-56 convention:
#' es = mean of saturation vp at Tmax and Tmin; ea = saturation vp at Tdew.
vpd_from_temps <- function(tmax_c, tmin_c, tdew_c) {
  es <- (svp_kpa(tmax_c) + svp_kpa(tmin_c)) / 2
  ea <- svp_kpa(tdew_c)
  pmax(0, es - ea) * 1000
}

# ---- combination -------------------------------------------------------------

#' Combine sub-indices into the daily index.
#'
#' "product" is the Jolly/NFDRS2016 form. The alternatives exist because the
#' product form is a strict bottleneck: any single index near zero zeroes the
#' whole index, which is why GSI collapses through the pre-monsoon Southwest
#' even where perennial vegetation is still functioning.
#'
#' @param idx named list of indicator vectors
#' @param method one of product, min, geometric, weighted
#' @param weights numeric vector, recycled over `idx`, used by "weighted"
gsi_combine <- function(idx, method = c("product", "min", "geometric", "weighted"),
                        weights = NULL) {
  method <- match.arg(method)
  k <- length(idx)
  switch(method,
    product   = Reduce(`*`, idx),
    min       = do.call(pmin, idx),
    geometric = Reduce(`*`, idx)^(1 / k),
    weighted  = {
      w <- if (is.null(weights)) rep(1, k) else rep_len(weights, k)
      Reduce(`+`, Map(function(v, wi) v * wi, idx, w)) / sum(w)
    }
  )
}

# ---- smoothing ---------------------------------------------------------------

#' Smooth the daily index.
#'
#' NFDRS2016 uses a 21-day running average. "ewma_asym" allows different rise
#' and fall time constants, because green-up in an annual grassland is far
#' faster than senescence and a single symmetric window cannot do both.
gsi_smooth <- function(x, method = c("box_trailing", "box_centered", "ewma_asym"),
                       window = 21, halflife_rise = 7, halflife_fall = 21) {
  method <- match.arg(method)
  switch(method,
    box_trailing = roll_mean_trailing(x, window),
    box_centered = roll_mean_centered(x, window),
    ewma_asym = {
      a_r <- 1 - 2^(-1 / max(halflife_rise, 1e-6))
      a_f <- 1 - 2^(-1 / max(halflife_fall, 1e-6))
      s <- numeric(length(x))
      s[1] <- x[1]
      for (t in seq_along(x)[-1]) {
        a <- if (x[t] > s[t - 1]) a_r else a_f
        s[t] <- s[t - 1] + a * (x[t] - s[t - 1])
      }
      s
    }
  )
}

#' One-way curing ratchet: after the seasonal maximum, disallow any increase.
#' Represents the irreversible senescence of annual herbaceous fuels -- a warm,
#' humid spell in October does not re-green cured cheatgrass, but the standard
#' index happily says it does.
apply_curing_ratchet <- function(gsi) {
  i_peak <- which.max(gsi)
  if (length(i_peak) == 0 || i_peak >= length(gsi)) return(gsi)
  out <- gsi
  tail_idx <- i_peak:length(gsi)
  out[tail_idx] <- cummin(gsi[tail_idx])
  out
}

# ---- phenophase state machine ------------------------------------------------

PHASES <- c("dormant", "greenup", "maintenance", "senescence", "cured")

#' Classify each day into a phenophase.
#'
#' NFDRS2016 has no explicit maintenance phase -- it is implied by GSI simply
#' staying high. Making the transitions explicit and tunable is the point: the
#' `maint` threshold and `persist` requirement become objects a local expert can
#' argue with.
#'
#' `seasons = "single"` treats "cured" as absorbing: one green-up per year, which
#' is what an operational single-season assumption amounts to. That assumption is
#' wrong across the Southwest, where cool-season ephemerals and warm-season C4
#' grasses give TWO distinct pulses in the same calendar year. Use
#' `seasons = "multiple"` to let the index green up again after curing, and
#' compare the two -- the difference is one of the clearest statements of what
#' the operational model cannot represent.
#'
#' @param gsi_rel GSI / GSImax, in [0, 1]
#' @param greenup relative GSI needed to start green-up (NFDRS default 0.5)
#' @param maint absolute relative-GSI level at or above which the canopy is
#'   "maintaining". Resolve any fraction-of-annual-max convention before calling.
#' @param persist consecutive days above `greenup` required to trigger green-up
#' @param seasons "single" (cured is absorbing) or "multiple" (re-green allowed)
gsi_phase <- function(gsi_rel, greenup = 0.5, maint = 0.8, persist = 3,
                      seasons = c("single", "multiple")) {
  seasons <- match.arg(seasons)
  n <- length(gsi_rel)
  st <- character(n)
  s <- "dormant"
  run <- 0L
  for (t in seq_len(n)) {
    g <- gsi_rel[t]
    if (s == "dormant") {
      run <- if (g >= greenup) run + 1L else 0L
      if (run >= persist) s <- "greenup"
    } else if (s == "greenup") {
      if (g >= maint) s <- "maintenance"
      else if (g < greenup) { s <- "dormant"; run <- 0L }
    } else if (s == "maintenance") {
      if (g < maint) s <- "senescence"
    } else if (s == "senescence") {
      if (g >= maint) s <- "maintenance"
      else if (g < greenup) { s <- "cured"; run <- 0L }
    } else if (s == "cured" && seasons == "multiple") {
      run <- if (g >= greenup) run + 1L else 0L
      if (run >= persist) s <- "greenup"
    }
    st[t] <- s
  }
  factor(st, levels = PHASES)
}

#' Start indices of each distinct green-up pulse in a phase series.
greenup_pulses <- function(phase) {
  is_gu <- phase == "greenup"
  which(is_gu & c(TRUE, !is_gu[-length(is_gu)]))
}

# ---- live fuel moisture ------------------------------------------------------

#' Map relative GSI to live fuel moisture, after NFDRS2016.
#'
#' Jolly, W.M., Freeborn, P.H., Bradshaw, L.S., Wallace, J. & Brittain, S.
#' (2024). Modernizing the US National Fire Danger Rating System (version 4).
#' Environmental Modelling & Software 181, 106181. Equation 10:
#'
#'     LFM = Min                                        if GSI' <  GU
#'           Min + (Max - Min) * (GSI' - GU) / (1 - GU)  if GSI' >= GU
#'
#' The ramp starts at the green-up threshold GU, NOT at zero. Earlier versions
#' of this function used the plain proportional form, which reports far too much
#' moisture through the lower half of the index -- at GSI' = 0.5 it gives 140 %
#' herbaceous where Eq. 10 with GU = 0.2 gives 112 %. Against the 120 % load
#' transfer threshold that is the difference between none of the fine fuel
#' counting as dead and some of it doing so. Setting gu = 0 recovers the old
#' behaviour exactly, which is why it stays a parameter rather than a constant.
#'
#' GU is a third undocumented lever, alongside the daylength and VPD
#' conventions. Three sources give three values: Jolly et al. (2024) Table 6
#' says 0.2; NWCG PMS 437 describes herbaceous moisture sitting at its minimum
#' below GSI 0.5; the FEMS operational defaults use a green-up threshold of 0.3.
#' At GSI' = 0.5 those give 112 %, 30 % and 96 % herbaceous moisture -- 8 %,
#' 100 % and 27 % of the fine fuel load transferred to dead. Establish which one
#' the operational code uses before tuning anything against it.
#'
#' Note what this mapping still means: herbaceous and woody moisture are the
#' SAME curve rescaled to different endpoints. Rooting depth, response lag, and
#' whether senescence is reversible are all represented by nothing more than a
#' different pair of numbers -- and now a second GU, which is the first place
#' the two fuel classes are allowed to differ in shape rather than in scale.
#'
#' @param gu green-up threshold below which moisture is held at its minimum
#'
#' @param gate if TRUE, hold at the minimum whenever the canopy is not in a
#'   green phase, i.e. both before green-up and after curing. Gating on curing
#'   matters: in single-season mode it is what makes an autumn GSI rise produce
#'   no moisture response, which is the honest depiction of a one-season-per-year
#'   assumption applied to a bimodal climate.
gsi_to_lfm <- function(gsi_rel, phase, lo, hi, gate = TRUE, gu = 0.2) {
  gu <- min(max(gu, 0), 0.999)          # gu = 1 would divide by zero
  v  <- ifelse(gsi_rel < gu, lo,
               lo + (hi - lo) * (gsi_rel - gu) / (1 - gu))
  if (gate) v[phase %in% c("dormant", "cured")] <- lo
  v
}

#' Fraction of live herbaceous load transferred to 1-h dead.
#' Approximate NFDRS2016 behaviour: transfer begins as herbaceous fuel moisture
#' falls below `none` and is complete at `full`. This is the step that actually
#' moves ERC and BI, so it belongs on screen next to the moisture curve.
herb_load_transfer <- function(lhfm, full = 30, none = 120) {
  pmin(1, pmax(0, (none - lhfm) / (none - full)))
}

# ---- parameters --------------------------------------------------------------

#' NFDRS2016 default parameter set.
gsi_defaults <- function() {
  list(
    tmin_lo = -2, tmin_hi = 5,
    vpd_lo = 900, vpd_hi = 4100,
    photo_lo = 36000, photo_hi = 39600, photo_method = "geometric",
    use_precip = FALSE, precip_window = 21, precip_lo = 0, precip_hi = 25,
    combine = "product", w_tmin = 1, w_vpd = 1, w_photo = 1, w_precip = 1,
    smooth = "box_trailing", window = 21,
    halflife_rise = 7, halflife_fall = 21,
    ratchet = FALSE,
    gsi_max = 1.0, greenup = 0.5, maint = 0.8, persist = 3,
    maint_ref = "gsimax", seasons = "single",
    lhfm_lo = 30, lhfm_hi = 250,
    lwfm_lo = 60, lwfm_hi = 200,
    gu_herb = 0.2, gu_woody = 0.2,   # Jolly et al. (2024) Table 6
    gate = TRUE
  )
}

# ---- main driver -------------------------------------------------------------

#' Run the GSI model over a daily met series.
#'
#' @param met data.frame with date, doy, tmin_c, vpd_pa, precip_mm
#' @param lat site latitude, decimal degrees
#' @param p parameter list, see gsi_defaults()
#' @return the input augmented with sub-indices, GSI, phase and fuel moistures
run_gsi <- function(met, lat, p = gsi_defaults()) {
  d <- met
  d$photo_s <- photoperiod_seconds(d$doy, lat, p$photo_method %||% "geometric")

  d$i_tmin  <- idx_tmin(d$tmin_c, p$tmin_lo, p$tmin_hi)
  d$i_vpd   <- idx_vpd(d$vpd_pa, p$vpd_lo, p$vpd_hi)
  d$i_photo <- idx_photo(d$photo_s, p$photo_lo, p$photo_hi)

  idx <- list(tmin = d$i_tmin, vpd = d$i_vpd, photo = d$i_photo)
  w   <- c(p$w_tmin, p$w_vpd, p$w_photo)

  if (isTRUE(p$use_precip)) {
    d$i_precip <- idx_precip(d$precip_mm, p$precip_window, p$precip_lo, p$precip_hi)
    idx$precip <- d$i_precip
    w <- c(w, p$w_precip)
  } else {
    d$i_precip <- NA_real_
  }

  d$igsi <- gsi_combine(idx, method = p$combine, weights = w)
  d$gsi  <- gsi_smooth(d$igsi, method = p$smooth, window = p$window,
                       halflife_rise = p$halflife_rise,
                       halflife_fall = p$halflife_fall)
  if (isTRUE(p$ratchet)) d$gsi <- apply_curing_ratchet(d$gsi)

  d$gsi_rel <- pmin(1, d$gsi / p$gsi_max)

  # Resolve the maintenance level. Held against GSImax (the NFDRS convention) a
  # site whose GSI never reaches `maint` never enters maintenance at all, and
  # every downstream phase silently fails to occur -- so the alternative
  # reference is offered explicitly rather than left as a trap.
  maint_abs <- if (identical(p$maint_ref %||% "gsimax", "annual_max")) {
    p$maint * max(d$gsi_rel, na.rm = TRUE)
  } else {
    p$maint
  }
  attr(d, "maint_abs") <- maint_abs
  attr(d, "maint_reachable") <- max(d$gsi_rel, na.rm = TRUE) >= maint_abs

  d$phase <- gsi_phase(d$gsi_rel, p$greenup, maint_abs, p$persist,
                       p$seasons %||% "single")

  d$lhfm <- gsi_to_lfm(d$gsi_rel, d$phase, p$lhfm_lo, p$lhfm_hi, p$gate,
                       p$gu_herb %||% 0.2)
  d$lwfm <- gsi_to_lfm(d$gsi_rel, d$phase, p$lwfm_lo, p$lwfm_hi, p$gate,
                       p$gu_woody %||% 0.2)
  d$herb_dead_frac <- herb_load_transfer(d$lhfm)

  # which sub-index is binding on each day
  M <- do.call(cbind, idx)
  colnames(M) <- names(idx)
  d$binding <- factor(colnames(M)[max.col(-M, ties.method = "first")],
                      levels = c("tmin", "vpd", "photo", "precip"))
  d$binding_value <- M[cbind(seq_len(nrow(M)), max.col(-M, ties.method = "first"))]

  d
}

# ---- season summary ----------------------------------------------------------

#' Key dates and durations a fire manager actually asks about.
gsi_summary <- function(d) {
  first_of <- function(ph) {
    i <- which(d$phase == ph)
    if (length(i)) d$date[i[1]] else as.Date(NA)
  }
  n_of <- function(ph) sum(d$phase == ph)
  green_days <- sum(d$phase %in% c("greenup", "maintenance", "senescence"))
  pulses <- greenup_pulses(d$phase)
  list(
    n_pulses        = length(pulses),
    pulse_dates     = d$date[pulses],
    maint_reachable = isTRUE(attr(d, "maint_reachable")),
    greenup_date    = first_of("greenup"),
    maint_date      = first_of("maintenance"),
    senescence_date = first_of("senescence"),
    cured_date      = first_of("cured"),
    season_length   = green_days,
    maint_days      = n_of("maintenance"),
    peak_gsi        = max(d$gsi, na.rm = TRUE),
    peak_date       = d$date[which.max(d$gsi)],
    min_lhfm        = min(d$lhfm, na.rm = TRUE),
    max_lhfm        = max(d$lhfm, na.rm = TRUE),
    days_full_cure  = sum(d$herb_dead_frac >= 0.999)
  )
}

# ---- sensitivity -------------------------------------------------------------

#' Numerical sensitivity of smoothed GSI to each threshold.
#'
#' Perturbs one threshold at a time and returns the resulting daily change in
#' GSI. This is what turns slider fiddling into reasoning: it shows which knob
#' can actually move the curve on a given date, and which are inert.
#'
#' Steps are deliberately COMMENSURATE -- each threshold is shifted by
#' `frac` of its own ramp width, so a degree of Tmin and a Pascal of VPD are
#' put on the same footing. Comparing raw unit steps across knobs with
#' different natural scales produces a ranking that is an artefact of units.
#'
#' @param frac fraction of each indicator's ramp width to perturb by
#' @return long data.frame: date, knob, step (in native units), d_gsi
gsi_sensitivity <- function(met, lat, p = gsi_defaults(), frac = 0.05) {
  base <- run_gsi(met, lat, p)$gsi

  span_tmin   <- p$tmin_hi - p$tmin_lo
  span_vpd    <- p$vpd_hi - p$vpd_lo
  span_photo  <- p$photo_hi - p$photo_lo
  steps <- list(
    tmin_lo  = frac * span_tmin,  tmin_hi  = frac * span_tmin,
    vpd_lo   = frac * span_vpd,   vpd_hi   = frac * span_vpd,
    photo_lo = frac * span_photo, photo_hi = frac * span_photo
  )
  if (isTRUE(p$use_precip)) {
    span_precip <- p$precip_hi - p$precip_lo
    steps$precip_lo <- frac * span_precip
    steps$precip_hi <- frac * span_precip
    steps$precip_window <- max(1, round(frac * p$precip_window))
  }
  steps$window <- max(1, round(frac * p$window))

  out <- lapply(names(steps), function(k) {
    p2 <- p
    p2[[k]] <- p[[k]] + steps[[k]]
    data.frame(
      date  = met$date,
      knob  = k,
      step  = steps[[k]],
      d_gsi = run_gsi(met, lat, p2)$gsi - base,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

## ---- from GrowingSeasonIndex/study/R/study_utils.R (unchanged) ----------------------------
#' Green-up and end-of-season dates from a GSI run, per calendar year.
#'
#' Onset = first day gsi_rel crosses `thr` upward and stays above for `persist`
#' days. End = last day, after the annual peak, that gsi_rel is >= `thr`.
#' These are deliberately simple threshold crossings, NOT the NFDRS phase
#' machine -- the point is to find which threshold (if any) matches the
#' observed phenology, so the threshold has to be the free variable.
#' Southwest sites can have two pulses, so a second onset (>= 60 days after the
#' first) is returned as `onset2_doy`.
gsi_dates <- function(run, thr = 0.5, persist = 3) {
  yrs <- sort(unique(run$year))
  do.call(rbind, lapply(yrs, function(y) {
    r <- run[run$year == y, ]
    above <- r$gsi_rel >= thr
    rl <- rle(above)
    ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1
    # a run already under way on Jan 1 is not an onset (evergreen/warm sites)
    ok <- rl$values & rl$lengths >= persist & starts > 1
    on_all <- r$doy[starts[ok]]
    ipk <- which.max(r$gsi_rel)
    after <- which(above & seq_along(above) >= ipk)
    later <- on_all[on_all >= (on_all[1] %||% Inf) + 60]
    data.frame(year = y,
               thr = thr,
               onset_doy  = if (length(on_all)) on_all[1] else NA_integer_,
               onset2_doy = if (length(later))  later[1]  else NA_integer_,
               end_doy   = if (length(after)) r$doy[max(after)] else NA_integer_,
               peak_doy  = r$doy[ipk],
               gsi_max_yr = max(r$gsi_rel, na.rm = TRUE))
  }))
}

## ---- NFDRS v4 parameters, as in GrowingSeasonIndex/study/05_build_tables.R ----------------
npn_gsi_params_v4 <- function() modifyList(gsi_defaults(), list(
  tmin_lo = -2, tmin_hi = 5, vpd_lo = 900, vpd_hi = 4100,
  photo_lo = 36000, photo_hi = 39600, photo_method = "geometric",
  use_precip = TRUE, precip_window = 28, precip_lo = 0, precip_hi = 10,
  combine = "product", smooth = "box_trailing", window = 28,
  gsi_max = 1, greenup = 0.2, gu_herb = 0.2, gu_woody = 0.2,
  lhfm_lo = 30, lhfm_hi = 250, lwfm_lo = 60, lwfm_hi = 200,
  gate = FALSE          # Eq. 10 alone; the phase machine is not part of v4
))
