# GSI model for gsi-pheno-regions (Block 3 and the GSI-per-region model test).
#
# Part 1 is copied from gsi-scout R/gsi.R (Mike Crimmins' GSI Scout app; file dated
# 2026-09-22, decided Oct 8 2026 to copy rather than reference) with two marked changes:
# compute_kbdi() takes an optional mean annual precipitation R, and run_gsi() passes
# p$kbdi_R_in to it. Part 1 is the point model, used as the reference in
# scripts/12b_gsi_check.R. Do not edit Part 1 except to re-sync it with gsi-scout.
#
# Part 2 is the grid version: the same calculations on matrices of cells x days (one row
# chunk), vectorized over cells, used by scripts/12_gsi.R. It must reproduce Part 1 exactly
# at any cell (checked by 12b).

# ==== Part 1: gsi-scout point model (copied) ===============================================

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

#' Minimum temperature indicator. FEMS default: -2 C to 5 C.
idx_tmin <- function(tmin_c, lo = -2, hi = 5) {
  gsi_ramp(tmin_c, lo, hi)
}

#' Vapour pressure deficit indicator. FEMS default: 1956 Pa to 3882 Pa.
#' Decreasing: unconstrained below `lo`, fully limiting above `hi`.
idx_vpd <- function(vpd_pa, lo = 1956, hi = 3882) {
  gsi_ramp(vpd_pa, lo, hi, decreasing = TRUE)
}

#' Photoperiod indicator. FEMS default: 39600 s (11 h) to 43200 s (12 h).
idx_photo <- function(photo_s, lo = 39600, hi = 43200) {
  gsi_ramp(photo_s, lo, hi)
}

#' Precipitation indicator (Daham et al. 2018), generalised.
#' FEMS default: 28-day running total, limiting at/below 0.4 in (10.16 mm),
#' unconstrained at/above 0.8 in (20.32 mm). Off by default -- see
#' scout_defaults()$use_precip.
idx_precip <- function(precip_mm, window = 28, lo = 10.16, hi = 20.32) {
  acc <- roll_sum_trailing(precip_mm, window)
  gsi_ramp(acc, lo, hi)
}

#' Soil moisture indicator -- GSI Scout's own addition, not in NFDRS2016 or
#' FEMS. Volumetric water content (m3/m3), from ERA5 via Open-Meteo at one of
#' three depths (see fetch.R): 0-7 cm, 7-28 cm, or 28-100 cm. Which depth
#' suits herbaceous vs. woody fuels is exactly the open question Cheryl and
#' Nick want to explore with this tool -- it is a UI control, not a constant.
#' Off by default -- see scout_defaults()$use_soilm.
idx_soilm <- function(sm_frac, lo = 0.10, hi = 0.25) {
  gsi_ramp(sm_frac, lo, hi)
}

#' KBDI indicator -- GSI Scout's own addition (2026-09-21, at Mike's request:
#' a land manager he's working with is already familiar with KBDI and wants
#' to see it used as a ramp into the GSI product), not in NFDRS2016 or FEMS.
#' Decreasing: unconstrained (index = 1) below `lo`, fully limiting (index =
#' 0) at/above `hi` -- high KBDI means deep drought, which should suppress
#' GSI, the same direction as VPD. lo/hi default to the drying/duff-burning
#' transition of the published TAMU/TICC interpretive bands (200 = "litter
#' and duff beginning to dry" to 600 = "duff burns actively and combines with
#' the deep litter to actively burn") -- see compute_kbdi()'s header and
#' docs/gsi-mechanics.qmd for the full band table and citation. Off by
#' default -- see scout_defaults()$use_kbdi.
idx_kbdi <- function(kbdi, lo = 200, hi = 600) {
  gsi_ramp(kbdi, lo, hi, decreasing = TRUE)
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

# ---- photoperiod -------------------------------------------------------------

#' Sun-elevation offsets defining the end of "day", in degrees below horizon.
PHOTO_CONVENTIONS <- c(
  geometric      =  0.000,  # sun centre on the horizon; FAO-56 / Allen et al.
  sunrise_sunset = -0.833,  # sun's upper limb + refraction: published daylength
  civil_twilight = -6.000   # usable light; arguably the biological cue
)

#' Daylength in seconds from day-of-year and latitude.
#'
#' `method = "geometric"` reproduces the FAO-56 relation and is the formulation
#' Jolly et al. (2005) and Daham et al. (2018) actually use -- GSI Scout's
#' default, chosen (2026-09-14) to match the literature over the alternative
#' of numerically matching FEMS's threshold band. That means the seeded FEMS
#' Day Length ramp (11-12 h) will NOT sit naturally against geometric
#' daylength values for many sites -- geometric daylength never reaches 11 h
#' in the desert Southwest in winter. That's intentional: it's the same
#' "here's where the national default doesn't fit the West" moment the FEMS
#' overlay exists to surface, not a bug. See claude/gsi-pointapp-scope.md.
#'
#' The convention stays a UI selector precisely because it is not a cosmetic
#' choice: at Tucson (32.2 N) the December solstice is 9.89 h geometric,
#' 10.04 h sunrise-to-sunset, 10.94 h civil twilight -- and against an 11 h
#' threshold floor that choice alone decides whether iPhoto forces GSI to
#' zero across the desert Southwest every winter.
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

# ---- Keetch-Byram Drought Index (KBDI) ----------------------------------------
#
# GSI Scout's own addition (2026-09-21, at Mike's request), not in NFDRS2016
# or FEMS -- brainstormed with Mike before any code was written, then built
# to four explicit decisions: stacked alongside precip/soil moisture as a
# fourth independent optional sub-index (not replacing either, despite the
# known double-counting risk -- KBDI's own composite nature already includes
# both temperature and rainfall); mean annual precipitation derived
# automatically from ERA5 (no user-override slider); net rainfall computed
# per-storm-event, not per-day (more code, more correct -- see
# kbdi_net_rainfall()); and idx_kbdi()'s lo/hi seeded from the published
# operational interpretive bands, not tuned.
#
# Keetch, J.J. & Byram, G.M. (1968). A Drought Index for Forest Fire Control.
# USDA Forest Service Research Paper SE-38. 0-800 scale, hundredths of an
# inch of soil moisture deficit in the top 8 inches of soil -- 0 is
# saturated, 800 is the theoretical maximum drought:
#
#     KBDI_t = Q + (800 - Q) * (0.968*e^(0.0486*T) - 8.30) * dt
#              / (1 + 10.88*e^(-0.0441*R)) * 1e-3
#
# where Q is yesterday's KBDI, T is daily MAXIMUM temperature in FAHRENHEIT
# (ERA5/fetch.R gives Celsius -- converted below), R is MEAN ANNUAL
# PRECIPITATION in inches (a fixed site constant, not a daily variable), and
# dt = 1 day. Order of operations, per MetricGate's KBDI documentation
# (TAMU's own reference PDF doesn't say explicitly): net rainfall is
# subtracted from yesterday's Q FIRST (floored at 0), and it's that
# already-rain-adjusted Q that feeds the (800 - Q) drying term.
#
# Published interpretive bands (Texas A&M Forest Service / Texas
# Interagency Coordination Center KBDI fact sheet), which idx_kbdi()'s
# default lo/hi are seeded from: 0-200 "high soil/fuel moisture, little
# contribution to fire intensity"; 200-400 "lower litter and duff layers
# beginning to dry, more fires starting to carry"; 400-600 "duff and litter
# layers actively contribute to fire intensity, deeper burning"; 600-800
# "intense, deep-burning fires with more extreme fire behaviour, isolated
# heavy spotting."

#' Net rainfall for KBDI: the standard 0.20 in canopy/litter interception
#' allowance applies ONCE PER RAIN EVENT (a run of one or more consecutive
#' wet days), not once per calendar day. A storm of 0.15 in, then 0.15 in,
#' then (after a dry day resets the event) 0.25 in nets 0.00 + 0.10 + 0.05 =
#' 0.15 in total; treating each day as its own event and subtracting 0.20
#' from each independently would net only 0.05 in for the same storm -- a
#' real under-count, not a rounding nuance, which is why this needs a loop
#' over running event totals rather than a vectorised pmax(0, precip -
#' interception).
#'
#' @param precip_in daily precipitation, INCHES (convert from precip_mm first)
#' @param interception in, subtracted once per rain event (KBDI standard: 0.20)
#' @return numeric vector, same length as precip_in, net rainfall in inches
kbdi_net_rainfall <- function(precip_in, interception = 0.20) {
  n   <- length(precip_in)
  net <- numeric(n)
  event_cum <- 0  # cumulative rainfall so far in the current consecutive-wet run
  for (t in seq_len(n)) {
    p <- precip_in[t]
    if (is.na(p) || p <= 0) {
      event_cum <- 0
      next
    }
    prev_cum  <- event_cum
    event_cum <- event_cum + p
    net[t] <- (event_cum - min(event_cum, interception)) -
             (prev_cum  - min(prev_cum, interception))
  }
  net
}

#' The daily KBDI recursion (Keetch & Byram 1968), 0-800.
#'
#' Computed over met's FULL date-ordered, multi-year record -- KBDI is
#' cumulative/recursive, so (like gsi_phase2()'s persistence counter) it
#' needs its own history to spin up, and run_gsi() always calls this on the
#' whole multi-year met, the same way it does everything else, before
#' gsi_sel() ever filters down to one selected year.
#'
#' Mean annual precipitation (R) is derived from met itself (2026-09-21, at
#' Mike's request: "derive from ERA5 automatically," no user-override
#' slider), averaged over full years within `band_years` -- deliberately the
#' SAME fixed 1991-2020 reference period every other "normal" in this app
#' uses (compute_climatology() in plots.R, called with CLIMATOLOGY_BAND in
#' app.R). This has to stay independent of whichever year happens to be
#' selected: sel_year only filters gsi_sel() AFTER run_gsi() has already run
#' over the whole met, so as long as met itself doesn't change, neither does
#' R, and neither does the recursive KBDI series it feeds -- changing only
#' the year dropdown must never silently shift KBDI's trajectory. If fewer
#' than the full band is actually present in met (e.g. a point fetched with
#' year_to earlier than band_years[2] -- see fetch_point()'s year_from =
#' CLIMATOLOGY_BAND[1] in app.R), R is simply averaged over whichever band
#' years ARE present.
#'
#' @param met data.frame with date, year, tmax_c, precip_mm -- any order;
#'   sorted internally by date before the recursion runs
#' @param band_years mean-annual-precipitation reference period, c(lo, hi)
#' @return numeric vector, same length and row order as met, KBDI in [0, 800]
compute_kbdi <- function(met, band_years = c(1991, 2020), R = NULL) {
  ord <- order(met$date)
  m   <- met[ord, ]

  precip_in <- m$precip_mm / 25.4
  if (is.null(R)) {   # [gsi-pheno-regions] R (inches) can be passed in, e.g. a PRISM-period mean
    ref_years <- m$year >= band_years[1] & m$year <= band_years[2]
    annual_totals <- tapply(precip_in[ref_years], m$year[ref_years], sum, na.rm = TRUE)
    R <- mean(annual_totals, na.rm = TRUE)
    if (!is.finite(R) || R <= 0) R <- mean(precip_in, na.rm = TRUE) * 365.25  # fallback
  }

  tmax_f <- m$tmax_c * 9 / 5 + 32
  net_in <- kbdi_net_rainfall(precip_in)
  denom  <- 1 + 10.88 * exp(-0.0441 * R)

  n <- nrow(m)
  kbdi <- numeric(n)
  Q <- 0  # start saturated; the first stretch of days is effectively spin-up
  for (t in seq_len(n)) {
    Q <- max(0, Q - net_in[t] * 100)  # net_in is inches; KBDI is hundredths of an inch
    Tf <- tmax_f[t]
    if (!is.na(Tf)) {
      dQ <- (800 - Q) * (0.968 * exp(0.0486 * Tf) - 8.30) / denom * 1e-3
      Q <- min(800, max(0, Q + dQ))
    }
    kbdi[t] <- Q
  }
  out <- numeric(n)
  out[ord] <- kbdi
  out
}

# ---- combination & smoothing --------------------------------------------------

#' Combine sub-indices into the daily index. "product" is the Jolly/NFDRS2016
#' form and FEMS's operational form -- kept as the only wired-up option in the
#' scout UI for v1, but the model supports the alternatives Glass Box explores.
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

#' Smooth the daily index. FEMS default: 28-day trailing mean.
gsi_smooth <- function(x, window = 28) roll_mean_trailing(x, window)

# ---- phenophase: 2-phase (greenup / dormant) ---------------------------------

PHASES2 <- c("dormant", "greenup")

#' Classify each day into greenup vs. dormant.
#'
#' GSI Scout's phase model, distinct from Glass Box's 5-state machine (see
#' PHASES there: dormant/greenup/maintenance/senescence/cured). Cheryl's ask,
#' directly: let drought -- which acts through the precip/soil-moisture
#' sub-indices when enabled -- push a site back into dormancy at any point,
#' rather than requiring it to pass through an explicit maintenance window
#' first. One threshold (`greenup`) governs both directions; `persist`
#' requires a run of consecutive days on the far side before the phase
#' actually flips, so single-day noise near the threshold doesn't flicker it.
#'
#' A side benefit of the 2-state design: unlike Glass Box's `cured` state,
#' `dormant` here is never absorbing, so a bimodal site (e.g. the Southwest's
#' cool-season/monsoon double green-up) is represented without needing an
#' explicit "multiple seasons" switch -- it falls out of the model for free.
#'
#' @param gsi_rel GSI / GSImax, in [0, 1]
#' @param greenup relative-GSI threshold separating the two phases (FEMS: 0.3)
#' @param persist consecutive days required on the far side of `greenup`
#'   before the phase actually switches
gsi_phase2 <- function(gsi_rel, greenup = 0.3, persist = 3) {
  n <- length(gsi_rel)
  st <- character(n)
  s <- "dormant"
  run <- 0L
  for (t in seq_len(n)) {
    want <- if (gsi_rel[t] >= greenup) "greenup" else "dormant"
    if (identical(want, s)) {
      run <- 0L
    } else {
      run <- run + 1L
      if (run >= persist) { s <- want; run <- 0L }
    }
    st[t] <- s
  }
  factor(st, levels = PHASES2)
}

#' Start indices of each distinct green-up pulse in a phase series.
greenup_pulses <- function(phase) {
  is_gu <- phase == "greenup"
  which(is_gu & c(TRUE, !is_gu[-length(is_gu)]))
}

#' Day-of-year "normal" phase -- gsi_phase2() run on a day-of-year GSI
#' climatology instead of one year's actual GSI, 2026-09-20 at Mike's
#' request, so the phase strip can show whether a typical year is greenup or
#' dormant on a given day, next to whichever year is actually selected.
#'
#' Deliberately the SAME classifier (gsi_phase2()) and the SAME live
#' greenup/persist/gsi_max settings the selected year's own phase strip
#' uses (from `p`, not scout_defaults()) -- an apples-to-apples read of
#' "given how I've set green-up sensitivity, is this year ahead of, behind,
#' or on pace with normal," not a comparison against a different threshold.
#'
#' @param clim compute_climatology(gsi_all_years, "gsi", band_years, sel_year)
#'   output (doy, lo, mean, hi) -- mean is on the same absolute GSI scale
#'   run_gsi() itself uses, not gsi_rel, hence dividing by p$gsi_max below
#' @param p params() list -- gsi_max/greenup/persist
#' @return data.frame(doy, phase), one row per day-of-year `clim` covers
climatology_phase <- function(clim, p) {
  clim <- clim[order(clim$doy), ]
  rel <- pmin(1, clim$mean / p$gsi_max)
  data.frame(doy = clim$doy, phase = gsi_phase2(rel, p$greenup, p$persist))
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
#' @param gate if TRUE, hold at the minimum whenever the site is dormant.
gsi_to_lfm <- function(gsi_rel, phase, lo, hi, gate = TRUE, gu = 0.2) {
  gu <- min(max(gu, 0), 0.999)          # gu = 1 would divide by zero
  v  <- ifelse(gsi_rel < gu, lo,
               lo + (hi - lo) * (gsi_rel - gu) / (1 - gu))
  if (gate) v[phase == "dormant"] <- lo
  v
}

#' Fraction of live herbaceous load transferred to 1-h dead.
herb_load_transfer <- function(lhfm, full = 30, none = 120) {
  pmin(1, pmax(0, (none - lhfm) / (none - full)))
}

# ---- parameters --------------------------------------------------------------

#' FEMS-seeded default parameter set for GSI Scout.
#'
#' Source: Cheryl's FEMS operational parameter reference (photographed table,
#' convo_on_gsi.docx), transcribed in claude/gsi-pointapp-scope.md. This is
#' NOT gsi_defaults() from Glass Box -- these are FEMS's tuned operational
#' values, which differ from the NFDRS2016 literature defaults Glass Box
#' seeds on (most visibly: VPD 1956-3882 Pa here vs. 900-4100 Pa there).
#'
#' Every value here is a UI seed, not a constant baked into the model --
#' run_gsi() takes the full parameter list and nothing defaults silently.
scout_defaults <- function() {
  list(
    tmin_lo = -2, tmin_hi = 5,
    vpd_lo = 1956, vpd_hi = 3882, vpd_driver = "native_max",  # FEMS: "VPD max"
    photo_lo = 39600, photo_hi = 43200, photo_method = "geometric",  # FAO-56
    use_precip = FALSE, precip_window = 28, precip_lo = 10.16, precip_hi = 20.32,
    use_soilm = FALSE, soilm_depth = "sm_0_7", soilm_lo = 0.10, soilm_hi = 0.25,
    use_kbdi = FALSE, kbdi_lo = 200, kbdi_hi = 600,  # TAMU/TICC bands, see gsi.R
    combine = "product",
    window = 28,
    gsi_max = 1.0, greenup = 0.3, persist = 3,
    lhfm_lo = 30, lhfm_hi = 250,
    lwfm_lo = 60, lwfm_hi = 200,
    gu_herb = 0.2, gu_woody = 0.2,   # Jolly et al. (2024) Table 6; FEMS is silent
    gate = TRUE
  )
}

# ---- main driver -------------------------------------------------------------

#' Run the GSI Scout model over a daily met series for one point.
#'
#' @param met data.frame with date, doy, tmin_c, vpd_pa, precip_mm, and (if
#'   soil moisture is toggled on) whichever of sm_0_7 / sm_7_28 / sm_28_100
#'   `p$soilm_depth` names.
#' @param lat site latitude, decimal degrees
#' @param p parameter list, see scout_defaults()
#' @return the input augmented with sub-indices, GSI, phase and fuel moistures
run_gsi <- function(met, lat, p = scout_defaults()) {
  d <- met
  d$photo_s <- photoperiod_seconds(d$doy, lat, p$photo_method %||% "geometric")

  d$i_tmin  <- idx_tmin(d$tmin_c, p$tmin_lo, p$tmin_hi)
  d$i_vpd   <- idx_vpd(d$vpd_pa, p$vpd_lo, p$vpd_hi)
  d$i_photo <- idx_photo(d$photo_s, p$photo_lo, p$photo_hi)

  idx <- list(tmin = d$i_tmin, vpd = d$i_vpd, photo = d$i_photo)

  if (isTRUE(p$use_precip)) {
    d$i_precip <- idx_precip(d$precip_mm, p$precip_window, p$precip_lo, p$precip_hi)
    idx$precip <- d$i_precip
  } else {
    d$i_precip <- NA_real_
  }

  if (isTRUE(p$use_soilm)) {
    depth_col <- p$soilm_depth %||% "sm_0_7"
    sm <- met[[depth_col]]
    if (is.null(sm)) stop("run_gsi: met has no column '", depth_col, "'")
    d$i_soilm <- idx_soilm(sm, p$soilm_lo, p$soilm_hi)
    idx$soilm <- d$i_soilm
  } else {
    d$i_soilm <- NA_real_
  }

  if (isTRUE(p$use_kbdi)) {
    d$kbdi   <- compute_kbdi(met, R = p$kbdi_R_in)   # [gsi-pheno-regions] NULL -> gsi-scout behaviour
    d$i_kbdi <- idx_kbdi(d$kbdi, p$kbdi_lo, p$kbdi_hi)
    idx$kbdi <- d$i_kbdi
  } else {
    d$kbdi   <- NA_real_
    d$i_kbdi <- NA_real_
  }

  d$igsi <- gsi_combine(idx, method = p$combine %||% "product")
  d$gsi  <- gsi_smooth(d$igsi, window = p$window %||% 28)
  d$gsi_rel <- pmin(1, d$gsi / p$gsi_max)

  d$phase <- gsi_phase2(d$gsi_rel, p$greenup, p$persist)

  d$lhfm <- gsi_to_lfm(d$gsi_rel, d$phase, p$lhfm_lo, p$lhfm_hi, p$gate,
                       p$gu_herb %||% 0.2)
  d$lwfm <- gsi_to_lfm(d$gsi_rel, d$phase, p$lwfm_lo, p$lwfm_hi, p$gate,
                       p$gu_woody %||% 0.2)
  d$herb_dead_frac <- herb_load_transfer(d$lhfm)

  # which sub-index is binding on each day
  M <- do.call(cbind, idx)
  colnames(M) <- names(idx)
  bind_levels <- c("tmin", "vpd", "photo", "precip", "soilm", "kbdi")
  d$binding <- factor(colnames(M)[max.col(-M, ties.method = "first")],
                      levels = bind_levels)
  d$binding_value <- M[cbind(seq_len(nrow(M)), max.col(-M, ties.method = "first"))]

  d
}


# ==== Part 2: grid version (matrices of cells x days) ======================================
# Same formulas as Part 1, vectorized over the rows (cells) of a matrix whose columns are
# consecutive days. Everything is row-wise; no apply() over cells.

gsi_metric_names <- c("gsi_mean", "gsi_max", "peak_doy", "sos20", "eos50", "gu_doy", "dorm_doy",
                      "green_days", "n_pulses", "lim_tmin", "lim_vpd", "lim_photo", "lim_moist")

# Variants: every param_set x moisture option, plus config extra_variants (e.g. no VPD).
gsi_variants <- function(g) {
  v <- list()
  for (ps in names(g$param_sets)) for (mo in unlist(g$moisture))
    v[[paste(ps, mo, sep = "_")]] <- list(param_set = ps, moisture = mo, vpd = TRUE)
  for (e in g$extra_variants) {
    nm <- paste(e$param_set, e$moisture, if (isFALSE(e$vpd)) "novpd" else "vpd", sep = "_")
    v[[nm]] <- list(param_set = e$param_set, moisture = e$moisture, vpd = !isFALSE(e$vpd))
  }
  v
}

# Trailing sum / mean along days, as roll_sum_trailing / roll_mean_trailing (NA -> 0;
# partial windows at the start, the mean dividing by the days available).
m_roll_sum <- function(X, n) {
  X[is.na(X)] <- 0
  C <- cbind(0, X); for (j in 2:ncol(C)) C[, j] <- C[, j - 1] + C[, j]
  idx <- seq_len(ncol(X)); lo <- pmax(0L, idx - n)
  C[, idx + 1L, drop = FALSE] - C[, lo + 1L, drop = FALSE]
}
m_roll_mean <- function(X, n) {
  idx <- seq_len(ncol(X)); lo <- pmax(0L, idx - n)
  sweep(m_roll_sum(X, n), 2, idx - lo, "/")
}

# Daylength (s), cells x days, from latitude per cell and day of year per column.
m_photoperiod <- function(lat, doy, method = "geometric") {
  h0 <- PHOTO_CONVENTIONS[[method]] * pi / 180
  phi <- lat * pi / 180
  delta <- 0.409 * sin(2 * pi / 365 * doy - 1.39)
  arg <- (sin(h0) - outer(sin(phi), sin(delta))) / outer(cos(phi), cos(delta))
  arg[arg > 1] <- 1; arg[arg < -1] <- -1                     # (keeps dim, unlike pmin(1, .))
  acos(arg) * (24 / pi) * 3600
}

# KBDI recursion (Keetch & Byram 1968) as compute_kbdi(): net rainfall per rain event, start
# saturated (Q = 0) at the first column. P mm, TX degC, R_in mean annual precipitation (in)
# per cell. Returns cells x days.
m_kbdi <- function(P, TX, R_in, interception = 0.2) {
  n <- ncol(P); m <- nrow(P)
  Pin <- P / 25.4
  denom <- 1 + 10.88 * exp(-0.0441 * R_in)
  K <- matrix(0, m, n); Q <- numeric(m); ev <- numeric(m)
  for (t in seq_len(n)) {
    p <- Pin[, t]; wet <- !is.na(p) & p > 0
    prev <- ev; ev <- ifelse(wet, ev + p, 0)
    net <- ifelse(wet, (ev - pmin(ev, interception)) - (prev - pmin(prev, interception)), 0)
    Q <- pmax(0, Q - net * 100)
    Tf <- TX[, t] * 9 / 5 + 32; ok <- !is.na(Tf)
    dQ <- (800 - Q) * (0.968 * exp(0.0486 * Tf) - 8.30) / denom * 1e-3
    Q[ok] <- pmin(800, pmax(0, Q[ok] + dQ[ok]))
    K[, t] <- Q
  }
  K
}

# Two-phase classifier (gsi_phase2) along days: TRUE = greenup. Starts dormant.
m_phase2 <- function(R, greenup, persist) {
  n <- ncol(R); m <- nrow(R)
  S <- matrix(FALSE, m, n); s <- logical(m); run <- integer(m)
  for (t in seq_len(n)) {
    want <- R[, t] >= greenup
    same <- want == s
    run <- ifelse(same, 0L, run + 1L)
    flip <- !same & run >= persist
    s[flip] <- want[flip]; run[flip] <- 0L
    S[, t] <- s
  }
  S
}

# gsi_ramp() for matrices: same values, keeps the dimensions (pmin/pmax take attributes
# from their first argument, which drops dim when that is a scalar).
m_ramp <- function(x, lo, hi, decreasing = FALSE) {
  y <- if (isTRUE(all.equal(lo, hi))) (x >= hi) + 0 else (x - lo) / (hi - lo)
  y[y < 0] <- 0; y[y > 1] <- 1
  if (decreasing) 1 - y else y
}

# Sub-index matrices shared by all variants of one parameter set.
gsi_subindices <- function(D, ps, g, photo_s, kbdi) {
  list(tmin = m_ramp(D$tmin, ps$tmin[[1]], ps$tmin[[2]]),
       vpd = m_ramp(D$vpd_pa, ps$vpd[[1]], ps$vpd[[2]], decreasing = TRUE),
       photo = m_ramp(photo_s, ps$photo_h[[1]] * 3600, ps$photo_h[[2]] * 3600),
       precip = m_ramp(m_roll_sum(D$ppt, ps$precip$window), ps$precip$lo, ps$precip$hi),
       kbdi = if (!is.null(kbdi)) m_ramp(kbdi, g$kbdi$lo, g$kbdi$hi, decreasing = TRUE) else NULL)
}

# Season of the smoothed GSI in year Y (peak in Y; relative 20 % / 50 % crossings within
# +/- season_halfwidth_days of the peak, as for the satellite curves). G: cells x days.
gsi_season_from_G <- function(G, dd, iy, g) {
  m <- nrow(G); n <- ncol(G); rows <- seq_len(m)
  out <- matrix(NA_real_, m, 5, dimnames = list(NULL, c("gsi_mean", "gsi_max", "peak_doy", "sos20", "eos50")))
  Gy <- G[, iy, drop = FALSE]
  out[, "gsi_mean"] <- rowMeans(Gy)
  pmx <- apply(Gy, 1, max)
  pk <- first_true(Gy >= pmx - 1e-6)             # first day at the maximum (plateaus at 1 are common)
  out[, "gsi_max"] <- pmx
  pcol <- iy[pk]; out[, "peak_doy"] <- dd[pcol]
  # relative crossings around the peak (as Cut B): base = min within +/- hw days
  hw <- g$season_halfwidth_days
  base <- pmx
  for (k in seq_len(hw)) {
    for (cc in list(pcol - k, pcol + k)) {
      ok <- cc >= 1 & cc <= n
      base[ok] <- pmin(base[ok], G[cbind(rows[ok], cc[ok])])
    }
  }
  amp <- pmx - base; seas <- amp >= g$min_amp
  thr20 <- base + 0.2 * amp; thr50 <- base + 0.5 * amp
  sos <- rep(NA_real_, m); eos <- rep(NA_real_, m)
  for (k in seq_len(hw)) {
    cb <- pcol - k; ok <- is.na(sos) & cb >= 1
    hit <- ok; hit[ok] <- G[cbind(rows[ok], cb[ok])] < thr20[ok]
    sos[hit] <- dd[cb[hit]] + 1
    ca <- pcol + k; ok <- is.na(eos) & ca <= n
    hit <- ok; hit[ok] <- G[cbind(rows[ok], ca[ok])] < thr50[ok]
    eos[hit] <- dd[ca[hit]]
  }
  out[, "sos20"] <- ifelse(seas, sos, NA); out[, "eos50"] <- ifelse(seas, eos, NA)
  out
}

# One variant: daily iGSI, smoothed GSI, phase, and the per-year metrics.
# dd: day of each column relative to Jan 1 of Y (1 = Jan 1); iy: columns within year Y.
gsi_variant_metrics <- function(SI, v, ps, g, dd, iy) {
  parts <- list(tmin = SI$tmin)
  if (isTRUE(v$vpd)) parts$vpd <- SI$vpd
  parts$photo <- SI$photo
  if (v$moisture == "precip") parts$moist <- SI$precip
  if (v$moisture == "kbdi") parts$moist <- SI$kbdi
  igsi <- Reduce(`*`, parts)
  G <- m_roll_mean(igsi, ps$smooth)
  G1 <- G; G1[G1 > 1] <- 1
  m <- nrow(G); n <- ncol(G); rows <- seq_len(m)
  out <- matrix(NA_real_, m, length(gsi_metric_names), dimnames = list(NULL, gsi_metric_names))
  ss <- gsi_season_from_G(G, dd, iy, g)
  out[, colnames(ss)] <- ss
  # phases (absolute threshold on GSI / gsi_max, gsi_max = 1 as in the operational default)
  S <- m_phase2(G1, ps$greenup, g$persist)
  on <- S[, iy, drop = FALSE] & !S[, iy - 1, drop = FALSE]
  out[, "n_pulses"] <- rowSums(on)
  out[, "green_days"] <- rowSums(S[, iy, drop = FALSE])
  f <- first_true(on); has <- f <= length(iy)
  out[has, "gu_doy"] <- dd[iy[f[has]]]
  # first switch back to dormant after that green-up (may fall after Dec 31 of Y)
  gcol <- ifelse(has, iy[pmin(f, length(iy))], n + 1)
  dorm <- rep(NA_real_, m)
  for (t in (min(iy) + 1):n) {
    ok <- is.na(dorm) & has & t > gcol & !S[, t] & S[, t - 1]
    dorm[ok] <- dd[t]
  }
  out[, "dorm_doy"] <- dorm
  # limiting factor: share of days in Y on which each sub-index is the smallest and < 1
  P <- lapply(parts, function(x) x[, iy, drop = FALSE])
  mn <- Reduce(pmin, P)
  lim <- mn < 1
  taken <- matrix(FALSE, m, length(iy))
  for (nm in names(P)) {
    b <- !taken & P[[nm]] == mn & lim
    out[, paste0("lim_", nm)] <- rowMeans(b)
    taken <- taken | b
  }
  list(metrics = out, gsi = G, igsi = igsi, phase = S)
}

# One year Y for all variants (runs in a worker). files: list(var -> list(year -> path)),
# with years Y-1, Y and (when available) Y+1. Writes one GeoTIFF per variant.
gsi_year <- function(Y, files, grid_file, map_file, out_files, g, variants, chunk_rows, memfrac) {
  tryCatch({
    terra::terraOptions(memfrac = memfrac, progress = 0)
    grid <- terra::rast(grid_file); mapr <- terra::rast(map_file)
    yrs <- names(files$ppt)
    R <- lapply(files, function(v) lapply(v, terra::rast))
    dates <- do.call(c, lapply(yrs, function(y) as.Date(names(R$ppt[[y]]))))
    last_day <- as.Date(sprintf("%d-03-31", Y + 1))
    keep <- lapply(yrs, function(y) { d <- as.Date(names(R$ppt[[y]])); which(d <= last_day) })
    names(keep) <- yrs
    dates <- dates[dates <= last_day]
    if (any(diff(dates) != 1)) stop("daily series is not contiguous")
    dd <- as.numeric(dates - as.Date(sprintf("%d-01-01", Y))) + 1
    iy <- which(format(dates, "%Y") == as.character(Y))
    if (min(iy) < 2) stop("need days before Jan 1 of ", Y)
    doy <- as.integer(format(dates, "%j"))
    res <- lapply(variants, function(v) matrix(NA_real_, terra::ncell(grid), length(gsi_metric_names),
                                               dimnames = list(NULL, gsi_metric_names)))
    nr_all <- terra::nrow(grid); nc <- terra::ncol(grid)
    for (row in seq(1, nr_all, by = chunk_rows)) {
      nr <- min(chunk_rows, nr_all - row + 1)
      cells <- (row - 1) * nc + seq_len(nr * nc)
      land <- !is.na(read_rows(grid, row, nr, 1)[, 1])
      if (!any(land)) next
      rd <- function(v) do.call(cbind, lapply(yrs, function(y) read_rows(R[[v]][[y]], row, nr, keep[[y]])))[land, , drop = FALSE]
      D <- list(ppt = rd("ppt"), tmin = rd("tmin"), tmax = rd("tmax"))
      D$vpd_pa <- rd("vpdmax") * 100                                   # hPa -> Pa
      bad <- rowSums(is.na(D$ppt) | is.na(D$tmin) | is.na(D$tmax) | is.na(D$vpd_pa)) > 0
      lat <- terra::yFromRow(grid, row + (seq_len(nr) - 1))[rep(seq_len(nr), each = nc)[land]]
      photo_s <- m_photoperiod(lat, doy, g$photo_method)
      need_kbdi <- any(vapply(variants, function(v) v$moisture == "kbdi", TRUE))
      kb <- if (need_kbdi) m_kbdi(D$ppt, D$tmax, read_rows(mapr, row, nr, 1)[land, 1] / 25.4,
                                  g$kbdi$interception_in) else NULL
      SIs <- list()
      for (nm in names(variants)) {
        v <- variants[[nm]]; ps <- g$param_sets[[v$param_set]]
        if (is.null(SIs[[v$param_set]])) SIs[[v$param_set]] <- gsi_subindices(D, ps, g, photo_s, kb)
        mt <- gsi_variant_metrics(SIs[[v$param_set]], v, ps, g, dd, iy)$metrics
        mt[bad, ] <- NA
        res[[nm]][cells[land], ] <- mt
      }
    }
    for (nm in names(variants)) write_feature_matrix(res[[nm]], grid, gsi_metric_names, out_files[[nm]])
    list(ok = TRUE, year = Y, n_cells = sum(!is.na(res[[1]][, 1])), to = format(max(dates)))
  }, error = function(e) list(ok = FALSE, year = Y, msg = conditionMessage(e)),
  finally = try(terra::tmpFiles(current = TRUE, orphan = FALSE, old = FALSE, remove = TRUE), silent = TRUE))
}
