# CLAUDE.md: gsi-pheno-regions

Fuel phenology regionalization for CONUS herbaceous fuels, part of the Growing Season
Index (GSI) for Wildland Fire project (Mike Crimmins, University of Arizona Cooperative
Extension). The full plan lives in the claude.ai project doc
`claude/gsi-pheno-regions-plan.md`; this file is the working summary for agents in the repo.

## Goal
Find regions with roughly common green-up, maintenance and curing patterns, so regional
phenology models can be parameterized within them. Regions are parameter domains by
functional group (herbaceous first) and can later roll up to the two NFDRS live classes.

## Design: two independent cuts, then compare
- **Cut A (climate/controls):** cluster on climate descriptors from **daily PRISM 4 km**.
  These regions are hypotheses about drivers, not tests of phenology.
- **Cut B (satellite curves):** cluster on observed seasonal greenness/curing curves
  (MODIS MOD13C1 0.05°, a water-sensitive index). This is the response and the test of A.
- Agreement = solid basis; disagreement = where to look for missing drivers.
- Stay in R at the 4 km PRISM grid. No Earth Engine.

### Cut A feature blocks (nested runs A1 ⊂ A2 ⊂ A3)
Every feature is computed **per year, then summarized across years** (median and spread).
Never compute features from mean-year climate. Daylength from latitude.
- **Block 1, seasonal climate:** seasonal T and P, fraction of P by season, aridity
  (P/PET, Hargreaves PET from Tmin/Tmax), freeze-free season length, GDD.
- **Block 2, daily variability and timing:** interannual variability; timing of season-start
  cues (last spring freeze, warm-season rain onset) and their year-to-year spread; rain-pulse
  frequency and dry-spell length; false starts; rule-based limiting-factor strata.
- **Block 3, GSI sub-indices:** daily Tmin, VPD, daylength and moisture ramps per year
  (reuse `gsi.R` from the existing GSI code; thresholds must be parameters); limiting-factor
  shares; onset/peak/decline dates from the 28-day smoothed index and their spread.
  No GSI-to-LFM mapping. Run with and without the VPD ramp (PRISM VPD is least reliable).
- **Static:** elevation (PRISM DEM), latitude, soil AWC/texture if easy.
- A1 = Block 1 + static; A2 = A1 + Block 2; A3 = A2 + Block 3.
- Standardize, reduce within each block (e.g. PCA), weight blocks equally.
- Precip event counts: thresholds above trace; test sensitivity (interpolation, day boundary).

### First-cut deliverables
Herbaceous mask (exclude cropland); Cut A clusters for A1–A3 (k-means; mclust on the P720)
with spatial smoothing; stability-based k; scoreboard vs EPA ecoregions, PSAs, NFDRS climate
classes; sampling-gap map (Globe-LFMC, NPN, iNaturalist). Then Cut B.
Deferred: fire data, gridMET, ERA5-Land.

## Open items: do NOT decide these; ask Mike
- Study period for Cut A (the `full` profile years are a placeholder, 1981–2025).
- Block 2 thresholds (rain-event size, onset rule, freeze definition).
- Moisture-ramp choice for Block 3 (precip index vs simple Hargreaves water balance).
- Copy vs reference `gsi.R`; P720 disk space and any existing local PRISM archive.

## Compute workflow
- Develop on the Windows laptop with the `dev` profile (Arizona box, 2023–2025, 2 workers);
  run at scale on the Ubuntu ThinkStation P720 (dual Xeon Gold, 128 GB) with `full`,
  via `Rscript` in `tmux` over SSH. Smoke test there with `dev` (or `test`) first.
- **Everything machine- or run-specific comes from `config.yml`** (data root, bbox, years,
  variables, workers, memfrac). Scripts are identical on both machines.
- Profiles: `test` (5 days, 2 vars, minutes), `dev`, `full`. Select with the first script
  argument, a `gsi_profile` variable in the global env (RStudio console / background job
  with `importEnv = TRUE`), or `R_CONFIG_ACTIVE`. Scripts must not call `quit()` when interactive. `GSI_DATA_ROOT` overrides the data root.
- Data never goes in git. Layout under `<data_root>/<run_name>/`:
  - `raw/prism_daily/<var>/<year>/` temporary daily tifs (deleted after conversion)
  - `prism_daily/<var>/prism_<var>_<year>.tif` one band per day, band names = ISO dates
  - later: `features/<block>/<year>/...`, `summaries/`, `clusters/`, `eval/`
- Logs go to `logs/` in the repo (git-ignored), one file per run.

## Coding conventions
- Scripts in `scripts/NN_name.R`, run from the project root; shared functions in `R/`,
  sourced with `here::here()`. No package structure.
- Resume-safe: skip existing outputs; write to a temp path and rename on success.
- Log with `log_msg()/log_warn()/log_err()` from `R/log.R`.
- Parallel: `future` + `furrr` with `plan(multisession)` only (Windows and Linux).
  Pass **file paths** to workers, never terra objects (SpatRasters don't serialize).
  Set `terraOptions(memfrac)` per worker from config; single-threaded BLAS.
- Work year by year; per-year outputs so blocks can be rerun without touching raw data.
- Raster outputs: GeoTIFF, FLT4S, `COMPRESS=DEFLATE, PREDICTOR=3, TILED=YES, INTERLEAVE=BAND`.
- Dependencies managed with `renv`; run `renv::snapshot()` after adding a package.
- Line endings LF (`.gitattributes`).

## PRISM download rules (checked Oct 2026)
- Web service: `https://services.nacse.org/prism/data/get/<region>/<res>/<var>/<YYYYMMDD>`
  returns a zip with a Cloud Optimized GeoTIFF. Docs: PRISM_downloads_web_service.pdf.
- A file may be downloaded **at most twice per 24 h**; excessive use can get the IP blocked.
  Downloads are sequential with a 2 s pause; days on disk are never re-fetched; one retry max.
- Data older than ~6 months are stable until a new time-series version is released.
- Free to use and redistribute; cite: "PRISM Group, Oregon State University,
  https://prism.oregonstate.edu, accessed <date>." PRISM advises against very long-term
  trend analysis with these grids.

## Pipeline
1. `scripts/01_prism_daily.R`: download + convert daily PRISM (ppt, tmin, tmax, vpdmax). **Built.**
2. Static layers, herbaceous mask. 3. Per-year features, one script per block.
4. Across-year summaries. 5. Clustering. 6. Evaluation.
