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
- **Static:** elevation (PRISM DEM), latitude. Soils deferred (decided Oct 4).
- A1 = Block 1 + static; A2 = A1 + Block 2; A3 = A2 + Block 3.
- Standardize, reduce within each block (e.g. PCA), weight blocks equally.
- Precip event counts: thresholds above trace; test sensitivity (interpolation, day boundary).

### Clustering and ML methods
Judge methods by how well regions serve as domains for regional phenology models, not by
clustering scores alone.
- **First cut:** Cut A with k-means (baseline), Gaussian mixture (`mclust`; soft membership
  maps transition zones) and Ward hierarchical (nested: fine for modeling, coarse for
  operational roll-up such as the two NFDRS live classes). Spatially constrained clustering
  (`ClustGeo`, or SKATER in `spdep`) to see what contiguity costs.
  Cut B: functional PCA of seasonal curves (`fda`, `fdapace`) + Gaussian mixture on the scores
  (`funFEM` as an alternative). Link A and B: cluster Cut B, then classify from Cut A
  features with `ranger` or `xgboost`; variable importance / SHAP tests whether Blocks 2–3 matter.
- **Worth trying:** SOMs (`kohonen`); gradient forests (`gradientForest`, maybe R-Forge) or
  GDM (`gdm`) to transform climate space by phenology response, then cluster.
- **Later:** autoencoders (`torch`, P720 GPU) benchmarked against the simpler methods;
  mixtures of regressions (`flexmix`) once LFMC / NPN / iNaturalist observations are in.
- **Evaluation for all:** bootstrap stability (`fpc::clusterboot`), adjusted Rand index between
  partitions (A1 vs A2 vs A3 vs B), spatial block cross-validation for any classifier.
- Ward, ClustGeo and SKATER need pairwise distances or graphs that won't fit at ~0.5 M CONUS
  cells; run them on k-means micro-clusters or an aggregated grid (two-stage).

### Scope and observations
Regions cover all of CONUS (decided Oct 5). NPN and iNaturalist phenology observations
(CONUS-wide) are part of preliminary prototyping: assign observations to grid cells and
regions, map sampling density per region, and use observed phenophase timing as an early,
independent check on the regions alongside Cut B.

### First-cut deliverables
Herbaceous mask (exclude cropland); Cut A clusters for A1–A3 (k-means, Gaussian mixture, Ward,
plus a spatially constrained version) with spatial smoothing; stability-based k; scoreboard vs
EPA ecoregions, PSAs, NFDRS climate classes; sampling-gap map (Globe-LFMC, NPN, iNaturalist).
Then Cut B (functional PCA + mixture) and the A-to-B classifier.
Deferred: fire data, gridMET, ERA5-Land.

## Open items: do NOT decide these; ask Mike
- Study period for Cut A. The `full` profile downloads 2001–2025 first (first full MODIS
  year onward); that is download order, not a study-period decision.
- Block 2 thresholds: decided Oct 7 (see pipeline step 5b); don't change without asking.
- Block 3 / GSI: decided Oct 8 (step 5c); don't change variants without asking.
- Working regions: decided Oct 9 (A2 k13 nested in k7; step 9).
- P720 disk space and any existing local PRISM archive.

## Compute workflow
- Develop on the Windows laptop with the `dev` profile (Arizona box, 2023–2025, 2 workers);
  run at scale on the Ubuntu ThinkStation P720 (dual Xeon Gold, 128 GB) with `full`,
  via `Rscript` in `tmux` over SSH. Smoke test there with `dev` (or `test`) first.
- **Everything machine- or run-specific comes from `config.yml`** (data root, bbox, years,
  variables, workers, memfrac). Scripts are identical on both machines.
- Profiles: `test` (5 days, 2 vars, minutes), `dev`, `full`. Select with the first script
  argument, a `gsi_profile` variable in the global env (RStudio console / background job
  with `importEnv = TRUE`), or `R_CONFIG_ACTIVE`. Worker counts: `step_workers(cfg, step)` (R/config.R)
  = env `GSI_WORKERS` (one run) > profile `workers_by_step: <step>` > profile `workers`; it
  also caps terra memfrac at 0.6 / workers. New parallel scripts should call it.
- **P720 hardware** (Oct 8): 2 x Intel Xeon Gold 6148 (20 cores each: 40 cores / 80
  threads), 125 GB RAM (+7 GB swap), NVIDIA RTX 3090 with 24 GB VRAM. Size per-year jobs to
  run in one round where memory allows (24 years -> 24 workers; 12_gsi.R ~2 GB per worker),
  keep disk-bound steps (many raster reads) near 8-16 workers. The GPU is unused so far: use
  it for heavy numeric work that fits in 24 GB, e.g. `torch` (CUDA) for autoencoders or for
  k-means / distance matrices on all ~480k cells x tens of features (fits easily), and check
  `nvidia-smi` first. Keep a CPU fallback so scripts still run on the laptop.
- Scripts must not call `quit()` when interactive. `GSI_DATA_ROOT` overrides the data root.
- `GSI_SCRATCH_ROOT` (or `scratch_root` in config) puts temporary raw downloads on a local
  disk. On the P720 data_root is the TrueNAS NFS share `/mnt/truenas_phenology` and scratch
  is local; set both in `~/.Renviron` there.
- Data never goes in git. Layout under `<data_root>/<run_name>/`:
  - `raw/prism_daily/<var>/<year>/` temporary daily tifs (deleted after conversion)
  - `prism_daily/<var>/prism_<var>_<year>.tif` one band per day, band names = ISO dates
  - later: `features/<block>/<year>/...`, `summaries/`, `clusters/`, `eval/`
- Logs go to `logs/` in the repo (git-ignored), one file per run.

## Coding conventions
- Scripts in `scripts/NN_name.R`, run from the project root; shared functions in `R/`,
  sourced with `here::here()`. No package structure.
- Resume-safe: skip existing outputs; write to a temp path and rename on success.
- Log with `log_msg()/log_warn()/log_err()` from `R/log.R`. `log_err()` also pushes an ntfy
  notification (capped per run). Long scripts call `notify_init(cfg)` and send start, progress,
  per-year summary and finish messages with `notify()` (`R/notify.R`). The topic comes from the
  `GSI_NTFY_TOPIC` env variable only; never put it in config or git. Notification failures
  must never stop a run.
- Parallel: `future` + `furrr` with `plan(multisession)` only (Windows and Linux).
  Pass **file paths** to workers, never terra objects (SpatRasters don't serialize).
  Set `terraOptions(memfrac)` per worker from config; single-threaded BLAS.
- Work year by year; per-year outputs so blocks can be rerun without touching raw data.
- Raster outputs: GeoTIFF, `COMPRESS=DEFLATE, TILED=YES, INTERLEAVE=BAND`. PRISM tmin, tmax,
  vpdmax are INT2S with GDAL scale 0.01 (set in `config.yml` `prism$scale`; terra returns
  physical units on read, max rounding error 0.005); ppt is FLT4S. Write stacks with
  `prism_write_stack()`, which QA's every band against the source before swapping the file in.
- Dates: use band names (ISO dates) as the date reference; `terra::time()` is not reliably
  stored in GeoTIFF across terra versions.
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
1. `scripts/01_prism_daily.R`: download + convert daily PRISM (ppt, tmin, tmax, vpdmax).
   `scripts/01b_prism_recode.R` re-encodes existing files to the configured storage.
   **Built; dev run verified Oct 3 2026** (12 variable-years, ~3.1 s/file, ~77 min per year
   of 4 variables).
2. `scripts/02_static.R`: `static/grid_mask.tif` (analysis domain = PRISM land cells),
   `elev.tif` (PRISM 4 km DEM, cropped if aligned with the daily grid, else resampled and
   logged), `lat.tif`, `lon.tif`. Daylength is computed from latitude where needed, not stored.
   Dev run Oct 4: DEM aligned with the daily grid (cropped, no resampling). Soils deferred.
   Grid helpers in `R/grid.R` (`to_grid()` matches any layer to the grid).
3. `scripts/03_herb_mask.R`: Annual NLCD C1.2 (30 m) -> per-year class shares on the grid
   (`static/nlcd/share_<class>_<year>.tif`, cached), averaged across years
   (`nlcd_shares.tif`), herb/crop share of land + across-year range (`herb_share.tif`), and
   `herb_mask.tif`. Decided Oct 4: herb = grassland (71) + shrub/scrub (52); years
   2001/2012/2024; herb >= 50 % of land and crops < 25 %. Pasture, wetlands, water stored but
   not herb. Rule lives in `config.yml` `herb_mask:`; changing it only re-runs the cheap
   combine step. Aggregation reads each year once (no 30 m temp files).
   CONUS run Oct 5: 155,427 mask cells (32 % of land), 34 min on the P720.
   **Scope = wall-to-wall CONUS** (decided Oct 5): Cut A clusters ALL PRISM land cells
   (`grid_mask.tif`) so NPN / iNaturalist observations anywhere fall in a region. The
   herbaceous layers describe and optionally weight cells; they are NOT the domain.
   `herb_mask.tif` (strict, 155k cells, West + Plains only) is for Cut B pixel selection
   and sensitivity runs. `open_herb_share` (grass + shrub + pasture) covers the East, where
   open herbaceous land is mostly pasture/hay. NLCD's grassland-vs-shrub split shows
   state-line artifacts (e.g. WY/NE, CO/KS, NM/TX): always use combined shares; never
   grassland-only or the grass fraction as a feature.
4. `scripts/04_cutb_download.R` (Cut B, decided Oct 5): MOD13A2 + MYD13A2 v061 (16-day VI,
   1 km; Terra + Aqua 8 days apart) via AppEEARS area requests (profile bbox, CONUS when
   null), one task per product-year, state in `cutb/tasks.csv`. AppEEARS has no 0.05 deg
   CMG products, so each year is aggregated to the 4 km grid as it is stacked: mean of
   1 km pixels with pixel reliability 0-1, plus `n_valid` and `snow_frac` per composite.
   Output `cutb/<product>/<var>/<product>_<var>_<year>.tif` on the PRISM grid, physical
   units. Layer names resolved from AppEEARS at run time (`config.yml` `cutb:`). Earthdata
   login from a netrc file only (`R/appeears.R`). Indices: NDVI (greenness) + NDII7 (NIR vs
   2.1 um MIR, curing). Upgrade path: MCD43A4 (500 m daily NBAR, in AppEEARS).
   Scoping notes: project doc `claude/gsi-cutb-scope.md`.
4b. `scripts/06_cutb_curves.R` (Cut B curves, built Oct 7; code `R/cutb_curves.R`, config
   `cutb_curves:`). Pass 1 per year: Terra - Aqua offset per cell (Terra interpolated to Aqua
   composite midpoints, median of valid pairs) -> `cutb/offsets/`, median across years. Pass 2
   per year Y: composites Oct (Y-1) - Mar (Y+1), year-boundary duplicates removed (each
   AppEEARS year file repeats the previous December composite), n_valid >= 5 and
   snow_frac <= 0.5, offset split between sensors, observation date from the composite
   `doy` layer (falls back to the midpoint when the 1 km doy average straddles the year
   boundary), weighted Whittaker (lambda 30 on a 4-day grid, one upper-envelope pass;
   vectorized banded Cholesky `whit_fit()`). Outputs `cutb/curves/{ndvi,ndii}_<Y>.tif` (DOY 1,
   9, ..., 361), `cutb/metrics/metrics_<Y>.tif` (18 bands: NDVI max/base/amp/peak, sos20/50,
   eos50/20, gsl20, n_green, NDII7 max/base/amp/peak, cure50, cure_days, n_obs, max_gap; DOY
   relative to Jan 1 of Y), summaries `cutb_metrics_*`, `cutb_{ndvi,ndii}_curve_*`. Verified on
   dev against an independent per-cell rebuild (dense solve): agreement to float32 precision.
   Dev offsets: NDII7 Terra - Aqua ~0.033 (5-95 %: 0.020-0.048); NDVI ~0. NDII7 curing dates
   can land in winter where NDII7 keeps falling with soil moisture (e.g. SE Arizona).
   `summarise_years()` now uses vectorized `row_quantiles()` (identical results, ~36x faster).
5. `scripts/10_block1_seasonal.R` (Block 1, decided Oct 6): per year label Y, water year
   Oct(Y-1)-Sep(Y) for temperature / precipitation / aridity, calendar year Y for freezes
   and GDD; seasons OND/JFM/AMJ/JAS; 21 bands (seasonal + annual T, T range, annual P,
   seasonal P fractions, P/PET annual and AMJ+JAS with Hargreaves PET, last spring / first
   fall freeze and freeze-free days for tmin <= 0 and <= -2.2 C, GDD base 5 and 10).
   `features/block1/block1_<Y>.tif`, then `summaries/block1_median.tif` and `_iqr.tif`.
   Shared helpers in `R/features.R` (row-chunked reads, Hargreaves, across-year summary);
   block code in `R/block1.R`. Freeze definition tmin <= 0 C (-2.2 C sensitivity) is also the
   Block 2 default. Verified against an independent per-cell calculation.
   Blocks 2-3: one script each (11_, 12_), same pattern.
5b. `scripts/11_block2_timing.R` (Block 2, decided Oct 7; code in `R/block2.R`, config
   `features: block2`): rain event >= 5 mm/day (2.5, 10 mm sensitivity); dry day < 1 mm;
   warm-season onset = first wet day from May 1 starting a 10-day total >= 20 mm with no
   dry spell >= 20 days in the next 30 (none by Sep 30 -> DOY 274, censored); hard freeze
   tmin <= -2.2 C. 12 per-year bands: onset_doy, ev_warm/ev_cool at each threshold,
   dsl_max_warm, n_dry20_warm, gdd5_lhf (GDD by the last spring hard freeze), cold_frac,
   dry_frac (not cold and 30-day P < 0.5 x 30-day Hargreaves PET). Summaries
   `block2_median`, `block2_iqr` and `block2_var` (cv_p_ann, cv_p_warm, sd_t_amj from the
   Block 1 per-year files; onset_frac). Verified against an independent per-cell
   calculation on the dev box (exact). Gotcha: `format()` on a vector pads ("5.0"); use
   `as.character()` per element for band-name tags.
5c. `scripts/12_gsi.R` (Block 3 / GSI, decided Oct 8; code `R/gsi.R`, config `features: gsi`):
   daily GSI per year Y on the grid, series Jan 1 (Y-1) .. Mar 31 (Y+1) (Dec 31 for the last
   year). Variants = param set {nfdrs: VPD 900-4100 Pa, daylength 10-11 h, precip 0-10 mm/28 d,
   21-d smoothing, green-up 0.5; fems: 1956-3882 Pa, 11-12 h, 10.16-20.32 mm/28 d, 28 d, 0.3}
   x moisture {none, precip, kbdi (ramp 200-600 decreasing, MAP from 2001-2025 ->
   static/kbdi_map.tif)} + no-VPD runs (nfdrs/fems with precip). Tmin -2..5 C, geometric
   daylength, gsi_max 1, persistence 3 d. 13 bands per year and variant (gsi_mean/max,
   peak_doy, sos20, eos50 (relative crossings as Cut B), gu_doy, dorm_doy, green_days,
   n_pulses, lim_tmin/vpd/photo/moist) -> `features/gsi_<variant>/`, summaries
   `gsi_<variant>_{median,iqr}`. `R/gsi.R` Part 1 = gsi-scout's point model copied (only
   change: compute_kbdi() takes R); Part 2 = grid version. `scripts/12b_gsi_check.R` checks
   Part 2 against run_gsi() at the check sites (dev Oct 8: daily GSI equal to ~1e-15, PASS).
   Gotchas: pmin/pmax(scalar, matrix) drop dim (use m_ramp); max.col() has a tolerance, and
   GSI plateaus at 1, so the peak is the first day within 1e-6 of the maximum.
6. `scripts/20_cluster_a1.R` (A1 = Block 1 + static; built and run on CONUS Oct 7):
   `Rscript scripts/20_cluster_a1.R <profile> [variant,...]` (RStudio: `gsi_profile`,
   `gsi_variant`). Variants in config `clustering: a1: variants`: `base` (one feature group;
   output `clusters/a1`), `groups` (Block 1 split into temperature / moisture amount /
   seasonality groups with equal variance) and `groups_herbw` (same, cells weighted
   0.1 + 0.9 x open_herb_share). Outputs `clusters/a1_<variant>/`: k-means for k = 4-30,
   GMM + two-stage Ward at 3 k values (stability-picked or config `k_detail`), spatial-block
   stability, agreement (incl. vs base run), per-cluster profiles, `check_points.csv`
   (region of 11 check sites per k), maps. Helpers in `R/cluster.R` (block / group PCA,
   weighted k-means, ARI, label matching). ~22 min per variant on the P720.
   **Working A1 = `groups` at k = 10 (decided Oct 7).** Lessons: weight feature groups equally
   within a block (one group = temperature-dominated); tallgrass prairie and eastern forest
   don't separate on Block 1 medians (needs Block 2); mclust VVV is overconfident (retune).
   `mclust::Mclust` needs `library(mclust)` attached.
7. `scripts/07_cutb_regions.R` (Cut B regions, Oct 7; config `cutb_regions:`): PCA of the
   across-year median NDVI + NDII7 curves (46 values each; on regular smooth curves this is
   functional PCA), each index total variance 1; variants `raw` (level + shape) and `shape`
   (each curve centred and divided by its amplitude: timing only). Spatial-block stability,
   k-means k = 4-20, two-stage Ward, ARI vs Cut A for every (B k, A k) on all cells and the
   strict herb mask, cluster-mean curves. Outputs `clusters/b_<variant>/` (`b_kmeans.tif`,
   `b_features.tif`, `ari_vs_A.csv`, `curves_k<K>.csv`, figs). ~30-50 min on CONUS.
8. Evaluation (config `evaluation:`; helpers `R/evaluate.R`): partitions = k-means label
   rasters under `clusters/` (A1_groups, A2_groups, B_raw, B_shape) + `geo` baseline (k-means
   on location only, cached in `clusters/geo/`).
   - `scripts/30_compare_ab.R`: share of the Cut B curve variance (and of single metrics:
     amplitude, green-up, peak, end, season length, curing) explained by each partition per k,
     all cells and strict herb mask -> `eval/ab/`. Cut B = upper reference (fitted to the
     curves); geo = what any compact regions give. Two curve feature sets (config
     `curve_features`): `level_shape` (b_raw) and `timing` (b_shape: curves scaled by their
     amplitude), column `features` in r2_curves.csv. Minutes.
   - `scripts/31_obs_regions.R`: ground check from the GSI study tables
     (`GrowingSeasonIndex/study/data/derived/*.rds`; config `evaluation: study_dir` or
     `GSI_STUDY_DIR`): NPN graminoid % green (green-up / curing / cured), NPN herbaceous
     onset / end, woody onset, Globe-LFMC herb / woody peak and 50 % decline, iNaturalist
     flowering onset. Site = place x plant type (NPN site x species, LFMC series, iNat point x
     genus), median date across years, then adjusted R^2 on region per partition and k, plus
     sites per region -> `eval/obs/`. Run on the laptop (study data live there).
     Spatial-block bootstrap at report_k (config `evaluation: boot`; 2-degree blocks, 500
     resamples, ~5 min): `obs_r2_boot.csv` (range per score), `obs_r2_diff.csv` (paired
     differences, e.g. A2 - A1). Gotcha: a YAML key named `n` reads as FALSE (YAML 1.1); use
     `n_boot`.
   - `scripts/33_synchrony.R` (P720; needs per-year `cutb/metrics`): share of cells' yearly
     anomalies (value - own median; jumps > 90 days dropped) in sos20, NDVI peak, eos20 and
     cure50 explained by their region's mean anomaly that year, per partition and k (k = 1 =
     CONUS), herb and all cells; per-region numbers and per-cell correlation rasters for the
     profiled sets -> `eval/sync/`. Page figure `eval_sync`.
   - `scripts/32_region_profiles.R` (laptop, ~1 min): per candidate region set (config
     `region_profiles: sets`; A1 k10, A2 k7, A2 k13) feature medians/IQR, Cut B season
     metrics and median curves (all cells and strict herb cells), ground-observation dates
     per region, fine->coarse crosswalks -> `eval/profiles/<partition>_k<KK>/`. Region
     numbers/colours from `region_display()` in `R/evaluate.R` (shared with the report maps:
     colours matched to `report: region_maps[1]`). Page figures `prof_*` via 90 (`profiles`).
     Keep config.yml ASCII (a degree sign broke YAML reading in a C locale).
   - `scripts/34_gsi_model.R` (P720, ~1-1.5 h with 24 forked workers; code `R/gsi_model.R`,
     config `gsi_model:`): samples 4,000 strict-herb cells, caches their daily PRISM series and
     yearly satellite sos20 / cure50 (`eval/gsi_model/sample_<profile>.rds`), fits GSI ramp
     thresholds (Tmin, VPD, daylength start, 28-day precip; FEMS smoothing) by Nelder-Mead to
     minimize the capped (90 d) mean absolute error of GSI sos20 / eos50 vs satellite sos20 /
     cure50: once for CONUS and per region for A1 / A2 / A2 herb-weighted at k 7/10/13; two
     folds (odd / even years). Compared with FEMS / NFDRS defaults and cell climatology on
     held-out years: MAE, median error, share within 30 d, anomaly correlation ->
     `scores.csv`, `region_scores.csv`, `params.csv`. Uses forked workers (parallel::mclapply)
     on Linux so the cached series are shared; sequential on Windows. Page figure `gsi_model`.
     v2 (Oct 8 evening): moisture options precip / kbdi / none fitted separately ("best" =
     lowest fitting-year loss per region and fold), maxit 250 + 1 restart, k 7 and 13, 36
     workers; the cache is extended with tmax and map_mm (no re-sampling). First run (v1,
     precip only, maxit 150, k 7/10/13): regional 28-31 d green-up vs climatology 27; k flat.
   - `scripts/35_gsi_curing.R` (P720, ~4-5 h with 24 forked workers, run alone; config `gsi_curing:`,
     unset keys from `gsi_model:`; needs 34's sample cache): GSI curing forms = moisture
     {kbdi, kbdi30 (30-day mean), precip60, precip90} x crossing {rel (20/50 % of amplitude,
     as 34), rel_split (fractions fitted), abs_nfdrs4 (absolute: green-up at GU, 50 % cured at
     (1+GU)/2 as NFDRS4 Cure()), abs_split (levels fitted)}, CONUS + A2 k7/k13, same folds and
     loss; "best" = lowest fitting-year loss per region and fold -> `eval/gsi_curing/`
     (scores, region_scores, params, best_form). Page figure `gsi_curing` (+ tables), shown on
     the page only when the full version exists. `R/gsi_model.R` gm_moist / gm_cross /
     gm_season_dates / gm_predict; cross = "rel" reproduces gsi_season_from_G exactly (checked).
     Memory (Oct 9): running it next to 20 with 36 workers and 32 CONUS fits at once ran the
     P720 out of memory and killed the tmux server; workers now copy only the 4 driver series
     their form needs and at most `max_big` (8) CONUS fits run at once.
     v1 results (Oct 10, eval/gsi_curing_v1 on the P720 after the v2 run): best regional curing
     ~48 d = climatology, but the fixed 60-day penalty let fits skip curing dates (coverage
     3-4 % in monsoon SW / mid-South). Honest forms: GSI beats the typical curing date in the
     northern Plains (16 vs 20 d, r 0.71), southern Great Basin, Snake River Plain; fails in the
     monsoon SW and California annual grasslands. 99 % of fits hit maxit 250.
     v2 (Oct 10): `missing: fallback` = a missing GSI date takes the cell's typical date from
     the fitting years (fits and scores); paired_* scores compare on the cell-years the GSI
     dates; maxit 600; moisture kbdi + precip90 only; each fit logs "fit i/n done".
     v2 results (Oct 10, ~5 h; no fit hit 1200 evaluations): with the fallback, best fits
     switch curing dates off (coverage < 1 %) in regions 2, 3, 4, 5, 11 at A2 k13, so national
     curing MAE 47.8 vs climatology 48.3 is mostly fallback; paired (years the GSI dates) 40 vs
     41. GSI beats the typical curing date in 10 (16 vs 20, r 0.71), 12 (20 vs 23), 6 (39 vs
     45). Green-up: precip90 forms 25 vs 27 d (paired 23 vs 25, r ~0.4). KBDI wins curing in
     25 of 34 region-folds. eval/gsi_curing_v1 keeps the first run (laptop and P720).
     Curing target option (Oct 10): `GSI_CURING_TARGET=eos50` (or config `target`) scores curing
     against NDVI eos50 instead of NDII7 cure50 (36 found eos matches ground curing with fewer
     season mismatches); the sample cache is extended with the yearly target once (needs
     cutb/metrics on the P720); config `target_forms` limits the forms (4 best at cure50);
     outputs eval/gsi_curing_eos50/; page figure `gsi_curing_eos50`.
     eos50 results (Oct 10, 144 fits, ~2.5 h): climatology 29.3 d (vs 48 for cure50); regional
     best 29.2-29.3, paired 30 vs 30; curing anomaly r 0.36-0.37 (vs 0.15-0.21). Regions: GSI
     beats typical in 10 (15 vs 17, r 0.72) and 6 (33 vs 44); tie 12; switched off or worse in
     2, 3, 4, 5, 11 (CA 37 vs 24). Failures are the GSI's form, not the target.
   - `scripts/36_curing_target_check.R` (laptop, < 1 min; config `target_check:`; ground events
     from `R/obs.R` obs_events(), same sets as 31): satellite sos20 / cure50 / eos50 / eos20 vs
     NPN grass green-up / curing / cured and Globe-LFMC herb decline, site medians and (with
     cutb/metrics/metrics_<Y>.tif copied from the P720 to the laptop) site-years; circular
     differences, |d| > 90 = season mismatch; anomaly correlation -> `eval/target_check/`.
     Oct 10 results (site-years): grass curing vs cure50 |d| 58, 26 % mismatch, anomaly r 0.73
     (0.83 without mismatches); vs eos50 |d| 50, 21 %, 0.70; cured vs eos20 |d| 43, 14 %, 0.76;
     green-up vs sos20 |d| 36, 7 %, 0.32; LFMC herb decline not comparable (71 % mismatch).
     Target good in CA annual (86 % within 30 d) and mid-South (72 %); monsoon SW 34 % mismatch
     (eos50 18 %); Snake River (6 sites) satellite ~4 months later than ground. NDVI end of
     season is at least as good a curing target as NDII7, with fewer season mismatches.
9. A3 (Block 3, GSI features; built Oct 9): `20_cluster_a1.R` variants `a3_groups` (A2 groups +
   Block 3 from GSI variant fems_precip, equal weight), `a3_groups_novpd` (fems_precip_novpd)
   and `a3_groups_kbdi` (fems_kbdi); config `block3: gsi_variant`, sources g3med / g3iqr
   (summaries/gsi_<variant>_{median,iqr}); groups timing (sos20, peak_doy, eos50), spread
   (IQR sos20, eos50), activity (gsi_mean, green_days, log n_pulses), limiting (lim_*).
   All-NA / constant features are dropped (lim_vpd without VPD). Cells with no GSI season
   (e.g. KBDI in the low desert, ~8 % of the dev box) get median dates; the activity group
   marks them. Evaluation partitions A3_groups / A3_novpd / A3_kbdi; boot pairs vs A2 / A3.
   Dev (Oct 9): stability lower than A2 (~0.55 at k 7-13). ~25-30 min per variant on the P720.
   **Working regions (decided Oct 9): A2 at k 13, nested into A2 k 7** for fitting where data
   are thin; A2 herb-weighted k 13 as western sensitivity. A1 k 10 was the Oct 7 working set.
   Note: `clusters/a2_groups` was made in the cloud workspace and lives on the laptop only;
   copy it to the P720 before running 30/33 there or A2 is silently skipped.

## Project page (report/)
`report/index.qmd` is a living Quarto page for the project team and Mike's website: what
the project is, decisions, open questions, status, changelog and preliminary figures. It
renders to one self-contained HTML (`quarto render report/index.qmd`; HTML is git-ignored).
Figures come from `scripts/90_report_figs.R [profile]`: dev writes small PNGs to
`report/figs/<figure>_dev.png`; full writes to `<run_dir>/report_figs/` on the P720 (config
`report: figs_in_run_dir`, keeps the P720 checkout clean) and 91 copies them to
`report/figs/`. PNGs + `manifest.csv` are committed from the laptop only; the page shows the
`full` version of a figure when present, else `dev`. Each figure skips if its inputs are
missing. Sample points for Cut B curves: config `report: sample_points`. Optional figure names
after the profile run only those (`Rscript scripts/90_report_figs.R full regions eval_obs`).
Region and evaluation figures (`regions`, `eval_curves`, `eval_obs`, `obs_sites`, plus
`eval_summary_<profile>.csv` shown as a table) are made on the laptop with the full profile
(`figs_in_run_dir` is `{windows: false, linux: true}`, so they go straight to report/figs),
because 31's output exists only there; 91 keeps laptop-made full rows when merging the
manifest. Workflow from Oct 7: Mike runs the analyses on his machines; agents update the
page (text + figure code) and read results from the laptop. When a decision
changes, also update the page's Decisions / Open questions / Status / Changelog sections.
The page is public-facing: no machine names, paths, credentials or ntfy topics.
`report/_worked_examples.qmd` (included in index.qmd as "How the numbers work", added Oct 10)
walks through each step on toy data with visible code, sourcing the pipeline's own functions
from `R/` (features, cluster, evaluate, gsi, gsi_model, cutb_curves; `sync_r2` is copied from
33). Numbers from random toy data are inline R, so the text follows the output. If a metric's
definition changes in `R/`, check the matching tab.
`scripts/91_sync_from_p720.R` (laptop) pulls `*_full.png` + merged manifest (and optionally
summaries/static rasters, logs) by scp; host from `GSI_P720_HOST`, remote paths in config `sync:`. It also renders the page
(`render`) and uploads `index.html` to S3 with the AWS CLI (`publish`, only when asked;
`GSI_S3_DEST`, optional `GSI_AWS_PROFILE`, `GSI_CF_DIST`). Never handle AWS credentials.
