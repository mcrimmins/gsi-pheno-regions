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
     New model forms (Oct 11; config `form_sets: newforms`, env GSI_CURING_FORMS=newforms ->
     eval/gsi_curing[_eos50]_newforms/, page figures gsi_curing[_eos50]_newforms):
     swb = soil-water bucket (`R/gsi_model.R` m_bucket; W += P - PET x W/wmax, Hargreaves PET
     from the cached tmin/tmax and latitude; index = ramp(relative fill) x ramp(Tmin), no VPD or
     daylength; fitted tmin_lo/w, wmax 10-400 mm, fill ramp; C++ via Rcpp when it compiles,
     identical R loop otherwise); rel_anchor = rel_split with the peak searched within +/-
     anchor_days (75) of the cell's usual mid-season from the fitting years (gm_anchor: midway
     between median green-up and curing), for two-season climates. Dev with 2 years can't judge
     rel_anchor (the "usual" season is one year).
     newforms results (Oct 11): bucket wins ~3/4 of region-folds; wmax median 180-205 mm. eos50
     target, k13: best 28.2 vs typical 29.3 (k7 27.7), paired 27.5 vs 29.0, r 0.36-0.44;
     cure50: 47.5 vs 48.3. Green-up 24.8-25.5 vs 27.2 (r ~0.5); gains in S Plains (30 vs 40)
     and low deserts (27 vs 35). Curing still switched off / worse in regions 2, 3, 5, 11;
     rel_anchor doesn't rescue the monsoon SW.
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
     Oct 11-12 fix: anomalies unwrap ground dates around the site median (cdiff) and leave out
     season mismatches (|d| > 90 put the satellite on the far side of the circle and made large
     shared departures). Revised anom r: curing vs cure50 0.77, vs eos50 0.70; cured 0.67 /
     0.66; green-up 0.49; mid-South (5) ~0.1 (ground varies little); Snake River 0.94-0.98;
     herb-mask cells 0.53-0.55 for curing. The Oct 10 anom r numbers above are superseded.
   - **Curing baseline decision (Oct 11):** the typical (median) curing date is the curing
     baseline in regions 3 (monsoon SW), 5 (mid-South), 2 (S Plains), 11 (CA annual), where no
     model beat it; GSI / bucket where they beat it (10, 6); bucket for green-up.
   - `scripts/37_curing_cues.R` (P720, minutes; config `curing_cues:`; `R/cues.R`): screens
     weather cues for early / late curing years in the hard regions before building
     region-specific rules. Uses sample_<profile>.rds (needs eos50, added by the eos50 run of
     35). Cues per cell-year (cue_year: onset, rain_end, spring_rain, fall_frz0/2, tmax_spring /
     summer, p_cool / spring / summer / fall, dry_spell; plus green-up departures for curing
     targets), cached to eval/curing_cues/cues_<profile>.rds. Targets eos50, cure50, sos20 as
     cell-demeaned departures, season jumps (> 90 d from the cell median) dropped. Per A2 k13
     region: pooled r per cue and cv_r2 of a forward-selected <= 3-cue linear model (odd / even
     years) -> cue_r.csv, cue_models.csv. Page figure `curing_cues` (heatmap of r, top row
     cv_r2; CSV table).
     Oct 12 results (eos50 target): cv_r2 CONUS 0.10, regions 0.08-0.33; CA annual (11) 0.33
     (tmax_spring r -0.53, p_cool 0.37, spring_rain 0.35); low deserts (4) 0.29; mid-South (5)
     0.25 (63 cells; target r ~0.1 there); Great Plains (10) 0.25 (GSI r 0.72 is better); S
     Plains (2) 0.20, mostly greenup; monsoon SW (3) 0.13. cure50 lower (CONUS 0.06). sos20
     CONUS 0.05. Regions 1, 8, 9 too few cells. Reading: a rule is worth trying only in CA
     annual grasslands (open; Mike to decide).
   - **REGION NUMBERING CORRECTION (Oct 12).** 34 and 35 label regions with the clustering's
     internal labels (eval_partitions); the page, 32, 33, 36, 37, 38 use the map numbers
     (region_display). The 35 region results above were read as map numbers and are WRONG
     in their region names. A2 k13 internal -> map (from cell counts; 35 and 38 typical-date
     MAE agree exactly): 1->4, 2->11, 3->12, 4->3, 5->6, 6->2, 9->10, 10->5, 11->7, 12->13
     (7, 8: < 25 cells). Corrected: GSI / bucket beat the typical curing date in 2 (S
     Plains) and 5 (mid-South), barely 3 (monsoon SW); nothing beats it in 11, 6, 12, 7, 10,
     13, 4. Bucket green-up gains are in 11 (CA) and 3 (monsoon SW). 35 now writes a
     `region_page` column (from region_display) next to `region`; any table of 34/35 region
     results must use region_page.
   - `scripts/38_curing_adjust.R` (P720, minutes; config `curing_adjust:`; R/cues.R adj_cues /
     adj_fit). Decision Oct 12: one general rule, typical date (cell median over fitting years)
     + linear adjustment of <= 3 lead-time cue departures (no intercept, clamp +/- 60 d), fitted
     per A2 k13 region on odd / even years, scored on the others by capped MAE as in 35; kept
     only where it beats the typical date. Cues end at typical date - lead (0 / 14 / 30 d; 14 =
     operational): p30 / p90 / p180, tx30 / tx60, dry60, since_rain, greenup (if before the
     cutoff, else typical). Outputs eval/curing_adjust/scores.csv, coefs.csv, cue caches. Page
     table `curing_adjust` (90 writes curing_adjust[_ndii]_<profile>.csv). Synthetic check:
     recovers a planted -5 d/C heat effect (-4.6 to -4.9). Last curing experiment before the
     write-up. Worked-examples tab 16.
     Oct 12 results (eos50, lead 14): typical 29.3 -> adjusted 27.2 nationally (lead 0 +2.3,
     30 +1.8); beats typical in every region: 2 +4.6, 10 +4.5, 4 +3.1, 13 +2.8, 11 +2.7
     (17 -> 14), 6 +1.6, 3 +1.2, 12 +1.1, 7 +0.6, 5 +0.3. Best approach in 8 of 10; bucket /
     GSI better in 2 (33 vs 39) and 5 (14 vs 17); tie in 3 (28). Cues: Tmax 60 d (-5 to -6
     d/C Great Plains, -2 to -3 CA), rain 30-90 d (+0.1-0.3 d/mm), dry days / days since rain
     (earlier), green-up (later). Proposed decision (pending Mike): typical + adjustment
     everywhere except bucket in 2 and 5; green-up bucket in 11, 3, 4, 2, typical elsewhere;
     stop curing model work and write up. DECIDED by Mike Oct 12 (yes).
   - `scripts/39_rule_ground_check.R` (P720; needs 38 coefs.csv, PRISM daily, cutb metrics, and
     npn_curing_metrics.rds via GSI_OBS_DIR): applies the eos50 lead-14 rule at NPN grass site
     cells (typical = cell median eos50 over fitting years; 38 region+fold coefs, CONUS where
     none) and compares held-out site-years with ground curing (g50 down) and cured dates:
     capped MAE without season mismatches, offset-removed error (site-demeaned ground vs
     site-demeaned adjustment, sites with >= 3 non-mismatch years), r. Outputs
     eval/curing_adjust/ground_site_years.csv, ground_scores.csv; page table
     `curing_adjust_ground` (90). The bucket rule (regions 2, 5) is not ground-checked.
     Oct 12 first run (plant rows): adjustment does not track ground curing (offset-removed
     error 25.5 -> 27.6 d, r 0.03) except CA (11: r 0.59, 5.8 -> 4.4) and monsoon SW (3:
     r 0.37). Check: adjustment vs satellite departures at the same site-years r 0.35 (code
     OK), but satellite vs ground at the site-year level only ~0.2-0.46 (0.70 on the page was
     plant rows, several plants per site-year, and conditioned on |sat - ground| <= 90).
   - **Plant rows fix (Oct 12 evening):** npn_curing_metrics has one row per PLANT; 36 and 39
     counted plant-years as site-years (4,928 rows = 645 site-years for curing vs eos50). R/obs.R
     obs_site_years() collapses to one row per set x site x year (circular median over
     plants); 36 and 39 call it after obs_events. Rerun both; target-check tables on the page
     must be refreshed from the rerun.
     Oct 13: also fixed 36's satellite departure (was ground + diff, which copied ground
     outliers into the satellite series and gave r 0.70): now unwrapped around the satellite's
     own site median. Corrected (computed from the Oct 13 site_years; rerun 36 to refresh its
     CSVs): curing vs eos50 r 0.35 (herb cells 0.56), vs cure50 0.27 (0.55), green-up 0.48;
     regions eos50: 2 0.62, 3 0.65, 10 0.55, 5 0.21, 11 0.48, 13 0.31, 1 -0.14.
     39 rerun (site-years): adjustment vs ground, offset removed, CONUS 19.5 -> 21.9 d, r 0.08;
     monsoon SW 21.5 -> 17.8 (r 0.43), CA 5.3 -> 4.2 (r 0.66, 10 site-years); worse in 10, 2,
     13, 5. Proposed amendment (pending Mike): typical date everywhere, adjustment only in
     monsoon SW and CA; bucket for 2 / 5 waits on a ground check. DECIDED by Mike Oct 13 (yes): typical date
     everywhere, adjustment in monsoon SW and (provisionally) CA; bucket rules kept only if the
     ground check supports them; then the write-up. Further curing rules should be fitted to
     ground data (the NPN curing study), not the satellite.
     39 extended (Oct 13): also the bucket rule (best swb form per internal region and fold from
     gsi_curing_eos50_newforms/params.csv, rebuilt with gm_predict at the site cells; PET via
     hargreaves; typical date where it gives none) for curing and green-up (grass_greenup);
     ground_scores.csv now has `rule` (adjustment / bucket) and mae_rule / dev_rule / r_rule.
     90 `curing_adjust_ground` writes curing_adjust_ground[_cured], curing_bucket_ground,
     greenup_bucket_ground CSVs.
     Oct 13 bucket ground results: worse than typical everywhere it matters (curing S Plains
     23 -> 43 d offset-removed, mid-South 8.6 -> 11.9, CONUS 19.5 -> 28.3; green-up S Plains
     21 -> 30, monsoon SW 23 -> 43; no CA green-up sites; Great Plains green-up r 0.42 but ~no
     gain). FINAL RULES (per the Oct 13 decision): typical date everywhere; curing adjustment
     (lead 14) in monsoon SW (3) and provisionally CA (11). Next: the write-up.
     Write-up drafted Oct 13 as a Docs artifact ("Herbaceous green-up and curing rules by
     region: findings and recommendations"; with a satellite-vs-ground gain chart). Page:
     "At a glance" has a findings summary; "How the analysis is built" has the four-step
     rule-testing method; worked-examples tab 17 = ground check of a rule (offset removed).
   - **Rules package (Oct 15):** `scripts/40_rules_build.R` (P720; config `rules:`) writes
     <run_dir>/rules/: typical_dates.tif (greenup_doy, curing_doy = cutb_metrics median sos20 /
     eos50; *_spread = IQR; herb; region = A2 k13 map number; adjust = herb cells in fitted
     adjust_regions [3, 11]), adjust_coefs.csv (refit on ALL sample years, adj_fit), 
     adjust_cue_means.tif (per-cell cue climatology at the cell's own cutoff; chunked, forked),
     README.txt. `scripts/41_rules_predict.R <year> [as-of] [profile]` writes predict_<Y>.tif
     (greenup_typical, curing_typical, curing_adjustment, curing_expected, adjusted, region)
     and, for past years, predict_<Y>_check.csv (in-sample vs satellite eos50). R/cues.R
     adj_apply(). 90 figures rules_typical (2x2 maps) and rules_adjust (example_years, west
     crop; needs predict files). Dev checked with regions 4/9 (dev regions 3/11 too small).
     Oct 15 full run: region 3 cues p30 (+0.26 d/mm), greenup (+0.16), p180 (-0.05), R2 0.11;
     region 11 p90 (+0.05), since_rain (-0.13), tx60 (-2.6 d/C), R2 0.33; 34,041 adjusted
     cells; 2021/2023 in-sample checks CA 23.1->19.7, 14.4->9.6. Wall to wall (Oct 15):
     typical_dates.tif gains `confidence` (1 herb in tested_regions [2,3,5,10,11,13], 2 other
     herb, 3 mixed cover); maps show every land cell; 90 figure rules_confidence.
   - **Current-year run (Oct 15):** `scripts/01_prism_daily.R [profile] [years] [vars]` (e.g.
     `full 2026 ppt,tmax`); a partial current year is fetched to yesterday, rebuilt when newer
     days exist, and keeps its raw days (drop_raw skipped) so updates download only new days.
     41 writes predict_<Y>_info.csv (as_of, weather_to, cells adjusted / pending). 90 figure
     rules_current (env GSI_RULES_YEAR, default this year; pending cells dark grey);
     rules_adjust / rules_current share draw_adjust().
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
from `R/` (features, cluster, evaluate, gsi, gsi_model, cutb_curves, block2, cues; `sync_r2`
is copied from 33). Numbers are inline R, so the text follows the output. Tabs 1-11: features
to GSI model test and profiles; 12 target check (36: circular differences, mismatch, unwrapped
anomaly r); 13 soil-water bucket; 14 season-aware reading (rel_anchor); 15 cue screening (37).
If a metric's definition changes in `R/` or a script, check the matching tab.
`scripts/91_sync_from_p720.R` (laptop) pulls `*_full.png` + merged manifest (and optionally
summaries/static rasters, logs) by scp; host from `GSI_P720_HOST`, remote paths in config `sync:`. It also renders the page
(`render`, both pages: index.qmd and npn_study.qmd) and uploads both HTML files to one S3 folder with the AWS CLI (`publish`, only when asked;
`GSI_S3_DEST`, optional `GSI_AWS_PROFILE`, `GSI_CF_DIST`). Never handle AWS credentials.

## NPN grass curing study (started Oct 15 2026) — separate companion study
**Boundary (Oct 16, Mike):** a separate, complementary study to the regions work, with its own
chat thread, page and paper. Code stays in this repo (scripts 50-59, R/npn.R) because it reuses
the PRISM, cue and region code and the map region numbers. Its page is report/npn_study.qmd
(-> npn_study.html, rendered and published next to index.html by 91); index.qmd keeps only a
short summary under "Ground dates: the NPN curing study" with a link. Put NPN methods, figures
and results on npn_study.qmd, not index.qmd. Handoff: project doc claude/gsi-npn-study-handoff.md.
Region-work changes that affect NPN (region numbers, typical-date product) should be noted on
both pages.
Scope: project doc claude/gsi-npn-curing-study-scope.md (Q1-Q6, model ladder L0-L4). Decided
Oct 15 (Mike): code lives here (scripts 50+); likely managed sites flagged and reported both
ways (needs site names, not in the derived tables yet); C3/C4 drafted in
`data/npn_species_pathway.csv` (pathway, life cycle, origin, crop) for Mike / Theresa to review;
config `npn_study:` defaults (min 12 readings, P-PET memory 15/30/60 d, day-length gates 12/11 h).
Record (Oct 15): 6,420 usable plant-years 2014-2024, 1,706 plants, 132 sites, 72 species (C4
4,265 plant-years, C3 2,155); g50 down 5,165, cured 3,607, g50 up 2,761; median 34 readings;
119 site x species with curing in >= 3 years. Strong in the West and Plains (FL, ND, KS, CO, UT,
OK, TX, NM, CA, AZ), not mostly East as the scope assumed.
Lessons carried from the rules work: NPN rows are plants (aggregate or model plants within
site-years explicitly); score year-to-year skill on departures from each site's own mean;
satellite dates are not a target here.
- `scripts/50_npn_prism.R` (P720): daily PRISM (ppt, tmin, tmax, vpdmax; 1990 to the last
  year) at every NPN grass site cell -> <run_dir>/npn_study/npn_site_daily.rds (sites x days),
  so metrics and the ladder (51+) run anywhere. PRISM daily on the P720 starts in 2001, so the
  study years are 2001-2025 and normals 2002-2025.
- `scripts/51_npn_metrics.R` (anywhere, seconds; GSI_OBS_DIR): series-years (site x species x
  year medians over plants, >= 12 readings) with pathway and region (from target_check
  site_years.csv); site-year metrics (thermal, freeze, seasonal rain, VPD, dry spells, onset,
  P-PET memory on fixed dates); site normals; green-up-anchored rain windows. Helpers R/npn.R.
  -> npn_study/npn_metrics.rds, npn_series_years.csv. Oct 16: 707 series-years, 181 series, 86 sites.
- `scripts/52_npn_ladder.R` (anywhere, < 1 min): L0-L4 with leave-one-year-out CV on
  series-years (crops out, >= 3 years, |dep| > 90 d dropped as mismatch); "lead" driver set
  (windows end by typical date - 14 d) vs "hindcast" (adds summer/fall windows). L2 curing =
  held downward crossing of the 30-d P-PET memory within typical +/- 60 d, or first fall
  freeze. -> npn_study/ladder/{predictions,scores,selected,importance}.csv.
  Result Oct 16 (lead set): no level beats the typical date. MAE L0 green-up 24.8 d, curing
  27.6, cured 23.7; best skill +0.01 (L1 cured, r 0.15); L2-L4 negative; same by pathway,
  by region and on site-year means. Species at one site-year share only 22-40 % of their
  departure variance (ICC: green-up 0.22, curing 0.24, cured 0.40), so most year-to-year
  scatter in NPN dates is plant / observer noise, which caps any weather model.
- `scripts/53_npn_variance.R` (anywhere, seconds; lme4): variance partition of plant-year dates,
  y ~ (1|site)+(1|species)+(1|site:species)+(1|year)+(1|region:year)+(1|site:year)+(1|site:species:year),
  all / C3 / C4 -> npn_study/variance/{components,summary}.csv. Oct 16: place 68-80 % of all
  variance; of year-to-year, shared by a site's grasses 39 / 26 / 38 % (green-up / curing /
  cured; ~18 d SD each), region-wide 28 / 10 / 3 %, single plant 39 / 53 / 39 % (C4 curing 62 %);
  reading gaps explain 1-3 % of the plant part. So a shared signal exists that the 52 drivers
  miss; curing's is site-local (weather, or management / observer change at the site).
- `scripts/54_npn_site_signal.R` (anywhere; lme4): BLUP shared signal per site-year (year +
  region:year + site:year) vs weather departures; LOYO forward selection; managed-site refit
  when data/npn_site_names.csv exists (`scripts/54a_npn_site_names.R`, laptop, rnpn; names kept
  out of git and off the page). -> npn_study/site_signal/. Oct 16: curing / cured signal
  follows moisture in the expected direction (within-year r 0.2-0.27: vpd_jja -, p_amj +,
  events_amjjas +, vpd30_lead -) but LOYO R2 <= 0.04 single driver, <= 0 selected. Green-up
  links are with summer temperature (after the event): year confounding.
  Site names (Oct 16, 54a): 96 % of usable plant-years (5,547 of 5,780) are NEON sites (site
  names "XXXX_NNN.phenology.phe", 36 sites; median 36 readings, 5 plants per species-year);
  managed sites ~4 %, too few to fit; non-NEON side 124-154 plant-years. 54 refits by flag
  (managed, neon; side needs >= 8 sites and 100 plant-years). Plant share is NEON technician
  plots, not volunteer noise. Next: PhenoCam at the NEON sites as an observer-free check
  of the site-year signal.
- Page figures (on npn_study.qmd; 90_report_figs.R, run on the laptop after 51-54): npn_sites (NEON filled; reads
  data/npn_site_names.csv only for the NEON flag, names never drawn), npn_dates, npn_ladder,
  npn_variance, npn_signal.

## Fire activity study (designed Oct 17 2026) — separate companion study
Scope and decisions: project doc claude/gsi-fire-study-scope.md (replaces
gsi-fire-activity-design.md where they differ). Own chat thread, page report/fire_study.qmd
(-> fire_study.html, rendered and published with the other pages by 91), code here as scripts
60-69 with helpers R/fire.R and config `fire_study:`; data under <run_dir>/fire_study/ (never
in git; FIRMS as Parquet). Put fire methods, figures and results on fire_study.qmd, not
index.qmd. Decided: PRISM is the only weather source (1000-h dead fuel moisture computed from
it; KBDI and VPD30 as checks; no gridMET); FPA-FOD 7th edition (1992-2024), MTBS, VIIRS S-NPP
(NOAA-20 check); windows 2001-2024 (starts) and 2012-2024 (starts + growth); primary outcome
wildfires >= 10 ac by discovery week, lightning vs human; pilot regions 3, 6, 12, 11, 10, 13;
NLCD herbaceous cover; region x week x year, pooled count model, leave-one-year-out.
Downloads run on the P720. Raw downloads go to <data_root>/_shared/raw/fire/ (shared by
profiles). Tables via R/fire.R fire_write_table(): Parquet with arrow, .rds without it.
Results come back to the laptop with `Rscript scripts/91_sync_from_p720.R fire` (every
fire_study/<step>/*.csv, plus the tables of the steps in config `sync: fire_tables`, default fod)
into <data_root>/full/fire_study/ on the laptop. Code still goes to the P720 by git.
- `scripts/60_fire_fod.R [profile] [force]` (P720, minutes; built Oct 17): FPA-FOD 7th edition
  GeoPackage (config `fire_study: fod: url`, downloaded once) -> <run_dir>/fire_study/fod/
  fod_fires.parquet, one row per fire in the study years on a grid land cell (AK/HI/PR out):
  date, doy, cause (lightning = NWCG_CAUSE_CLASSIFICATION "Natural", human, unknown),
  size_ac, c_plus (>= 10 ac), debris (general cause "Debris and open burning"), rx_name
  (prescribed-looking name; flag only), mtbs_id, ics209_id, cell, region (A2 k13 map number
  from rules/typical_dates.tif, else region_display), herb_mask, herb/open_herb/crop share of
  the 4 km cell, pilot. Also fod_fields.csv (field check; stops if a required field is
  missing), fod_values.csv (cause and size-class counts), fod_filters.csv,
  fod_region_summary.csv, fod_region_year.csv. Discovery week is set in 65, not here.
  FPA-FOD 7th edition checked Oct 17 (metadata): field names as in the 6th edition; cause
  classes Human / Natural / Missing data/not specified/undetermined; NAD83; planned ignitions
  excluded except escapes needing suppression (Short 2014; the 7th edition metadata is
  silent). Cite as "7th Edition" (the abstract says "sixth update"). Tested on a synthetic
  GeoPackage in the cloud workspace (filters, flags, region lookup checked independently).
  Full run Oct 17 (63 s): 1,941,893 fires on the grid; 193,405 C+ wildfires (no debris), 52,924
  in pilot regions. Found: (1) reporting change 2020-21 (all INTERAGCY from 2021): all-cause
  totals continuous, but unknown cause rises to 10-30 % in pilot regions and 50-86 % in CA (11);
  unknown fires' seasonal profile matches human (r 0.90-1.00); lightning all-size counts drop
  ~40-47 % in 3, 6, 12, 11 after 2019. (2) Region 13: 5 % of C+ fires on strict herb cells, but
  fire density rises ~5x from the least to the most herbaceous cells (grass + shrub + pasture +
  wetland herb), seasonal shape identical across cover terciles; fire cells 5 % wetland herb vs
  1.4 % landscape; landscape median 52 % crops (satellite dates track crops).
  DECIDED Oct 17 (Mike, config fire_study: outcome / ignition_term / era_check / fuel_classes /
  pheno_source / cover_sensitivity): primary outcome = all-cause C+ wildfires (debris removed);
  lightning vs human only 2001-2019, unknown left out, sensitivity unknown = human; region 11
  totals only; no lightning-count ignition term (region x year>=2020 indicator as a check only);
  fuel share = grass + shrub + pasture + wetland_herb; region 13 keeps all fires, phenology
  weighted by fuel share; phenology predictor in 13 and 10 from cells with fuel share >= 0.35,
  applied regionwide; sensitivity = fires in the top third of 13's cells by fuel share.
  R/fire.R fire_fuel_share(cfg): that fuel share per grid cell from static/nlcd_shares.tif.
- `scripts/61_fire_mtbs.R [profile] [force]` (P720, minutes; built Oct 17): MTBS burned area
  boundaries (v12.0 April 2025, 1984-2024; 2023-24 still being mapped per the MTBS schedule;
  config `fire_study: mtbs:`) -> fire_study/mtbs/: mtbs_fires (event_id, type, wildfire =
  Wildfire or Wildland Fire Use, ig_date, acres, area_km2 on the grid, majority region,
  herb-mask and fuel-weighted km2, linked fod_id) and mtbs_cells (fire x cell: exact covered
  fraction x cell area, region, fuel_share), plus completeness by year (share of FPA-FOD
  wildfires >= 1,000 ac whose MTBS_ID is in MTBS: flags unfinished recent years), date check
  vs FPA-FOD discovery, region x year burned area. State from Event_ID's first two letters.
  Ig_Date parsed whatever its storage (Date, text, epoch days / s / ms). Tested on synthetic
  Albers perimeters (incl. an invalid one): grid area = polygon area to 1e-7, regions right.
  v12.0 file has lower-case field names: 61 matches fields ignoring case and renames. Full run
  Oct 17 (2 min): 30,932 perimeters (types Wildfire 16,779, Prescribed Fire 9,639, Other 4,311,
  Wildland Fire Use 203; "Other" 2 % linked to FPA-FOD, mostly SE: not wildfire); 21,759 in
  2001-2024 CONUS, 10,428 wildfires; 96.5 % linked to FPA-FOD; ignition vs discovery |d| <= 1 d
  96 %; completeness (FPA-FOD wildfires >= 1,000 ac found) 0.71-0.91 through 2023, 0.62 in 2024
  (not finished). Pilot wildfires: 3 1,198 / 6 1,653 / 10 867 / 11 595 / 12 928 / 13 ~186;
  burned km2 3 48,935 / 6 92,827 / 10 32,442 / 11 38,568 / 12 44,568 / 13 2,959; fuel-weighted
  share 62-83 % (13: 50 %).
- `scripts/62_fire_firms.R [profile] [sensors] [years] [force]` (P720; built Oct 17): VIIRS 375 m
  active fire from the FIRMS area API (config `fire_study: firms:`): key ONLY from env
  FIRMS_MAP_KEY (~/.Renviron; masked in logs; a key error stops the run). Source dates from
  FIRMS data_availability (config fallback_dates if unreachable): SP where it exists, NRT after.
  5-day pieces (config window_days; FIRMS rejects more since Oct 2026; HTTP 400 stops the run)
  over the CONUS box cached as csv.gz in _shared/raw/fire/firms/<SOURCE>/<year>/
  (never re-fetched; resumes), key transaction count checked before each request (fail-soft).
  Per sensor-year table fire_study/firms/<SENSOR>/firms_<SENSOR>_<year>.parquet: SP type 0 only
  (volcano, static land sources, offshore dropped), all confidence classes kept (conf column),
  grid land cells only; date_utc, time_utc, date_solar (UTC + lon/15 h), frp, cell, region,
  fuel_share, crop_share, herb_mask, pilot. firms_summary.csv, firms_gaps.csv (months with days
  of no CONUS detections = outages), firms_requests.csv. Tables are not synced (size); CSVs are.
  Tested on synthetic pieces (cells, regions, type / confidence handling, gap detection).
  Full run Oct 17-18 (~1,520 pieces; after a gap-check crash on NOAA-20's mid-year 2018 start,
  fixed; summaries now written per year): kept 259k-593k CONUS detections / yr, pilot 57k-246k;
  static etc. 13-26 % of raw; low confidence 4-5 %; NOAA-20 / S-NPP pilot ratio 0.95-1.04 in
  2019-23, 1.22 in 2024. Outages (>= 3 days without CONUS data): S-NPP 2022-07 (5), 2022-08
  (10), 2024-06 (3), 2024-07 (11), 2024-11 (3); NOAA-20 2023-07 (4).
  DECIDED Oct 17 (Mike; config fire_study: burned_area / viirs_fill): burned area = FPA-FOD final
  size by discovery date (log / seasonal-yearly reading); MTBS = linking detections (63) and
  fuel-burned check 2001-2023 (2024 unfinished); S-NPP primary, whole CONUS days without S-NPP
  data filled from NOAA-20 (flagged; sensitivity drops filled weeks). Implemented in 65.
- `scripts/63_fire_link.R [profile] [sensors] [years] [force]` (P720, forked by sensor-year,
  workers_by_step fire_link 11; built Oct 18): flags every 62 detection (config `fire_study:
  link:`). (1) MTBS: point in perimeter buffered 750 m (EPSG:5070, terra::relate) with solar date
  in ig - 1 .. ig + 180 d, nearest ignition; all MTBS types (mtbs_type kept). (2) else FPA-FOD
  (all sizes, debris included): within 1.5 km + equivalent-circle radius of the final size, date
  in discovery - 1 .. end + 1 (end = containment date if within a year, else 90th-pct duration by
  size class: A-B 1, C 2, D 6, E 15, F 47, G 84 d); candidates from 0.1-degree buckets, nearest
  wins. wildfire = MTBS Wildfire / WFU or FPA-FOD non-debris; crop = cell crop share >= 0.5;
  keep = wildfire & !crop. -> fire_study/link/<SENSOR>/link_<SENSOR>_<Y>.parquet (62 columns +
  link, mtbs_id, mtbs_type, fod_id, dist_km, wildfire, crop, keep), link_summary.csv (sensor x
  year x region incl. all / pilot: counts per link class and screen). County-week burn-season
  screen (sensitivity) not built (needs counties). Tested with planted detections (FPA-FOD 39/40
  linked, the 40th inside an MTBS perimeter; MTBS 5/5; far / late plants unlinked). Bug caught by
  the test: bucket keys must repeat each fire's position `each` = offsets.
  Full run Oct 18 (S-NPP 2012-2025, pilot): MTBS wildfire 56 %, FPA-FOD wildfire 3 %, unlinked
  40 %, crop 13 %, kept 59 %. By region kept / unlinked / crop: 12 79/20/6, 6 78/21/2, 11 60/39/22,
  3 57/40/1, 10 32/67/16, 13 3/95/66. 2025 kept 0 % (FPA-FOD and MTBS end 2024). Low-fire years
  (2019, 2023) ~30 % kept. Without arguments on existing tables 63 only writes link_months.csv
  (pilot region x month x class kept / unlinked_noncrop / crop / other) and logs monthly shares.
  DECIDED Oct 18 (Mike; config fire_study: viirs_outcome): VIIRS growth outcome 2012-2024 only;
  none in region 13 (FPA-FOD counts / acres only there); region 10 kept as a secondary check.
  Month check (link_months.csv, S-NPP <= 2024): in 3, 6, 11, 12 unlinked non-crop detections are
  spread over the year (Jun-Sep 27-35 % vs kept 64-94 %; profile r with kept -0.47..0.33): screen
  OK. Region 10: unlinked share the summer peak (profile r 0.86, 53 % Jun-Sep) but yearly totals
  barely track kept (r 0.18). DECIDED Oct 18 (Mike): screen unchanged; 66 adds a region 10 run
  counting all non-crop detections (config viirs_unlinked_sensitivity_regions); discovery week =
  fixed 7-day blocks from Jan 1, days 365-366 in week 52 (config fire_study: week; R/fire.R
  fire_week()).
- `scripts/64_fire_dfm.R [profile] [force]` (P720, forked by row chunk, workers_by_step fire_dfm
  16; built Oct 18; equations in R/fire_dfm.R, config `fire_study: dfm`): per cell and day for
  all PRISM years, then region x day means (weights all / fuel share): dfm1000 (NFDRS 1978
  1000-h: Simard EMC, day/night EMC weighted by daylength, boundary with rain duration, 7-day
  mean boundary, 7-day lagged update fr = 1 - 0.82 exp(-0.168); equations checked Oct 18 against
  Bradshaw et al. 1984 INT-GTR-169 eqs 38-42 (eq 42's bracket sign is a typo), Fosberg et al.
  1981 and PNW-RP-336), ppet60 (60-d P - Hargreaves PET, mm), vpd30 (kPa), kbdi (m_kbdi, MAP from
  the PRISM years), degree-day swe / snow_share and yearly snow-off DOY per cell. RH from PRISM:
  ea = es(Tmax) - VPDmax; rain hours = round((P_in + 0.02) / rate) <= 8, rate 0.25 in/h where
  annual P/PET < 0.8 or Apr-Sep P/PET < 0.5 (NFDRS classes 1-2), else 0.05. Spin-up: first year
  run once first. Outputs fire_study/dfm/: region_daily (csv + parquet), region_snowoff.csv,
  rain_class.csv, dfm_climatology.csv (week quantiles), dfm_drivers.csv (weekly anomaly r of
  dfm1000 with ppet60 / vpd30 / kbdi, R2 ppet60 + vpd30, by region and season),
  snowoff_cells.tif, rain_rate.tif; chunks_<profile>/ for resume. Synthetic test: region means
  equal an independent per-cell day loop to 5e-5 % moisture.
  DECIDED Oct 18 (Mike): ppet60 is the simple baseline fitted first; dfm1000 swapped in after;
  report the simple one if dfm1000 doesn't improve held-out years (config baseline: simple).
