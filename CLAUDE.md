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
5. `scripts/10_block1_seasonal.R` (Block 1, decided Oct 6): per year label Y, water year
   Oct(Y-1)-Sep(Y) for temperature / precipitation / aridity, calendar year Y for freezes
   and GDD; seasons OND/JFM/AMJ/JAS; 21 bands (seasonal + annual T, T range, annual P,
   seasonal P fractions, P/PET annual and AMJ+JAS with Hargreaves PET, last spring / first
   fall freeze and freeze-free days for tmin <= 0 and <= -2.2 C, GDD base 5 and 10).
   `features/block1/block1_<Y>.tif`, then `summaries/block1_median.tif` and `_iqr.tif`.
   Shared helpers in `R/features.R` (row-chunked reads, Hargreaves, across-year summary);
   block code in `R/block1.R`. Freeze definition tmin <= 0 C (-2.2 C sensitivity) is also the
   Block 2 default. Verified against an independent per-cell calculation.
   Blocks 2-3: one script each (11_, 12_), same pattern. 6. Cut B curves and summaries. 7. Across-year summaries. 8. Clustering. 9. Evaluation.

## Project page (report/)
`report/index.qmd` is a living Quarto page for the project team and Mike's website: what
the project is, decisions, open questions, status, changelog and preliminary figures. It
renders to one self-contained HTML (`quarto render report/index.qmd`; HTML is git-ignored).
Figures come from `scripts/90_report_figs.R [profile]`, which writes small PNGs to
`report/figs/<figure>_<profile>.png` plus `manifest.csv` (committed); the page shows the
`full` version of a figure when present, else `dev`. Each figure skips if its inputs are
missing. Sample points for Cut B curves: config `report: sample_points`. When a decision
changes, also update the page's Decisions / Open questions / Status / Changelog sections.
The page is public-facing: no machine names, paths, credentials or ntfy topics.
`scripts/91_sync_from_p720.R` (laptop) pulls `*_full.png` + merged manifest (and optionally
summaries/static rasters, logs) by scp; host from `GSI_P720_HOST`, remote paths in config `sync:`.
