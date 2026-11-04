# npn_fire — USA-NPN plant phenology and fire activity (CONUS)

Self-contained study folder inside `gsi-pheno-regions` (like `npn_study/`). Scope and decisions:
project doc `claude/gsi-npn-fire-scope.md` (31 Oct 2026).

**Question.** Across EPA Level II ecoregions, do ground-observed plant phenophases mark when fire
seasons open and close, and where the record allows, do early or late years shift fire activity?
Where do NPN observations and fire activity line up (ecoregion x plant group x season x year)?

**Data.** USA-NPN plant status and intensity records (all plant groups incl. conifers, no animals,
2009-2024); fire outcomes from the parent pipeline's data outputs (FPA-FOD, MTBS, VIIRS / FIRMS);
PRISM ppet60 as the dryness baseline. No satellite phenology.

## Layout

- `config.yml` — own config (profiles `full`, `test`); never reads the parent's config.
- `R/nf_setup.R` (paths, config, logging), `R/nf.R` (helpers).
- `scripts/NN_nf_*.R` — run from the repo root: `Rscript npn_fire/scripts/<script> [profile] ...`
- `data/` — small tables that go in git (none yet). `logs/`, `report/*.html` stay out of git.
- Data: `<data_root>/<run_name>/npn_fire/`; raw NPN downloads `<data_root>/_shared/raw/npn/`.
  `data_root` = env `NF_DATA_ROOT`, else `GSI_DATA_ROOT`, else config.

## Scripts

| Script | Runs on | Does | Status |
|---|---|---|---|
| 01_nf_npn_download.R | P720 (internet) | NPN status + intensity, plants, CONUS, state x year, resumable; lookups | run 1 Nov: 29.4 M records |
| 02a_nf_species_groups.R | laptop (internet) | plant group per species (NPN type + PLANTS growth habit) -> data/nf_species_groups.csv | built 1 Nov |
| 02b_nf_npn_clean.R | P720 | onset / end / intensity-threshold dates per plant x phenophase x year; sites, Level II, effort | built 1 Nov |
| 03_nf_fire.R | P720 | FPA-FOD (all sizes, C+), MTBS, VIIRS -> Level II x year x week x LANDFIRE fuel group | built 1 Nov |
| 04_nf_prism.R | P720 | ppet60 (60-day rain - Hargreaves PET) and tmean by Level II x day, all cells and by fuel weight | built 1 Nov; tmean added 2 Nov (rerun with `force`) |
| 05_nf_inventory.R | laptop | Tier 1 (36 events: all 13 plant phenophase classes): usable site-year dates per event, fire in the matched fuel around the typical date, lead fire record, which tests are possible | built 1 Nov |
| 06_nf_seasonal.R | laptop | Tier 2: fire rate after vs before the typical event date (NB, year + dryness), placebo dates, bootstrap; ~1.5 h | built 1 Nov |
| 06b_nf_calendar.R | laptop | fixed-calendar check on 06: own date vs national median and other ecoregions' dates (~10 min) | built 2 Nov |
| 07_nf_yeartoyear.R | laptop | Tier 3: step at the year's NPN date vs the calendar date and a temperature-predicted date, leave-one-year-out + bootstrap (minutes) | built 2 Nov |
| 07b_nf_thinning.R | laptop | step 3 rerun on thinned records (k sites a year, 20 draws): how many sites the signal needs (~30 min) | built 3 Nov |
| 08_nf_explore.R | laptop | EXPLORATORY: season size and in-season progress vs 8 fire records, FDR; summarizes 06 / 07 reruns (`outcome=<x>`) (~40 min) | built 3 Nov |
| 09_nf_confirm.R | laptop | CONFIRMATORY: eastern deciduous-tree leaf-out status vs fire given spring temperature (minutes) | built 3 Nov |
| 10_nf_weekly.R | P720 | roadmap phase 1: weekly status share and mean intensity by ecoregion x group x pheno class x year x week (and by site, not published) from the raw records (~1 h) | built 3 Nov |

### 01_nf_npn_download.R

    Rscript npn_fire/scripts/01_nf_npn_download.R test            # smoke test: AZ 2023
    Rscript npn_fire/scripts/01_nf_npn_download.R full            # all states, 2009-2024
    Rscript npn_fire/scripts/01_nf_npn_download.R full 2015:2016 AZ,NM force

One `rnpn::npn_download_status_data()` request per state x year, plant species ids in the request
(rnpn >= 1.4), animals dropped again after download, site names never requested. Each state-year
lands in `_shared/raw/npn/status/<Y>/npn_status_<S>_<Y>.parquet` (`.empty` when no records);
existing files are skipped, so rerunning retries only failed requests. Lookups (species,
phenophases, pheno classes, datasets) in `_shared/raw/npn/lookups/`. Per-request log
`npn_fire/download/requests.csv`; by-year `download_summary.csv`. The log lists which
`additional_fields` came back; fix names in config if any are reported missing.

### 02a / 02b

    Rscript npn_fire/scripts/02a_nf_species_groups.R full          # laptop: PLANTS lookups, ~10 min first time
    Rscript npn_fire/scripts/02b_nf_npn_clean.R full               # P720: all years, then combine + ecoregions

02b writes one row per plant x phenophase x year with dates and the gaps to the bracketing "no"
readings (`onset_gap`, `end_gap`, `up50_gap`, `down50_gap`); the 14-day rule is applied when the
tables are read, so 7 / 30-day sensitivity needs no rerun. Results:
`<data_root>/<run>/npn_fire/clean/` (`plant_years`, `site_years.csv`, `sites.csv`, `effort_eco.csv`,
`phenophase_map.csv`, `intensity_map.csv`, `clean_summary.csv`); to the laptop with
`Rscript scripts/91_sync_from_p720.R npn_fire`. `sites.csv` has coordinates: never publish it.

### 03

    Rscript npn_fire/scripts/03_nf_fire.R full                     # P720, minutes

Reads the fire study's tables (60 FPA-FOD, 61 MTBS, 63 VIIRS link flags) and the ecofuel study's
LANDFIRE fuel-group shares and Level II raster, read-only. Each fire is split over its 4 km cell's
fuel shares (GR, GS, SH, TU, TL, SB, NB; `unmapped` where LANDFIRE maps < half the cell). Output
`npn_fire/fire/fire_weekly.parquet`: outcome x Level II x year x week x fuel group (outcomes
`fod_all`, by cause, `fod_cplus`, `fod_cplus_acres`, `mtbs_fires`, `mtbs_km2`, `viirs_px`), plus
`fire_inventory.csv` (coverage by ecoregion x fuel x outcome) and `fire_year.csv`. Fuel groups are
matched to NPN plant groups in 05.

## Study page

`report/npn_fire.qmd` (-> `npn_fire.html`), rendered and published with the other project pages by
`scripts/91_sync_from_p720.R render publish`. Figures and tables come from `scripts/90_nf_figs.R`
(laptop, after `91 ... npn_fire`): `report/figs/<figure>_full.png` + `manifest.csv` and
`report/tables/*.csv`, both in git. Everything on the page is summarised to Level II ecoregions:
no site names or site locations. Missing figures show as "Figure pending".

    Rscript npn_fire/scripts/90_nf_figs.R full
    Rscript scripts/91_sync_from_p720.R render publish
