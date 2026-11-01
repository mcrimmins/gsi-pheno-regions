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
| 01_nf_npn_download.R | P720 (internet) | NPN status + intensity, plants, CONUS, state x year, resumable; lookups | built 31 Oct |
| 02_nf_npn_clean.R | laptop / P720 | plant groups, Level II / III per site, onsets / ends per plant-year, effort | next |
| 03_nf_fire.R | P720 | FPA-FOD, MTBS, FIRMS -> Level II x fuel class x week x year | planned |
| 04_nf_prism.R | P720 | ppet60 by Level II x day | planned |
| 05_nf_inventory.R | laptop | Tier 1 alignment; outcome choice | planned |
| 06 / 07 | laptop | Tier 2 seasonal, Tier 3 year to year | planned |

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
