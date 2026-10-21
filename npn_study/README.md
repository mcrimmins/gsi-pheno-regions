# NPN grass curing study

Ground observations of grass green-up and curing across the contiguous United States from
USA National Phenology Network (NPN) percent-green records, most of them from NEON field sites,
with daily PRISM climate. Questions: the record (Q1), typical dates by ecoregion and grass type
(Q2), place versus year (Q3), climate across space (Q4), weather from year to year (Q5) and
whether simple rules are enough (Q6).

This folder is self-contained: it reads nothing from the rest of `gsi-pheno-regions`, so it can
be exported as its own repository at publication. Scripts find this folder whether it sits
inside that repo (`npn_study/`) or is the repository root.

## Layout

| Path | What |
|---|---|
| `config.yml` | data paths (per OS), climate thresholds, study settings, sensitivity runs |
| `R/npn_setup.R` | paths, config loader, logger, PET / freeze / rain-onset helpers |
| `R/npn.R` | study helpers (series-years, site-year metrics, normals, regions, sensitivity) |
| `scripts/` | the analysis, in run order below |
| `data/` | `npn_species_pathway.csv` (C3 / C4, life cycle, origin, crop), `npn_site_ecoregions.csv` (EPA Level II per site), `npn_species_fire.csv` (invasive grasses that raise fire occurrence, Fusco et al. 2019), `npn_species_traits.csv` (58d), `npn_site_fuelmodels.csv` (59b); `npn_site_names.csv` is local only (not in git) |
| `report/` | `npn_study.qmd` (study page), `figs/` (PNGs + `manifest.csv`) and `tables/` (`npn_species.csv`, from 58c) |
| `logs/` | run logs (not in git) |

## Inputs

- NPN graminoid percent-green tables `npn_curing_metrics.rds` and `npn_green_obs.rds`, built from
  NPN status and intensity data (`rnpn`) by `GrowingSeasonIndex/study/02c_get_npn_intensity.R`
  (config `study_dir`; `GSI_OBS_DIR` overrides).
- Daily PRISM 4 km stacks `prism_<var>_<year>.tif` (ppt, tmin, tmax, vpdmax; config
  `prism_daily_dir`), read only by script 50.
- EPA / CEC Level II ecoregions (downloaded by script 59).
- USDA PLANTS checklist and web service (58d, internet); LANDFIRE FBFM40 CONUS GeoTIFF, downloaded by 59a from landfire.gov (read by 59b).

Outputs go to `<data_root>/full/npn_study/`.

## Run order (from the repository root)

```
Rscript npn_study/scripts/50_npn_prism.R full            # daily PRISM at every site (where PRISM lives)
Rscript npn_study/scripts/54a_npn_site_names.R full      # NPN site names, local only (rnpn, internet)
Rscript npn_study/scripts/59_npn_ecoregions.R full       # EPA Level II ecoregion per site (sf, internet)
Rscript npn_study/scripts/51_npn_metrics.R full          # series-years, site-year metrics, normals
Rscript npn_study/scripts/52_npn_ladder.R full           # model ladder L0-L4 (Q5, Q6)
Rscript npn_study/scripts/53_npn_variance.R full         # variance partition (Q3)
Rscript npn_study/scripts/54_npn_site_signal.R full      # shared site-year signal vs weather
Rscript npn_study/scripts/55_npn_filter_inventory.R full # what each filter removes
Rscript npn_study/scripts/58_npn_record_climatology.R full  # record, climatology, climate across space (Q1, Q2, Q4)
Rscript npn_study/scripts/58b_npn_c3c4_gap.R full         # C3 -> C4 gap by region (optional: region CSV + tag)
Rscript npn_study/scripts/58d_npn_species_traits.R full    # USDA PLANTS height, growth form (internet) -> data/npn_species_traits.csv
Rscript npn_study/scripts/59a_npn_landfire_download.R full # LANDFIRE FBFM40 CONUS GeoTIFF (3.27 GB zip, internet) -> <data_root>/landfire
Rscript npn_study/scripts/59b_npn_fuelmodels.R full        # LANDFIRE fuel model per site (local FBFM40 GeoTIFF) -> data/npn_site_fuelmodels.csv
Rscript npn_study/scripts/58c_npn_species.R full           # species list for the page (report/tables/npn_species.csv)
Rscript npn_study/scripts/90_npn_figs.R full             # figures for the page
```

Sensitivity without 2016-2017: add `no1617` as a second argument to 52, 53 and 54.
Exploratory PhenoCam comparison: `56_phenocam_fetch.R` (phenocamr, internet), `57_phenocam_compare.R`.
The page is rendered and published with the parent repo's `scripts/91_sync_from_p720.R render`.

Packages: config, here, lme4, ranger, sf, terra, maps, rnpn, phenocamr.
