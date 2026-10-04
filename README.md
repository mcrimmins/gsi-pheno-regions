# gsi-pheno-regions

Fuel phenology regionalization for CONUS herbaceous fuels (GSI for Wildland Fire project).
See `CLAUDE.md` for the design and conventions.

## Setup
```r
renv::restore()
```
Data are stored outside the repo under `data_root` in `config.yml`
(override with the `GSI_DATA_ROOT` environment variable).

## Run
```sh
Rscript scripts/01_prism_daily.R test   # minutes: 5 days, 2 variables
Rscript scripts/01_prism_daily.R dev    # Arizona, 2023-2025
Rscript scripts/01_prism_daily.R full   # CONUS (P720, inside tmux)
```

From RStudio with the project open (background job keeps the console free):
```r
gsi_profile <- "test"   # or "dev"
rstudioapi::jobRunScript("scripts/01_prism_daily.R", workingDir = getwd(), importEnv = TRUE)
```

Data: PRISM Group, Oregon State University, https://prism.oregonstate.edu.
