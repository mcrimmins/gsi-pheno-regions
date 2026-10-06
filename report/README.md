# Project page

`index.qmd` is the living project page: what the project is, the analysis choices so far
and preliminary figures. It renders to one self-contained HTML file (figures embedded)
that can be copied to a website or shared directly.

## Update the figures

Figures come from `scripts/90_report_figs.R`, run wherever the data for a profile lives.
Each figure is skipped (with a log line) if its inputs don't exist yet.

```r
# Laptop, RStudio (Arizona test area)
gsi_profile <- "dev"; source("scripts/90_report_figs.R")
```

```bash
# P720 (CONUS)
Rscript scripts/90_report_figs.R full
```

Then pull the CONUS figures to the laptop (merges `manifest.csv`, keeps the dev rows):

```r
source("scripts/91_sync_from_p720.R")                         # figures only
gsi_sync <- c("figs", "summaries"); source("scripts/91_sync_from_p720.R")  # + rasters
```

`91_sync_from_p720.R` uses scp with key login. One-time setup on the laptop:

1. In `~/.Renviron` (`usethis::edit_r_environ()`), add `GSI_P720_HOST=user@p720-address`
   and restart R.
2. If `ssh user@p720-address` asks for a password, set up a key (PowerShell):
   ```powershell
   ssh-keygen -t ed25519            # accept defaults; skip if ~/.ssh/id_ed25519 exists
   type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh user@p720-address "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
   ```

Outputs are `report/figs/<figure>_<profile>.png` plus `report/figs/manifest.csv`
(extent, years and date of each figure). PNGs are small and go in git. The page shows the
CONUS (`full`) version of a figure when it exists, otherwise the `dev` one.

## Render

In RStudio, open `report/index.qmd` and click **Render**, or from a terminal:

```bash
quarto render report/index.qmd
```

This writes `report/index.html` (ignored by git). Copy it to the website.

## When something changes

- New results: rerun the figure script, render, copy the HTML.
- New decision: add a row to **Decisions so far** and remove it from **Open questions**.
- Every update: add a line to the **Changelog** and update the **Status** table.
- New figure: add a block to `90_report_figs.R` (use `fig()` for maps) and a
  `show_fig("<name>", "<caption>")` chunk in `index.qmd`.
