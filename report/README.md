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

1. Open `~/.Renviron` with `file.edit("~/.Renviron")` (no package needed; on Windows R's
   `~` is the Documents folder, which is where R looks), add
   `GSI_P720_HOST=user@p720-address`, save and restart R.
2. If `ssh user@p720-address` asks for a password, set up a key (PowerShell):
   ```powershell
   ssh-keygen -t ed25519            # accept defaults; skip if ~/.ssh/id_ed25519 exists
   type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh user@p720-address "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys"
   ```

Dev figures go to `report/figs/<figure>_dev.png`. CONUS (`full`) figures are written on
the P720 next to the data (`<data_root>/full/report_figs/`), so the P720's git checkout
stays clean, and `91_sync_from_p720.R` copies them into `report/figs/` on the laptop.
`manifest.csv` records the extent, years and date of each figure. PNGs are small and are
committed from the laptop only. The page shows the
CONUS (`full`) version of a figure when it exists, otherwise the `dev` one.

## Sync, render and publish (laptop)

`scripts/91_sync_from_p720.R` does the whole update in one step:

```r
source("scripts/91_sync_from_p720.R")                 # pull CONUS figures + render
gsi_sync <- "all"; source("scripts/91_sync_from_p720.R")   # ... + upload to S3
gsi_sync <- "render"; source("scripts/91_sync_from_p720.R")  # just render
```

Steps: `figs`, `summaries`, `logs`, `render`, `publish`, `all` (= figs + render +
publish). Default is figs + render; nothing is uploaded unless `publish` or `all` is given.

Render uses Quarto (on PATH, `QUARTO_PATH`, or RStudio's bundled copy) and writes
`report/index.html` (git-ignored). The RStudio **Render** button works too.

Publish uses the AWS CLI, already logged in (`aws configure` or `aws sso login`; check
with `aws sts get-caller-identity`). Add to `~/.Renviron`:

```
GSI_S3_DEST=s3://your-bucket/path/      # folder (ends in /) or full key of the page
GSI_AWS_PROFILE=your-profile            # optional
GSI_CF_DIST=E1234567890                 # optional: CloudFront distribution to invalidate
```

The page is uploaded as `index.html` (text/html, 5-minute cache).

## Candidate-region pages

A stand-alone page for one candidate set of regions (any partition and number of regions
that was clustered), with the same map, profile table, heatmap, satellite curves and ground
dates as the profile tabs on `index.qmd`, plus crosswalks to the working 13 and 7 regions
and the evaluation scores at that number of regions. Laptop, full profile:

```bash
Rscript scripts/95_region_page.R full A2_groups 5           # -> report/regions_A2_groups_k05.html
Rscript scripts/95_region_page.R full A2_groups 5 tables render   # rebuild tables + page only
Rscript scripts/95_region_page.R full A2_groups 5 publish         # upload to the pages' S3 folder
```

`publish` runs only when named and uses the same `~/.Renviron` settings as 91
(`GSI_S3_DEST`, optional `GSI_AWS_PROFILE`, `GSI_CF_DIST`). The page goes next to
`index.html` as `regions_<partition>_k<KK>.html`; it is not linked from `index.qmd`.

The template is `regions.qmd` (parameters `set_id`, `label`); comparison sets and region
names are in config `region_pages:`. The working regions are not changed.

## When something changes

- New results: rerun the figure script, render, copy the HTML.
- New decision: add a row to **Decisions so far** and remove it from **Open questions**.
- Every update: add a line to the **Changelog** and update the **Status** table.
- New figure: add a block to `90_report_figs.R` (use `fig()` for maps) and a
  `show_fig("<name>", "<caption>")` chunk in `index.qmd`.
