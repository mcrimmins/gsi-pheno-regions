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

## When something changes

- New results: rerun the figure script, render, copy the HTML.
- New decision: add a row to **Decisions so far** and remove it from **Open questions**.
- Every update: add a line to the **Changelog** and update the **Status** table.
- New figure: add a block to `90_report_figs.R` (use `fig()` for maps) and a
  `show_fig("<name>", "<caption>")` chunk in `index.qmd`.
