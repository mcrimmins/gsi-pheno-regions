#!/usr/bin/env Rscript
# Step 95 (project page, laptop): a stand-alone page for one candidate region set, e.g. climate
# A2 at 5 regions, built from results that already exist. No new clustering: k-means was fitted
# for every k 4-30 by 20_cluster_a1.R, and 30 / 31 / 33 scored k 4-20.
#
# Steps (default: all, in this order):
#   profiles  32_region_profiles.R for this set only (eval/profiles/<id>/: regions, features,
#             season metrics, curves, ground dates). Laptop: needs 31's eval/obs/sites.csv.
#   figs      90_report_figs.R "profiles" for this set only: prof_map_<id>, prof_heat_<id>,
#             prof_curves_<id>, prof_obs_<id> (.png) and prof_table_<id>_<profile>.csv in
#             report/figs (same figures as the profile tabs on index.qmd).
#   tables    page tables in report/figs, all regionpage_<id>_*_<profile>.csv:
#             xw_<compare id>  crosswalk with each config `region_pages: compare` set (cells,
#                              share of the candidate region, share of the comparison region,
#                              comparison region name from `region_pages: region_names`)
#             scores           scores at this k, the compare k's and the location-only baseline,
#                              from eval/ab/r2_curves.csv (30), eval/obs/obs_r2.csv (31),
#                              eval/sync/sync_r2.csv (33) and clusters/<variant>/kmeans_fit.csv
#             obs              ground sites per region and data set (from step profiles)
#   render    quarto render report/regions.qmd -> report/regions_<id>.html (self-contained)
#   publish   (only when named) upload report/regions_<id>.html to the S3 folder of the project
#             pages, as 91_sync_from_p720.R does: GSI_S3_DEST, optional GSI_AWS_PROFILE and
#             GSI_CF_DIST in ~/.Renviron; AWS CLI already logged in. Never handles credentials.
#
# The set is passed to 32 and 90 with the GSI_REGION_SETS environment variable
# (R/config.R region_sets_override()); config.yml and the other region sets are untouched.
# Region numbers follow the project-page maps (region_display(), reference = first entry of
# `report: region_maps`), so a region number means the same area here and on index.qmd.
#
# Usage: Rscript scripts/95_region_page.R [profile] <partition> <k> [step ...]
#   e.g. Rscript scripts/95_region_page.R full A2_groups 5
#        Rscript scripts/95_region_page.R full A2_groups 5 tables render   (page only)
#        Rscript scripts/95_region_page.R full A2_groups 5 publish         (upload the page)
#   RStudio: gsi_profile <- "full"; gsi_page <- list(partition = "A2_groups", k = 5)
#            source("scripts/95_region_page.R")

`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
gp <-if (exists("gsi_page", envir = globalenv())) base::get("gsi_page", envir = globalenv()) else list()
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "full")
partition <- if (length(args) >= 2) args[2] else gp$partition
k <- as.integer(if (length(args) >= 3) args[3] else gp$k)
all_steps <- c("profiles", "figs", "tables", "render", "publish")
steps <- if (length(args) >= 4) args[-(1:3)] else gp$steps %||% setdiff(all_steps, "publish")
if (is.null(partition) || !length(k) || is.na(k)) stop("usage: Rscript scripts/95_region_page.R <profile> <partition> <k> [steps]")
bad <- setdiff(steps, all_steps)
if (length(bad)) stop("unknown step(s): ", paste(bad, collapse = ", "), "; use: ", paste(all_steps, collapse = ", "))

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("95_region_page", profile)
rp <- cfg$region_pages %||% list()
target <- list(partition = partition, k = k)
id <- region_set_id(target)
os <- if (.Platform$OS.type == "windows") "windows" else "linux"
t0 <- Sys.time()
log_msg("=== 95_region_page | profile: ", profile, " | set: ", id, " | steps: ", paste(steps, collapse = ", "))

# ---- Checks ---------------------------------------------------------------------------------
pf <- cfg$evaluation$partitions[[partition]]
if (is.null(pf)) stop("partition ", partition, " is not in config evaluation: partitions (",
                      paste(names(cfg$evaluation$partitions), collapse = ", "), ")")
pf <- file.path(cfg$run_dir, "clusters", pf)
if (!file.exists(pf)) stop("no cluster raster for ", partition, ": ", pf)
kb <- names(terra::rast(pf))
if (!sprintf("k%02d", k) %in% kb) stop(partition, " has no k = ", k, " (bands: ", paste(kb, collapse = " "), ")")
in_run <- cfg$report$figs_in_run_dir
if (is.list(in_run)) in_run <- in_run[[os]]
fig_dir <- if (isTRUE(in_run)) file.path(cfg$run_dir, "report_figs") else here::here("report", "figs")
if (isTRUE(in_run) && "render" %in% steps)
  stop("figures for this profile go to ", fig_dir, " (report: figs_in_run_dir), not report/figs: run 95 on the laptop")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
compare <- Filter(function(s) !(identical(s$partition, partition) && as.integer(s$k) == k), rp$compare %||% list())
label <- sprintf("%s, %d regions", rp$partition_labels[[partition]] %||% partition, k)

# Run another pipeline script in a fresh R (from the repo root, so renv and here() work) with
# this set as the only region_profiles set.
run_step <- function(script, extra = character()) {
  old <- Sys.getenv("GSI_REGION_SETS", unset = NA)
  Sys.setenv(GSI_REGION_SETS = sprintf("%s:%d", partition, k))
  on.exit(if (is.na(old)) Sys.unsetenv("GSI_REGION_SETS") else Sys.setenv(GSI_REGION_SETS = old))
  owd <- setwd(here::here()); on.exit(setwd(owd), add = TRUE)
  rs <- file.path(R.home("bin"), if (os == "windows") "Rscript.exe" else "Rscript")
  log_msg("running ", script, " ", profile, " ", paste(extra, collapse = " "))
  st <- system2(rs, c(file.path("scripts", script), profile, extra))
  if (!identical(as.integer(st), 0L)) stop(script, " failed (exit ", st, "); see its log in logs/")
}

# ---- profiles / figs ------------------------------------------------------------------------
if ("profiles" %in% steps) run_step("32_region_profiles.R")
if ("figs" %in% steps) run_step("90_report_figs.R", "profiles")

# ---- tables ---------------------------------------------------------------------------------
out_csv <- function(x, what) {
  f <- file.path(fig_dir, sprintf("regionpage_%s_%s_%s.csv", id, what, profile))
  utils::write.csv(x, f, row.names = FALSE); log_msg("wrote ", basename(f))
}
if ("tables" %in% steps) {
  grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
  # crosswalks with the comparison sets (same numbering as the maps)
  D <- region_display(cfg, grid, c(list(target), compare))
  for (s in compare) {
    cid <- region_set_id(s); Dc <- D[[cid]]
    if (is.null(Dc)) { log_warn("compare set ", cid, ": not found (skipped)"); next }
    tab <- table(D[[id]]$num, Dc$num)
    xw <- as.data.frame(tab, stringsAsFactors = FALSE); names(xw) <- c("region", "compare_region", "cells")
    xw <- xw[xw$cells > 0, ]
    xw$region <- as.integer(xw$region); xw$compare_region <- as.integer(xw$compare_region)
    xw$share_of_region <- xw$cells / rowSums(tab)[as.character(xw$region)]
    xw$share_of_compare <- xw$cells / colSums(tab)[as.character(xw$compare_region)]
    nm <- rp$region_names[[cid]]
    xw$compare_name <- if (is.null(nm)) "" else vapply(as.character(xw$compare_region), function(j) nm[[j]] %||% "", "")
    xw$compare_set <- cid
    xw <- xw[order(xw$region, -xw$share_of_region), ]
    out_csv(xw, paste0("xw_", cid))
  }

  # scores at this k, the compare k's, and the location-only baseline
  ks <- sort(unique(c(k, vapply(compare, function(s) as.integer(s$k), 1L))))
  parts <- unique(c(partition, vapply(compare, `[[`, "", "partition"), rp$baseline %||% "geo"))
  keep_pk <- function(d) d[d$partition %in% parts & d$k %in% ks, , drop = FALSE]
  sc <- data.frame()
  add <- function(test, measure, d, col) if (nrow(d))
    sc <<- rbind(sc, data.frame(test = test, measure = measure, partition = d$partition, k = d$k, value = d[[col]]))
  ev_dir <- file.path(cfg$run_dir, "eval")
  f <- file.path(ev_dir, "ab", "r2_curves.csv")
  if (file.exists(f)) {
    rc <- keep_pk(utils::read.csv(f, stringsAsFactors = FALSE))
    if (is.null(rc$features)) rc$features <- "level_shape"
    lab <- c(timing = "Timing of the satellite season", level_shape = "Satellite curves as observed")
    for (fs in names(lab)) for (dm in c("herb", "all"))
      add("satellite", sprintf("%s, %s", lab[[fs]], if (dm == "herb") "herbaceous cells" else "all cells"),
          rc[rc$features == fs & rc$domain == dm, ], "r2")
  } else log_warn("no ", f, " (satellite scores skipped)")
  f <- file.path(ev_dir, "obs", "obs_r2.csv")
  if (file.exists(f)) {
    orr <- keep_pk(utils::read.csv(f, stringsAsFactors = FALSE))
    olab <- c(grass_greenup = "NPN grasses: green-up", grass_curing = "NPN grasses: curing",
              grass_cured = "NPN grasses: cured", lfmc_herb_peak = "Live fuel moisture, herbaceous: peak",
              lfmc_woody_peak = "Live fuel moisture, woody: peak", inat_flower = "iNaturalist: flowering onset")
    for (st in intersect(unlist(cfg$report$eval_sets), unique(orr$set)))
      add("ground", if (st %in% names(olab)) olab[[st]] else st, orr[orr$set == st, ], "adj_r2")
  } else log_warn("no ", f, " (ground scores skipped)")
  f <- file.path(ev_dir, "sync", "sync_r2.csv")
  if (file.exists(f)) {
    sy <- keep_pk(utils::read.csv(f, stringsAsFactors = FALSE))
    sy <- sy[sy$domain == (cfg$synchrony$region_domain %||% "herb"), ]
    mlab <- c(sos20 = "Green-up", ndvi_peak_doy = "NDVI peak", eos20 = "End of season", cure50 = "Curing (NDII7)")
    for (m in intersect(names(mlab), unique(sy$metric))) add("synchrony", mlab[[m]], sy[sy$metric == m, ], "sync_r2")
  } else log_warn("no ", f, " (synchrony skipped)")
  for (p in setdiff(parts, rp$baseline %||% "geo")) {
    f <- file.path(cfg$run_dir, "clusters", dirname(cfg$evaluation$partitions[[p]] %||% "none"), "kmeans_fit.csv")
    if (!file.exists(f)) { log_warn("no ", f, " (stability skipped for ", p, ")"); next }
    kf <- utils::read.csv(f); kf <- kf[kf$k %in% ks, ]; kf$partition <- p
    if (!is.null(kf$mean)) add("stability", "Stability (ARI between spatial half-samples)", kf, "mean")
    if (!is.null(kf$min_size)) add("stability", "Smallest region (cells in the fit sample)", kf, "min_size")
  }
  if (nrow(sc)) out_csv(sc, "scores") else log_warn("no scores found")

  # ground sites per region
  f <- file.path(cfg$run_dir, "eval", "profiles", id, "obs.csv")
  if (file.exists(f)) out_csv(utils::read.csv(f, stringsAsFactors = FALSE), "obs") else
    log_warn("no ", f, " (run the profiles step)")
}

# ---- render ---------------------------------------------------------------------------------
if ("render" %in% steps) {
  find_quarto <- function() {
    q1 <- Sys.getenv("QUARTO_PATH", ""); if (nzchar(q1) && file.exists(q1)) return(q1)
    q2 <- Sys.which("quarto"); if (nzchar(q2)) return(unname(q2))
    pd <- Sys.getenv("RSTUDIO_PANDOC", "")          # RStudio's bundled Quarto
    if (nzchar(pd)) for (cand in file.path(dirname(pd), c("quarto.exe", "quarto"))) if (file.exists(cand)) return(cand)
    ""
  }
  qbin <- find_quarto()
  if (!nzchar(qbin)) stop("Quarto not found: install it, or set QUARTO_PATH in ~/.Renviron")
  qs <- function(x) shQuote(x, type = if (os == "windows") "cmd" else "sh")
  old <- Sys.getenv("R_LIBS", unset = NA)             # Quarto's R does not source the renv .Rprofile
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  # Render a per-set copy of the template (title and parameters filled in) from inside report/:
  # Quarto's --output breaks embed-resources when the render starts in another folder.
  src <- readLines(here::here("report", "regions.qmd"), warn = FALSE)
  hdr <- seq_len(which(src == "---")[2])                # YAML header lines only
  src[hdr] <- sub('^title: .*$', sprintf('title: "%s"', label), src[hdr])
  src[hdr] <- sub('^  set_id: .*$', sprintf('  set_id: "%s"', id), src[hdr])
  src[hdr] <- sub('^  label: .*$', sprintf('  label: "%s"', label), src[hdr])
  base <- sprintf("regions_%s", id); out <- paste0(base, ".html")
  qmd <- here::here("report", paste0(base, ".qmd"))
  writeLines(src, qmd)
  owd <- setwd(here::here("report"))
  st <- system2(qbin, c("render", qs(basename(qmd))))
  setwd(owd)
  unlink(qmd); unlink(here::here("report", paste0(base, "_files")), recursive = TRUE)
  if (is.na(old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old)
  html <- here::here("report", out)
  if (!identical(as.integer(st), 0L) || !file.exists(html)) stop("quarto render failed (exit ", st, ")")
  log_msg(sprintf("rendered report/%s (%.1f MB)", out, file.size(html) / 1e6))
}

# ---- publish (S3, only when asked) ----------------------------------------------------------
if ("publish" %in% steps) {
  html <- here::here("report", sprintf("regions_%s.html", id))
  if (!file.exists(html)) stop(html, " missing: run the render step first")
  dest <- Sys.getenv("GSI_S3_DEST", "")
  if (!grepl("^s3://", dest)) stop("Set GSI_S3_DEST=s3://bucket/path/ in ~/.Renviron")
  folder <- if (endsWith(dest, "/")) dest else paste0(dirname(dest), "/")   # same folder as the other pages
  aws <- Sys.which("aws")
  if (!nzchar(aws)) stop("AWS CLI not found (install it and run aws configure / aws sso login)")
  prof <- Sys.getenv("GSI_AWS_PROFILE", "")
  pa <- if (nzchar(prof)) c("--profile", prof) else character()
  qs <- function(x) shQuote(x, type = if (os == "windows") "cmd" else "sh")
  d <- paste0(folder, basename(html))
  st <- system2(aws, c("s3", "cp", qs(html), qs(d), "--content-type", qs("text/html; charset=utf-8"),
                       "--cache-control", qs("max-age=300"), "--only-show-errors", pa))
  if (!identical(as.integer(st), 0L))
    stop("aws s3 cp failed (exit ", st, "). Logged in? Try: aws sts get-caller-identity", if (nzchar(prof)) paste(" --profile", prof))
  log_msg("published -> ", d)
  cf <- Sys.getenv("GSI_CF_DIST", "")
  if (nzchar(cf)) {
    p_ <- sub("^s3://[^/]+", "", d)
    st <- system2(aws, c("cloudfront", "create-invalidation", "--distribution-id", cf, "--paths", qs(p_), pa))
    if (identical(as.integer(st), 0L)) log_msg("CloudFront invalidation: ", p_) else log_warn("CloudFront invalidation failed (exit ", st, ")")
  }
}

log_msg(sprintf("=== done | %s | %.1f min | log: %s", id, as.numeric(difftime(Sys.time(), t0, units = "mins")), log_file))
