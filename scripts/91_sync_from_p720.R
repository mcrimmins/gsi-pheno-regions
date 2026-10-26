#!/usr/bin/env Rscript
# Step 91 (laptop): copy results from the P720 with scp.
#   figs       <remote run dir>/report_figs/*_full.png (written there by 90_report_figs.R
#              full) into report/figs, merged into the local report/figs/manifest.csv
#              (rows for other profiles, e.g. dev, are kept)            [default]
#   summaries  full-run across-year summaries and static layers, into the local data
#              root (<data_root>/full/summaries, <data_root>/full/static) for plotting
#   logs       the P720 logs/ folder, into logs/p720/
#   fire       fire study results: every <remote run dir>/fire_study/<step>/*.csv, plus the
#              per-record tables (*.parquet / *.rds) of the steps in config `sync: fire_tables`
#              (default fod), into <data_root>/full/fire_study/<step>/ on the laptop
#   ecofuel    ecoregion x fuel type work (70-79): <remote run dir>/ecofuel/*.tif, *.csv and
#              ecofuel/inventory/* into <data_root>/full/ecofuel/ on the laptop (rasters are 4 km, small)
#   render     quarto render every page: report/index.qmd, npn_study/report/npn_study.qmd (NPN
#              curing study, its own subfolder), report/fire_study.qmd and report/ecofuel.qmd
#              -> <page>.html next to its .qmd
#   publish    upload report/<page>.html (all pages, same folder) to S3 with the AWS CLI (optional CloudFront
#              invalidation). Never runs unless asked for ("publish" or "all").
#   all        figs + render + publish
#   Default: figs + render.
#
# Needs: OpenSSH client (built into Windows 10/11) and key-based SSH login to the P720
# (scp can't ask for a password from RStudio). Host comes from the GSI_P720_HOST
# environment variable (e.g. user@p720.local or an alias from ~/.ssh/config), set in
# ~/.Renviron so it never goes in git. Remote paths: config `sync:`.
#
# Render needs Quarto: on PATH, QUARTO_PATH, or the copy bundled with RStudio. The R
# session's library paths (renv) are passed to Quarto's R process.
#
# Publish needs the AWS CLI already logged in (aws configure / aws sso login); this script
# never handles credentials. Settings in ~/.Renviron:
#   GSI_S3_DEST=s3://bucket/path/         folder (ends in /) or full key of the HTML file
#   GSI_AWS_PROFILE=name                  optional AWS CLI profile
#   GSI_CF_DIST=E123ABC                   optional CloudFront distribution to invalidate
#
# Usage: Rscript scripts/91_sync_from_p720.R [figs] [summaries] [logs] [fire] [ecofuel] [render] [publish] [all]
#   RStudio: gsi_sync <- "all"; source("scripts/91_sync_from_p720.R")

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args else if (exists("gsi_sync", envir = globalenv())) {
  base::get("gsi_sync", envir = globalenv())
} else c("figs", "render")
what <- match.arg(what, c("figs", "summaries", "logs", "fire", "ecofuel", "render", "publish", "all"),
                  several.ok = TRUE)
if ("all" %in% what) what <- unique(c(setdiff(what, "all"), "figs", "render", "publish"))
need_ssh <- any(c("figs", "summaries", "logs", "fire", "ecofuel") %in% what)

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

host <- Sys.getenv("GSI_P720_HOST", "")
if (need_ssh && !nzchar(host))
  stop("Set GSI_P720_HOST (e.g. user@p720.local) in ~/.Renviron, then restart R")
cfg <- load_cfg("full")                                     # local paths for the full run
sync <- cfg$sync %||% list()
remote_repo <- sync$remote_repo %||% "~/RProjects/gsi-pheno-regions"
remote_run <- file.path(sync$remote_data_root %||% "~/data/gsi-pheno-regions", cfg$run_name)
log_file <- log_init("91_sync_from_p720", "full")
log_msg("=== 91_sync_from_p720 | ", if (need_ssh) host else "local", " | ",
        paste(what, collapse = ", "))

win <- .Platform$OS.type == "windows"
q <- function(x) shQuote(x, type = if (win) "cmd" else "sh")
run <- function(cmd, args) {
  out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE))
  list(ok = is.null(attr(out, "status")) || attr(out, "status") == 0, out = out)
}
# Copy remote path(s) (globs allowed; expanded on the P720) into a local directory.
# Returns TRUE on success, FALSE (with a log warning) if nothing matched or scp failed.
scp_get <- function(remote, local_dir) {
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  r <- run("scp", c("-p", "-q", "-o", "BatchMode=yes",
                    q(paste0(host, ":", remote)), q(normalizePath(local_dir, "/"))))
  if (!r$ok) { log_warn("scp ", remote, ": ", paste(r$out, collapse = " ")); return(FALSE) }
  TRUE
}

# Connection check first: clear message instead of a hang or a password prompt.
if (need_ssh) {
  chk <- run("ssh", c("-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, "echo ok"))
  if (!chk$ok || !any(grepl("^ok$", chk$out))) {
    log_err("can't reach ", host, " with key login: ", paste(chk$out, collapse = " "))
    stop("SSH to ", host, " failed. Check the host name, that the P720 is on, and that your ",
         "SSH key is in ~/.ssh/authorized_keys there (see report/README.md).")
  }
}

# ---- figs -------------------------------------------------------------------------------
if ("figs" %in% what) {
  fig_dir <- here::here("report", "figs")
  if (scp_get(file.path(remote_run, "report_figs/*_full.png"), fig_dir)) {
    log_msg("figures: ", paste(basename(list.files(fig_dir, "_full\\.png$")), collapse = ", "))
  }
  tmp <- tempfile(); dir.create(tmp)
  if (scp_get(file.path(remote_run, "report_figs/manifest.csv"), tmp)) {
    rem <- utils::read.csv(file.path(tmp, "manifest.csv"), stringsAsFactors = FALSE)
    mf <- file.path(fig_dir, "manifest.csv")
    loc <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else rem[0, ]
    # keep local rows except full figures the P720 also made (full figures made on the
    # laptop, e.g. the evaluation figures, stay)
    remf <- rem[rem$profile == "full", , drop = FALSE]
    m <- rbind(loc[!(loc$profile == "full" & loc$figure %in% remf$figure), , drop = FALSE], remf)
    utils::write.csv(m[order(m$figure, m$profile), ], mf, row.names = FALSE)
    log_msg("manifest merged: ", sum(m$profile == "full"), " full + ",
            sum(m$profile != "full"), " other rows")
  }
  unlink(tmp, recursive = TRUE)
}

# ---- summaries and static layers ----------------------------------------------------------
if ("summaries" %in% what) {
  if (scp_get(file.path(remote_run, "summaries/*.tif"), file.path(cfg$run_dir, "summaries")))
    log_msg("summaries -> ", file.path(cfg$run_dir, "summaries"))
  for (f in c("grid_mask.tif", "elev.tif", "herb_share.tif", "herb_mask.tif")) {
    scp_get(file.path(remote_run, "static", f), file.path(cfg$run_dir, "static"))
  }
  log_msg("static -> ", file.path(cfg$run_dir, "static"))
}

# ---- logs ---------------------------------------------------------------------------------
if ("logs" %in% what) {
  if (scp_get(file.path(remote_repo, "logs/*.log"), here::here("logs", "p720")))
    log_msg("logs -> logs/p720")
}

# ---- fire study results ----------------------------------------------------------------------
if ("fire" %in% what) {
  rfire <- file.path(remote_run, "fire_study")
  ls <- run("ssh", c("-o", "BatchMode=yes", host, q(paste0("ls -1d ", rfire, "/*/"))))
  steps <- if (ls$ok) basename(sub("/+$", "", ls$out[nzchar(ls$out)])) else character()
  if (!length(steps)) log_warn("no fire_study results on the P720 (", rfire, ")")
  tables <- unlist(sync$fire_tables %||% list("fod"))
  for (st in steps) {
    ldir <- file.path(cfg$run_dir, "fire_study", st)
    if (scp_get(file.path(rfire, st, "*.csv"), ldir)) log_msg("fire_study/", st, ": ", paste(basename(list.files(ldir, "\\.csv$")), collapse = ", "))
    if (st %in% tables) {
      tb <- run("ssh", c("-o", "BatchMode=yes", host, q(paste0("ls -1 ", rfire, "/", st, "/*.parquet ", rfire, "/", st, "/*.rds 2>/dev/null"))))
      for (f in tb$out[grepl("\\.(parquet|rds)$", tb$out)]) {
        if (scp_get(f, ldir)) log_msg(sprintf("fire_study/%s/%s (%.0f MB)", st, basename(f), file.size(file.path(ldir, basename(f))) / 2^20))
      }
    }
  }
}

# ---- ecoregion x fuel type (70-79) ------------------------------------------------------------
if ("ecofuel" %in% what) {
  reco <- file.path(remote_run, "ecofuel"); leco <- file.path(cfg$run_dir, "ecofuel")
  got <- c(scp_get(file.path(reco, "*.tif"), leco), scp_get(file.path(reco, "*.csv"), leco))
  inv <- scp_get(file.path(reco, "inventory", "*"), file.path(leco, "inventory"))
  if (any(got)) log_msg("ecofuel: ", paste(basename(list.files(leco, "\\.(tif|csv)$")), collapse = ", "))
  if (inv) log_msg("ecofuel/inventory: ", paste(basename(list.files(file.path(leco, "inventory"))), collapse = ", "))
  if (!any(got) && !inv) log_warn("no ecofuel results on the P720 (", reco, ")")
}

# ---- render -------------------------------------------------------------------------------
pages <- c("index", "npn_study", "fire_study", "ecofuel")          # <page>.qmd -> <page>.html
# Folder of each page: report/, except the NPN study, which lives in npn_study/report/ (moved Oct 18
# 2026 so it can become its own repo). All pages still publish to the same S3 folder.
page_dir <- function(pg) if (pg == "npn_study") here::here("npn_study", "report") else here::here("report")
pages <- pages[file.exists(vapply(pages, function(pg) file.path(page_dir(pg), paste0(pg, ".qmd")), ""))]
html <- here::here("report", "index.html")
find_quarto <- function() {
  q1 <- Sys.getenv("QUARTO_PATH", ""); if (nzchar(q1) && file.exists(q1)) return(q1)
  q2 <- Sys.which("quarto"); if (nzchar(q2)) return(unname(q2))
  # RStudio bundles Quarto next to its pandoc (RSTUDIO_PANDOC = .../quarto/bin/tools)
  pd <- Sys.getenv("RSTUDIO_PANDOC", "")
  if (nzchar(pd)) for (cand in file.path(dirname(pd), c("quarto.exe", "quarto")))
    if (file.exists(cand)) return(cand)
  ""
}
if ("render" %in% what) {
  qbin <- find_quarto()
  if (!nzchar(qbin)) stop("Quarto not found: install it, or set QUARTO_PATH in ~/.Renviron")
  # Quarto starts its own R in report/, where the project .Rprofile (renv) isn't sourced,
  # so hand it this session's library paths.
  old <- Sys.getenv("R_LIBS", unset = NA)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  for (pg in pages) {
    out <- file.path(page_dir(pg), paste0(pg, ".html"))
    r <- run(qbin, c("render", q(file.path(page_dir(pg), paste0(pg, ".qmd")))))
    if (!r$ok || !file.exists(out)) {
      if (is.na(old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old)
      log_err("render failed (", pg, "):\n", paste(utils::tail(r$out, 15), collapse = "\n"))
      stop("quarto render failed for ", pg, ".qmd (see log)")
    }
    log_msg("rendered ", out, " (", round(file.size(out) / 1e6, 1), " MB)")
  }
  if (is.na(old)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = old)
}

# ---- publish (S3) ---------------------------------------------------------------------------
if ("publish" %in% what) {
  dest <- Sys.getenv("GSI_S3_DEST", "")
  if (!grepl("^s3://", dest)) stop("Set GSI_S3_DEST=s3://bucket/path/ in ~/.Renviron")
  folder <- if (endsWith(dest, "/")) dest else paste0(dirname(dest), "/")   # all pages go in one folder
  aws <- Sys.which("aws")
  if (!nzchar(aws)) stop("AWS CLI not found (install it and run aws configure / aws sso login)")
  prof <- Sys.getenv("GSI_AWS_PROFILE", "")
  pa <- if (nzchar(prof)) c("--profile", prof) else character()
  paths <- character()
  for (pg in pages) {
    h <- file.path(page_dir(pg), paste0(pg, ".html"))
    if (!file.exists(h)) stop(h, " missing: run the render step first")
    newest_input <- max(file.mtime(c(file.path(page_dir(pg), paste0(pg, ".qmd")),
                                     list.files(file.path(page_dir(pg), "figs"), full.names = TRUE))))
    if (file.mtime(h) < newest_input)
      log_warn(pg, ".html is older than its .qmd or a figure; publishing it anyway")
    d <- if (pg == "index" && !endsWith(dest, "/")) dest else paste0(folder, pg, ".html")   # a full key names the main page
    r <- run(aws, c("s3", "cp", q(h), q(d), "--content-type", q("text/html; charset=utf-8"),
                    "--cache-control", q("max-age=300"), "--only-show-errors", pa))
    if (!r$ok) {
      log_err("S3 upload failed: ", paste(r$out, collapse = " "))
      stop("aws s3 cp failed (see log). Logged in? Try: aws sts get-caller-identity", if (nzchar(prof)) paste(" --profile", prof))
    }
    log_msg("published -> ", d)
    paths <- c(paths, sub("^s3://[^/]+", "", d))
  }
  cf <- Sys.getenv("GSI_CF_DIST", "")
  if (nzchar(cf)) {
    paths <- unique(c(paths, sub("index\\.html$", "", paths[grepl("index\\.html$", paths)])))  # files and the folder URL
    r <- run(aws, c("cloudfront", "create-invalidation", "--distribution-id", cf,
                    "--paths", vapply(paths, q, ""), pa))
    if (r$ok) log_msg("CloudFront invalidation: ", paste(paths, collapse = " "))
    else log_warn("CloudFront invalidation failed: ", paste(r$out, collapse = " "))
  }
}

log_msg("=== done")
if ("figs" %in% what) message("\nRemember to commit report/figs (PNGs + manifest.csv).")
