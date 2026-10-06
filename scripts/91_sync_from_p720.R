#!/usr/bin/env Rscript
# Step 91 (laptop): copy results from the P720 with scp.
#   figs       report/figs/*_full.png, merged into the local report/figs/manifest.csv
#              (rows for other profiles, e.g. dev, are kept)            [default]
#   summaries  full-run across-year summaries and static layers, into the local data
#              root (<data_root>/full/summaries, <data_root>/full/static) for plotting
#   logs       the P720 logs/ folder, into logs/p720/
#
# Needs: OpenSSH client (built into Windows 10/11) and key-based SSH login to the P720
# (scp can't ask for a password from RStudio). Host comes from the GSI_P720_HOST
# environment variable (e.g. user@p720.local or an alias from ~/.ssh/config), set in
# ~/.Renviron so it never goes in git. Remote paths: config `sync:`.
#
# Usage: Rscript scripts/91_sync_from_p720.R [figs] [summaries] [logs]
#   RStudio: gsi_sync <- c("figs", "summaries"); source("scripts/91_sync_from_p720.R")

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args)) args else if (exists("gsi_sync", envir = globalenv())) {
  base::get("gsi_sync", envir = globalenv())
} else "figs"
what <- match.arg(what, c("figs", "summaries", "logs"), several.ok = TRUE)

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

host <- Sys.getenv("GSI_P720_HOST", "")
if (!nzchar(host)) stop("Set GSI_P720_HOST (e.g. user@p720.local) in ~/.Renviron, then restart R")
cfg <- load_cfg("full")                                     # local paths for the full run
sync <- cfg$sync %||% list()
remote_repo <- sync$remote_repo %||% "~/RProjects/gsi-pheno-regions"
remote_run <- file.path(sync$remote_data_root %||% "~/data/gsi-pheno-regions", cfg$run_name)
log_file <- log_init("91_sync_from_p720", "full")
log_msg("=== 91_sync_from_p720 | ", host, " | ", paste(what, collapse = ", "))

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
chk <- run("ssh", c("-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, "echo ok"))
if (!chk$ok || !any(grepl("^ok$", chk$out))) {
  log_err("can't reach ", host, " with key login: ", paste(chk$out, collapse = " "))
  stop("SSH to ", host, " failed. Check the host name, that the P720 is on, and that your ",
       "SSH key is in ~/.ssh/authorized_keys there (see report/README.md).")
}

# ---- figs -------------------------------------------------------------------------------
if ("figs" %in% what) {
  fig_dir <- here::here("report", "figs")
  if (scp_get(file.path(remote_repo, "report/figs/*_full.png"), fig_dir)) {
    log_msg("figures: ", paste(basename(list.files(fig_dir, "_full\\.png$")), collapse = ", "))
  }
  tmp <- tempfile(); dir.create(tmp)
  if (scp_get(file.path(remote_repo, "report/figs/manifest.csv"), tmp)) {
    rem <- utils::read.csv(file.path(tmp, "manifest.csv"), stringsAsFactors = FALSE)
    mf <- file.path(fig_dir, "manifest.csv")
    loc <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else rem[0, ]
    m <- rbind(loc[loc$profile != "full", , drop = FALSE], rem[rem$profile == "full", , drop = FALSE])
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

log_msg("=== done")
if ("figs" %in% what) message("\nNext: render report/index.qmd, then commit report/figs.")
