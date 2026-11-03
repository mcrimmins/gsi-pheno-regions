# NPN x fire study: paths, config and logging. Self-contained: nothing here sources the parent
# gsi-pheno-regions code. Works with npn_fire/ as a subfolder of that repo or, after export, as the
# root of its own repo. Each script starts with:
#   nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
#   source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))

if (!exists("nf_root")) nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
nf_path <- function(...) file.path(nf_root, ...)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- config ----------------------------------------------------------------------------------
# Profile from the first script argument (else R_CONFIG_ACTIVE, else full); npn_fire/config.yml.
# Adds data_root (resolved per OS; NF_DATA_ROOT > GSI_DATA_ROOT > config), run_dir, out_dir
# (<run_dir>/npn_fire) and raw_dir (<data_root>/_shared/raw).
load_nf_cfg <- function(profile = NULL) {
  if (is.null(profile) || !nzchar(profile)) profile <- Sys.getenv("R_CONFIG_ACTIVE", "full")
  cfg <- config::get(config = profile, file = nf_path("config.yml"), use_parent = FALSE)
  os <- if (.Platform$OS.type == "windows") "windows" else "linux"
  pick <- function(x) if (is.list(x)) x[[os]] else x
  root <- Sys.getenv("NF_DATA_ROOT", "")
  if (!nzchar(root)) root <- Sys.getenv("GSI_DATA_ROOT", "")
  if (!nzchar(root)) root <- pick(cfg$data_root)
  cfg$data_root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)
  cfg$run_dir <- file.path(cfg$data_root, cfg$run_name)
  cfg$out_dir <- file.path(cfg$run_dir, "npn_fire")
  cfg$raw_dir <- file.path(cfg$data_root, "_shared", "raw")
  cfg$profile <- profile
  cfg
}
nf_out <- function(cfg, ...) { d <- file.path(cfg$out_dir, ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
nf_raw <- function(cfg, ...) { d <- file.path(cfg$raw_dir, ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }

# ---- logging: timestamped lines to the console and npn_fire/logs/<script>_<profile>_<time>.log ---
log_init <- function(script, profile) {
  dir.create(nf_path("logs"), showWarnings = FALSE)
  ov <- sub("^outcome=", "", grep("^outcome=", commandArgs(trailingOnly = TRUE), value = TRUE))   # 06 / 07 exploratory reruns
  tag <- if (length(ov)) paste0(profile, "_", ov[1]) else profile
  f <- nf_path("logs", sprintf("%s_%s_%s.log", script, tag, format(Sys.time(), "%Y%m%d-%H%M%S")))
  options(nf.log_file = f)
  f
}
log_msg <- function(..., level = "INFO") {
  line <- sprintf("%s [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, paste0(...))
  message(line)
  f <- getOption("nf.log_file")
  if (!is.null(f)) cat(line, "\n", file = f, append = TRUE, sep = "")
  invisible(line)
}
log_warn <- function(...) log_msg(..., level = "WARN")
log_err <- function(...) log_msg(..., level = "ERROR")
