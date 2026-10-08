# Load a config profile and resolve machine-specific values.

load_cfg <- function(profile = NULL) {
  if (is.null(profile) || !nzchar(profile)) {
    profile <- Sys.getenv("R_CONFIG_ACTIVE", "dev")
  }
  cfg <- config::get(config = profile, file = here::here("config.yml"),
                     use_parent = FALSE)

  os <- if (.Platform$OS.type == "windows") "windows" else "linux"
  root <- Sys.getenv("GSI_DATA_ROOT", "")
  if (!nzchar(root)) root <- cfg$data_root[[os]]
  if (is.null(root)) stop("No data_root for OS '", os, "' in profile '", profile, "'")
  cfg$data_root <- normalizePath(path.expand(root), winslash = "/", mustWork = FALSE)

  cfg$profile <- profile
  cfg$run_dir <- file.path(cfg$data_root, cfg$run_name)

  # Scratch for temporary raw downloads: GSI_SCRATCH_ROOT env > config scratch_root > run_dir.
  # Point it at a local disk when data_root is on a network share.
  scratch <- Sys.getenv("GSI_SCRATCH_ROOT", "")
  if (!nzchar(scratch) && !is.null(cfg$scratch_root)) scratch <- cfg$scratch_root[[os]]
  cfg$scratch_dir <- if (is.null(scratch) || !nzchar(scratch)) cfg$run_dir else
    file.path(normalizePath(path.expand(scratch), winslash = "/", mustWork = FALSE), cfg$run_name)
  cfg$years <- seq(as.integer(cfg$years$start), as.integer(cfg$years$end))


  if (!is.null(cfg$bbox)) {
    b <- unlist(cfg$bbox)
    need <- c("xmin", "ymin", "xmax", "ymax")
    if (!all(need %in% names(b))) stop("bbox needs ", paste(need, collapse = ", "))
    cfg$bbox <- b[need]
  }
  cfg
}

# Workers for one step: GSI_WORKERS env (one run) > profile `workers_by_step: <step>` >
# profile `workers`. Also caps terra's memfrac so all workers together stay under ~60 % of RAM.
step_workers <- function(cfg, step) {
  w <- suppressWarnings(as.integer(Sys.getenv("GSI_WORKERS", "")))
  if (is.na(w) || w < 1) w <- cfg$workers_by_step[[step]] %||% cfg$workers
  w <- as.integer(w)
  list(n = w, memfrac = min(cfg$memfrac, 0.6 / w))
}
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (is.null(a)) b else a

# Paths used by the PRISM step, all under the run directory.
prism_paths <- function(cfg) {
  list(
    out_root = file.path(cfg$run_dir, "prism_daily"),          # prism_<var>_<year>.tif
    raw_root = file.path(cfg$scratch_dir, "raw", "prism_daily") # temporary daily tifs
  )
}

prism_out_file <- function(cfg, var, year) {
  file.path(prism_paths(cfg)$out_root, var, sprintf("prism_%s_%d.tif", var, year))
}

prism_raw_dir <- function(cfg, var, year) {
  file.path(prism_paths(cfg)$raw_root, var, year)
}
