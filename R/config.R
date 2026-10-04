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
  cfg$years <- seq(as.integer(cfg$years$start), as.integer(cfg$years$end))

  if (!is.null(cfg$bbox)) {
    b <- unlist(cfg$bbox)
    need <- c("xmin", "ymin", "xmax", "ymax")
    if (!all(need %in% names(b))) stop("bbox needs ", paste(need, collapse = ", "))
    cfg$bbox <- b[need]
  }
  cfg
}

# Paths used by the PRISM step, all under the run directory.
prism_paths <- function(cfg) {
  list(
    out_root = file.path(cfg$run_dir, "prism_daily"),          # prism_<var>_<year>.tif
    raw_root = file.path(cfg$run_dir, "raw", "prism_daily")    # temporary daily tifs
  )
}

prism_out_file <- function(cfg, var, year) {
  file.path(prism_paths(cfg)$out_root, var, sprintf("prism_%s_%d.tif", var, year))
}

prism_raw_dir <- function(cfg, var, year) {
  file.path(prism_paths(cfg)$raw_root, var, year)
}
