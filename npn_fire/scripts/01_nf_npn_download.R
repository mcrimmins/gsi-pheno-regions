#!/usr/bin/env Rscript
# Step 01 (NPN x fire study; P720, internet, rnpn): USA-NPN plant status and intensity records for
# CONUS, one request per state x year, cached raw and never re-fetched.
# Scope: project doc claude/gsi-npn-fire-scope.md. Decided 31 Oct 2026 (Mike): status AND intensity
# (one record type), all plant groups incl. conifers, no animals, 2009-2024, Level II ecoregions
# downstream (02). NPN onsets / ends are computed by our own rules in 02, not taken from NPN's
# phenometric products.
#
#   - lookups (once; `force` refreshes): species (kingdom, functional_type, species_type, taxonomy),
#     phenophases (pheno_class_id, category), pheno classes, datasets, plant types
#     -> <data_root>/_shared/raw/npn/lookups/*.csv, copied to <run_dir>/npn_fire/download/lookups/
#   - status records: rnpn::npn_download_status_data(years = Y, states = S, plant species ids)
#     -> <data_root>/_shared/raw/npn/status/<Y>/npn_status_<S>_<Y>.parquet (plants only; site names
#     never requested; drop_columns removed). A request with no records leaves an .empty marker.
#     Existing files / markers are skipped, so a rerun resumes and retries only failed requests.
#   - <run_dir>/npn_fire/download/requests.csv: one row per request ever made (state, year, status,
#     rows, sites, individuals, species, phenophases, share of yes, seconds), and
#     download_summary.csv (by year: rows, sites, states with data, failed requests).
#
# Usage: Rscript npn_fire/scripts/01_nf_npn_download.R [profile] [years] [states] [force]
#   years: 2015 | 2015:2018 | 2015,2017 ; states: AZ | AZ,NM ; force = refetch lookups and the
#   selected state-years. Smoke test first: Rscript npn_fire/scripts/01_nf_npn_download.R test
# Needs: rnpn (>= 1.4 for the plant-species filter in the request), arrow (else csv.gz).

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && !grepl("^[0-9]|^force$|^[A-Z]{2}(,|$)", args[1])) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
rest <- if (length(args) && identical(args[1], profile)) args[-1] else args
force <- "force" %in% rest; rest <- setdiff(rest, "force")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("rnpn", quietly = TRUE)) stop("needs the rnpn package: renv::install('rnpn')")

cfg <- load_nf_cfg(profile)
log_file <- log_init("01_nf_npn_download", profile)
nf <- cfg$npn_fire; dl <- nf$download
years <- nf_parse_years(rest[grepl("^[0-9]", rest)][1]) %||% seq(nf$years$start, nf$years$end)
states <- nf_parse_states(rest[grepl("^[A-Za-z]", rest)][1]) %||% unlist(nf$states)
bad <- setdiff(states, unlist(nf$states))
if (length(bad)) stop("states not in config npn_fire: states: ", paste(bad, collapse = ", "))
rv <- utils::packageVersion("rnpn")
log_msg("=== 01_nf_npn_download | profile: ", profile, " | rnpn ", as.character(rv), " | years ",
        paste(range(years), collapse = "-"), " | ", length(states), " states | force: ", force)
log_msg("data_root: ", cfg$data_root)

raw <- nf_raw(cfg, "npn"); lk_dir <- nf_raw(cfg, "npn", "lookups"); out <- nf_out(cfg, "download")

# ---- lookups -----------------------------------------------------------------------------------
get_lookup <- function(name, fun) {
  f <- file.path(lk_dir, paste0(name, ".csv"))
  if (force || !file.exists(f)) {
    x <- tryCatch(as.data.frame(fun()), error = function(e) { log_err(name, " lookup failed: ", conditionMessage(e)); NULL })
    if (is.null(x)) { if (file.exists(f)) return(utils::read.csv(f, stringsAsFactors = FALSE)); return(NULL) }
    # list columns (if any) to text so the CSV is flat
    for (k in names(x)) if (is.list(x[[k]])) x[[k]] <- vapply(x[[k]], function(v) paste(v, collapse = ";"), "")
    utils::write.csv(x, f, row.names = FALSE)
    log_msg(sprintf("lookup %-14s %6d rows (fetched)", name, nrow(x)))
  } else {
    x <- utils::read.csv(f, stringsAsFactors = FALSE)
    log_msg(sprintf("lookup %-14s %6d rows (cached %s)", name, nrow(x), format(file.mtime(f), "%Y-%m-%d")))
  }
  dir.create(file.path(out, "lookups"), showWarnings = FALSE)
  file.copy(f, file.path(out, "lookups", basename(f)), overwrite = TRUE)
  x
}
sp <- get_lookup("species", rnpn::npn_species)
if (is.null(sp)) stop("no species lookup: NPN services unreachable?")
get_lookup("phenophases", rnpn::npn_phenophases)
get_lookup("pheno_classes", rnpn::npn_pheno_classes)
get_lookup("datasets", rnpn::npn_datasets)
get_lookup("plant_types", function() rnpn::npn_species_types("Plantae"))

plants <- sp[sp$kingdom == "Plantae" & !(sp$functional_type %in% unlist(dl$drop_functional_types)), ]
tb <- sort(table(plants$functional_type), decreasing = TRUE)
log_msg(sprintf("plant species: %d (of %d in NPN) | %s", nrow(plants), nrow(sp),
                paste(sprintf("%s %d", names(tb), tb), collapse = ", ")))
pf <- dl$plant_filter %||% "auto"
if (pf == "auto") pf <- if (rv >= "1.4.0") "species_ids" else "after"
if (pf == "species_ids" && rv < "1.4.0") log_warn("rnpn < 1.4 sends GET requests; a list of ", nrow(plants), " species ids may be too long. Set plant_filter: after if requests fail.")
log_msg("plant filter: ", pf, if (pf == "after") " (all records requested, animals dropped after download)" else " (plant species ids in each request)")
plant_ids <- if (pf == "species_ids") plants$species_id else NULL

# ---- requests ----------------------------------------------------------------------------------
req_file <- file.path(out, "requests.csv")
log_req <- function(row) utils::write.table(row, req_file, sep = ",", row.names = FALSE,
                                            col.names = !file.exists(req_file), append = file.exists(req_file))
af <- unlist(dl$additional_fields); drop_cols <- tolower(unlist(dl$drop_columns))
checked_fields <- FALSE; n_err <- 0; n_done <- 0; n_skip <- 0
jobs <- expand.grid(state = states, year = years, stringsAsFactors = FALSE)
jobs <- jobs[order(jobs$year, jobs$state), ]
t_all <- Sys.time()

for (j in seq_len(nrow(jobs))) {
  S <- jobs$state[j]; Y <- jobs$year[j]
  ydir <- file.path(raw, "status", Y); dir.create(ydir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(ydir, sprintf("npn_status_%s_%d", S, Y)); empty <- paste0(base, ".empty")
  have <- nf_table_path(base)
  if (!force && (!is.na(have) || file.exists(empty))) { n_skip <- n_skip + 1; next }
  tmp <- tempfile(sprintf("npn_%s_%d_", S, Y), fileext = ".csv")
  t0 <- Sys.time()
  res <- tryCatch({
    suppressMessages(rnpn::npn_download_status_data(
      request_source = dl$request_source, years = Y, states = S, species_ids = plant_ids,
      additional_fields = af, climate_data = FALSE, email = dl$email, download_path = tmp))
    "ok"
  }, error = function(e) conditionMessage(e))
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  row <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), state = S, year = Y, status = "error",
                    rows = NA, plant_rows = NA, sites = NA, individuals = NA, species = NA,
                    phenophases = NA, yes_share = NA, with_intensity = NA, seconds = secs, message = "")
  if (!identical(res, "ok")) {
    n_err <- n_err + 1
    row$message <- substr(gsub("[\r\n,]", " ", res), 1, 200)
    log_warn(sprintf("%s %d: request failed (%s s): %s", S, Y, secs, row$message))
    log_req(row)
    if (n_err >= (dl$max_errors %||% 25)) { log_err("too many failed requests (", n_err, "); stopping. Rerun to resume."); break }
    Sys.sleep(dl$sleep_sec %||% 2); next
  }
  if (!file.exists(tmp) || file.size(tmp) == 0) {
    writeLines(format(Sys.time()), empty)
    row$status <- "empty"; row$rows <- 0; row$plant_rows <- 0
    log_msg(sprintf("%s %d: no records (%s s)", S, Y, secs))
    log_req(row); n_done <- n_done + 1; Sys.sleep(dl$sleep_sec %||% 2); next
  }
  d <- nf_read_npn_csv(tmp); unlink(tmp)
  n_raw <- nrow(d)
  if ("kingdom" %in% names(d)) d <- d[d$kingdom %in% "Plantae", , drop = FALSE]
  if ("species_id" %in% names(d)) d <- d[d$species_id %in% plants$species_id, , drop = FALSE]
  d <- d[, setdiff(names(d), drop_cols), drop = FALSE]
  if (!checked_fields && n_raw > 0) {
    got <- tolower(af) %in% names(d)
    log_msg("returned columns: ", paste(names(d), collapse = ", "))
    if (any(!got)) log_warn("additional fields not returned (check names in config): ", paste(af[!got], collapse = ", "))
    checked_fields <- TRUE
  }
  if (!nrow(d)) {
    writeLines(format(Sys.time()), empty)
    row$status <- "empty"; row$rows <- n_raw; row$plant_rows <- 0
    log_msg(sprintf("%s %d: %d records, none plants (%s s)", S, Y, n_raw, secs))
    log_req(row); n_done <- n_done + 1; Sys.sleep(dl$sleep_sec %||% 2); next
  }
  w <- nf_write_table(d, base)
  if (file.exists(empty)) file.remove(empty)
  row$status <- "ok"; row$rows <- n_raw; row$plant_rows <- nrow(d)
  row$sites <- length(unique(d$site_id)); row$individuals <- length(unique(d$individual_id))
  row$species <- length(unique(d$species_id)); row$phenophases <- length(unique(d$phenophase_id))
  row$yes_share <- round(mean(d$phenophase_status == 1, na.rm = TRUE), 3)
  row$with_intensity <- if ("intensity_value" %in% names(d)) round(mean(!is.na(d$intensity_value)), 3) else NA
  log_msg(sprintf("%s %d: %8d plant records | %4d sites | %5d plants | %3d species | yes %.2f | %s s | %s",
                  S, Y, nrow(d), row$sites, row$individuals, row$species, row$yes_share, secs, basename(w)))
  log_req(row); n_done <- n_done + 1
  Sys.sleep(dl$sleep_sec %||% 2)
}

# ---- summary (latest row per state-year) -------------------------------------------------------
if (file.exists(req_file)) {
  R <- utils::read.csv(req_file, stringsAsFactors = FALSE)
  R <- R[order(R$time), ]; R <- R[!duplicated(paste(R$state, R$year), fromLast = TRUE), ]
  R <- R[R$year %in% seq(nf$years$start, nf$years$end) & R$state %in% unlist(nf$states), ]
  yrs <- sort(unique(R$year))
  Ssum <- do.call(rbind, lapply(yrs, function(y) {
    r <- R[R$year == y, ]
    data.frame(year = y, requests = nrow(r), ok = sum(r$status == "ok"), empty = sum(r$status == "empty"),
               failed = sum(r$status == "error"), plant_rows = sum(r$plant_rows, na.rm = TRUE),
               sites = sum(r$sites, na.rm = TRUE), individuals = sum(r$individuals, na.rm = TRUE))
  }))
  if (!is.null(Ssum)) {
    utils::write.csv(Ssum, file.path(out, "download_summary.csv"), row.names = FALSE)
    for (i in seq_len(nrow(Ssum))) with(Ssum[i, ], log_msg(sprintf(
      "  %d: %2d states with data, %2d empty, %2d failed | %9d plant records | %5d sites | %6d plants",
      year, ok, empty, failed, plant_rows, sites, individuals)))
  }
}
mins <- round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1)
log_msg(sprintf("=== done | %d requested, %d skipped (on disk), %d failed | %s min | raw: %s | log: %s",
                n_done + n_err, n_skip, n_err, mins, file.path(raw, "status"), log_file))
if (n_err > 0) log_warn("some requests failed: rerun the same command to retry only those")
