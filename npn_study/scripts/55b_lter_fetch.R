#!/usr/bin/env Rscript
# Step 55b (NPN grass curing study; laptop with internet): fetch LTER plant phenology data
# packages from EDI for the independent cross-checks, and write a structure report so the
# comparison scripts can be built against the real columns and codes.
# Packages (config npn_study: lter_packages, each with scope, id and revision; 23 Oct 2026: EDI
# refuses the public "list revisions" call, so revisions are pinned in the config):
#   knb-lter-jrn 210287001  Plant phenology observations at 15 NPP sites, Jornada Basin LTER (1992-)
#   knb-lter-sev 137        Core Site Phenology Study, Sevilleta LTER (2000-2015)
# For each package: every data entity and the EML metadata are saved to
# <run_dir>/npn_study/lter/<scope>.<id>.<rev>/, and a report goes to
# npn_study/logs/55b_lter_structure_<profile>.txt: entities, columns, first rows, codes of
# categorical columns, date ranges, EML attribute definitions and code lists, coordinates.
# Fetching: EML from https://pasta.lternet.edu/package/metadata/eml/<scope>/<id>/<rev>, then each
# entity from the download URL in its EML. If EDI refuses that too, download the package zip by
# hand from the EDI portal ("Download Data Package") into <run_dir>/npn_study/lter/; the script
# unzips any <scope>.<id>*.zip there and reads the files from it.
# Re-running skips files already downloaded.
#
# Usage: Rscript npn_study/scripts/55b_lter_fetch.R [profile]
# Needs: xml2

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
for (p in c("xml2")) if (!requireNamespace(p, quietly = TRUE)) stop("needs ", p, ": install.packages('", p, "')")

cfg <- load_npn_cfg(profile)
log_file <- log_init("55b_lter_fetch", profile)
ns <- cfg$npn_study %||% list()
pk <- ns$lter_packages %||% list(list(scope = "knb-lter-jrn", id = 210287001, rev = 60), list(scope = "knb-lter-sev", id = 137, rev = 208261))
pasta <- ns$pasta_url %||% "https://pasta.lternet.edu/package"
get <- function(u, f) tryCatch({ utils::download.file(u, f, mode = "wb", quiet = TRUE); file.exists(f) && file.size(f) > 0 },
                               error = function(e) { log_warn(basename(f), ": ", conditionMessage(e)); unlink(f); FALSE },
                               warning = function(w) { log_warn(basename(f), ": ", conditionMessage(w)); unlink(f); FALSE })
root <- file.path(cfg$run_dir, "npn_study", "lter"); dir.create(root, showWarnings = FALSE, recursive = TRUE)
rep_file <- npn_path("logs", sprintf("55b_lter_structure_%s.txt", profile))
rep <- function(...) cat(paste0(...), "\n", file = rep_file, append = TRUE, sep = "")
cat("", file = rep_file)
log_msg("=== 55b_lter_fetch | ", length(pk), " packages -> ", root)

safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
for (p in pk) {
  pid <- paste(p$scope, p$id, p$rev, sep = ".")
  pdir <- file.path(root, pid); dir.create(pdir, showWarnings = FALSE)
  log_msg("package ", pid)
  rep("\n################ ", pid, " ################")
  # a hand-downloaded package zip, if any
  zz <- list.files(root, "\\.zip$", full.names = TRUE); zz <- zz[startsWith(basename(zz), paste0(p$scope, ".", p$id))]
  if (length(zz)) { utils::unzip(zz[1], exdir = pdir); log_msg("  unzipped ", basename(zz[1])) }
  # metadata: from the zip (<pid>.xml) or from PASTA
  ef <- file.path(pdir, "eml.xml")
  zx <- list.files(pdir, "\\.xml$", full.names = TRUE); zx <- zx[!grepl("report\\.xml$|^eml\\.xml$", basename(zx))]   # the EML, not EDI's quality report
  if (length(zx)) file.copy(zx[1], ef, overwrite = TRUE) else if (!file.exists(ef)) get(sprintf("%s/metadata/eml/%s/%s/%s", pasta, p$scope, p$id, p$rev), ef)
  if (!file.exists(ef)) {
    log_err(pid, ": no metadata. Download the package zip from https://portal.edirepository.org/nis/mapbrowse?scope=", p$scope,
            "&identifier=", p$id, " into ", root, " and run again")
    rep("NO METADATA: download the zip by hand"); next
  }
  eml <- xml2::read_xml(ef); xml2::xml_ns_strip(eml)
  rep("TITLE: ", xml2::xml_text(xml2::xml_find_first(eml, ".//dataset/title")))
  rep("TEMPORAL: ", paste(xml2::xml_text(xml2::xml_find_all(eml, ".//temporalCoverage//calendarDate")), collapse = " .. "))
  meth <- xml2::xml_text(xml2::xml_find_all(eml, ".//methods//para"))
  if (length(meth)) rep("METHODS (first 6000 chars):\n", substr(paste(meth, collapse = "\n"), 1, 6000))
  gc <- xml2::xml_find_all(eml, ".//geographicCoverage")
  rep("GEOGRAPHIC COVERAGE (", length(gc), "):")
  for (g in gc) rep("  ", gsub("\\s+", " ", xml2::xml_text(xml2::xml_find_first(g, "geographicDescription"))), " | W ",
                    xml2::xml_text(xml2::xml_find_first(g, ".//westBoundingCoordinate")), " E ", xml2::xml_text(xml2::xml_find_first(g, ".//eastBoundingCoordinate")),
                    " N ", xml2::xml_text(xml2::xml_find_first(g, ".//northBoundingCoordinate")), " S ", xml2::xml_text(xml2::xml_find_first(g, ".//southBoundingCoordinate")))
  # attribute definitions and code lists
  for (dt in xml2::xml_find_all(eml, ".//dataTable|.//otherEntity")) {
    rep("\nENTITY (EML): ", xml2::xml_text(xml2::xml_find_first(dt, "entityName")), " | ", xml2::xml_text(xml2::xml_find_first(dt, ".//objectName")))
    for (a in xml2::xml_find_all(dt, ".//attribute")) {
      rep("  ", xml2::xml_text(xml2::xml_find_first(a, "attributeName")), ": ", gsub("\\s+", " ", xml2::xml_text(xml2::xml_find_first(a, "attributeDefinition"))))
      cd <- xml2::xml_find_all(a, ".//codeDefinition")
      for (cc in cd) rep("      ", xml2::xml_text(xml2::xml_find_first(cc, "code")), " = ", gsub("\\s+", " ", xml2::xml_text(xml2::xml_find_first(cc, "definition"))))
    }
  }
  # data entities: object name and download URL from the EML
  for (dt in xml2::xml_find_all(eml, ".//dataTable|.//otherEntity")) {
    nm <- xml2::xml_text(xml2::xml_find_first(dt, "entityName"))
    on <- xml2::xml_text(xml2::xml_find_first(dt, ".//physical/objectName")); if (is.na(on)) on <- paste0(safe(nm), ".csv")
    u <- xml2::xml_text(xml2::xml_find_first(dt, ".//distribution/online/url"))
    f <- file.path(pdir, on)
    if (!file.exists(f)) {                                                   # maybe in a subfolder of the zip
      hit <- list.files(pdir, recursive = TRUE, full.names = TRUE); hit <- hit[basename(hit) == on]
      if (length(hit)) f <- hit[1] else if (!is.na(u)) get(u, f)
    }
    if (!file.exists(f)) { log_err("  ", on, ": not downloaded (", u, ")"); rep("\n=== DATA ENTITY: ", nm, " -> NOT DOWNLOADED"); next }
    log_msg(sprintf("  %s (%.1f MB)", basename(f), file.size(f) / 1e6))
    rep("\n=== DATA ENTITY: ", nm, " -> ", basename(f), sprintf(" (%.1f MB)", file.size(f) / 1e6))
    if (!grepl("\\.(csv|txt|tsv)$", f, ignore.case = TRUE)) { rep("  (not a table; skipped)"); next }
    sep <- if (grepl("\\t", readLines(f, n = 1, warn = FALSE))) "\t" else ","
    d <- tryCatch(utils::read.csv(f, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA", ".", "-9999", "NaN")),
                  error = function(e) { rep("  read error: ", conditionMessage(e)); NULL })
    if (is.null(d)) next
    rep(sprintf("  %d rows x %d columns", nrow(d), ncol(d)))
    rep("  HEAD:"); rep(paste(utils::capture.output(print(utils::head(d, 6))), collapse = "\n"))
    for (k in names(d)) {
      x <- d[[k]]; ux <- unique(x[!is.na(x)])
      if (is.numeric(x)) rep(sprintf("  [num] %-28s range %s .. %s, %d unique, %d NA", k, format(min(x, na.rm = TRUE)), format(max(x, na.rm = TRUE)), length(ux), sum(is.na(x))))
      else if (length(ux) <= 80) {
        tb <- sort(table(x), decreasing = TRUE)
        rep(sprintf("  [chr] %-28s %d values: %s", k, length(ux), paste(sprintf("%s (%d)", names(tb), as.integer(tb)), collapse = ", ")))
      } else rep(sprintf("  [chr] %-28s %d unique; e.g. %s", k, length(ux), paste(utils::head(ux, 12), collapse = " | ")))
    }
  }
}
log_msg("=== done | report: ", rep_file, " | log: ", log_file)
