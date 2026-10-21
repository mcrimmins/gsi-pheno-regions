#!/usr/bin/env Rscript
# Step 58d (NPN grass curing study; laptop with internet): fuel-relevant plant traits per species
# from the USDA PLANTS database: accepted symbol and name, mature height, growth form (bunch, sod,
# rhizomatous ...), growth rate, active growth period, fire tolerance and fire resistance.
# Fire tolerance / resistance in PLANTS rate how well the plant survives fire, not how much fuel or
# fire risk it adds: reported as is, not used as a fuel rating.
# Steps
#   1. PLANTS complete checklist (config npn_study: plants_list_url) -> accepted symbol for each
#      species (synonyms followed to the accepted name, e.g. Pennisetum ciliare -> Cenchrus ciliaris)
#   2. PLANTS web service (plants_api): <api>/PlantProfile?symbol=<sym> -> plant Id, durations,
#      growth habits; <api>/PlantCharacteristics/<Id> -> characteristics (height, growth form ...).
#      PLANTS publishes characteristics for about 2,200 species; others get status "no characteristics".
#   3. Manual values in data/npn_species_traits_manual.csv (if present; same columns) override.
# Output data/npn_species_traits.csv (in git). The first profile and characteristics replies are
# saved to logs/ so field names can be checked if PLANTS changes its service.
#
# Usage: Rscript npn_study/scripts/58d_npn_species_traits.R [profile]
# Needs: jsonlite

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("needs jsonlite: renv::install('jsonlite')")

cfg <- load_npn_cfg(profile)
log_file <- log_init("58d_npn_species_traits", profile)
ns <- cfg$npn_study %||% list()
list_url <- unique(c(unlist(ns$plants_list_url), "https://plants.sc.egov.usda.gov/DocumentLibrary/Txt/plantlst.txt",
                     "https://plants.usda.gov/DocumentLibrary/Txt/plantlst.txt"))
api <- ns$plants_api %||% "https://plantsservices.sc.egov.usda.gov/api"
log_msg("=== 58d_npn_species_traits | ", api)

sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)

# ---- 1. accepted symbols from the checklist ------------------------------------------------------
cache <- file.path(cfg$run_dir, "npn_study", "plants"); dir.create(cache, showWarnings = FALSE, recursive = TRUE)
lf <- file.path(cache, "plantlst.txt")
is_list <- function(f) file.exists(f) && grepl("Symbol", readLines(f, n = 1, warn = FALSE)[1])   # not an HTML page
if (!is_list(lf)) {
  unlink(lf)
  for (u in list_url) {
    ok <- tryCatch({ utils::download.file(u, lf, mode = "wb", quiet = TRUE); is_list(lf) }, error = function(e) FALSE)
    if (ok) { log_msg("checklist: ", u); break }
    log_warn("not the checklist: ", u); unlink(lf)
  }
  if (!is_list(lf)) stop("could not download the PLANTS checklist; save it by hand as ", lf, " (from plants.usda.gov, Downloads)")
}
PL <- utils::read.csv(lf, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
names(PL) <- gsub("[^a-z]", "_", tolower(names(PL)))                     # symbol, synonym_symbol, scientific_name_with_author ...
sci <- names(PL)[grepl("^scientific", names(PL))][1]; syn <- names(PL)[grepl("^synonym", names(PL))][1]
PL$binom <- vapply(strsplit(PL[[sci]], " "), function(z) paste(z[1:min(2, length(z))], collapse = " "), "")
PL$trinom <- sub("^(\\S+ \\S+)(?: .*? (?:ssp\\.|var\\.|subsp\\.) (\\S+)).*$", "\\1 \\2", PL[[sci]], perl = TRUE)
acc_name <- stats::setNames(PL[[sci]][!nzchar(PL[[syn]])], PL$symbol[!nzchar(PL[[syn]])])
lookup <- function(g, s) {
  nm <- if (grepl("-", s)) paste(g, sub("-", " ", s)) else paste(g, s)      # e.g. australis-americanus -> trinomial
  hit <- PL[(if (grepl("-", s)) PL$trinom else PL$binom) == nm, ]
  if (!nrow(hit)) return(c(symbol = NA, how = "not in checklist"))
  if (!grepl("-", s)) { sp_only <- !grepl(" (ssp\\.|subsp\\.|var\\.) ", hit[[sci]]); if (any(sp_only)) hit <- hit[sp_only, ] }
  a <- hit[!nzchar(hit[[syn]]), ]
  if (nrow(a)) return(c(symbol = a$symbol[1], how = "accepted"))
  c(symbol = hit$symbol[1], how = paste("synonym", hit[[syn]][1]))           # synonym rows: Symbol = accepted symbol
}
L <- t(mapply(lookup, sp$genus, sp$species)); sp$usda_symbol <- L[, "symbol"]; sp$symbol_match <- L[, "how"]
sp$accepted_name <- unname(acc_name[sp$usda_symbol])
log_msg(sprintf("symbols: %d of %d species matched (%d through a synonym)", sum(!is.na(sp$usda_symbol)), nrow(sp), sum(grepl("^synonym", sp$symbol_match))))

# ---- 2. profiles and characteristics -------------------------------------------------------------
get_json <- function(u) tryCatch(jsonlite::fromJSON(u, simplifyVector = TRUE), error = function(e) { log_warn(u, ": ", conditionMessage(e)); NULL })
first <- TRUE
want <- c(height_ft = "Height, Mature \\(feet\\)", growth_form = "^Growth Form$", growth_rate = "^Growth Rate$",
          active_growth = "^Active Growth Period$", lifespan = "^Lifespan$", fire_tolerance = "^Fire Tolerance$",
          fire_resistance = "^Fire Resistance$", fall_conspicuous = "^Fall Conspicuous$")
out <- lapply(seq_len(nrow(sp)), function(i) {
  r <- c(list(species_id = sp$species_id[i], usda_symbol = sp$usda_symbol[i], accepted_name = sp$accepted_name[i],
              symbol_match = sp$symbol_match[i], growth_habit = NA, plants_duration = NA, status = "no symbol"),
         stats::setNames(as.list(rep(NA, length(want))), names(want)))
  if (is.na(sp$usda_symbol[i])) return(r)
  pr <- get_json(sprintf("%s/PlantProfile?symbol=%s", api, sp$usda_symbol[i]))
  if (first && !is.null(pr)) writeLines(jsonlite::toJSON(pr, auto_unbox = TRUE, pretty = TRUE), npn_path("logs", "58d_example_profile.json"))
  if (is.null(pr) || is.null(pr$Id)) { r$status <- "no profile"; return(r) }
  r$growth_habit <- paste(unlist(pr$GrowthHabits), collapse = ", "); r$plants_duration <- paste(unlist(pr$Durations), collapse = ", ")
  ch <- get_json(sprintf("%s/PlantCharacteristics/%s", api, pr$Id))
  if (first && !is.null(ch)) { writeLines(jsonlite::toJSON(ch, auto_unbox = TRUE, pretty = TRUE), npn_path("logs", "58d_example_characteristics.json")); first <<- FALSE }
  if (is.null(ch) || !length(ch) || !"PlantCharacteristicName" %in% names(ch)) { r$status <- "no characteristics"; return(r) }
  for (k in names(want)) { v <- ch$PlantCharacteristicValue[grepl(want[[k]], ch$PlantCharacteristicName)]; if (length(v)) r[[k]] <- v[1] }
  r$status <- "characteristics"; Sys.sleep(0.3); r
})
TR <- do.call(rbind, lapply(out, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
TR$height_ft <- suppressWarnings(as.numeric(TR$height_ft))
TR$source <- ifelse(TR$status == "characteristics", "USDA PLANTS", NA)

# ---- 3. manual overrides ------------------------------------------------------------------------
mf <- npn_path("data", "npn_species_traits_manual.csv")
if (file.exists(mf)) {
  M <- utils::read.csv(mf, stringsAsFactors = FALSE)
  for (j in seq_len(nrow(M))) { i <- match(M$species_id[j], TR$species_id); if (is.na(i)) next
    for (k in intersect(setdiff(names(M), "species_id"), names(TR))) if (!is.na(M[[k]][j]) && nzchar(as.character(M[[k]][j]))) TR[[k]][i] <- M[[k]][j] }
  log_msg(sprintf("manual values for %d species from %s", nrow(M), basename(mf)))
}
utils::write.csv(TR, npn_path("data", "npn_species_traits.csv"), row.names = FALSE)
log_msg(sprintf("status: %s", paste(names(table(TR$status)), table(TR$status), collapse = ", ")))
log_msg(sprintf("height for %d species, growth form for %d", sum(!is.na(TR$height_ft)), sum(!is.na(TR$growth_form))))
log_msg("=== done | data/npn_species_traits.csv | log: ", log_file)
