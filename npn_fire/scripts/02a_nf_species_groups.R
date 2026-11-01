#!/usr/bin/env Rscript
# Step 02a (NPN x fire study; laptop with internet): plant group for every NPN plant species.
# NPN functional types don't separate trees from shrubs (e.g. "Deciduous broadleaf" holds oaks and
# sumacs), so woody broadleaf species take their growth habit from USDA PLANTS.
#   1. NPN species lookup from 01 (<data_root>/_shared/raw/npn/lookups/species.csv), plants only
#   2. PLANTS checklist (accepted symbol; synonyms followed) and PLANTS profile -> growth habits,
#      duration. Cached in <data_root>/_shared/raw/plants/ (profiles as one .rds, resumable).
#   3. Group rule (config npn_fire: groups):
#        Graminoid -> grass; Forb, Semi-evergreen forb, Evergreen forb -> forb;
#        Evergreen conifer, Pine, Deciduous conifer -> conifer; Cactus -> cactus;
#        broadleaf types -> <leaf>_<form>, leaf = deciduous (Deciduous broadleaf, Drought deciduous
#        broadleaf) or evergreen (Evergreen broadleaf, Semi-evergreen broadleaf); form from PLANTS
#        habits: Tree only -> tree; Shrub / Subshrub only -> shrub; both -> tree_shrub;
#        Vine only -> vine; no PLANTS habit -> unknown
#   3b. "Shrub, Tree" broadleaf species: PLANTS mature height (PlantCharacteristics) >= height_tree_ft
#      (20) -> tree, below -> shrub, none -> stays tree_shrub (cached in heights.rds)
#   4. Manual values in data/nf_species_groups_manual.csv (species_id, group[, note]) override.
# Output data/nf_species_groups.csv (in git): species_id, genus, species, common_name,
# functional_type, usda_symbol, symbol_match, growth_habit, plants_duration, group, group_rule.
#
# Usage: Rscript npn_fire/scripts/02a_nf_species_groups.R [profile] [force]
# Needs: jsonlite. ~1,500 PLANTS profile calls the first time (~10 min); cached after.

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1 && args[1] != "force") args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
force <- "force" %in% args
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("needs jsonlite: renv::install('jsonlite')")

cfg <- load_nf_cfg(profile)
log_file <- log_init("02a_nf_species_groups", profile)
gc <- cfg$npn_fire$groups
api <- gc$plants_api
log_msg("=== 02a_nf_species_groups | profile: ", profile, " | ", api)

sf <- file.path(cfg$raw_dir, "npn", "lookups", "species.csv")
if (!file.exists(sf)) stop("no NPN species lookup at ", sf, ": run 01 first")
sp <- utils::read.csv(sf, stringsAsFactors = FALSE)
sp <- sp[sp$kingdom == "Plantae" & !(sp$functional_type %in% unlist(cfg$npn_fire$download$drop_functional_types)),
         c("species_id", "genus", "species", "common_name", "functional_type")]
log_msg(sprintf("NPN plant species: %d", nrow(sp)))

# ---- PLANTS checklist -> accepted symbols ------------------------------------------------------
pdir <- nf_raw(cfg, "plants"); lf <- file.path(pdir, "plantlst.txt")
is_list <- function(f) file.exists(f) && grepl("Symbol", readLines(f, n = 1, warn = FALSE)[1])
if (!is_list(lf)) {
  alt <- file.path(cfg$run_dir, "npn_study", "plants", "plantlst.txt")        # the NPN study's copy, if present
  if (is_list(alt)) { file.copy(alt, lf, overwrite = TRUE); log_msg("checklist: copied from ", alt) }
}
if (!is_list(lf)) {
  unlink(lf)
  for (u in unlist(gc$plants_list_url)) {
    ok <- tryCatch({ utils::download.file(u, lf, mode = "wb", quiet = TRUE); is_list(lf) }, error = function(e) FALSE)
    if (ok) { log_msg("checklist: ", u); break }
    log_warn("not the checklist: ", u); unlink(lf)
  }
  if (!is_list(lf)) stop("could not download the PLANTS checklist; save it by hand as ", lf)
}
PL <- utils::read.csv(lf, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
names(PL) <- gsub("[^a-z]", "_", tolower(names(PL)))
sci <- names(PL)[grepl("^scientific", names(PL))][1]; syn <- names(PL)[grepl("^synonym", names(PL))][1]
PL$binom <- vapply(strsplit(PL[[sci]], " "), function(z) paste(z[1:min(2, length(z))], collapse = " "), "")
PL$trinom <- sub("^(\\S+ \\S+)(?: .*? (?:ssp\\.|var\\.|subsp\\.) (\\S+)).*$", "\\1 \\2", PL[[sci]], perl = TRUE)
lookup <- function(g, s) {
  if (is.na(s) || !nzchar(s) || grepl("^sp\\.?$|^spp\\.?$", s)) return(c(symbol = NA, how = "genus only"))
  tri <- grepl("-", s)
  nm <- if (tri) paste(g, sub("-", " ", s)) else paste(g, s)                 # e.g. australis-americanus -> trinomial
  hit <- PL[(if (tri) PL$trinom else PL$binom) == nm, ]
  if (!nrow(hit)) return(c(symbol = NA, how = "not in checklist"))
  if (!tri) { sp_only <- !grepl(" (ssp\\.|subsp\\.|var\\.) ", hit[[sci]]); if (any(sp_only)) hit <- hit[sp_only, ] }
  a <- hit[!nzchar(hit[[syn]]), ]
  if (nrow(a)) return(c(symbol = a$symbol[1], how = "accepted"))
  c(symbol = hit$symbol[1], how = "synonym")                                  # synonym rows: Symbol = accepted symbol
}
L <- t(mapply(lookup, sp$genus, sp$species))
sp$usda_symbol <- unname(L[, "symbol"]); sp$symbol_match <- unname(L[, "how"])
log_msg(sprintf("PLANTS symbols: %d of %d (%d through a synonym)", sum(!is.na(sp$usda_symbol)), nrow(sp),
                sum(sp$symbol_match == "synonym")))

# ---- PLANTS profiles (growth habit), cached -----------------------------------------------------
cache <- file.path(pdir, "profiles.rds")
prof <- if (file.exists(cache) && !force) readRDS(cache) else list()
syms <- setdiff(unique(stats::na.omit(sp$usda_symbol)), names(prof))
log_msg(sprintf("PLANTS profiles: %d cached, %d to fetch", length(prof), length(syms)))
n_fail <- 0
for (i in seq_along(syms)) {
  s <- syms[i]
  pr <- tryCatch(jsonlite::fromJSON(sprintf("%s/PlantProfile?symbol=%s", api, s), simplifyVector = TRUE),
                 error = function(e) { n_fail <<- n_fail + 1; NULL })
  if (!is.null(pr)) prof[[s]] <- list(habit = paste(unlist(pr$GrowthHabits), collapse = ", "),
                                      duration = paste(unlist(pr$Durations), collapse = ", "),
                                      id = if (length(pr$Id)) pr$Id[1] else NA)
  if (i %% 100 == 0) { saveRDS(prof, cache); log_msg(sprintf("  %d / %d profiles (%d failed)", i, length(syms), n_fail)) }
  if (n_fail >= 25 && n_fail == i) stop("PLANTS service not answering (", n_fail, " failures); try later")
  Sys.sleep(gc$plants_sleep_sec %||% 0.3)
}
saveRDS(prof, cache)
if (n_fail) log_warn(n_fail, " profile requests failed (rerun to retry)")
sp$growth_habit <- vapply(sp$usda_symbol, function(s) if (!is.na(s) && !is.null(prof[[s]])) prof[[s]]$habit else NA_character_, "")
sp$plants_duration <- vapply(sp$usda_symbol, function(s) if (!is.na(s) && !is.null(prof[[s]])) prof[[s]]$duration else NA_character_, "")

# ---- group rule ----------------------------------------------------------------------------------
ft <- sp$functional_type; h <- tolower(ifelse(is.na(sp$growth_habit), "", sp$growth_habit))
has <- function(x) grepl(x, h)
form <- ifelse(has("tree") & (has("shrub") | has("subshrub")), "tree_shrub",
        ifelse(has("tree"), "tree",
        ifelse(has("shrub") | has("subshrub"), "shrub",
        ifelse(has("vine"), "vine", "unknown"))))
leaf <- ifelse(ft %in% unlist(gc$deciduous_broadleaf), "deciduous",
        ifelse(ft %in% unlist(gc$evergreen_broadleaf), "evergreen", NA))

# ---- mature height for broadleaf species PLANTS lists as both tree and shrub (decided 1 Nov) -------
# PLANTS "Shrub, Tree" covers sugar maple as well as manzanita, so the habit alone can't split them:
# height >= height_tree_ft -> tree, below -> shrub, no height -> stays tree_shrub (manual list).
# PlantCharacteristics needs the PLANTS record id; profiles cached before ids were kept are refetched.
cut_ft <- gc$height_tree_ft %||% 20
need <- !is.na(leaf) & form == "tree_shrub" & !is.na(sp$usda_symbol)
hsyms <- unique(sp$usda_symbol[need])
ccache <- file.path(pdir, "heights.rds")
hts <- if (file.exists(ccache) && !force) readRDS(ccache) else list()
todo <- setdiff(hsyms, names(hts))
log_msg(sprintf("mature height: %d tree_shrub symbols, %d cached, %d to fetch (cutoff %s ft)", length(hsyms),
                length(hsyms) - length(todo), length(todo), cut_ft))
h_fail <- 0
for (i in seq_along(todo)) {
  s <- todo[i]
  id <- prof[[s]]$id %||% NA
  if (is.na(id)) {
    pr <- tryCatch(jsonlite::fromJSON(sprintf("%s/PlantProfile?symbol=%s", api, s), simplifyVector = TRUE), error = function(e) NULL)
    if (!is.null(pr) && length(pr$Id)) { id <- pr$Id[1]; prof[[s]]$id <- id }
    Sys.sleep(gc$plants_sleep_sec %||% 0.3)
  }
  if (is.na(id)) { h_fail <- h_fail + 1; next }                      # not cached: retried on the next run
  ch <- tryCatch(jsonlite::fromJSON(sprintf("%s/PlantCharacteristics/%s", api, id), simplifyVector = TRUE),
                 error = function(e) { h_fail <<- h_fail + 1; NULL })
  if (!is.null(ch)) {
    v <- if (is.data.frame(ch) && all(c("PlantCharacteristicName", "PlantCharacteristicValue") %in% names(ch)))
      ch$PlantCharacteristicValue[grepl("^Height, Mature \\(feet\\)", ch$PlantCharacteristicName)] else character()
    hts[[s]] <- if (length(v)) suppressWarnings(as.numeric(v[1])) else NA_real_   # NA = PLANTS has no height
  }
  if (i %% 50 == 0) { saveRDS(hts, ccache); saveRDS(prof, cache); log_msg(sprintf("  %d / %d heights", i, length(todo))) }
  Sys.sleep(gc$plants_sleep_sec %||% 0.3)
}
saveRDS(hts, ccache); saveRDS(prof, cache)
if (h_fail) log_warn(h_fail, " height lookups failed (rerun to retry)")
sp$height_ft <- vapply(sp$usda_symbol, function(s) if (!is.na(s) && !is.null(hts[[s]])) hts[[s]] else NA_real_, 0)
by_h <- need & !is.na(sp$height_ft)
form[by_h] <- ifelse(sp$height_ft[by_h] >= cut_ft, "tree", "shrub")
log_msg(sprintf("tree_shrub split by height: %d tree, %d shrub, %d without a PLANTS height (stay tree_shrub)",
                sum(by_h & form == "tree"), sum(by_h & form == "shrub"), sum(need & !by_h)))

sp$group <- ifelse(ft %in% unlist(gc$grass), "grass",
            ifelse(ft %in% unlist(gc$forb), "forb",
            ifelse(ft %in% unlist(gc$conifer), "conifer",
            ifelse(ft %in% unlist(gc$cactus), "cactus",
            ifelse(!is.na(leaf), paste(leaf, form, sep = "_"), "other")))))
sp$group_rule <- ifelse(!is.na(leaf), paste0("NPN ", ft, " + PLANTS habit (", ifelse(nzchar(h), sp$growth_habit, "none"), ")",
                                              ifelse(by_h, sprintf(" + height %g ft", sp$height_ft), "")),
                        paste0("NPN ", ft))

mf <- nf_path("data", "nf_species_groups_manual.csv")
if (file.exists(mf)) {
  M <- utils::read.csv(mf, stringsAsFactors = FALSE)
  i <- match(M$species_id, sp$species_id); ok <- !is.na(i) & nzchar(M$group)
  sp$group[i[ok]] <- M$group[ok]; sp$group_rule[i[ok]] <- "manual"
  log_msg(sprintf("manual groups for %d species (%s)", sum(ok), basename(mf)))
}
sp <- sp[order(sp$group, sp$genus, sp$species), ]
dir.create(nf_path("data"), showWarnings = FALSE)
utils::write.csv(sp, nf_path("data", "nf_species_groups.csv"), row.names = FALSE)
# species still needing a hand decision (tree_shrub without a height; unknown): fill a group into
# data/nf_species_groups_manual.csv (species_id, group, note) and rerun
rv <- sp[grepl("_(tree_shrub|unknown)$", sp$group), c("species_id", "genus", "species", "common_name", "functional_type",
                                                      "symbol_match", "growth_habit", "group")]
utils::write.csv(rv, nf_path("data", "nf_species_groups_review.csv"), row.names = FALSE)
log_msg(sprintf("for review: %d species -> npn_fire/data/nf_species_groups_review.csv", nrow(rv)))
tb <- table(sp$group)
for (g in names(tb)) log_msg(sprintf("  %-22s %4d species", g, tb[[g]]))
log_msg("=== done | npn_fire/data/nf_species_groups.csv | log: ", log_file)
