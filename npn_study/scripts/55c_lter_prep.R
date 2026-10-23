#!/usr/bin/env Rscript
# Step 55c (NPN grass curing study; laptop): tidy the LTER phenology records (55b) into monthly
# grass tables on a common footing with the NPN record, for the independent cross-checks (55e).
# Grasses only (family Poaceae in the USDA PLANTS checklist cached by 58d); grass type and life
# cycle from data/lter_grass_pathway.csv (species rows override genus rows).
#   Jornada (knb-lter-jrn.210287001): monthly counts per NPP site and species of dormant, non-
#     reproductive, bud, flower and fruit plants. green = not dormant. Rows with a missing (M)
#     quality flag are dropped.
#   Sevilleta (knb-lter-sev.137): monthly, the first 10 plants per species per transect (web);
#     foliage N (new green), O (old green), B (brown), Z (no leaves). green = N or O,
#     new = N, brown = B or Z. Webs are pooled per core site.
# Both give, per site x species x month, the SHARE OF PLANTS green, not the percent green of a
# plant. The NPN equivalent is built from the NPN readings (usable plant-years, relaxed set):
# per NPN site x species x month, the share of readings with any green (>= 5 % class) and with
# >= 50 % green, and the mean percent green.
# Outputs <run_dir>/npn_study/lter/
#   lter_grass_monthly.csv   source, site, species, name, pathway, life_cycle, year, month, date
#                            (median sampling date), doy, n, green, new, brown (shares)
#   lter_species.csv         grass species per source: name, pathway, life cycle, sites, years, n
#   lter_sites.csv           source, site, label, lat, lon (EML bounding-box centres; the 15
#                            Jornada NPP sites share the basin box) -> PRISM extraction (55d)
#   npn_monthly.csv          NPN site x species x year x month: n, share_any, share50, mean_pct,
#                            lat, lon, pathway, neon
#
# Usage: Rscript npn_study/scripts/55c_lter_prep.R [profile]   (after 55b and 58d)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("xml2", quietly = TRUE)) stop("needs xml2")

cfg <- load_npn_cfg(profile)
log_file <- log_init("55c_lter_prep", profile)
ns <- cfg$npn_study %||% list()
ldir <- file.path(cfg$run_dir, "npn_study", "lter")
log_msg("=== 55c_lter_prep | ", ldir)
pk <- ns$lter_packages %||% list(list(scope = "knb-lter-jrn", id = 210287001, rev = 60), list(scope = "knb-lter-sev", id = 137, rev = 208261))
pdir <- function(scope) { d <- list.dirs(ldir, recursive = FALSE); d <- d[startsWith(basename(d), scope)]; if (!length(d)) stop("no ", scope, " folder in ", ldir, " (run 55b)"); d[1] }
find1 <- function(d, pat) { f <- list.files(d, pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE); if (!length(f)) stop("no ", pat, " in ", d); f[1] }

# ---- grass names and types -----------------------------------------------------------------------
pl_f <- file.path(cfg$run_dir, "npn_study", "plants", "plantlst.txt")
if (!file.exists(pl_f)) stop("no PLANTS checklist ", pl_f, " (run 58d first)")
PL <- utils::read.csv(pl_f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
names(PL) <- c("symbol", "synonym", "sci", "common", "family")[seq_len(ncol(PL))]
acc <- PL[!nzchar(PL$synonym), ]
syn <- PL[nzchar(PL$synonym), ]
sym_acc <- function(code) { a <- code; i <- !(a %in% acc$symbol) & a %in% syn$synonym; a[i] <- syn$symbol[match(a[i], syn$synonym)]; a }
binom <- function(s) vapply(strsplit(s, " "), function(z) paste(z[1:min(2, length(z))], collapse = " "), "")
info <- function(code) {
  a <- sym_acc(code); m <- match(a, acc$symbol)
  data.frame(code = code, accepted = a, family = acc$family[m], name = binom(acc$sci[m]), common = acc$common[m], stringsAsFactors = FALSE)
}
GP <- utils::read.csv(npn_path("data", "lter_grass_pathway.csv"), stringsAsFactors = FALSE)
gtype <- function(name) {
  g <- sub(" .*", "", name)
  sp <- GP[GP$level == "species", ]; ge <- GP[GP$level == "genus", ]
  pw <- ifelse(name %in% sp$name, sp$pathway[match(name, sp$name)], ge$pathway[match(g, ge$name)])
  lc <- ifelse(name %in% sp$name, sp$life_cycle[match(name, sp$name)], ge$life_cycle[match(g, ge$name)])
  data.frame(pathway = pw, life_cycle = lc, stringsAsFactors = FALSE)
}

# ---- Sevilleta ------------------------------------------------------------------------------------
sd <- pdir("knb-lter-sev.137")
S <- utils::read.csv(find1(sd, "sev137.*phenology.*\\.csv$"), stringsAsFactors = FALSE)
S$date <- as.Date(S$date); S <- S[!is.na(S$date) & !is.na(S$month) & S$foliage %in% c("N", "O", "B", "Z"), ]
dm <- as.integer(format(S$date, "%m")); S$year <- as.integer(format(S$date, "%Y"))
S$year[S$month - dm > 6] <- S$year[S$month - dm > 6] - 1                   # e.g. a December round sampled on Jan 2
S$year[dm - S$month > 6] <- S$year[dm - S$month > 6] + 1
IS <- info(unique(S$species)); IS <- IS[IS$family %in% "Poaceae" & !is.na(IS$name), ]
S <- S[S$species %in% IS$code, ]
S$key <- paste(S$site, S$species, S$year, S$month)
agg <- function(x, f) tapply(x, S$key, f)
SM <- data.frame(key = names(agg(S$foliage, length)), n = as.integer(agg(S$foliage, length)),
                 green = as.numeric(agg(S$foliage %in% c("N", "O"), mean)), new = as.numeric(agg(S$foliage == "N", mean)),
                 brown = as.numeric(agg(S$foliage %in% c("B", "Z"), mean)),
                 date = as.Date(as.numeric(agg(as.numeric(S$date), stats::median)), origin = "1970-01-01"), stringsAsFactors = FALSE)
kk <- do.call(rbind, strsplit(SM$key, " ")); SM$site <- kk[, 1]; SM$species <- kk[, 2]; SM$year <- as.integer(kk[, 3]); SM$month <- as.integer(kk[, 4])
SM$source <- "SEV"
log_msg(sprintf("Sevilleta: %d grass observations, %d species, %d site-species-months", nrow(S), length(unique(S$species)), nrow(SM)))

# ---- Jornada --------------------------------------------------------------------------------------
jd <- pdir("knb-lter-jrn.210287001")
J <- utils::read.csv(find1(jd, "NPP_Phenology\\.csv$"), stringsAsFactors = FALSE, na.strings = c("", "NA", "."))
J$date <- as.Date(J$date)
qf <- J[, grep("^QF_", names(J))]; bad <- apply(qf, 1, function(z) any(z %in% "M"))
J <- J[!bad & !is.na(J$date) & !is.na(J$total_obs) & J$total_obs > 0 & !is.na(J$dormant), ]
IJ <- info(unique(J$USDA_code)); IJ <- IJ[IJ$family %in% "Poaceae" & !is.na(IJ$name), ]
J <- J[J$USDA_code %in% IJ$code, ]
J$year <- as.integer(format(J$date, "%Y")); J$month <- as.integer(format(J$date, "%m"))
JM <- stats::aggregate(cbind(total_obs, dormant) ~ zone + site + USDA_code + year + month, J, sum)
JD <- stats::aggregate(as.numeric(date) ~ zone + site + USDA_code + year + month, J, stats::median); names(JD)[6] <- "d"
JM <- merge(JM, JD)
JM <- data.frame(source = "JRN", site = paste(JM$zone, JM$site, sep = "-"), species = JM$USDA_code, year = JM$year, month = JM$month,
                 date = as.Date(JM$d, origin = "1970-01-01"), n = JM$total_obs, green = 1 - JM$dormant / JM$total_obs,
                 new = NA_real_, brown = JM$dormant / JM$total_obs, stringsAsFactors = FALSE)
log_msg(sprintf("Jornada: %d grass rows, %d species, %d site-species-months", nrow(J), length(unique(J$USDA_code)), nrow(JM)))

# ---- combine --------------------------------------------------------------------------------------
I <- unique(rbind(IS, IJ)); I <- cbind(I, gtype(I$name))
if (any(is.na(I$pathway))) log_warn("grass type unknown for: ", paste(unique(I$name[is.na(I$pathway)]), collapse = ", "), " (add to data/lter_grass_pathway.csv)")
M <- rbind(SM[, c("source", "site", "species", "year", "month", "date", "n", "green", "new", "brown")], JM)
M <- cbind(M, I[match(M$species, I$code), c("name", "common", "pathway", "life_cycle")])
M$doy <- as.integer(format(M$date, "%j"))
M <- M[order(M$source, M$site, M$species, M$year, M$month), ]
utils::write.csv(M, file.path(ldir, "lter_grass_monthly.csv"), row.names = FALSE)
SPc <- do.call(rbind, lapply(split(M, paste(M$source, M$species)), function(d) data.frame(source = d$source[1], species = d$species[1], name = d$name[1],
         common = d$common[1], pathway = d$pathway[1], life_cycle = d$life_cycle[1], sites = paste(sort(unique(d$site)), collapse = " "),
         first_year = min(d$year), last_year = max(d$year), months = nrow(d), plants_obs = sum(d$n))))
SPc <- SPc[order(SPc$source, -SPc$plants_obs), ]
utils::write.csv(SPc, file.path(ldir, "lter_species.csv"), row.names = FALSE)
for (s in unique(SPc$source)) { x <- SPc[SPc$source == s, ]
  log_msg(sprintf("  %s: %d grasses (C3 %d, C4 %d); most observed: %s", s, nrow(x), sum(x$pathway %in% "C3"), sum(x$pathway %in% "C4"),
                  paste(utils::head(sprintf("%s (%s)", x$name, x$pathway), 6), collapse = ", "))) }

# ---- site coordinates (bounding-box centres) -----------------------------------------------------
box <- function(d) {
  e <- xml2::read_xml(list.files(d, "\\.xml$", full.names = TRUE)[!grepl("report\\.xml$|^eml\\.xml$", list.files(d, "\\.xml$"))][1]); xml2::xml_ns_strip(e)
  g <- xml2::xml_find_all(e, ".//geographicCoverage")
  v <- function(x, p) as.numeric(xml2::xml_text(xml2::xml_find_first(x, p)))
  data.frame(desc = vapply(g, function(x) gsub("\\s+", " ", xml2::xml_text(xml2::xml_find_first(x, "geographicDescription"))), ""),
             lon = vapply(g, function(x) mean(c(v(x, ".//westBoundingCoordinate"), v(x, ".//eastBoundingCoordinate"))), 0),
             lat = vapply(g, function(x) mean(c(v(x, ".//northBoundingCoordinate"), v(x, ".//southBoundingCoordinate"))), 0), stringsAsFactors = FALSE)
}
bs <- box(sd); bj <- box(jd)
sev_lab <- c(B = "Blue grama", G = "Black grama", C = "Creosote", J = "Juniper", P = "Pinon", ISO = "Seth's webs", JSAV = "Juniper")
sev_xy <- function(code) { k <- grep(paste0("^", sev_lab[[code]]), bs$desc, ignore.case = TRUE)[1]; if (is.na(k)) c(NA, NA) else c(bs$lat[k], bs$lon[k]) }
ST <- rbind(
  do.call(rbind, lapply(sort(unique(SM$site)), function(s) { xy <- sev_xy(s); data.frame(source = "SEV", site = s, label = sev_lab[[s]] %||% s, lat = xy[1], lon = xy[2]) })),
  do.call(rbind, lapply(sort(unique(JM$site)), function(s) data.frame(source = "JRN", site = s, label = s, lat = bj$lat[1], lon = bj$lon[1]))))
ST$site_id <- paste(ST$source, ST$site, sep = "_")
utils::write.csv(ST, file.path(ldir, "lter_sites.csv"), row.names = FALSE)
log_msg(sprintf("sites: %d (SEV %d, JRN %d); missing coordinates: %d", nrow(ST), sum(ST$source == "SEV"), sum(ST$source == "JRN"), sum(is.na(ST$lat))))

# ---- the NPN record as monthly shares --------------------------------------------------------------
go <- readRDS(file.path(npn_obs_dir(cfg), "npn_green_obs.rds"))
cm <- readRDS(file.path(npn_obs_dir(cfg), "npn_curing_metrics.rds"))
sp <- utils::read.csv(npn_path("data", "npn_species_pathway.csv"), stringsAsFactors = FALSE)
ok <- paste(cm$individual_id, cm$year)[cm$usable %in% TRUE & !(sp$crop[match(cm$species_id, sp$species_id)] %in% "yes")]
go <- go[paste(go$individual_id, go$year) %in% ok & !is.na(go$pct_green), ]
go$month <- as.integer(substr(go$observation_date, 6, 7))
go$key <- paste(go$site_id, go$species_id, go$year, go$month)
ag <- function(x, f) tapply(x, go$key, f)
NM <- data.frame(key = names(ag(go$pct_green, length)), n = as.integer(ag(go$pct_green, length)),
                 share_any = as.numeric(ag(go$pct_green >= 5, mean)), share50 = as.numeric(ag(go$pct_green >= 50, mean)),
                 mean_pct = as.numeric(ag(go$pct_green, mean)), stringsAsFactors = FALSE)
kk <- do.call(rbind, strsplit(NM$key, " ")); NM$site_id <- as.integer(kk[, 1]); NM$species_id <- as.integer(kk[, 2]); NM$year <- as.integer(kk[, 3]); NM$month <- as.integer(kk[, 4])
xy <- stats::aggregate(cbind(latitude, longitude) ~ site_id, go, stats::median)
NM$lat <- xy$latitude[match(NM$site_id, xy$site_id)]; NM$lon <- xy$longitude[match(NM$site_id, xy$site_id)]
NM$genus <- sp$genus[match(NM$species_id, sp$species_id)]; NM$species <- sp$species[match(NM$species_id, sp$species_id)]
NM$name <- paste(NM$genus, NM$species); NM$pathway <- sp$pathway[match(NM$species_id, sp$species_id)]
ty_f <- file.path(cfg$run_dir, "npn_study", "record", "series_typical.csv")
if (file.exists(ty_f)) { ty <- utils::read.csv(ty_f, stringsAsFactors = FALSE); NM$neon <- NM$site_id %in% ty$site_id[ty$group == "NEON"] } else NM$neon <- NA
NM$key <- NULL
utils::write.csv(NM, file.path(ldir, "npn_monthly.csv"), row.names = FALSE)
log_msg(sprintf("NPN monthly: %d site-species-months from %d readings; shared grasses with the LTER records: %s", nrow(NM), nrow(go),
                paste(intersect(unique(NM$name), unique(M$name)), collapse = ", ")))
log_msg("=== done | ", ldir, " | log: ", log_file)
