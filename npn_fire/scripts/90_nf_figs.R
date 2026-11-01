#!/usr/bin/env Rscript
# Step 90 (NPN x fire study; laptop): figures and tables for the study page (npn_fire/report/npn_fire.qmd).
# Reads the synced results in <run_dir>/npn_fire/ (91_sync_from_p720.R npn_fire) and the current
# data/nf_species_groups.csv. Each figure is skipped, with a log line, when its inputs aren't there yet.
# Public page: everything is summarised to EPA Level II ecoregions; no site names or site locations.
#
# Writes npn_fire/report/figs/<figure>_<profile>.png + manifest.csv (figure, profile, extent, years, made)
# and npn_fire/report/tables/*.csv (small, ecoregion-level, in git).
#   nf_record      the record by year: status records, plants and sites, plant-years with a usable onset (01, 02b)
#   nf_groups      usable onsets by plant group and kind of phenophase (leaves, flowers, fruits) (02b, 02a)
#   nf_effort_map  Level II ecoregions shaded by sites with a usable leaf-type onset, per plant group (02b)
#   nf_effort_years ecoregion x year: sites with any usable onset (02b)
#   nf_fire_fuel   fire by LANDFIRE fuel group per ecoregion: FPA-FOD all sizes and VIIRS pixel-days (03)
# Tables: nf_record_years.csv, nf_inventory_sites.csv, nf_fire_inventory.csv (when 03 is synced).
#
# Usage: Rscript npn_fire/scripts/90_nf_figs.R [profile] [figure ...]   (default: all)
# Needs: maps, sf (ecoregion layer), arrow (plant_years.parquet).

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
nf_root <- if (dir.exists(here::here("npn_fire", "R"))) here::here("npn_fire") else here::here()
source(file.path(nf_root, "R", "nf_setup.R")); source(file.path(nf_root, "R", "nf.R"))
if (!requireNamespace("maps", quietly = TRUE)) stop("needs the maps package: renv::install('maps')")

all_figs <- c("nf_record", "nf_groups", "nf_effort_map", "nf_effort_years", "nf_fire_fuel")
want <- if (length(args) >= 2) args[-1] else all_figs
bad <- setdiff(want, all_figs)
if (length(bad)) stop("unknown figure(s): ", paste(bad, collapse = ", "), "; use: ", paste(all_figs, collapse = ", "))
do <- function(name) name %in% want

cfg <- load_nf_cfg(profile)
log_file <- log_init("90_nf_figs", profile)
gap_ok <- cfg$npn_fire$clean$gap_days %||% 14
extent_label <- "CONUS"
fig_dir <- nf_path("report", "figs"); tab_dir <- nf_path("report", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE); dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
log_msg("=== 90_nf_figs | profile: ", profile, " | ", paste(want, collapse = ", "), " -> ", fig_dir)

# ---- style (same as the other project pages) ------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; bg <- "#fcfcfb"; landgray <- "#e4e3df"
# plant groups: categorical slots 1-7 (validated on the page background, light mode); the rest gray
grp_col <- c(deciduous_tree = "#2a78d6", deciduous_shrub = "#eb6834", conifer = "#1baf7a", evergreen_shrub = "#eda100",
             grass = "#e87ba4", forb = "#008300", evergreen_tree = "#6250d6", other = "#b8b7b2")
grp_lab <- c(deciduous_tree = "Deciduous trees", deciduous_shrub = "Deciduous shrubs", conifer = "Conifers",
             evergreen_shrub = "Evergreen shrubs", grass = "Grasses", forb = "Forbs", evergreen_tree = "Evergreen broadleaf trees",
             other = "Other (cacti, vines, unresolved)")
grp_main <- names(grp_col)[names(grp_col) != "other"]
grp_fold <- function(g) ifelse(g %in% grp_main, g, "other")
# LANDFIRE fuel groups: the same slots in fuel order; nonburnable / unmapped in grays
fuel_col <- c(GR = "#e87ba4", GS = "#eda100", SH = "#eb6834", TU = "#1baf7a", TL = "#2a78d6", SB = "#6250d6", NB = "#b8b7b2", unmapped = "#dddcd8")
fuel_lab <- c(GR = "Grass", GS = "Grass-shrub", SH = "Shrub", TU = "Timber-understory", TL = "Timber litter", SB = "Slash-blowdown",
              NB = "Nonburnable (ag, urban, barren)", unmapped = "Not mapped")

record <- function(name, years) {
  mf <- file.path(fig_dir, "manifest.csv")
  m <- if (file.exists(mf)) utils::read.csv(mf, stringsAsFactors = FALSE) else
    data.frame(figure = character(), profile = character(), extent = character(), years = character(), made = character())
  m <- m[!(m$figure == name & m$profile == profile), , drop = FALSE]
  m <- rbind(m, data.frame(figure = name, profile = profile, extent = extent_label, years = years, made = format(Sys.Date())))
  utils::write.csv(m[order(m$figure, m$profile), ], mf, row.names = FALSE)
  log_msg("wrote ", name, "_", profile, ".png")
}
plot_png <- function(name, years, w, h, code) {
  f <- file.path(fig_dir, sprintf("%s_%s.png", name, profile))
  grDevices::png(f, width = w, height = h, res = 170, bg = bg)
  ok <- tryCatch({ par(bg = bg, family = "sans", col.axis = ink2, fg = ink2, las = 1); code(); TRUE },
                 error = function(e) { log_err(name, ": ", conditionMessage(e)); FALSE }, finally = grDevices::dev.off())
  if (ok) record(name, years) else unlink(f)
  invisible(ok)
}
have <- function(...) { f <- c(...); ok <- all(!is.na(f) & file.exists(f))
  if (!ok) log_msg("skip: missing ", paste(basename(f[is.na(f) | !file.exists(f)]), collapse = ", "))
  ok }
yrs_label <- function(y) if (length(y) == 1) as.character(y) else paste0(min(y), "-", max(y))
title2 <- function(t, s, cex_t = 1.1, cex_s = 0.8) {
  mtext(t, side = 3, line = 1.9, adj = 0, cex = cex_t, font = 2, col = ink)
  mtext(s, side = 3, line = 0.6, adj = 0, cex = cex_s, col = ink2)
}
st_map <- function(sub, title) {
  s <- maps::map("state", fill = TRUE, plot = FALSE)
  plot(NA, xlim = c(-124.8, -66.9), ylim = c(24.4, 49.4), asp = 1 / cos(37 * pi / 180), axes = FALSE, xlab = "", ylab = "")
  graphics::polygon(s$x, s$y, col = landgray, border = "#ffffff", lwd = 0.7)
  title2(title, sub, 0.95, 0.72)
}

# ---- inputs ------------------------------------------------------------------------------------------------
cd <- file.path(cfg$out_dir, "clean"); fd <- file.path(cfg$out_dir, "fire")
py_f <- nf_table_path(file.path(cd, "plant_years"))
SG <- if (file.exists(nf_path("data", "nf_species_groups.csv"))) utils::read.csv(nf_path("data", "nf_species_groups.csv"), stringsAsFactors = FALSE) else NULL
P <- NULL
if (any(c("nf_record", "nf_groups", "nf_effort_map", "nf_effort_years") %in% want) && have(py_f)) {
  P <- nf_read_table(py_f)
  P$eco_l2_code <- as.character(P$eco_l2_code)
  if (!is.null(SG)) P$group <- SG$group[match(P$species_id, SG$species_id)]
  P$g <- grp_fold(P$group)
  P$ok <- !is.na(P$onset_gap) & P$onset_gap <= gap_ok
  log_msg(sprintf("plant_years: %d rows, %d with a usable onset (gap <= %d d)", nrow(P), sum(P$ok), gap_ok))
}
S <- if (file.exists(file.path(cd, "sites.csv"))) utils::read.csv(file.path(cd, "sites.csv"), stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character")) else NULL
eco_name <- if (!is.null(S)) { u <- S[!duplicated(S$eco_l2_code) & !is.na(S$eco_l2_code), ]
  stats::setNames(gsub("\\bUsa\\b", "USA", tools::toTitleCase(tolower(u$eco_l2_name))), u$eco_l2_code) } else character()
eco_order <- function(codes) { codes <- unique(codes[!is.na(codes)]); codes[order(as.numeric(sub("^(\\d+)\\.(\\d+)$", "\\1", codes)), as.numeric(sub("^\\d+\\.", "", codes)))] }
eco_lab <- function(codes, n = 34) { nm <- unname(eco_name[codes]); nm[is.na(nm)] <- ""; sprintf("%s %s", codes, substr(nm, 1, n)) }

# Level II polygons for maps (CONUS, simplified): the NPN study's cached layer if present, else built once
eco_cache <- NULL
eco_l2 <- function() {
  if (!is.null(eco_cache)) return(eco_cache)
  if (!requireNamespace("sf", quietly = TRUE)) { log_msg("skip ecoregion layer: needs sf"); return(NULL) }
  cands <- c(file.path(cfg$out_dir, "ecoregions", "l2_conus_plot.rds"), file.path(cfg$run_dir, "npn_study", "ecoregions", "l2_conus_plot.rds"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) { eco_cache <<- readRDS(hit[1]); return(eco_cache) }
  dirs <- c(file.path(cfg$raw_dir, "ecoregions", "l2"), file.path(cfg$run_dir, "npn_study", "ecoregions"))
  shp <- unlist(lapply(dirs[dir.exists(dirs)], list.files, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE))
  if (!length(shp)) { log_msg("skip ecoregion layer: no Level II shapefile in ", paste(dirs, collapse = " or ")); return(NULL) }
  shp <- shp[grepl("l2|level2", basename(shp), ignore.case = TRUE)][1] %||% shp[1]
  log_msg("building the ecoregion map layer from ", basename(shp), " (once; cached)")
  sf::sf_use_s2(FALSE)
  E <- sf::st_read(shp, quiet = TRUE); gcl <- attr(E, "sf_column"); names(E)[names(E) != gcl] <- toupper(names(E)[names(E) != gcl])
  E <- E[!as.character(E$NA_L2CODE) %in% c("0.0", "0"), ]
  E <- sf::st_simplify(sf::st_make_valid(sf::st_transform(E, 5070)), dTolerance = 1500, preserveTopology = TRUE)
  us <- sf::st_union(sf::st_make_valid(sf::st_transform(sf::st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE)), 5070)))
  E <- suppressWarnings(sf::st_intersection(sf::st_make_valid(E), us)); E <- E[!sf::st_is_empty(E), ]
  codes <- as.character(E$NA_L2CODE); g <- sf::st_geometry(E); sp <- split(seq_along(codes), codes)
  out <- sf::st_sf(code = names(sp), geometry = do.call(c, lapply(sp, function(i) sf::st_union(g[i]))))
  pcs <- suppressWarnings(sf::st_cast(sf::st_cast(out, "MULTIPOLYGON"), "POLYGON"))
  pcs$a <- as.numeric(sf::st_area(pcs))
  big <- do.call(rbind, lapply(split(pcs, pcs$code), function(d) d[which.max(d$a), ]))
  xy <- sf::st_coordinates(sf::st_transform(suppressWarnings(sf::st_point_on_surface(big)), 4326))
  out <- sf::st_transform(out, 4326)
  out$lab_x <- xy[match(out$code, big$code), 1]; out$lab_y <- xy[match(out$code, big$code), 2]
  out$area_km2 <- as.numeric(sf::st_area(sf::st_transform(out, 5070))) / 1e6
  dir.create(file.path(cfg$out_dir, "ecoregions"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, cands[1]); eco_cache <<- out; out
}
eco_draw <- function(E, fill = NA, border = "#ffffff", lwd = 0.5) plot(sf::st_geometry(E), col = fill, border = border, lwd = lwd, add = TRUE)

# ---- the record by year ----------------------------------------------------------------------------------
if (do("nf_record") && !is.null(P) && have(file.path(cd, "site_years.csv"))) {
  SY <- utils::read.csv(file.path(cd, "site_years.csv"))
  R <- data.frame(year = sort(unique(P$year)))
  R$records <- as.numeric(tapply(SY$n_records, SY$year, sum)[as.character(R$year)])
  R$plants <- as.integer(tapply(P$individual_id, P$year, function(z) length(unique(z)))[as.character(R$year)])
  R$sites <- as.integer(tapply(P$site_id, P$year, function(z) length(unique(z)))[as.character(R$year)])
  R$plant_years <- as.integer(table(P$year)[as.character(R$year)])
  R$usable_onsets <- as.integer(tapply(P$ok, P$year, sum)[as.character(R$year)])
  R$neon_share <- round(as.numeric(tapply(SY$neon_share * SY$n_records, SY$year, sum) / tapply(SY$n_records, SY$year, sum))[match(R$year, sort(unique(SY$year)))], 3)
  utils::write.csv(R, file.path(tab_dir, "nf_record_years.csv"), row.names = FALSE)
  plot_png("nf_record", yrs_label(R$year), 2600, 1050, function() {
    par(mfrow = c(1, 3), mar = c(3.2, 4.6, 4.2, 1))
    panel <- function(y, ttl, sub, col, fmt) {
      plot(R$year, y, type = "n", ylim = c(0, max(y, na.rm = TRUE) * 1.12), xlab = "", ylab = "", axes = FALSE)
      abline(h = pretty(c(0, max(y, na.rm = TRUE))), col = "#ecebe7", lwd = 0.8)
      rect(R$year - 0.38, 0, R$year + 0.38, y, col = col, border = bg, lwd = 2)
      axis(1, at = R$year[R$year %% 3 == 0], cex.axis = 0.8, col = NA, col.ticks = ink2); axis(2, at = pretty(c(0, max(y, na.rm = TRUE))), labels = fmt(pretty(c(0, max(y, na.rm = TRUE)))), cex.axis = 0.8, col = NA)
      i <- c(1, length(y)); text(R$year[i], y[i], fmt(y[i]), pos = 3, cex = 0.72, col = ink, xpd = NA)
      title2(ttl, sub, 0.95, 0.7)
    }
    panel(R$records / 1e6, "Status records (millions)", "plant phenophase yes / no reports", "#2a78d6", function(x) format(round(x, 1), nsmall = 1))
    panel(R$plants, "Plants observed", sprintf("at %s to %s sites a year", format(min(R$sites), big.mark = ","), format(max(R$sites), big.mark = ",")),
          "#2a78d6", function(x) format(x, big.mark = ",", trim = TRUE))
    panel(R$usable_onsets, "Usable onsets", sprintf("plant x phenophase x year, first yes within %d d of a no", gap_ok), "#2a78d6",
          function(x) format(x, big.mark = ",", trim = TRUE))
  })
}

# ---- usable onsets by plant group and kind of phenophase ----------------------------------------------------
if (do("nf_groups") && !is.null(P)) {
  kind <- ifelse(P$pheno_class_id <= 5, "Leaves", ifelse(P$pheno_class_id <= 9, "Flowers", "Fruits"))
  Q <- P[P$ok, ]; kq <- kind[P$ok]
  M <- table(factor(Q$g, levels = rev(names(grp_col))), factor(kq, levels = c("Leaves", "Flowers", "Fruits")))
  utils::write.csv(data.frame(group = rownames(M), label = grp_lab[rownames(M)], unclass(M), total = rowSums(M), check.names = FALSE),
                   file.path(tab_dir, "nf_groups.csv"), row.names = FALSE)
  kcol <- c(Leaves = "#2a78d6", Flowers = "#e87ba4", Fruits = "#eda100")
  plot_png("nf_groups", yrs_label(P$year), 2400, 1150, function() {
    par(mar = c(3.6, 13.5, 4.2, 6))
    tot <- rowSums(M); xm <- max(tot) * 1.05
    plot(NA, xlim = c(0, xm), ylim = c(0.4, nrow(M) + 0.6), axes = FALSE, xlab = "", ylab = "")
    abline(v = pretty(c(0, xm)), col = "#ecebe7", lwd = 0.8)
    for (i in seq_len(nrow(M))) { x0 <- 0
      for (k in colnames(M)) { v <- M[i, k]; if (v > 0) rect(x0, i - 0.36, x0 + v, i + 0.36, col = kcol[[k]], border = bg, lwd = 2); x0 <- x0 + v } }
    axis(1, at = pretty(c(0, xm)), labels = format(pretty(c(0, xm)), big.mark = ",", scientific = FALSE, trim = TRUE), cex.axis = 0.8, col = NA)
    axis(2, at = seq_len(nrow(M)), labels = grp_lab[rownames(M)], tick = FALSE, cex.axis = 0.85, col.axis = ink)
    text(tot, seq_len(nrow(M)), format(tot, big.mark = ","), pos = 4, cex = 0.75, col = ink, xpd = NA)
    title2("Usable onsets by plant group", sprintf("plant x phenophase x year with a first yes within %d d of a no, %s; groups from NPN functional type and USDA PLANTS growth habit", gap_ok, yrs_label(P$year)), 1.0, 0.72)
    legend("bottomright", inset = c(0, 0.02), fill = kcol, border = NA, legend = names(kcol), bty = "n", cex = 0.85, text.col = ink2, title = "phenophase", title.col = ink2)
  })
}

# ---- sites with a usable leaf-type onset by ecoregion, per group (maps) --------------------------------------
INV <- NULL
if (!is.null(P)) {
  L <- P[P$ok & P$pheno_class_id %in% 1:3 & !is.na(P$eco_l2_code), ]       # initial growth, young leaves, leaves
  key <- paste(L$eco_l2_code, L$g)
  INV <- data.frame(eco_l2_code = sub(" .*", "", unique(key)), group = sub("^\\S+ ", "", unique(key)), stringsAsFactors = FALSE)
  INV$sites <- as.integer(tapply(L$site_id, key, function(z) length(unique(z)))[unique(key)])
  INV$plants <- as.integer(tapply(L$individual_id, key, function(z) length(unique(z)))[unique(key)])
  INV$plant_years <- as.integer(table(key)[unique(key)])
  sy <- tapply(paste(L$site_id, L$year), key, function(z) length(unique(z)))
  INV$site_years <- as.integer(sy[unique(key)])
  yr10 <- tapply(seq_len(nrow(L)), key, function(i) { t <- tapply(L$site_id[i], L$year[i], function(z) length(unique(z))); sum(t >= 10) })
  INV$years_10_sites <- as.integer(yr10[unique(key)])
  INV$eco_l2_name <- unname(eco_name[INV$eco_l2_code])
  INV <- INV[order(match(INV$eco_l2_code, eco_order(INV$eco_l2_code)), match(INV$group, names(grp_col))), ]
  utils::write.csv(INV, file.path(tab_dir, "nf_inventory_sites.csv"), row.names = FALSE)
}
if (do("nf_effort_map") && !is.null(INV) && !is.null(E <- eco_l2())) {
  brk <- c(0, 1, 5, 20, 100, Inf); blab <- c("none", "1-4", "5-19", "20-99", "100+")
  seqc <- c("#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab")
  pan <- c("deciduous_tree", "deciduous_shrub", "conifer", "evergreen_shrub", "grass", "forb")
  plot_png("nf_effort_map", yrs_label(P$year), 2900, 1700, function() {
    par(mfrow = c(2, 3), mar = c(0.6, 0.4, 3.4, 0.4), oma = c(2.6, 0, 0, 0))
    for (g in pan) {
      d <- INV[INV$group == g, ]
      n <- d$sites[match(E$code, d$eco_l2_code)]; n[is.na(n)] <- 0
      k <- findInterval(n, brk, left.open = FALSE, rightmost.closed = TRUE)
      st_map(sprintf("%s sites in %d ecoregions", format(sum(d$sites), big.mark = ","), sum(d$sites > 0)), grp_lab[[g]])
      eco_draw(E, seqc[k], "#ffffff", 0.6)
      s <- maps::map("state", plot = FALSE); lines(s$x, s$y, col = grDevices::adjustcolor("#ffffff", 0.55), lwd = 0.35)
      lab <- n >= 5
      text(E$lab_x[lab], E$lab_y[lab], n[lab], cex = 0.62, font = 2, col = ifelse(k[lab] >= 4, "#ffffff", ink))
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", fill = seqc, border = "#d6d5d1", legend = blab, cex = 0.9, text.col = ink2,
           title = sprintf("sites with a usable leaf-out-type onset (initial growth, young leaves or leaves), any year %s", yrs_label(P$year)), title.col = ink2)
  })
}

# ---- ecoregion x year: sites with any usable onset ---------------------------------------------------------
if (do("nf_effort_years") && !is.null(P)) {
  Q <- P[P$ok & !is.na(P$eco_l2_code), ]
  codes <- eco_order(Q$eco_l2_code); yy <- sort(unique(P$year))
  M <- tapply(Q$site_id, list(factor(Q$eco_l2_code, levels = codes), factor(Q$year, levels = yy)), function(z) length(unique(z)))
  M[is.na(M)] <- 0
  utils::write.csv(data.frame(eco_l2_code = codes, eco_l2_name = unname(eco_name[codes]), unclass(M), check.names = FALSE),
                   file.path(tab_dir, "nf_effort_years.csv"), row.names = FALSE)
  brk <- c(0, 1, 5, 10, 25, 100, Inf); seqc <- c("#f0efec", "#cde2fb", "#9ec5f4", "#5598e7", "#256abf", "#0d366b")
  plot_png("nf_effort_years", yrs_label(yy), 2500, 1500, function() {
    par(mar = c(3.4, 17, 4.4, 1))
    plot(NA, xlim = c(0.5, length(yy) + 0.5), ylim = c(length(codes) + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_along(codes)) for (j in seq_along(yy)) {
      v <- M[i, j]; k <- findInterval(v, brk, rightmost.closed = TRUE)
      rect(j - 0.47, i - 0.45, j + 0.47, i + 0.45, col = seqc[k], border = NA)
      if (v >= 10) text(j, i, v, cex = 0.55, col = if (k >= 4) "#ffffff" else ink)
    }
    axis(1, at = seq_along(yy), labels = yy, cex.axis = 0.72, col = NA, las = 2)
    axis(2, at = seq_along(codes), labels = eco_lab(codes), tick = FALSE, cex.axis = 0.72, col.axis = ink)
    title2("Sites with a usable onset, by ecoregion and year", sprintf("any plant group or phenophase; numbers shown where 10 or more; the NPN study found ~10 sites give a reliable regional year mean"), 1.0, 0.72)
    legend("topright", inset = c(0, -0.085), horiz = TRUE, xpd = NA, bty = "n", fill = seqc, border = "#d6d5d1",
           legend = c("0", "1-4", "5-9", "10-24", "25-99", "100+"), cex = 0.72, text.col = ink2)
  })
}

# ---- fire by LANDFIRE fuel group per ecoregion (03) ----------------------------------------------------------
fi_f <- file.path(fd, "fire_inventory.csv")
if (do("nf_fire_fuel") && have(fi_f)) {
  FI <- utils::read.csv(fi_f, stringsAsFactors = FALSE, colClasses = c(eco_l2_code = "character"))
  keep_o <- c("fod_all", "fod_cplus", "mtbs_km2", "viirs_px")
  FT <- FI[FI$outcome %in% keep_o, c("outcome", "eco_l2_code", "eco_l2_name", "fuel", "total", "per_year", "years_with")]
  nm_f <- FI[!duplicated(FI$eco_l2_code) & !is.na(FI$eco_l2_name), ]                 # names for ecoregions without NPN sites
  miss_nm <- setdiff(nm_f$eco_l2_code, names(eco_name))
  if (length(miss_nm)) eco_name <- c(eco_name, stats::setNames(gsub("\\bUsa\\b", "USA", tools::toTitleCase(tolower(nm_f$eco_l2_name[match(miss_nm, nm_f$eco_l2_code)]))), miss_nm))
  utils::write.csv(FT, file.path(tab_dir, "nf_fire_inventory.csv"), row.names = FALSE)
  pan <- list(fod_all = "FPA-FOD wildfires, all sizes", viirs_px = "VIIRS active-fire pixel-days")
  codes <- eco_order(FT$eco_l2_code[FT$outcome == "fod_all"])
  fuels <- names(fuel_col)
  plot_png("nf_fire_fuel", "FPA-FOD 2001-2024; VIIRS 2012-2024", 2900, 1600, function() {
    layout(matrix(1:2, 1), widths = c(1.45, 1)); par(oma = c(3.2, 0, 0, 0))
    for (o in names(pan)) {
      par(mar = c(3.6, if (o == names(pan)[1]) 15.5 else 1, 4.4, 1.6))
      d <- FT[FT$outcome == o, ]
      M <- matrix(0, length(codes), length(fuels), dimnames = list(codes, fuels))
      i <- match(d$eco_l2_code, codes); j <- match(d$fuel, fuels); ok <- !is.na(i) & !is.na(j)
      M[cbind(i[ok], j[ok])] <- d$per_year[ok]
      tot <- rowSums(M); M <- M / ifelse(tot > 0, tot, NA); M[is.na(M)] <- 0; xm <- 1.2   # shares; totals as labels
      plot(NA, xlim = c(0, xm), ylim = c(length(codes) + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "")
      abline(v = seq(0, 1, 0.25), col = "#ecebe7", lwd = 0.8)
      for (r in seq_along(codes)) { x0 <- 0
        for (f in fuels) { v <- M[r, f]; if (v > 0) rect(x0, r - 0.38, x0 + v, r + 0.38, col = fuel_col[[f]], border = bg, lwd = 1); x0 <- x0 + v } }
      text(1, seq_along(codes), sprintf("%s / yr", format(round(tot), big.mark = ",")), pos = 4, cex = 0.62, col = ink, xpd = NA)
      axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 100, 25), " %"), cex.axis = 0.75, col = NA)
      if (o == names(pan)[1]) axis(2, at = seq_along(codes), labels = eco_lab(codes), tick = FALSE, cex.axis = 0.72, col.axis = ink)
      title2(pan[[o]], "share by LANDFIRE fuel group of each fire's 4 km cell; number = per year", 0.95, 0.72)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", horiz = TRUE, bty = "n", fill = fuel_col, border = NA, legend = fuel_lab, cex = 0.72, text.col = ink2)
  })
}
log_msg("=== done | ", fig_dir, " | log: ", log_file)
