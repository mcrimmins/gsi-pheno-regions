#!/usr/bin/env Rscript
# Step 58e (NPN grass curing study; laptop, lme4): robustness of the "climate across space" result
# (58, Q4) to nearby sites. 58 scores each model by leaving one site out; but many NPN "sites" are
# plots a few km apart (41 of 130 sites have another within 2 km), so a left-out plot's neighbour
# stays in the training data. Here the same four models (species only, ecoregion, site climate,
# ecoregion + climate; species random effect; weights = years observed, capped at 5) are scored by
# leaving out whole clusters of sites: single-linkage clusters at 0 km (= one site, as 58), 10 km
# and 50 km.
# Output <run_dir>/npn_study/record/spatial_cv.csv (set, target, pathway, model, radius_km, folds,
# r2, mae) and report/tables/npn_spatial_cv.csv; the log prints climate vs ecoregion R2 by radius.
#
# Usage: Rscript npn_study/scripts/58e_npn_spatial_cv.R [profile]   (after 58)

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else Sys.getenv("R_CONFIG_ACTIVE", "full")
npn_root <- if (dir.exists(here::here("npn_study", "R"))) here::here("npn_study") else here::here()
source(file.path(npn_root, "R", "npn_setup.R")); source(file.path(npn_root, "R", "npn.R"))
if (!requireNamespace("lme4", quietly = TRUE)) stop("needs lme4")

cfg <- load_npn_cfg(profile)
log_file <- log_init("58e_npn_spatial_cv", profile)
base <- file.path(cfg$run_dir, "npn_study"); rdir <- file.path(base, "record")
t0 <- Sys.time()
log_msg("=== 58e_npn_spatial_cv | profile: ", profile)
radii <- c(0, 10, 50)

cmed <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); a <- 2 * pi * x / 365.25; (atan2(mean(sin(a)), mean(cos(a))) * 365.25 / (2 * pi)) %% 365.25 }
cdiff <- function(a, b) { d <- (a - b) %% 365.25; ifelse(d > 182.625, d - 365.25, d) }
TY <- utils::read.csv(file.path(rdir, "series_typical.csv"), stringsAsFactors = FALSE, colClasses = c(region = "character"))
N <- readRDS(file.path(base, "npn_metrics.rds"))$normals
TY$wt <- pmin(TY$years, 5)
TN <- merge(TY, N[, c("site_id", setdiff(names(N), c("site_id", "lon", "lat")))], by = "site_id", all.x = TRUE)

# single-linkage clusters of sites at a radius (km)
S <- unique(TY[, c("site_id", "lat", "lon")]); S <- S[!duplicated(S$site_id), ]
la <- S$lat * pi / 180; lo <- S$lon * pi / 180
cz <- outer(sin(la), sin(la)) + outer(cos(la), cos(la)) * cos(outer(lo, lo, "-")); cz[cz > 1] <- 1; cz[cz < -1] <- -1
D <- 6371 * acos(cz)
clus <- function(r) { if (r <= 0) return(stats::setNames(seq_len(nrow(S)), S$site_id))
  g <- stats::cutree(stats::hclust(stats::as.dist(D), method = "single"), h = r); stats::setNames(g, S$site_id) }
CL <- lapply(radii, clus); names(CL) <- radii
for (r in radii) log_msg(sprintf("radius %3d km: %d clusters of %d sites (largest %d sites)", r, length(unique(CL[[as.character(r)]])), nrow(S), max(table(CL[[as.character(r)]]))))

preds <- list(greenup = c("t_marapr", "aridity", "p_jas_share", "lat"), curing = c("t_julaug", "aridity", "p_jas_share", "lat"),
              cured = c("t_julaug", "aridity", "p_jas_share", "lat"))
ctl <- lme4::lmerControl(calc.derivs = FALSE)
OUT <- list()
for (setname in c("relaxed", "yty")) for (tg in names(preds)) for (pw in c("all", "C3", "C4")) {
  x <- if (setname == "yty") TN[TN$yty %in% TRUE, ] else TN
  x <- x[!is.na(x[[tg]]) & stats::complete.cases(x[, preds[[tg]]]) & x$region != "none" & x$pathway %in% c("C3", "C4"), ]
  if (pw != "all") x <- x[x$pathway %in% pw, ]
  if (nrow(x) < 20 || length(unique(x$site_id)) < 10) next
  x$y <- cmed(x[[tg]]) + cdiff(x[[tg]], cmed(x[[tg]]))
  for (v in preds[[tg]]) x[[paste0("z_", v)]] <- (x[[v]] - mean(x[[v]])) / stats::sd(x[[v]])
  x$species <- factor(x$species_id); x$regionf <- factor(x$region)
  zc <- paste0("z_", preds[[tg]]); pwt <- if (pw == "all" && length(unique(x$pathway)) > 1) " + pathway" else ""
  forms <- list(species = paste0("y ~ 1", pwt), region = paste0("y ~ regionf", pwt), climate = paste0("y ~ ", paste(zc, collapse = " + "), pwt),
                both = paste0("y ~ regionf + ", paste(zc, collapse = " + "), pwt))
  for (mn in names(forms)) {
    fm <- stats::as.formula(paste(forms[[mn]], "+ (1 | species)"))
    for (r in radii) {
      fold <- CL[[as.character(r)]][as.character(x$site_id)]
      pr <- rep(NA_real_, nrow(x))
      for (k in unique(fold)) {
        te <- fold == k; tr <- x[!te, ]
        if (mn %in% c("region", "both")) te <- te & x$regionf %in% unique(tr$regionf)
        if (!any(te) || nrow(tr) < 10) next
        g <- tryCatch(suppressMessages(lme4::lmer(fm, data = droplevels(tr), weights = wt, REML = TRUE, control = ctl)), error = function(e) NULL)
        if (!is.null(g)) pr[te] <- tryCatch(stats::predict(g, newdata = x[te, ], allow.new.levels = TRUE), error = function(e) NA)
      }
      ok <- !is.na(pr)
      OUT[[length(OUT) + 1]] <- data.frame(set = setname, target = tg, pathway = pw, model = mn, radius_km = r, series = nrow(x), sites = length(unique(x$site_id)),
        folds = length(unique(fold)), n = sum(ok), r2 = 1 - sum((x$y[ok] - pr[ok])^2) / sum((x$y[ok] - mean(x$y[ok]))^2), mae = mean(abs(x$y[ok] - pr[ok])))
    }
  }
  o <- do.call(rbind, utils::tail(OUT, length(forms) * length(radii)))
  g <- function(m, r) o$r2[o$model == m & o$radius_km == r]
  log_msg(sprintf("%-7s %-8s %-3s %3d series %3d sites | R2 site / 10 km / 50 km: species %s | region %s | climate %s | both %s", setname, tg, pw, o$series[1], o$sites[1],
                  paste(sprintf("%.2f", sapply(radii, function(r) g("species", r))), collapse = "/"), paste(sprintf("%.2f", sapply(radii, function(r) g("region", r))), collapse = "/"),
                  paste(sprintf("%.2f", sapply(radii, function(r) g("climate", r))), collapse = "/"), paste(sprintf("%.2f", sapply(radii, function(r) g("both", r))), collapse = "/")))
}
O <- do.call(rbind, OUT)
utils::write.csv(O, file.path(rdir, "spatial_cv.csv"), row.names = FALSE)
tdir <- npn_path("report", "tables"); dir.create(tdir, showWarnings = FALSE); utils::write.csv(O, file.path(tdir, "npn_spatial_cv.csv"), row.names = FALSE)
log_msg(sprintf("=== done | %.1f min | %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), file.path(rdir, "spatial_cv.csv"), log_file))
