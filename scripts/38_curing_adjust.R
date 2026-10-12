#!/usr/bin/env Rscript
# Step 38 (evaluation, curing rule): does "each cell's typical curing date plus a small weather
# adjustment" beat the typical date alone, using only weather known some days before the date?
# Follows 37 (cue screening), which used whole seasonal windows and so could see weather after
# curing. Here every cue ends at a cutoff = the cell's typical date minus the lead time
# (R/cues.R adj_cues: rain in the last 30 / 90 / 180 days, mean Tmax in the last 30 / 60 days,
# dry days in the last 60, days since the last 10 mm rain, and this year's green-up date if it
# came before the cutoff).
# Same held-out design and scores as 35: the typical date is the cell's median over the
# fitting years (odd or even years); the adjustment is a linear model of up to three cue
# departures (no intercept; cues picked on the fitting years, pooled over the region's cells);
# both are scored on the other years by mean absolute error, errors capped at 90 days.
# Season jumps (> 90 days from the typical date) are left out of fitting but scored, as in 35.
# Outputs <run_dir>/eval/curing_adjust/: scores.csv (per target, lead, region), coefs.csv
# (cues and days per unit per fold), cues_<target>_<lead>_<fold>_<profile>.rds (cache).
# Runs on the P720 (needs eval/gsi_model/sample_<profile>.rds with eos50; minutes).
#
# Usage: Rscript scripts/38_curing_adjust.R [profile]

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1) args[1] else if (exists("gsi_profile", envir = globalenv())) {
  base::get("gsi_profile", envir = globalenv())
} else Sys.getenv("R_CONFIG_ACTIVE", "dev")

source(here::here("R", "config.R"))
source(here::here("R", "log.R"))
source(here::here("R", "grid.R"))
source(here::here("R", "features.R"))
source(here::here("R", "cluster.R"))
source(here::here("R", "evaluate.R"))
source(here::here("R", "block2.R"))
source(here::here("R", "cues.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- load_cfg(profile)
log_file <- log_init("38_curing_adjust", profile)
ca <- cfg$curing_adjust %||% list()
out_dir <- file.path(cfg$run_dir, "eval", "curing_adjust"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
log_msg("=== 38_curing_adjust | profile: ", profile)

cache <- file.path(cfg$run_dir, "eval", "gsi_model", sprintf("sample_%s.rds", profile))
if (!file.exists(cache)) stop("no ", cache, ": run 34 (and 35 with GSI_CURING_TARGET=eos50) first")
S <- readRDS(cache)
if (is.null(S$obs$eos50)) stop("sample cache lacks eos50: run 35 with GSI_CURING_TARGET=eos50 once to add it")
years <- S$years; nY <- length(years); n <- nrow(S$cells)

grid <- terra::rast(file.path(static_dir(cfg), "grid_mask.tif"))
rs <- list(partition = ca$partition %||% "A2_groups", k = ca$k %||% 13)
D <- region_display(cfg, grid, list(rs))[[1]]
if (is.null(D)) stop("region set ", rs$partition, " k", rs$k, " not found")
reg <- D$num[S$cells$row]
folds <- list(odd = which(years %% 2 == 1), even = which(years %% 2 == 0))
leads <- unlist(ca$lead_days %||% list(0, 14, 30))
cap <- ca$cap_days %||% 90; jump <- ca$max_abs_anomaly %||% 90
tg_all <- list(eos50 = S$obs$eos50, cure50 = S$obs$cure)
targets <- tg_all[intersect(unlist(ca$targets %||% list("eos50", "cure50")), names(tg_all))]
regions <- c("CONUS", sort(unique(reg[!is.na(reg)])))
regions <- regions[vapply(regions, function(rn) (if (rn == "CONUS") n else sum(reg == as.integer(rn), na.rm = TRUE)) >= (ca$min_cells %||% 25), TRUE)]

rows_s <- list(); rows_c <- list()
for (tn in names(targets)) for (lead in leads) {
  Y <- targets[[tn]]
  # per fold: typical date, cues up to the cutoff, departures; then fit per region
  FX <- list()
  for (f in names(folds)) {
    tr <- folds[[f]]; te <- setdiff(seq_len(nY), tr)
    typ <- apply(Y[, tr, drop = FALSE], 1, stats::median, na.rm = TRUE)
    cfile <- file.path(out_dir, sprintf("cues_%s_%d_%s_%s.rds", tn, lead, f, profile))
    if (file.exists(cfile)) CU <- readRDS(cfile) else {
      CU <- lapply(stats::setNames(adj_cue_names, adj_cue_names), function(x) matrix(NA_real_, n, nY))
      for (j in seq_len(nY)) {
        v <- adj_cues(S$drv$ppt, S$drv$tmax, S$dates, years[j], typ - lead, S$obs$sos[, j])
        for (cn in adj_cue_names) CU[[cn]][, j] <- v[, cn]
      }
      saveRDS(CU, cfile)
    }
    # cue departures from the fitting-year cell means; a missing cue (e.g. green-up after the
    # cutoff) counts as typical (0) in adj_fit
    X <- lapply(CU, function(x) x - rowMeans(x[, tr, drop = FALSE], na.rm = TRUE))
    dep <- Y - typ
    dtr <- dep; dtr[, te] <- NA; dtr[abs(dtr) > jump] <- NA
    for (rn in regions) {
      i <- if (rn == "CONUS") seq_len(n) else which(reg == as.integer(rn))
      fx <- adj_fit(dtr[i, tr, drop = FALSE], lapply(X, function(x) x[i, tr, drop = FALSE]),
                    lapply(X, function(x) x[i, te, drop = FALSE]),
                    ca$max_cues %||% 3, ca$min_gain %||% 0.01, ca$max_adj_days %||% 60)
      FX[[paste(rn, f)]] <- list(i = i, te = te, typ = typ[i], fx = fx)
      if (length(fx$cues)) rows_c[[length(rows_c) + 1]] <- data.frame(target = tn, lead = lead, region = rn, fold = f,
                                                                       cue = fx$cues, b = unname(fx$b), r2_fit = fx$r2_fit)
    }
  }
  for (rn in regions) {
    e_t <- e_a <- d_obs <- a_pred <- c(); cues <- c(odd = "", even = "")
    for (f in names(folds)) {
      z <- FX[[paste(rn, f)]]; obs <- Y[z$i, z$te, drop = FALSE]
      pt <- matrix(z$typ, nrow(obs), ncol(obs)); pa <- pt + z$fx$adj
      u <- !is.na(obs) & !is.na(pt)
      e_t <- c(e_t, pmin(abs(pt - obs)[u], cap)); e_a <- c(e_a, pmin(abs(pa - obs)[u], cap))
      d_obs <- c(d_obs, (obs - pt)[u]); a_pred <- c(a_pred, z$fx$adj[u])
      cues[[f]] <- paste(z$fx$cues, collapse = "+")
    }
    nj <- abs(d_obs) <= jump
    rows_s[[length(rows_s) + 1]] <- data.frame(
      target = tn, lead = lead, region = rn, cells = length(FX[[paste(rn, "odd")]]$i), n = length(e_t),
      mae_typical = mean(e_t), mae_adjusted = mean(e_a), gain = mean(e_t) - mean(e_a),
      better_share = mean(e_a < e_t - 0.5), worse_share = mean(e_a > e_t + 0.5),
      r_adj = if (sum(nj) > 30 && stats::sd(a_pred[nj]) > 0) stats::cor(a_pred[nj], d_obs[nj]) else NA,
      cues_odd = cues[["odd"]], cues_even = cues[["even"]])
  }
  log_msg(sprintf("%s lead %2d d done (%.1f min)", tn, lead, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
Sc <- do.call(rbind, rows_s); Co <- if (length(rows_c)) do.call(rbind, rows_c) else data.frame()
utils::write.csv(Sc, file.path(out_dir, "scores.csv"), row.names = FALSE)
utils::write.csv(Co, file.path(out_dir, "coefs.csv"), row.names = FALSE)
for (i in seq_len(nrow(Sc))) with(Sc[i, ], log_msg(sprintf(
  "%-6s lead %2d  region %-5s %4d cells | MAE typical %5.1f  adjusted %5.1f  gain %+5.1f d | better %3.0f %% worse %3.0f %% | r %5.2f | %s / %s",
  target, lead, region, cells, mae_typical, mae_adjusted, gain, 100 * better_share, 100 * worse_share, r_adj, cues_odd, cues_even)))
log_msg(sprintf("=== done | %.1f min | outputs: %s | log: %s", as.numeric(difftime(Sys.time(), t0, units = "mins")), out_dir, log_file))
