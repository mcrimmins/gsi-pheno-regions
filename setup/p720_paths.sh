# setup/p720_paths.sh: where gsi-pheno-regions data live on the P720.
# Sourced by setup/p720_inventory.sh and setup/p720_backup.sh (not run on its own).
#
# Roots come from the same places the R scripts use (R/config.R, npn_study/R/npn_setup.R):
#   environment variable  >  ~/.Renviron  >  config.yml default
#   GSI_DATA_ROOT    data_root     (CLAUDE.md: the TrueNAS NFS share /mnt/truenas_phenology)
#   GSI_SCRATCH_ROOT scratch_root  (local disk; transient raw downloads only)
#   GSI_STUDY_DIR    GrowingSeasonIndex/study (NPN / Globe-LFMC / iNat tables read by 31, 39, 50-59)
#   GSI_OBS_DIR      optional folder with npn_curing_metrics.rds (39, NPN study)
#   GSI_REPO         this repo's checkout (logs/ live here, git-ignored)

gsi_renv_get() {   # value of KEY in ~/.Renviron (last match; quotes stripped), never sourced
  [ -f "$HOME/.Renviron" ] || return 0
  grep -E "^[[:space:]]*$1[[:space:]]*=" "$HOME/.Renviron" | tail -n 1 |
    sed -E 's/^[^=]*=[[:space:]]*//; s/[[:space:]]+$//; s/^"(.*)"$/\1/; s/^'"'"'(.*)'"'"'$/\1/'
}
gsi_expand() { local p="$1"; p="${p/#\~/$HOME}"; printf '%s' "${p%/}"; }
gsi_root() {       # gsi_root VAR default
  local v="${!1:-}"; [ -n "$v" ] || v="$(gsi_renv_get "$1")"; [ -n "$v" ] || v="$2"
  gsi_expand "$v"
}

GSI_REPO="$(gsi_root GSI_REPO "$HOME/RProjects/gsi-pheno-regions")"
DATA_ROOT="$(gsi_root GSI_DATA_ROOT "$HOME/data/gsi-pheno-regions")"
SCRATCH_ROOT="$(gsi_root GSI_SCRATCH_ROOT "$DATA_ROOT")"
STUDY_DIR="$(gsi_root GSI_STUDY_DIR "$HOME/RProjects/GrowingSeasonIndex/study")"
OBS_DIR="$(gsi_root GSI_OBS_DIR "")"
RUN_NAME="${GSI_RUN_NAME:-full}"
FIRST_YEAR=2001; LAST_YEAR=2025      # config.yml `full: years`

# Known locations. Columns (tab-separated):
#   tier   1 = slow or impossible to get again (rate-limited / queued / hand downloads, outside inputs)
#          2 = derived results (rebuildable, but hours to days of P720 compute; decisions rest on them)
#          3 = one-shot source downloads (re-downloadable in minutes to an hour)
#          0 = transient (deleted by the scripts or safe to drop; never backed up)
#   base   DATA | SCRATCH | STUDY | REPO | OBS
#   path   relative to base ({run} = run name)
#   made_by, what
gsi_manifest() {
  local R="$RUN_NAME"
  cat <<EOF
1	DATA	$R/prism_daily	01_prism_daily.R	Daily PRISM 4 km ppt/tmin/tmax/vpdmax, one stack per var-year (2001-2025 + partial 2026); ~32 h of rate-limited downloads
1	DATA	$R/cutb	04_cutb_download.R + 06 + 04b	MODIS MOD13A2/MYD13A2 16-day VI on the 4 km grid (AppEEARS, days of queue) + offsets/, curves/, metrics/ from 06
1	DATA	$R/viirs_daily	05_viirs_daily.R	VNP09GA daily reflectance for regions 3, 10, 2012-2024 (AppEEARS, ~12 h queue per task)
1	DATA	_shared/raw/fire/firms	62_fire_firms.R	NASA FIRMS VIIRS S-NPP / NOAA-20 pieces (rate-limited API; pieces never re-fetched)
1	DATA	$R/npn_study	50-59, 90 (npn_study/)	NPN curing study outputs, sens/ runs, ecoregions, plants list, and LTER packages downloaded by hand (lter/)
1	STUDY	data/derived	GrowingSeasonIndex/study	NPN / Globe-LFMC / iNaturalist tables (lfmc_us_all.rds, npn_curing_metrics.rds ...) read by 31, 36, 39, 50-59
1	DATA	obs	GSI_OBS_DIR copy (39, 51)	NPN curing tables copied from the laptop's GrowingSeasonIndex study (npn_curing_metrics.rds ...)
2	DATA	$R/static	02_static.R, 03_herb_mask.R, 12_gsi.R	grid_mask, elev, lat/lon, NLCD shares, herb_share/mask, kbdi_map
2	DATA	$R/features	10, 11, 12_gsi.R	Per-year Block 1 / Block 2 / GSI-variant feature rasters
2	DATA	$R/summaries	10, 11, 12, 06	Across-year medians / IQR / variability rasters
2	DATA	$R/clusters	20_cluster_a1.R, 07_cutb_regions.R	A1/A2/A3 + Cut B partitions (a2_groups copied from the laptop), geo baseline
2	DATA	$R/eval	30-39, 12b	Evaluation tables: ab, obs, sync, profiles, gsi_model (sample cache), gsi_curing*, target_check, curing_cues, curing_adjust
2	DATA	$R/rules	40, 41_rules_predict.R	typical_dates.tif, adjust coefs / cue means, predict_<year>
2	DATA	$R/fire_study	60-69	FPA-FOD / MTBS tables, FIRMS Parquet, link, dfm (region_daily), weekly, realtime, test, forecast, c3c4
2	DATA	$R/report_figs	90, 69	Full-profile page figures (copied to the laptop by 91)
2	REPO	logs	all scripts	Run logs (git-ignored, one per run)
2	REPO	npn_study/logs	npn_study scripts	NPN study run logs (git-ignored)
2	OBS	.	GSI_OBS_DIR	npn_curing_metrics.rds folder, if set separately
3	DATA	_shared/raw/PRISM_us_dem_4km_bil.zip	02_static.R	PRISM 4 km DEM zip (elev.tif source)
3	DATA	_shared/raw/PRISM_us_dem_4km	02_static.R	Unzipped PRISM DEM
3	DATA	$R/viirs_daily_test	05 (GSI_VD_TEST=1)	VIIRS daily smoke test, Oct 21
3	DATA	_shared/raw/nlcd	03_herb_mask.R	Annual NLCD C1.2 land cover zips, 2001 / 2012 / 2024
3	DATA	_shared/raw/fire/fpa_fod	60_fire_fod.R	FPA-FOD 7th edition GeoPackage
3	DATA	_shared/raw/fire/mtbs	61_fire_mtbs.R	MTBS perimeters v12
3	DATA	landfire	59a_npn_landfire_download.R	LANDFIRE LF2023 FBFM40 CONUS (3.3 GB zip)
3	DATA	dev	dev profile runs	Arizona test box outputs (laptop-style dev runs on the P720)
3	DATA	test	test profile runs	Smoke-test outputs
0	SCRATCH	$R/raw/prism_daily	01	Daily raw tifs (deleted after conversion; a partial current year keeps its days)
0	SCRATCH	$R/raw/cutb	04	AppEEARS 1 km raw (deleted once a year is stacked)
0	SCRATCH	$R/raw/viirs_daily	05	AppEEARS VNP09GA raw (deleted after reduction)
0	DATA	$R/fire_study/dfm/chunks_$R	64_fire_dfm.R	Per-chunk intermediates (resume only; region_daily is the product)
EOF
}

gsi_base_dir() {
  case "$1" in
    DATA) printf '%s' "$DATA_ROOT" ;; SCRATCH) printf '%s' "$SCRATCH_ROOT" ;;
    STUDY) printf '%s' "$STUDY_DIR" ;; REPO) printf '%s' "$GSI_REPO" ;; OBS) printf '%s' "$OBS_DIR" ;;
  esac
}

# rsync / du exclusions for transient files (terra temp files, partial transfers, raw scratch)
GSI_EXCLUDES=(
  --exclude='.tmp/' --exclude='*.tmp' --exclude='spat_*' --exclude='.rsync-partial/'
  --exclude='raw/prism_daily/' --exclude='raw/cutb/' --exclude='raw/viirs_daily/'
  --exclude='.Rproj.user/'
)
