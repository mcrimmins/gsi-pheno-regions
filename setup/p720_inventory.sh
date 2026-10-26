#!/usr/bin/env bash
# setup/p720_inventory.sh (P720): what gsi-pheno-regions data are on this machine, where, how big,
# and on which disk (local vs the NAS). Read-only: it never changes or deletes anything.
#
#   cd ~/RProjects/gsi-pheno-regions && bash setup/p720_inventory.sh && echo "=== ALL DONE ==="
#
# Writes ~/gsi_inventory/inventory_<date>.txt (report) and inventory_<date>.tsv (one row per
# location: tier, path, size, files, newest file, filesystem). du over NFS can take several
# minutes for the large trees (PRISM, cutb, FIRMS); run it in tmux. QUICK=1 skips file counts.
# Roots: GSI_DATA_ROOT / GSI_SCRATCH_ROOT / GSI_STUDY_DIR / GSI_OBS_DIR as in ~/.Renviron
# (see setup/p720_paths.sh).

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=p720_paths.sh
source "$HERE/p720_paths.sh"

OUT_DIR="${GSI_INVENTORY_DIR:-$HOME/gsi_inventory}"; mkdir -p "$OUT_DIR"
STAMP="$(date +%Y%m%d_%H%M)"
TXT="$OUT_DIR/inventory_$STAMP.txt"; TSV="$OUT_DIR/inventory_$STAMP.tsv"
exec > >(tee "$TXT") 2>&1

hr()   { printf '\n==== %s ====\n' "$1"; }
human(){ numfmt --to=iec --suffix=B "${1:-0}" 2>/dev/null || echo "${1:-0}B"; }
fs_of(){ # filesystem source, type and target holding path $1 (nearest existing parent)
  local p="$1"; while [ ! -e "$p" ] && [ "$p" != "/" ]; do p="$(dirname "$p")"; done
  findmnt -no SOURCE,FSTYPE,TARGET -T "$p" 2>/dev/null | head -n 1 | awk '{print $1"\t"$2"\t"$3}'
}
size_b(){ du -sb "${GSI_EXCLUDES_DU[@]}" "$1" 2>/dev/null | awk '{print $1}'; }
nfiles(){ [ -n "${QUICK:-}" ] && { echo "-"; return; }; find "$1" -type f 2>/dev/null | wc -l; }
newest(){ find "$1" -type f -printf '%TY-%Tm-%Td\n' 2>/dev/null | sort | tail -n 1; }
GSI_EXCLUDES_DU=(--exclude='.tmp' --exclude='*.tmp' --exclude='spat_*')

hr "Machine and roots"
echo "host:        $(hostname)    $(date '+%Y-%m-%d %H:%M %Z')"
echo "repo:        $GSI_REPO  $( [ -d "$GSI_REPO/.git" ] && git -C "$GSI_REPO" log -1 --format='(%h %cs)' 2>/dev/null)"
for v in DATA_ROOT SCRATCH_ROOT STUDY_DIR OBS_DIR; do
  p="${!v}"; [ -n "$p" ] || { printf '%-12s (not set)\n' "$v:"; continue; }
  IFS=$'\t' read -r src fst tgt <<<"$(fs_of "$p")"
  printf '%-12s %s\n             on %s (%s, mounted at %s)%s\n' "$v:" "$p" "$src" "$fst" "$tgt" \
    "$( [ -e "$p" ] || echo '  [MISSING]')"
done
if [ "$DATA_ROOT" = "$SCRATCH_ROOT" ]; then
  echo "NOTE: GSI_SCRATCH_ROOT is not set, so raw downloads land under data_root."
fi

hr "Free space"
for p in "$DATA_ROOT" "$SCRATCH_ROOT" "$HOME" /tmp; do [ -e "$p" ] && df -h "$p" | tail -n 1; done | sort -u
echo; echo "Mounted network shares (candidates for the backup destination):"
findmnt -rn -t nfs,nfs4,cifs,smb3 -o TARGET,SOURCE,FSTYPE,SIZE,AVAIL 2>/dev/null | sed 's/^/  /' || echo "  (none)"

hr "Known locations (from the scripts)"
printf 'tier\tbase\tpath\tsize_bytes\tsize\tfiles\tnewest\tfs_source\tfs_type\tmade_by\twhat\n' > "$TSV"
printf '%-4s %-48s %9s %9s %-10s %s\n' tier path size files newest filesystem
declare -A SEEN
total_by_tier=( 0 0 0 0 )
while IFS=$'\t' read -r tier base rel made what; do
  [ -n "$tier" ] || continue
  b="$(gsi_base_dir "$base")"; [ -n "$b" ] || continue
  p="$b/$rel"; p="${p%/.}"; SEEN["$p"]=1
  if [ ! -e "$p" ]; then
    printf '%-4s %-48s %9s\n' "$tier" "$base:$rel" "(none)"
    printf '%s\t%s\t%s\t\t\t\t\t\t\t%s\t%s\n' "$tier" "$base" "$p" "$made" "$what" >> "$TSV"; continue
  fi
  s="$(size_b "$p")"; n="$(nfiles "$p")"; d="$(newest "$p")"
  IFS=$'\t' read -r src fst _ <<<"$(fs_of "$p")"
  total_by_tier[$tier]=$(( ${total_by_tier[$tier]} + ${s:-0} ))
  printf '%-4s %-48s %9s %9s %-10s %s\n' "$tier" "$base:$rel" "$(human "$s")" "$n" "$d" "$src"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tier" "$base" "$p" "$s" "$(human "$s")" \
    "$n" "$d" "$src" "$fst" "$made" "$what" >> "$TSV"
done < <(gsi_manifest)
echo
for t in 1 2 3 0; do
  label=( "transient (not backed up)" "tier 1, slow or impossible to get again" "tier 2, derived results" "tier 3, re-downloadable" )
  printf '  %-45s %s\n' "${label[$t]}:" "$(human "${total_by_tier[$t]}")"
done
echo "  (tier 1 includes cutb/curves, metrics, offsets and npn_study results that sit in the same trees)"

hr "Anything else under data_root (not made by a known step)"
found=0
for d in "$DATA_ROOT"/* "$DATA_ROOT/$RUN_NAME"/* "$DATA_ROOT/_shared/raw"/* "$DATA_ROOT/_shared/raw/fire"/*; do
  [ -e "$d" ] || continue
  case "$(basename "$d")" in spat_*|.tmp|*.tmp) continue ;; esac
  known=0
  for k in "${!SEEN[@]}"; do case "$k" in "$d"|"$d"/*) known=1; break ;; esac; done
  # parents of known paths are containers, not unknowns
  [ "$d" = "$DATA_ROOT/$RUN_NAME" ] || [ "$d" = "$DATA_ROOT/_shared" ] || [ "$d" = "$DATA_ROOT/_shared/raw" ] ||
    [ "$d" = "$DATA_ROOT/_shared/raw/fire" ] && known=1
  if [ "$known" = 0 ]; then
    found=1; printf '  %-60s %9s  %s\n' "${d#"$DATA_ROOT"/}" "$(human "$(size_b "$d")")" "$(newest "$d")"
    printf '?\tDATA\t%s\t%s\t%s\t\t%s\t\t\tunknown\tnot in the manifest\n' "$d" "$(size_b "$d")" \
      "$(human "$(size_b "$d")")" "$(newest "$d")" >> "$TSV"
  fi
done
[ "$found" = 1 ] || echo "  (nothing)"

hr "Completeness checks"
R="$DATA_ROOT/$RUN_NAME"
yrs() { # yrs DIR GLOBPREFIX: list the 4-digit years found in file names, compressed into ranges
  ls "$1" 2>/dev/null | grep -oE '(19|20)[0-9]{2}' | sort -un | awk '
    NR==1{s=$1;p=$1;next} $1==p+1{p=$1;next} {printf "%s%s ", s, (s==p?"":"-"p); s=$1;p=$1}
    END{if(NR) printf "%s%s", s, (s==p?"":"-"p)}'
}
echo "PRISM daily (expect ${FIRST_YEAR}-${LAST_YEAR} for each variable, plus a partial current year):"
for v in ppt tmin tmax vpdmax; do printf '  %-7s %s\n' "$v" "$(yrs "$R/prism_daily/$v")"; done
echo "Cut B (AppEEARS stacks, per product / variable):"
for pd in "$R"/cutb/*/; do
  b="$(basename "$pd")"; case "$b" in curves|metrics|offsets) continue ;; esac
  for vd in "$pd"*/; do [ -d "$vd" ] && printf '  %-28s %s\n' "$b/$(basename "$vd")" "$(yrs "$vd")"; done
done
for s in curves metrics offsets; do printf '  %-28s %s\n' "$s" "$(yrs "$R/cutb/$s")"; done
echo "Feature blocks (per-year files):"
for fd in "$R"/features/*/; do [ -d "$fd" ] && printf '  %-28s %s\n' "$(basename "$fd")" "$(yrs "$fd")"; done
echo "Clusters:      $(ls "$R/clusters" 2>/dev/null | tr '\n' ' ')"
echo "Eval:          $(ls "$R/eval" 2>/dev/null | tr '\n' ' ')"
echo "Fire study:    $(ls "$R/fire_study" 2>/dev/null | tr '\n' ' ')"
echo "Rules:         $(ls "$R/rules" 2>/dev/null | tr '\n' ' ')"
echo "FIRMS pieces (source / years):"
for sd in "$DATA_ROOT"/_shared/raw/fire/firms/*/; do
  [ -d "$sd" ] && printf '  %-20s %s  (%s pieces)\n' "$(basename "$sd")" "$(yrs "$sd")" \
    "$(find "$sd" -type f -name '*.csv*' 2>/dev/null | wc -l)"
done
echo "VIIRS daily:   $(find "$R/viirs_daily" -name 'daily_*.rds' 2>/dev/null | wc -l) region-year-chunk files"
[ -f "$R/viirs_daily/tasks.csv" ] && echo "               tasks.csv status: $(awk -F, 'NR==1{for(i=1;i<=NF;i++) if($i~/status/) c=i; next} c{n[$c]++} END{for(k in n) printf "%s %d  ", k, n[k]}' "$R/viirs_daily/tasks.csv")"

hr "Leftovers that are safe to clean (not deleted here)"
for d in "$SCRATCH_ROOT/$RUN_NAME/raw" "$DATA_ROOT/$RUN_NAME/raw"; do
  [ -d "$d" ] && printf '  %-60s %9s\n' "$d" "$(human "$(size_b "$d")")"
done
tmpb="$(find /tmp "${TMPDIR:-/tmp}" -maxdepth 2 -name 'spat_*' -type f -printf '%s\n' 2>/dev/null | awk '{s+=$1} END{print s+0}')"
echo "  terra temp files in /tmp:                                  $(human "$tmpb")"
for d in "$R"/fire_study/dfm/chunks_*; do [ -d "$d" ] && printf '  %-60s %9s\n' "${d#"$DATA_ROOT"/}" "$(human "$(size_b "$d")")"; done

hr "Outside data_root"
echo "repo logs:     $(ls "$GSI_REPO/logs" 2>/dev/null | wc -l) files, $(human "$(size_b "$GSI_REPO/logs")")"
if [ -d "$GSI_REPO/.git" ]; then
  ch="$(git -C "$GSI_REPO" status --porcelain 2>/dev/null | wc -l)"
  echo "repo git:      $ch uncommitted / untracked paths (code is backed up by GitHub only once pushed)"
fi
echo "~/.Renviron:   $( [ -f "$HOME/.Renviron" ] && echo present || echo missing ) (holds keys and topics: not copied by the backup script)"

echo; echo "Report: $TXT"; echo "Table:  $TSV"
