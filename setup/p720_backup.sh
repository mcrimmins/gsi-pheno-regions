#!/usr/bin/env bash
# setup/p720_backup.sh (P720): rsync gsi-pheno-regions data to a backup location on the NAS.
#
# DRY RUN BY DEFAULT: lists what would be copied and how much; nothing is written until --run.
#   cd ~/RProjects/gsi-pheno-regions && bash setup/p720_backup.sh && echo "=== ALL DONE ==="
#   cd ~/RProjects/gsi-pheno-regions && bash setup/p720_backup.sh --run && echo "=== ALL DONE ==="
#
# Destination: GSI_BACKUP_DEST (environment or ~/.Renviron), either
#   a mounted folder        GSI_BACKUP_DEST=/mnt/truenas_backup/p720
#   or rsync over ssh       GSI_BACKUP_DEST=mike@truenas.local:/mnt/tank/backups/p720
# Everything lands in <dest>/gsi-pheno-regions/{data_root,study,repo,obs}/ with the same tree.
#
# Options
#   --run            copy (default is a dry run)
#   --tiers 1,2,3    which tiers (setup/p720_paths.sh); default 1,2,3, run in that order, so the
#                    slow-to-replace data go first. Tier 3 also sweeps the rest of data_root, so a
#                    full run misses nothing outside the transient exclusions.
#   --mirror         also delete files at the destination that are gone from the P720 (at most
#                    GSI_BACKUP_MAX_DELETE per item, default 1000; an empty source is refused).
#                    Without --mirror nothing at the destination is deleted. Either way, a file
#                    that is replaced or deleted at the destination is moved to
#                    <dest>/gsi-pheno-regions/_replaced/<date>/, never thrown away.
#   --with-chunks    include fire_study/dfm/chunks_* (64's resume intermediates; skipped by default)
#   --allow-same     copy even when a source is on the same share as the destination
#   --bwlimit N      rsync --bwlimit (KB/s), to leave room for runs reading data over NFS
#
# Never copied: ~/.Renviron (keys, ntfy topic), transient raw downloads, terra temp files,
# rsync partials. Sends an ntfy message at the end when GSI_NTFY_TOPIC is set.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=p720_paths.sh
source "$HERE/p720_paths.sh"

RUN=0; TIERS="1,2,3"; MIRROR=0; MAX_DELETE="${GSI_BACKUP_MAX_DELETE:-1000}"; CHUNKS=0; ALLOW_SAME=0; BWLIMIT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --run) RUN=1 ;; --tiers) TIERS="$2"; shift ;; --mirror) MIRROR=1 ;;
    --with-chunks) CHUNKS=1 ;; --allow-same) ALLOW_SAME=1 ;; --bwlimit) BWLIMIT="$2"; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 2 ;;
  esac; shift
done

DEST="$(gsi_root GSI_BACKUP_DEST "")"
if [ -z "$DEST" ]; then
  echo "Set GSI_BACKUP_DEST (a mounted NAS folder or user@host:/path), e.g. in ~/.Renviron:"
  echo "  GSI_BACKUP_DEST=/mnt/truenas_backup/p720"
  exit 2
fi
REMOTE=0; case "$DEST" in /*) ;; *:*) REMOTE=1 ;; esac
TOP="$DEST/gsi-pheno-regions"
STAMP="$(date +%Y%m%d_%H%M)"
LOG_DIR="${GSI_BACKUP_LOG_DIR:-$HOME/gsi_backup_logs}"; mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/backup_${STAMP}$( [ $RUN = 1 ] || echo _dryrun).log"
exec > >(tee "$LOG") 2>&1

human(){ numfmt --to=iec --suffix=B "${1:-0}" 2>/dev/null || echo "${1:-0}B"; }
fs_src(){ local p="$1"; while [ ! -e "$p" ] && [ "$p" != "/" ]; do p="$(dirname "$p")"; done
          findmnt -no SOURCE -T "$p" 2>/dev/null | head -n 1; }
nfs_host(){ case "$1" in *:*) printf '%s' "${1%%:*}" ;; esac; }
notify(){ local t; t="$(gsi_root GSI_NTFY_TOPIC "")"; [ -n "$t" ] || return 0
  local s; s="$(gsi_root GSI_NTFY_SERVER "https://ntfy.sh")"
  curl -fsS -m 20 -H "Title: P720 backup" -d "$1" "$s/$t" >/dev/null 2>&1 || true; }

echo "=== p720_backup | $(hostname) | $(date '+%Y-%m-%d %H:%M') | $( [ $RUN = 1 ] && echo RUN || echo 'DRY RUN') ==="
echo "destination: $TOP"
echo "tiers:       $TIERS$( [ $MIRROR = 1 ] && echo '   (mirror: deletions mirrored, old copies kept in _replaced/)')"

# ---- destination checks -----------------------------------------------------------------------
DEST_SRC=""
if [ $REMOTE = 0 ]; then
  [ -d "$DEST" ] || { echo "Destination folder $DEST does not exist (is the share mounted?)"; exit 2; }
  touch "$DEST/.gsi_write_test" 2>/dev/null && rm -f "$DEST/.gsi_write_test" ||
    { echo "Cannot write to $DEST"; exit 2; }
  DEST_SRC="$(fs_src "$DEST")"
  echo "dest share:  $DEST_SRC ($(df -h --output=avail "$DEST" | tail -n 1 | tr -d ' ') free)"
else
  rhost="${DEST%%:*}"; rpath="${DEST#*:}"
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$rhost" "mkdir -p '$rpath' && df -h '$rpath' | tail -n 1" ||
    { echo "Cannot reach $rhost over ssh with key login (BatchMode)."; exit 2; }
fi

# Same-share guard: copying data_root onto the share it already lives on is not a backup.
check_same() { # check_same SRC_DIR -> 0 ok, 1 skip
  local s; s="$(fs_src "$1")"
  if [ $REMOTE = 0 ] && [ -n "$DEST_SRC" ] && [ "$s" = "$DEST_SRC" ]; then
    echo "  !! $1 is on the destination share itself ($s): a copy there protects nothing."
    [ $ALLOW_SAME = 1 ] && { echo "     copying anyway (--allow-same)"; return 0; }
    echo "     skipped (use --allow-same to force)"; return 1
  fi
  local sh dh; sh="$(nfs_host "$s")"
  if [ $REMOTE = 0 ]; then dh="$(nfs_host "$DEST_SRC")"; else dh="${DEST%%:*}"; dh="${dh#*@}"; fi
  if [ -n "$sh" ] && [ -n "$dh" ] && [ "$sh" = "$dh" ]; then
    echo "  !  $1 already lives on NAS $sh; this copy guards against deleting or overwriting"
    echo "     files, not against losing the NAS. Keep a second copy elsewhere for that."
  fi
  return 0
}

# ---- rsync ------------------------------------------------------------------------------------
OPTS=(-rlt --human-readable --partial-dir=.rsync-partial --info=stats1,progress2 "${GSI_EXCLUDES[@]}")
[ $RUN = 1 ] || OPTS+=(--dry-run --stats --info=stats2)
[ $CHUNKS = 1 ] || OPTS+=(--exclude='chunks_*/')
[ -n "$BWLIMIT" ] && OPTS+=(--bwlimit="$BWLIMIT")
# -rlt: no owner / group / permission copying (NFS root squash and TrueNAS ACLs reject them)

FAILED=0; DONE=0
sync_one() { # sync_one BASE_DIR REL DEST_SUBDIR LABEL
  local base="$1" rel="$2" sub="$3" label="$4" src
  src="$base/$rel"; src="${src%/.}"
  [ -e "$src" ] || { echo "  -- $label: not on this machine ($src)"; return 0; }
  check_same "$src" || return 0
  # Changed files: the previous copy is always moved to _replaced/<date>/ (a bad file on the
  # P720 never silently overwrites a good backup). Deletions only with --mirror.
  local bdir="$TOP/_replaced/$STAMP/$sub"; [ $REMOTE = 1 ] && bdir="${TOP#*:}/_replaced/$STAMP/$sub"
  local extra=(--backup --backup-dir="$bdir")
  if [ $MIRROR = 1 ]; then
    # an empty source usually means an unmounted share: never mirror that onto the backup
    if [ -d "$src" ] && [ -z "$(find "$src" -type f -print -quit 2>/dev/null)" ]; then
      echo "  !! $label: source has no files (share not mounted?); mirroring skipped for it"
      FAILED=$((FAILED + 1)); return 0
    fi
    extra+=(--delete --max-delete="$MAX_DELETE")
  fi
  echo; echo "---- $label"
  echo "     $src  ->  $TOP/$sub/$rel"
  [ $REMOTE = 1 ] || mkdir -p "$TOP/$sub"
  [ $REMOTE = 1 ] && [ $RUN = 1 ] && ssh -o BatchMode=yes "${DEST%%:*}" "mkdir -p '${TOP#*:}/$sub'"
  # -R with the /./ anchor keeps the tree below the base folder
  rsync "${OPTS[@]}" "${extra[@]}" -R "$base/./$rel" "$TOP/$sub/" 2>&1 |
    grep -vE '^\s*$' | grep -E 'Number of (regular )?files|Number of deleted|deletions|Total file size|Total transferred|transferred file size|sent |total size|rsync error|failed|No such' | sed 's/^/     /'
  local rc=${PIPESTATUS[0]}
  if [ "$rc" = 25 ]; then echo "     !! stopped at --max-delete $MAX_DELETE: check before raising GSI_BACKUP_MAX_DELETE"; fi
  if [ "$rc" -ne 0 ] && [ "$rc" -ne 24 ]; then   # 24 = files vanished during transfer (a run is writing)
    echo "     !! rsync exit $rc"; FAILED=$((FAILED + 1))
  else DONE=$((DONE + 1)); fi
}

sub_of(){ case "$1" in DATA) echo data_root ;; STUDY) echo study ;; REPO) echo repo ;; OBS) echo obs ;; *) echo other ;; esac; }

IFS=',' read -ra TIER_LIST <<<"$TIERS"
for t in "${TIER_LIST[@]}"; do
  echo; echo "======== tier $t ========"
  while IFS=$'\t' read -r tier base rel made what; do
    [ "$tier" = "$t" ] || continue
    b="$(gsi_base_dir "$base")"; [ -n "$b" ] || continue
    sync_one "$b" "$rel" "$(sub_of "$base")" "$base:$rel  [$made]"
  done < <(gsi_manifest)
  if [ "$t" = 3 ]; then   # sweep: the rest of data_root (other run dirs, DEM, anything new)
    sync_one "$DATA_ROOT" "." data_root "rest of data_root (sweep; transient files excluded)"
  fi
done

echo
if [ $RUN = 1 ]; then
  echo "=== backup finished: $DONE items OK, $FAILED failed | log $LOG ==="
  notify "P720 backup finished: $DONE items OK, $FAILED failed"
else
  echo "=== dry run finished: nothing was copied. Re-run with --run to copy. | log $LOG ==="
fi
[ $FAILED = 0 ]
