#!/usr/bin/env bash
# Regenerate missing block*.eigen.bin from existing block*.ldm.bin using GCTB --make-ldm-eigen.
#
# Usage:
#   GCTB=/path/to/gctb THREADS=4 ./scripts/make_missing_eigen_from_ldm.sh /path/to/ldm_folder
#
# Block IDs are read from ldm.info (first column). For each ID, if block<ID>.ldm.bin exists
# and block<ID>.eigen.bin is missing or empty (size 0), runs:
#   gctb --ldm <dir> --make-ldm-eigen --block <id> --out <dir> --thread <n>
#
# If your ldm.info uses consecutive IDs 1..964 only, you can instead run manually:
#   for i in $(seq 1 964); do
#     [[ -f "$ldm/block${i}.ldm.bin" && ! -s "$ldm/block${i}.eigen.bin" ]] || continue
#     gctb --ldm "$ldm" --make-ldm-eigen --block "$i" --out "$ldm" --thread 4
#   done

set -euo pipefail

LDM_DIR="${1:-${LDM:-}}"
GCTB="${GCTB:-gctb}"
THREADS="${THREADS:-4}"

if [[ -z "${LDM_DIR}" || ! -d "${LDM_DIR}" ]]; then
  echo "Usage: $0 <ldm_directory>" >&2
  echo "  Environment: GCTB (default: gctb), THREADS (default: 4), LDM_DIR as first arg." >&2
  exit 1
fi

LDM_INFO="${LDM_DIR}/ldm.info"
if [[ ! -f "${LDM_INFO}" ]]; then
  echo "Missing ${LDM_INFO}" >&2
  exit 1
fi

# First column = block ID (skip header line containing "Block")
mapfile -t BLOCK_IDS < <(awk 'NR>1 && $1 ~ /^[0-9]+$/ { print $1 }' "${LDM_INFO}")

if [[ ${#BLOCK_IDS[@]} -eq 0 ]]; then
  echo "No block IDs parsed from ${LDM_INFO}" >&2
  exit 1
fi

echo "Parsed ${#BLOCK_IDS[@]} block ID(s) from ldm.info."

n_run=0
n_skip_ldm=0
n_skip_have_eig=0

for id in "${BLOCK_IDS[@]}"; do
  ldmf="${LDM_DIR}/block${id}.ldm.bin"
  eigf="${LDM_DIR}/block${id}.eigen.bin"
  if [[ ! -f "${ldmf}" ]]; then
    n_skip_ldm=$((n_skip_ldm + 1))
    continue
  fi
  # -s: file exists and size > 0; ! -s catches missing or zero-byte .eigen.bin
  if [[ -s "${eigf}" ]]; then
    n_skip_have_eig=$((n_skip_have_eig + 1))
    continue
  fi
  echo "=== block ${id}: LDM present, eigen missing or empty — running make-ldm-eigen ==="
  "${GCTB}" --ldm "${LDM_DIR}" --make-ldm-eigen --block "${id}" --out "${LDM_DIR}" --thread "${THREADS}"
  n_run=$((n_run + 1))
done

echo "Done. Ran make-ldm-eigen for ${n_run} block(s); skipped ${n_skip_ldm} without LDM; skipped ${n_skip_have_eig} with non-empty eigen already."
