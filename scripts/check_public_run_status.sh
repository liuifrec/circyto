#!/usr/bin/env bash
set -euo pipefail

IMR90_ROOT="${IMR90_ROOT:-/user/ifrec/liuyuchen/circyto_redo/scrr_imr90}"
HAP1_ROOT="${HAP1_ROOT:-/user/ifrec/liuyuchen/circyto_redo/scrr_hap1}"
IMR90_WORKDIR="${IMR90_WORKDIR:-${IMR90_ROOT}/work_hg38_fullstack}"
HAP1_WORKDIR="${HAP1_WORKDIR:-${HAP1_ROOT}/work_hg38_fullstack}"
CIRCYTO_BIN="${CIRCYTO_BIN:-circyto}"

check_one() {
  local label="$1"
  local root="$2"
  local workdir="$3"

  echo "== ${label} =="
  echo "root: ${root}"
  echo "workdir: ${workdir}"

  if [[ ! -d "${root}" ]]; then
    echo "[WARN] root missing: ${root}"
    return 0
  fi
  if [[ ! -d "${workdir}" ]]; then
    echo "[WARN] workdir missing: ${workdir}"
    return 0
  fi

  if [[ -f "${workdir}/workflow_summary.json" ]]; then
    echo "[INFO] workflow_summary.json present"
  else
    echo "[WARN] workflow_summary.json missing"
  fi

  if [[ -f "${workdir}/mudata/full_length.h5mu" ]]; then
    echo "[INFO] h5mu: ${workdir}/mudata/full_length.h5mu"
    du -sh "${workdir}/mudata/full_length.h5mu" || true
  else
    echo "[INFO] h5mu not present yet"
  fi

  du -sh "${workdir}" || true

  if [[ -f "${workdir}/qc/cleanup_alignments_dryrun.json" ]]; then
    echo "[INFO] cleanup dry-run JSON present"
  else
    echo "[INFO] cleanup dry-run JSON not present"
  fi

  if command -v "${CIRCYTO_BIN}" >/dev/null 2>&1; then
    "${CIRCYTO_BIN}" inspect-workdir --workdir "${workdir}" --json || true
    "${CIRCYTO_BIN}" validate-workdir --workdir "${workdir}" --json || true
    if [[ -f "${workdir}/mudata/full_length.h5mu" ]]; then
      "${CIRCYTO_BIN}" inspect-mudata --input "${workdir}/mudata/full_length.h5mu" --json || true
      "${CIRCYTO_BIN}" summarize-mudata-qc --input "${workdir}/mudata/full_length.h5mu" --json || true
    fi
  else
    echo "[WARN] circyto command not found on PATH"
  fi
}

check_one "IMR90" "${IMR90_ROOT}" "${IMR90_WORKDIR}"
echo
check_one "HAP1" "${HAP1_ROOT}" "${HAP1_WORKDIR}"
