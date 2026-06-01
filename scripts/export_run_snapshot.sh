#!/usr/bin/env bash
set -euo pipefail

WORKDIR=""
DATASET_NAME=""
CIRCYTO_BIN="${CIRCYTO_BIN:-circyto}"
SNAPSHOT_BASE="${SNAPSHOT_BASE:-snapshots}"
DATE_TAG="$(date +%Y%m%d)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --dataset-name)
      DATASET_NAME="$2"
      shift 2
      ;;
    --circyto-bin)
      CIRCYTO_BIN="$2"
      shift 2
      ;;
    -h|--help)
      echo "usage: $0 --workdir WORKDIR [--dataset-name NAME] [--circyto-bin circyto]" >&2
      exit 0
      ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${WORKDIR}" ]]; then
  echo "[ERROR] --workdir is required" >&2
  exit 1
fi

if [[ ! -d "${WORKDIR}" ]]; then
  echo "[ERROR] workdir not found: ${WORKDIR}" >&2
  exit 1
fi

if [[ -z "${DATASET_NAME}" ]]; then
  DATASET_NAME="$(basename "${WORKDIR}")"
fi

OUTDIR="${SNAPSHOT_BASE}/${DATASET_NAME}_${DATE_TAG}"
mkdir -p "${OUTDIR}"

if ! command -v "${CIRCYTO_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] circyto command not found: ${CIRCYTO_BIN}" >&2
  exit 1
fi

"${CIRCYTO_BIN}" inspect-workdir --workdir "${WORKDIR}" --json > "${OUTDIR}/inspect_workdir.json"
"${CIRCYTO_BIN}" validate-workdir --workdir "${WORKDIR}" --json > "${OUTDIR}/validate_workdir.json"

if [[ -f "${WORKDIR}/mudata/full_length.h5mu" ]]; then
  "${CIRCYTO_BIN}" inspect-mudata --input "${WORKDIR}/mudata/full_length.h5mu" --json > "${OUTDIR}/inspect_mudata.json"
  "${CIRCYTO_BIN}" summarize-mudata-qc --input "${WORKDIR}/mudata/full_length.h5mu" --json > "${OUTDIR}/summarize_mudata_qc.json"
fi

du -sh "${WORKDIR}" > "${OUTDIR}/disk_usage.txt"

cat <<EOF
snapshot_dir=${OUTDIR}
inspect_workdir_json=${OUTDIR}/inspect_workdir.json
validate_workdir_json=${OUTDIR}/validate_workdir.json
disk_usage_txt=${OUTDIR}/disk_usage.txt
EOF
