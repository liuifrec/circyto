#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/user/ifrec/liuyuchen/circyto_redo/scrr_hap1"
RAW_DIR="${ROOT_DIR}/raw"
MANIFEST_ALL="${ROOT_DIR}/manifest_all.tsv"
TMP_DIR="${ROOT_DIR}/tmp_downloads"
DATE_TAG="$(date +%Y%m%d)"
LOG_DIR="${ROOT_DIR}/logs"
LOG_PATH="${LOG_DIR}/prepare_scrr_hap1_remaining_downloads_${DATE_TAG}.log"

mkdir -p "${RAW_DIR}" "${TMP_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_PATH}") 2>&1

REQUIRED_FREE_GB=80
AVAILABLE_GB="$(df -BG "${ROOT_DIR}" | awk 'NR==2 {gsub(/G/, "", $4); print $4}')"

echo "[INFO] $(date -Is) preparing remaining HAP1 downloads"
echo "[INFO] root_dir=${ROOT_DIR}"
echo "[INFO] raw_dir=${RAW_DIR}"
echo "[INFO] required_free_gb=${REQUIRED_FREE_GB}"
echo "[INFO] available_free_gb=${AVAILABLE_GB}"

if [[ "${AVAILABLE_GB}" -lt "${REQUIRED_FREE_GB}" ]]; then
  echo "[ERROR] insufficient free space for staged HAP1 download prep: have ${AVAILABLE_GB}G, require at least ${REQUIRED_FREE_GB}G" >&2
  exit 1
fi

download_if_missing() {
  local srr="$1"
  local fq1="${RAW_DIR}/${srr}_1.fastq.gz"
  local fq2="${RAW_DIR}/${srr}_2.fastq.gz"
  if [[ -f "${fq1}" && -f "${fq2}" ]]; then
    echo "[INFO] ${srr} already present, skipping download"
    return 0
  fi

  echo "[INFO] downloading ${srr} with prefetch + fasterq-dump"
  prefetch --output-directory "${TMP_DIR}" "${srr}"
  fasterq-dump \
    --split-files \
    --threads 8 \
    --outdir "${RAW_DIR}" \
    "${TMP_DIR}/${srr}/${srr}.sra"
  gzip -f "${RAW_DIR}/${srr}_1.fastq" "${RAW_DIR}/${srr}_2.fastq"
}

# Existing pilot run is intentionally not redownloaded:
echo "[INFO] preserving existing pilot run SRR30911454"

download_if_missing "SRR30911453"
download_if_missing "SRR30911559"

cat > "${MANIFEST_ALL}" <<EOF
sample_id	fastq_1	fastq_2	protocol	strandedness	read_layout
SRR30911454_HAP1_scRR	${RAW_DIR}/SRR30911454_1.fastq.gz	${RAW_DIR}/SRR30911454_2.fastq.gz	ramda	unstranded	paired
SRR30911453_HAP1_scRR	${RAW_DIR}/SRR30911453_1.fastq.gz	${RAW_DIR}/SRR30911453_2.fastq.gz	ramda	unstranded	paired
SRR30911559_HAP1_scRR	${RAW_DIR}/SRR30911559_1.fastq.gz	${RAW_DIR}/SRR30911559_2.fastq.gz	ramda	unstranded	paired
EOF

echo "[INFO] wrote ${MANIFEST_ALL}"
echo "[INFO] $(date -Is) completed HAP1 remaining-download preparation"
