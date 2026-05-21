#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="/user/ifrec/liuyuchen/circyto_redo/scrr_imr90"
RAW_DIR="${ROOT_DIR}/raw"
MANIFEST_ALL="${ROOT_DIR}/manifest_all.tsv"
TMP_DIR="${ROOT_DIR}/tmp_downloads"
DATE_TAG="$(date +%Y%m%d)"
LOG_DIR="${ROOT_DIR}/logs"
LOG_PATH="${LOG_DIR}/prepare_scrr_imr90_remaining_downloads_${DATE_TAG}.log"

mkdir -p "${RAW_DIR}" "${TMP_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_PATH}") 2>&1

REQUIRED_FREE_GB=20
AVAILABLE_GB="$(df -BG "${ROOT_DIR}" | awk 'NR==2 {gsub(/G/, "", $4); print $4}')"

echo "[INFO] $(date -Is) preparing remaining IMR90 downloads"
echo "[INFO] root_dir=${ROOT_DIR}"
echo "[INFO] raw_dir=${RAW_DIR}"
echo "[INFO] required_free_gb=${REQUIRED_FREE_GB}"
echo "[INFO] available_free_gb=${AVAILABLE_GB}"
echo "[INFO] validated_imr90_srrs=SRR30918117 SRR30918126"

if [[ "${AVAILABLE_GB}" -lt "${REQUIRED_FREE_GB}" ]]; then
  echo "[ERROR] insufficient free space for staged IMR90 download prep: have ${AVAILABLE_GB}G, require at least ${REQUIRED_FREE_GB}G" >&2
  exit 1
fi

download_if_missing() {
  local srr="$1"
  local fq="${RAW_DIR}/${srr}.fastq.gz"
  if [[ -f "${fq}" ]]; then
    echo "[INFO] ${srr} already present, skipping download"
    return 0
  fi

  echo "[INFO] downloading ${srr} with prefetch + fasterq-dump"
  prefetch --output-directory "${TMP_DIR}" "${srr}"
  fasterq-dump \
    --threads 8 \
    --outdir "${RAW_DIR}" \
    "${TMP_DIR}/${srr}/${srr}.sra"
  gzip -f "${RAW_DIR}/${srr}.fastq"
}

# Current validated public IMR90 RNA-side set in the repo metadata is exactly two runs.
download_if_missing "SRR30918117"
download_if_missing "SRR30918126"

cat > "${MANIFEST_ALL}" <<EOF
sample_id	fastq_1	fastq_2	protocol	strandedness	read_layout
SRR30918117_IMR90_scRR	${RAW_DIR}/SRR30918117.fastq.gz		ramda	unstranded	single
SRR30918126_IMR90_scRR	${RAW_DIR}/SRR30918126.fastq.gz		ramda	unstranded	single
EOF

echo "[INFO] wrote ${MANIFEST_ALL}"
echo "[INFO] current metadata contains no additional validated public IMR90 RNA-side SRRs beyond the two pilot runs"
echo "[INFO] $(date -Is) completed IMR90 remaining-download preparation"
