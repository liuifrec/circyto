#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <manifest.tsv> <dataset_root>" >&2
  exit 1
fi

MANIFEST="$1"
DATASET_ROOT="$2"
RAW_DIR="${DATASET_ROOT}/raw"
TMP_DIR="${DATASET_ROOT}/tmp_downloads"
LOG_DIR="${DATASET_ROOT}/logs"
DATE_TAG="$(date +%Y%m%d)"
LOG_PATH="${LOG_DIR}/download_fastqs_from_manifest_${DATE_TAG}.log"

mkdir -p "${RAW_DIR}" "${TMP_DIR}" "${LOG_DIR}"
exec > >(tee -a "${LOG_PATH}") 2>&1

if [[ ! -f "${MANIFEST}" ]]; then
  echo "[ERROR] missing manifest: ${MANIFEST}" >&2
  exit 1
fi

echo "[INFO] $(date -Is) starting manifest-driven FASTQ download prep"
echo "[INFO] manifest=${MANIFEST}"
echo "[INFO] dataset_root=${DATASET_ROOT}"
echo "[INFO] raw_dir=${RAW_DIR}"
echo "[INFO] tmp_dir=${TMP_DIR}"
echo "[INFO] log_path=${LOG_PATH}"

trim_cr() {
  local value="$1"
  value="${value%$'\r'}"
  printf '%s' "${value}"
}

check_disk_space() {
  local root_dir="$1"
  local required_gb="$2"
  local available_gb
  available_gb="$(df -BG "${root_dir}" | awk 'NR==2 {gsub(/G/, "", $4); print $4}')"
  echo "[INFO] disk_check required_free_gb=${required_gb} available_free_gb=${available_gb}"
  if [[ "${available_gb}" -lt "${required_gb}" ]]; then
    echo "[ERROR] insufficient free space under ${root_dir}: have ${available_gb}G, require at least ${required_gb}G" >&2
    exit 1
  fi
}

download_single_end_if_missing() {
  local srr="$1"
  local expected_fastq="${RAW_DIR}/${srr}.fastq.gz"
  if [[ -f "${expected_fastq}" ]]; then
    echo "[INFO] ${srr} single-end FASTQ already present, skipping"
    return 0
  fi

  check_disk_space "${DATASET_ROOT}" 20
  echo "[INFO] downloading single-end ${srr} with prefetch + fasterq-dump"
  prefetch --output-directory "${TMP_DIR}" "${srr}"
  fasterq-dump \
    --threads 8 \
    --outdir "${RAW_DIR}" \
    "${TMP_DIR}/${srr}/${srr}.sra"
  gzip -f "${RAW_DIR}/${srr}.fastq"
}

download_paired_end_if_missing() {
  local srr="$1"
  local expected_fastq_1="${RAW_DIR}/${srr}_1.fastq.gz"
  local expected_fastq_2="${RAW_DIR}/${srr}_2.fastq.gz"
  if [[ -f "${expected_fastq_1}" && -f "${expected_fastq_2}" ]]; then
    echo "[INFO] ${srr} paired-end FASTQs already present, skipping"
    return 0
  fi

  check_disk_space "${DATASET_ROOT}" 80
  echo "[INFO] downloading paired-end ${srr} with prefetch + fasterq-dump"
  prefetch --output-directory "${TMP_DIR}" "${srr}"
  fasterq-dump \
    --split-files \
    --threads 8 \
    --outdir "${RAW_DIR}" \
    "${TMP_DIR}/${srr}/${srr}.sra"
  gzip -f "${RAW_DIR}/${srr}_1.fastq" "${RAW_DIR}/${srr}_2.fastq"
}

header_line="$(head -n 1 "${MANIFEST}")"
if [[ -z "${header_line}" ]]; then
  echo "[ERROR] manifest is empty: ${MANIFEST}" >&2
  exit 1
fi

IFS=$'\t' read -r -a header_fields <<< "${header_line}"
sample_id_idx=-1
fastq_1_idx=-1
fastq_2_idx=-1
read_layout_idx=-1

for i in "${!header_fields[@]}"; do
  column_name="$(trim_cr "${header_fields[$i]}")"
  case "${column_name}" in
    sample_id)
      sample_id_idx="${i}"
      ;;
    fastq_1)
      fastq_1_idx="${i}"
      ;;
    fastq_2)
      fastq_2_idx="${i}"
      ;;
    read_layout)
      read_layout_idx="${i}"
      ;;
  esac
done

if [[ "${sample_id_idx}" -lt 0 || "${fastq_1_idx}" -lt 0 || "${fastq_2_idx}" -lt 0 || "${read_layout_idx}" -lt 0 ]]; then
  echo "[ERROR] manifest must contain header columns: sample_id, fastq_1, fastq_2, read_layout" >&2
  exit 1
fi

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r -a fields; do
  sample_id="$(trim_cr "${fields[$sample_id_idx]:-}")"
  if [[ -z "${sample_id}" ]]; then
    continue
  fi
  fastq_1="$(trim_cr "${fields[$fastq_1_idx]:-}")"
  fastq_2="$(trim_cr "${fields[$fastq_2_idx]:-}")"
  read_layout="$(trim_cr "${fields[$read_layout_idx]:-}")"

  layout_normalized="$(printf '%s' "${read_layout}" | tr '[:upper:]' '[:lower:]')"
  if [[ "${layout_normalized}" == "single" || "${layout_normalized}" == "single-end" ]]; then
    if [[ "${fastq_1}" =~ raw/(SRR[0-9]+)\.fastq\.gz$ ]]; then
      srr="${BASH_REMATCH[1]}"
    else
      echo "[ERROR] could not infer single-end SRR from fastq_1 path: ${fastq_1}" >&2
      exit 1
    fi
    download_single_end_if_missing "${srr}"
  elif [[ "${layout_normalized}" == "paired" || "${layout_normalized}" == "paired-end" ]]; then
    if [[ "${fastq_1}" =~ raw/(SRR[0-9]+)_1\.fastq\.gz$ ]]; then
      srr="${BASH_REMATCH[1]}"
    else
      echo "[ERROR] could not infer paired-end SRR from fastq_1 path: ${fastq_1}" >&2
      exit 1
    fi
    if [[ ! "${fastq_2}" =~ raw/${srr}_2\.fastq\.gz$ ]]; then
      echo "[ERROR] fastq_2 does not match inferred paired-end SRR ${srr}: ${fastq_2}" >&2
      exit 1
    fi
    download_paired_end_if_missing "${srr}"
  else
    echo "[ERROR] unsupported read_layout for sample ${sample_id}: ${read_layout}" >&2
    exit 1
  fi
done

echo "[INFO] $(date -Is) completed manifest-driven FASTQ download prep"
