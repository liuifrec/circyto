#!/usr/bin/env bash
set -euo pipefail

DATE_TAG="$(date +%Y%m%d)"
ROOT_DIR="/user/ifrec/liuyuchen/circyto_redo/scrr_imr90"
MANIFEST_DEFAULT="${ROOT_DIR}/manifest.tsv"
MANIFEST_ALL="${ROOT_DIR}/manifest_all.tsv"
if [[ -f "${MANIFEST_ALL}" ]]; then
  MANIFEST="${MANIFEST_ALL}"
else
  MANIFEST="${MANIFEST_DEFAULT}"
fi
RAW_DIR="${ROOT_DIR}/raw"
OUTDIR="${ROOT_DIR}/work_hg38_fullstack"
LOG_DIR="${ROOT_DIR}/logs"
LOG_PATH="${LOG_DIR}/run_scrr_imr90_full_stack_${DATE_TAG}.log"
MUDATA_OUT="${OUTDIR}/mudata/full_length.h5mu"

HG38_FA="/user/ifrec/liuyuchen/ref_clean/hg38/hg38.fa"
GTF="/user/ifrec/liuyuchen/ref_clean/hg38/gtf/gencode.v45.annotation.gtf"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_PATH}") 2>&1

echo "[INFO] $(date -Is) starting IMR90 full-stack run"
echo "[INFO] root_dir=${ROOT_DIR}"
echo "[INFO] raw_dir=${RAW_DIR}"
echo "[INFO] manifest=${MANIFEST}"
echo "[INFO] manifest_selection_rule=use manifest_all.tsv if present, else manifest.tsv"
echo "[INFO] outdir=${OUTDIR}"
echo "[INFO] expected_route=single-end FASTQ -> BWA-MEM -> direct SAM -> CIRI3 -> matrix -> h5ad"
echo "[INFO] current available runs:"
printf '  - %s\n' \
  "${RAW_DIR}/SRR30918117.fastq.gz" \
  "${RAW_DIR}/SRR30918126.fastq.gz"

echo "[INFO] manifest preview"
sed -n '1,10p' "${MANIFEST}"

circyto workflow full-length-circrna \
  --manifest "${MANIFEST}" \
  --outdir "${OUTDIR}" \
  --protocol ramda \
  --genome-fasta "${HG38_FA}" \
  --gtf "${GTF}" \
  --detector ciri3 \
  --export-h5ad \
  --gene-expression-method simple-overlap

circyto refresh-rna-qc \
  --workdir "${OUTDIR}"

circyto summarize-rna-circ \
  --workdir "${OUTDIR}" \
  --write-summary

circyto export-mudata \
  --workdir "${OUTDIR}" \
  --output "${MUDATA_OUT}" \
  --overwrite

circyto inspect-mudata \
  --input "${MUDATA_OUT}" \
  --json > "${OUTDIR}/qc/mudata_inspect.json"

circyto summarize-mudata-qc \
  --input "${MUDATA_OUT}" \
  --json > "${OUTDIR}/qc/mudata_qc_summary.json"

circyto cleanup-workflow \
  --workdir "${OUTDIR}" \
  --scope alignments \
  --dry-run > "${OUTDIR}/qc/cleanup_alignments_dryrun.json"

# Storage note:
# IMR90 is lighter than HAP1, but align/ SAM outputs can still become much larger than FASTQs.
# Review qc/cleanup_alignments_dryrun.json before any real cleanup.

# Real cleanup is intentionally not executed automatically.
# Uncomment only after verifying all final outputs:
# circyto cleanup-workflow \
#   --workdir "${OUTDIR}" \
#   --scope alignments

# Optional exploratory Scanpy QC, not part of the default production stack:
# circyto scanpy-qc-report \
#   --input "${MUDATA_OUT}" \
#   --output-dir "${OUTDIR}/scanpy_qc"

echo "[INFO] $(date -Is) completed IMR90 full-stack run"
