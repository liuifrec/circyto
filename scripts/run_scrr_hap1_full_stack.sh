#!/usr/bin/env bash
set -euo pipefail

DATE_TAG="$(date +%Y%m%d)"
ROOT_DIR="/user/ifrec/liuyuchen/circyto_redo/scrr_hap1"
RAW_DIR="${ROOT_DIR}/raw"
MANIFEST_DEFAULT="${ROOT_DIR}/manifest.tsv"
MANIFEST_ALL="${ROOT_DIR}/manifest_all.tsv"
if [[ -f "${MANIFEST_ALL}" ]]; then
  MANIFEST="${MANIFEST_ALL}"
else
  MANIFEST="${MANIFEST_DEFAULT}"
fi
OUTDIR="${ROOT_DIR}/work_hg38_fullstack"
LOG_DIR="${ROOT_DIR}/logs"
LOG_PATH="${LOG_DIR}/run_scrr_hap1_full_stack_${DATE_TAG}.log"
MUDATA_OUT="${OUTDIR}/mudata/full_length.h5mu"

HG38_FA="/user/ifrec/liuyuchen/ref_clean/hg38/hg38.fa"
GTF="/user/ifrec/liuyuchen/ref_clean/hg38/gtf/gencode.v45.annotation.gtf"
STAR_INDEX="/user/ifrec/liuyuchen/ref_clean/hg38/star_index"

mkdir -p "${LOG_DIR}"
exec > >(tee -a "${LOG_PATH}") 2>&1

echo "[INFO] $(date -Is) starting HAP1 full-stack run"
echo "[INFO] root_dir=${ROOT_DIR}"
echo "[INFO] raw_dir=${RAW_DIR}"
echo "[INFO] manifest=${MANIFEST}"
echo "[INFO] outdir=${OUTDIR}"
echo "[INFO] expected_route=paired-end FASTQ -> STAR -> BWA rescue -> CIRI3 STAR tuple mode -> matrix -> h5ad"
echo "[INFO] current available runs:"
printf '  - %s\n' \
  "${RAW_DIR}/SRR30911454_1.fastq.gz" \
  "${RAW_DIR}/SRR30911454_2.fastq.gz" \
  "${RAW_DIR}/SRR30911453_1.fastq.gz" \
  "${RAW_DIR}/SRR30911453_2.fastq.gz" \
  "${RAW_DIR}/SRR30911559_1.fastq.gz" \
  "${RAW_DIR}/SRR30911559_2.fastq.gz"

echo "[INFO] manifest preview"
sed -n '1,10p' "${MANIFEST}"

circyto workflow full-length-circrna \
  --manifest "${MANIFEST}" \
  --outdir "${OUTDIR}" \
  --protocol ramda \
  --genome-fasta "${HG38_FA}" \
  --gtf "${GTF}" \
  --star-index "${STAR_INDEX}" \
  --allow-paired-ramda \
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
# HAP1 paired-end STAR/BWA rescue intermediates can become very large.
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

echo "[INFO] $(date -Is) completed HAP1 full-stack run"
