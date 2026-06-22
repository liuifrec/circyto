#!/usr/bin/env bash
set -euo pipefail

# Replace these paths with local manuscript objects and an output directory.
SMART_H5MU="/path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu"
HAP1_H5MU="/path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu"
IMR90_H5MU="/path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu"
OUTDIR="/path/to/results/manuscript"

mkdir -p "${OUTDIR}"

python scripts/manuscript/summarize_mudata_inventory.py \
  "${SMART_H5MU}" \
  "${HAP1_H5MU}" \
  "${IMR90_H5MU}" \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --out "${OUTDIR}/mudata_inventory.tsv"

python scripts/manuscript/hap1_rt_circ_analysis.py \
  --input "${HAP1_H5MU}" \
  --correlations-out "${OUTDIR}/hap1_rt_circ_correlations.tsv" \
  --ols-out "${OUTDIR}/hap1_rt_circ_ols.tsv" \
  --scatter-out "${OUTDIR}/hap1_rt_circ_cell_metrics.tsv"

python scripts/manuscript/imr90_cnv_circ_analysis.py \
  --input "${IMR90_H5MU}" \
  --correlations-out "${OUTDIR}/imr90_cnv_circ_correlations.tsv" \
  --enrichment-out "${OUTDIR}/imr90_cnv_high_host_genes.tsv" \
  --cell-metrics-out "${OUTDIR}/imr90_cnv_circ_cell_metrics.tsv" \
  --local-cnv-out "${OUTDIR}/imr90_local_cnv_at_circ_loci.tsv"

python scripts/manuscript/cross_dataset_host_overlap.py \
  "${SMART_H5MU}" \
  "${HAP1_H5MU}" \
  "${IMR90_H5MU}" \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --outdir "${OUTDIR}/host_overlap" \
  --enrichment-table HAP1="${OUTDIR}/hap1_rt_high_host_genes.tsv" \
  --enrichment-table IMR90="${OUTDIR}/imr90_cnv_high_host_genes.tsv"

python scripts/manuscript/known_novel_circ_summary.py \
  "${SMART_H5MU}" \
  "${HAP1_H5MU}" \
  "${IMR90_H5MU}" \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --out "${OUTDIR}/known_novel_circ_summary.tsv"
