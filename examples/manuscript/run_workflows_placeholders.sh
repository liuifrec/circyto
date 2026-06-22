#!/usr/bin/env bash
set -euo pipefail

# This file documents command shape only. It does not assume private data are present.

REF_FASTA="/path/to/ref/hg38.fa"
GTF="/path/to/ref/gencode.v38.annotation.gtf"
STAR_INDEX="/path/to/ref/star_index_hg38"

circyto workflow smartseq3-ciri3 \
  --read1 /path/to/data/E-MTAB-8735/Smartseq3.R1.fastq.gz \
  --read2 /path/to/data/E-MTAB-8735/Smartseq3.R2.fastq.gz \
  --index1 /path/to/data/E-MTAB-8735/Smartseq3.I1.fastq.gz \
  --index2 /path/to/data/E-MTAB-8735/Smartseq3.I2.fastq.gz \
  --annotation /path/to/data/E-MTAB-8735/sample_annotation.tsv \
  --cell-id-column cell_id \
  --index1-column index1 \
  --index2-column index2 \
  --ref-fa "${REF_FASTA}" \
  --gtf "${GTF}" \
  --star-index "${STAR_INDEX}" \
  --gene-counts /path/to/data/E-MTAB-8735/gene_counts.tsv \
  --gene-counts-format tsv \
  --outdir /path/to/work/mtab8735_smartseq3 \
  --threads 24 \
  --resume \
  --export-h5ad \
  --export-mudata

circyto workflow full-length-circrna \
  --manifest /path/to/data/HAP1/hap1_rna_manifest.tsv \
  --outdir /path/to/work/scrr_hap1/full_length \
  --protocol ramda \
  --genome-fasta "${REF_FASTA}" \
  --gtf "${GTF}" \
  --gene-counts /path/to/data/HAP1/gene_counts.tsv \
  --export-h5ad

circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278952_family.soft.gz \
  --out /path/to/work/scrr_hap1/scrr_cell_map.tsv

circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_hap1/scrr_cell_map.tsv \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu

circyto import-scrr-rt \
  --rt-table /path/to/data/GSE278952/HAP1_binarized_rt_state_hg38.txt.gz \
  --avg-rt-bedgraph /path/to/data/GSE278952/HAP1_midS_Avg_RT_hg38.bedGraph.gz \
  --outdir /path/to/work/scrr_hap1/rt

circyto merge-scrr-rt \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu \
  --rt /path/to/work/scrr_hap1/rt/rt.h5ad \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.h5mu

circyto import-scrr-cnv \
  --cnv-states /path/to/data/IMR90/summary_CNV_states_bin50kb.txt.gz \
  --cnv-mappabilitynorm /path/to/data/IMR90/summary_CNV_mappabilitynorm_bin50kb.txt.gz \
  --outdir /path/to/work/scrr_imr90/cnv

circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278958_family.soft.gz \
  --out /path/to/work/scrr_imr90/scrr_cell_map.tsv

circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_imr90/scrr_cell_map.tsv \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu

circyto merge-scrr-cnv \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu \
  --cnv /path/to/work/scrr_imr90/cnv/cnv.h5ad \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.h5mu
