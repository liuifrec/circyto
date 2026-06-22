# Manuscript Reproducibility

This page documents the expected inputs, processed outputs, and command shapes for the `circyto` manuscript:

> Genome-state-associated circular RNA programs revealed by multimodal full-length single-cell sequencing with circyto

The repository does not contain large FASTQ files, reference genomes, processed GEO DNA tables, circAtlas databases, or completed manuscript-scale `.h5mu` objects. Replace placeholder paths such as `/path/to/data/...` with local or institutional paths.

## Expected Inputs

- Full-length single-cell RNA FASTQs or workflow manifests.
- Reference genome FASTA, for example `/path/to/ref/hg38.fa`.
- Gene annotation GTF/GFF, for example `/path/to/ref/gencode.v38.annotation.gtf`.
- STAR index for Smart-seq3 or paired-end full-length workflows.
- Optional gene-count tables for RNA modality export.
- Optional circAtlas or compatible known-circRNA annotation table.
- scRR GEO SOFT metadata for RNA/DNA cell mapping.
- Processed scRR CNV state and mappability-normalized tables for IMR90.
- Processed scRR RT/state table and average RT bedGraph for HAP1.
- Optional SComatic outputs for exploratory RNA-derived candidate variant signal layers.

## Expected Processed Objects

- `mtab8735_smartseq3/full_length.hostgene_fixed.h5mu`
- `scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu`
- `scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu`

These paths are local manuscript targets and are not committed data.

## Regenerate Smart-seq3 MuData

```bash
circyto workflow smartseq3-ciri3 \
  --read1 /path/to/data/E-MTAB-8735/Smartseq3.R1.fastq.gz \
  --read2 /path/to/data/E-MTAB-8735/Smartseq3.R2.fastq.gz \
  --index1 /path/to/data/E-MTAB-8735/Smartseq3.I1.fastq.gz \
  --index2 /path/to/data/E-MTAB-8735/Smartseq3.I2.fastq.gz \
  --annotation /path/to/data/E-MTAB-8735/sample_annotation.tsv \
  --cell-id-column cell_id \
  --index1-column index1 \
  --index2-column index2 \
  --ref-fa /path/to/ref/hg38.fa \
  --gtf /path/to/ref/gencode.v38.annotation.gtf \
  --star-index /path/to/ref/star_index_hg38 \
  --gene-counts /path/to/data/E-MTAB-8735/gene_counts.tsv \
  --gene-counts-format tsv \
  --outdir /path/to/work/mtab8735_smartseq3 \
  --threads 24 \
  --resume \
  --export-h5ad \
  --export-mudata
```

## Regenerate HAP1 RNA+circRNA+RT MuData

```bash
circyto workflow full-length-circrna \
  --manifest /path/to/data/HAP1/hap1_rna_manifest.tsv \
  --outdir /path/to/work/scrr_hap1/full_length \
  --protocol ramda \
  --genome-fasta /path/to/ref/hg38.fa \
  --gtf /path/to/ref/gencode.v38.annotation.gtf \
  --gene-counts /path/to/data/HAP1/gene_counts.tsv \
  --export-h5ad
```

```bash
circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278952_family.soft.gz \
  --out /path/to/work/scrr_hap1/scrr_cell_map.tsv
```

```bash
circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_hap1/scrr_cell_map.tsv \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu
```

```bash
circyto import-scrr-rt \
  --rt-table /path/to/data/GSE278952/HAP1_binarized_rt_state_hg38.txt.gz \
  --avg-rt-bedgraph /path/to/data/GSE278952/HAP1_midS_Avg_RT_hg38.bedGraph.gz \
  --outdir /path/to/work/scrr_hap1/rt
```

```bash
circyto merge-scrr-rt \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ.remapped.h5mu \
  --rt /path/to/work/scrr_hap1/rt/rt.h5ad \
  --output /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.h5mu
```

## Regenerate IMR90 RNA+circRNA+CNV MuData

```bash
circyto build-scrr-cell-map \
  --soft /path/to/data/GSE278958_family.soft.gz \
  --out /path/to/work/scrr_imr90/scrr_cell_map.tsv
```

```bash
circyto import-scrr-cnv \
  --cnv-states /path/to/data/IMR90/summary_CNV_states_bin50kb.txt.gz \
  --cnv-mappabilitynorm /path/to/data/IMR90/summary_CNV_mappabilitynorm_bin50kb.txt.gz \
  --outdir /path/to/work/scrr_imr90/cnv
```

```bash
circyto remap-scrr-mudata-obs \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.h5mu \
  --cell-map /path/to/work/scrr_imr90/scrr_cell_map.tsv \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu
```

```bash
circyto merge-scrr-cnv \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ.remapped.h5mu \
  --cnv /path/to/work/scrr_imr90/cnv/cnv.h5ad \
  --output /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.h5mu
```

## Repair Host-Gene Annotations

Run this when older objects lack the complete host-gene provenance fields or when GTF/GFF coordinate overlap should be applied.

```bash
circyto repair-host-genes \
  --input /path/to/work/input.h5mu \
  --output /path/to/work/output.hostgene_fixed.h5mu \
  --circ-mod circ \
  --gtf /path/to/ref/gencode.v38.annotation.gtf
```

Expected circRNA feature fields after repair:

- `host_gene`
- `host_gene_source`
- `host_gene_from_gtf`
- `host_gene_from_circatlas`
- `host_gene_from_circatlas_id`
- `host_gene_id`
- `host_genes_multi`
- `host_gene_ids_multi`
- `host_gene_n`

## Manuscript Summary Scripts

Inventory and host-gene recovery:

```bash
python scripts/manuscript/summarize_mudata_inventory.py \
  /path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --out /path/to/results/manuscript/mudata_inventory.tsv
```

HAP1 RT/circRNA analysis:

```bash
python scripts/manuscript/hap1_rt_circ_analysis.py \
  --input /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  --correlations-out /path/to/results/manuscript/hap1_rt_circ_correlations.tsv \
  --ols-out /path/to/results/manuscript/hap1_rt_circ_ols.tsv \
  --scatter-out /path/to/results/manuscript/hap1_rt_circ_cell_metrics.tsv
```

IMR90 CNV/circRNA analysis:

```bash
python scripts/manuscript/imr90_cnv_circ_analysis.py \
  --input /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --correlations-out /path/to/results/manuscript/imr90_cnv_circ_correlations.tsv \
  --enrichment-out /path/to/results/manuscript/imr90_cnv_high_host_genes.tsv \
  --cell-metrics-out /path/to/results/manuscript/imr90_cnv_circ_cell_metrics.tsv \
  --local-cnv-out /path/to/results/manuscript/imr90_local_cnv_at_circ_loci.tsv
```

Cross-dataset host-gene overlap:

```bash
python scripts/manuscript/cross_dataset_host_overlap.py \
  /path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --outdir /path/to/results/manuscript/host_overlap
```

Known versus novel circRNA summary:

```bash
python scripts/manuscript/known_novel_circ_summary.py \
  /path/to/work/mtab8735_smartseq3/full_length.hostgene_fixed.h5mu \
  /path/to/work/scrr_hap1/mudata/full_length_rna_circ_rt.hostgene_fixed.h5mu \
  /path/to/work/scrr_imr90/mudata/full_length_rna_circ_cnv.hostgene_fixed.h5mu \
  --dataset-name Smart-seq3_E-MTAB-8735 \
  --dataset-name HAP1_scRR_RT \
  --dataset-name IMR90_scRR_CNV \
  --out /path/to/results/manuscript/known_novel_circ_summary.tsv
```

## Producing Figures and Tables

The scripts above write plot-ready TSV files. Figure rendering can be done in notebooks or plotting scripts outside the repo, but manuscript tables should be regenerated from the TSV outputs rather than manually edited.
