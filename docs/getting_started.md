# Getting started

## Installation

Clone and install in editable mode:

```bash
git clone https://github.com/liuifrec/circyto.git
cd circyto
pip install -e .
```

Optional extras for detector-specific helpers:

```bash
pip install -e .[detectors]
```

Ensure external dependencies are installed for your chosen detector(s), e.g.:

- CIRI-full: Java, BWA, samtools
- Future detectors: minimap2, STAR, etc.

## Basic workflow

The typical `circyto` workflow is:

1. Prepare a manifest of cells and their FASTQs/BAMs
2. Run a detector via `circyto run` or `circyto run-manifest`
3. Normalize detector output into TSV per cell
4. Use `circyto collect` to build a circRNA × cell matrix
5. Load the matrix into Scanpy / Seurat / ML pipelines

## Manifest format

A minimal manifest (tab-separated) looks like:

```tsv
cell_id    r1    r2
cellA      fastq/cellA_R1.fastq.gz  fastq/cellA_R2.fastq.gz
cellB      fastq/cellB_R1.fastq.gz  fastq/cellB_R2.fastq.gz
```

For single-end data, you can leave the `r2` column empty or omit it, depending on the detector adapter.

## Running a detector (CIRI-full example)

```bash
circyto run-manifest   --manifest manifest.tsv   --outdir work/ciri_full   --ref-fa ref/hg38.fa   --gtf ref/hg38.gtf   --cmd-template 'bash -lc "
    export R1={r1} R2={r2} REF_FA={ref_fa} GTF={gtf} OUT_TSV={out_tsv} THREADS=8 ;
    tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh
  "'
```

This produces one normalized TSV per cell under `work/ciri_full/`.

## Building a matrix

```bash
circyto collect   --cirifull-dir work/ciri_full   --matrix circ.mtx   --circ-index circ_ids.txt   --cell-index cell_ids.txt
```

You can then load this into Scanpy:

```python
import scanpy as sc

adata = sc.read_mtx("circ.mtx").T
adata.obs_names = open("cell_ids.txt").read().split()
adata.var_names = open("circ_ids.txt").read().split()

## Multimodal export (mRNA + circRNA)

If you have a standard Scanpy `.h5ad` with gene expression, and a circRNA × cell
matrix built by `circyto collect`, you can combine them:

```bash
circyto export-multimodal \
  --genes-h5ad genes.h5ad \
  --circ-matrix work_smartseq2/circ_chr21_all.mtx \
  --circ-index work_smartseq2/circ_chr21_all_ids.txt \
  --cell-index work_smartseq2/cell_chr21_all_ids.txt \
  --out multimodal_circ_genes.h5ad
