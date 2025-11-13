# CIRI-full integration (chr21 Smart-seq2 example)

This page documents a concrete, tested example of using `circyto` with **CIRI-full v2.x** on a chr21-only Smart-seq2 subset (E-MTAB-6072).

## Overview

We demonstrate:

1. Preparing a manifest of paired-end FASTQs
2. Running CIRI-full via `circyto run-manifest`
3. Normalizing detector outputs into TSV per cell
4. Building a 12 × 16 circRNA × cell matrix
5. Guarding the pipeline with an integration test

## File layout

Example layout for the chr21 subset:

```text
ref/
  chr21.fa
  chr21.fa.{amb,ann,bwt,pac,sa}
  chr21.gtf

fastq/E-MTAB-6072/
  ERR2139486_1.fastq.gz
  ERR2139486_2.fastq.gz
  ERR2139559_1.fastq.gz
  ERR2139559_2.fastq.gz
  ... (other cells)

tools/CIRI-full_v2.0/
  CIRI-full.jar
  bin/ciri_full_adapter.sh

manifest.tsv
manifest_2.tsv
```

## Manifest

Full manifest (`manifest.tsv`):

```tsv
cell_id    r1                                               r2
ERR2139486 fastq/E-MTAB-6072/ERR2139486_1.fastq.gz         fastq/E-MTAB-6072/ERR2139486_2.fastq.gz
ERR2139559 fastq/E-MTAB-6072/ERR2139559_1.fastq.gz         fastq/E-MTAB-6072/ERR2139559_2.fastq.gz
...
```

2-cell test manifest (`manifest_2.tsv`):

```bash
(head -n1 manifest.tsv && tail -n +2 manifest.tsv | head -n2) > manifest_2.tsv
```

## Running CIRI-full via run-manifest

```bash
circyto run-manifest   --manifest manifest.tsv   --outdir work_smartseq2/ciri_full_chr21_all   --ref-fa ref/chr21.fa   --gtf ref/chr21.gtf   --cmd-template 'bash -lc "
    export R1={r1} R2={r2} REF_FA={ref_fa} GTF={gtf} OUT_TSV={out_tsv} THREADS=8 ;
    tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh
  "'
```

This creates normalized per-cell TSVs such as:

```text
work_smartseq2/ciri_full_chr21_all/ERR2139486.tsv
work_smartseq2/ciri_full_chr21_all/ERR2139559.tsv
...
```

Each TSV follows the common circRNA schema:

```text
circ_id    chr    start    end    strand    support
```

## Building the matrix

```bash
circyto collect   --cirifull-dir work_smartseq2/ciri_full_chr21_all   --matrix work_smartseq2/circ_chr21_all.mtx   --circ-index work_smartseq2/circ_chr21_all_ids.txt   --cell-index work_smartseq2/cell_chr21_all_ids.txt
```

For the 16-cell chr21 subset, the matrix header looks like:

```bash
grep -v '^%' work_smartseq2/circ_chr21_all.mtx | head -n 5
# 12 16 13
# 1 1 1
# 2 2 1
# 3 3 1
# 4 4 1
```

This indicates:

- 12 circRNAs
- 16 cells
- 13 non-zero entries

Cell index consistency check:

```bash
wc -l work_smartseq2/cell_chr21_all_ids.txt
# 16

wc -l manifest.tsv
# 17   (1 header + 16 data rows)
```

## Integration test

An end-to-end integration test is provided in:

```text
tests/test_cirifull_chr21_integration.py
```

It runs `run-manifest` on `manifest_2.tsv`, calls `collect`, and checks that:

- The resulting matrix is non-empty (nnz > 0)
- The number of matrix columns matches the number of cells in the manifest

Run:

```bash
pytest tests/test_cirifull_chr21_integration.py -vv
```

This test protects the CIRI-full integration against future refactoring.

## Notes / TODO
- MatrixMarket header uses `general`, so the matrix can be read directly with `scipy.io.mmread` and Scanpy.
- circRNA IDs will be normalized to a consistent format such as `chr:start|end|strand`.
- A dedicated `circyto export-h5ad` command is planned for direct AnnData export.
