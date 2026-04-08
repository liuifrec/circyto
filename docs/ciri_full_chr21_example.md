# CIRI-full integration (chr21 Smart-seq2 example)

This page documents a concrete paired-end chr21 example of using `circyto` with the bundled **CIRI-full v2.x** assets on a Smart-seq2 subset (E-MTAB-6072).

This is intentionally a **paired-end** example, because upstream CIRI-full Pipeline mode is paired-end only.

## Overview

We demonstrate:

1. Preparing a manifest of paired-end FASTQs
2. Running `ciri-full` via `circyto run-detector`
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

## Running CIRI-full via `run-detector`

```bash
circyto run-detector ciri-full \
  --manifest manifest.tsv \
  --outdir work_smartseq2/ciri_full_chr21_all \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 8 \
  --parallel 1
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
circyto collect-matrix \
  --detector ciri-full \
  --indir work_smartseq2/ciri_full_chr21_all \
  --matrix work_smartseq2/circ_chr21_all.mtx \
  --circ-index work_smartseq2/circ_chr21_all_ids.txt \
  --cell-index work_smartseq2/cell_chr21_all_ids.txt
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

It runs `run-detector` on `manifest_2.tsv`, calls the matrix collector, and checks that:

- The resulting matrix is non-empty (nnz > 0)
- The number of matrix columns matches the number of cells in the manifest

Run:

```bash
pytest tests/test_cirifull_chr21_integration.py -vv
```

This test protects the paired-end `ciri-full` integration against future refactoring.

## Notes
- Public `ciri-full` semantics are layout-dependent:
  - paired-end rows use the upstream bundled CIRI-full Java Pipeline
  - single-end rows use a bundled CIRI2-based fallback path
- Both layouts normalize to the same TSV schema, but only the paired-end path is true upstream CIRI-full Pipeline behavior.
- MatrixMarket header uses `general`, so the matrix can be read directly with `scipy.io.mmread` and Scanpy.
- circRNA IDs will be normalized to a consistent format such as `chr:start|end|strand`.
