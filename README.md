![CI](https://github.com/liuifrec/circyto/actions/workflows/test.yml/badge.svg)

<p align="center">
  <img src="assets/circCyto_logo.png" alt="circCyto logo" width="500">
</p>

# circyto

![GitHub release (latest by date)](https://img.shields.io/github/v/release/liuifrec/circyto)
![CI](https://github.com/liuifrec/circyto/actions/workflows/ci.yml/badge.svg)

A Python CLI toolkit for single-cell circRNA detection, merging, and conversion.
Supports modular detectors (CIRI-full, CIRI-long, find_circ, CIRCexplorer2) and
outputs sparse matrices for downstream analysis (Scanpy, scVelo, etc.).

A CLI to run **CIRI-full** on:
- **Plate/full-length** scRNA-seq (Smart-seq/2, RamDA-seq, Quartz-seq, SUPeR-seq, Tang) via a **manifest** of per-cell FASTQs.
- **10x Genomics** 3' and 5' chemistries from a **BAM** (experimental; best fidelity comes from manifest).

## Commands
- `circyto prepare` — 10x BAM → batched FASTQs (`--chemistry tenx-3p|tenx-5p`)
- `circyto run` — run CIRI-full on batches from `prepare`
- `circyto run-manifest` — run CIRI-full per-cell via a manifest TSV (recommended for plate data)
- `circyto collect` — build circ×cell matrix from outputs (auto plate/batch mode)
- `circyto convert` — write .loom / .h5ad
- `circyto make` — all-in-one; accepts either `--bam` (10x) or `--manifest` (plate)

## Manifest format (plate/full-length)
TSV with headers:

cell_id r1 r2
A01 /data/A01_R1.fastq.gz /data/A01_R2.fastq.gz
A02 /data/A02_R1.fastq.gz /data/A02_R2.fastq.gz

For single-end protocols, omit `r2` or leave it blank.

## Examples
Plate:
```bash
circyto make \
  --manifest cells.tsv \
  --ref-fa ref.fa --gtf genes.gtf \
  --outdir work --threads 32 \
  --cmd-template "CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv}"


10x 5':

circyto make \
  --bam possorted_genome_bam.bam \
  --chemistry tenx-5p \
  --ref-fa ref.fa --gtf genes.gtf \
  --outdir work --threads 16 \
  --cmd-template 'CIRI-full -r {ref_fa} -a {gtf} -1 {r1} -2 {r2} -o {out_tsv}'


Note: For 10x BAM, R1/R2 reconstruction is heuristic. For best results, use plate/manifest mode (true per-cell FASTQs).
