
# Multi-detector circRNA workflows

`circyto` v0.7.0 adds a small **multi-detector framework** on top of the
per-detector engines (currently **CIRI-full** and experimental **CIRI2.pl**).

The goal is to make it easy to:

- run multiple circRNA detectors on the **same cells**,
- normalize their outputs to a single circRNA schema,
- **merge** results into union / long-form tables,
- **compare** detectors with simple metrics (e.g. Jaccard overlap).

This document walks through the standard workflow:

1. `run-multidetector` – run multiple detectors on a manifest
2. `merge-detectors` – build union and long-form circRNA tables
3. `compare-detectors` – compute overlap / summary metrics

---

## 1. Prerequisites

You should already have:

- A working `circyto` installation:

  ```
  git clone https://github.com/liuifrec/circyto.git
  cd circyto
  pip install -e .
  ```

- External tools installed and on `PATH` for the detectors you want to use:

  - **CIRI-full** (Java pipeline, BWA, samtools, etc.)
  - **CIRI2.pl** (Perl script; typically bundled under `tools/CIRI-full_v2.0/bin/`)

- A **manifest TSV** describing per-cell FASTQs. Minimal columns:

  ```
  cell_id    r1                       r2
  ERR2139486 fastq/.../ERR2139486_1.fastq.gz fastq/.../ERR2139486_2.fastq.gz
  ERR2139559 fastq/.../ERR2139559_1.fastq.gz fastq/.../ERR2139559_2.fastq.gz
  ```

- A reference genome and annotation:

  - `ref/genome.fa` (or `ref/chr21.fa` for a toy test)
  - `ref/genes.gtf` (or `ref/chr21.gtf`)

---

## 2. Run multiple detectors on the same cells

Use `run-multidetector` to execute multiple engines against the same manifest.
The example below uses **CIRI-full** and **CIRI2.pl** on a small Smart-seq2
chr21 subset:

```
circyto run-multidetector ciri-full ciri2 \
  --manifest manifest.tsv \
  --outdir work/multi_chr21 \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 8 \
  --parallel 1
```

Notes:

- The detector names (`ciri-full`, `ciri2`) correspond to the internal
  detector engines exposed by `circyto`.
- For some environments, CIRI-full is more stable with `--parallel 1`
  (especially inside constrained Codespaces).
- Each detector writes its outputs into a dedicated subdirectory under
  `--outdir`.

Resulting layout:

```
work/multi_chr21/
  ├── ciri-full/
  ├── ciri2/
  └── summary.json
```

---

## 3. Merge detector outputs into union / long-form tables

```
circyto merge-detectors \
  work/multi_chr21 \
  work/multi_chr21/merged
```

Produces:

```
merged/
  ├── circ_union.tsv
  ├── circ_by_detector.tsv
  └── metadata.json
```

---

## 4. Compare detectors (overlap & summary metrics)

```
circyto compare-detectors \
  work/multi_chr21/merged \
  work/multi_chr21/compare
```

Produces:

```
compare/
  ├── jaccard.tsv
  ├── detector_summary.tsv
  └── compare_metadata.json
```

---

## 5. Practical notes

- CIRI2 may still filter extremely low-support circRNAs even with `-0`.
- CIRI-full may require `--parallel 1` in constrained environments.
- The union / long-form tables are designed to integrate with Scanpy, Seurat,
  and ML-based workflows.
- Future versions will add more detectors and consensus calling.

---

## 6. End-to-end example

```
circyto run-multidetector ciri-full ciri2 --manifest manifest.tsv --outdir work/multi --ref-fa ref.fa --gtf genes.gtf
circyto merge-detectors work/multi work/multi/merged
circyto compare-detectors work/multi/merged work/multi/compare
```
