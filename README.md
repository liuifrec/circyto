<p align="center">
  <img src="assets/cirCyto_logo.png" width="440" alt="circyto logo">
</p>

<p align="center">
  <a href="https://github.com/liuifrec/circyto/actions"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/liuifrec/circyto/tests.yml?branch=main"></a>
  <a href="https://github.com/liuifrec/circyto/releases"><img alt="Release" src="https://img.shields.io/github/v/release/liuifrec/circyto?display_name=tag&sort=semver"></a>
  <a href="https://github.com/liuifrec/circyto/blob/main/LICENSE"><img alt="License" src="https://img.shields.io/github/license/liuifrec/circyto"></a>
  <a href="https://pypi.org/project/circyto/"><img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-blue"></a>
</p>

# circyto

**circyto** is a unified Python CLI framework for *single-cell circRNA detection*, *detector orchestration*, *matrix generation*, and *multimodal export*.

It provides a reproducible interface over multiple circRNA detectors and produces standardized outputs suitable for downstream single-cell analysis (Scanpy, Seurat, scVI, etc.).

---

## ✨ Features

### Core
- Unified CLI for single-cell circRNA workflows
- Detector runners for:
  - **CIRI-full** — full-length circRNA detection (Smart-seq2 / full-length short reads)
  - **find-circ3** — modernized Python 3 rewrite of `find_circ` (junction/anchor-based)
- Manifest-driven batch execution (many cells / samples)
- Sparse MatrixMarket circRNA × cell export
- Host-gene annotation using GTF
- Multimodal export to `.h5ad`:
  - `.X` = mRNA (from an existing AnnData)
  - `obsm["X_circ"]` = circRNA counts
  - `uns["circ"]` = circRNA metadata (feature table, indices, host-gene map)

### Pipeline overview

```
FASTQ
  └─> circyto run-batch / run-detector / run-multidetector
         └─ per-cell circRNA calls
              └─> circyto collect-matrix
                    └─ circ_counts.mtx  (+ circ_index / cell_index)
                          └─> circyto annotate-host-genes
                                  └─> circyto export-multimodal
                                          └─ AnnData with mRNA + circRNA
```

---

## ✅ Command-line “truth table” (current interface)

The most common source of confusion is **where the detector name and outdir go**. Here’s the canonical map:

### Recommended entry point: `run-batch` (single detector)
- **Detector is a flag** (`--detector`)
- **Output directory is a flag** (`--outdir`)

```bash
circyto run-batch \
  --detector <DETECTOR> \
  --manifest manifest.tsv \
  --outdir work/<run_name> \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 4
```

### `run-detector` (single detector)
- **Detector is positional** (comes immediately after `run-detector`)
- **Output directory is a flag** (`--outdir`)

```bash
circyto run-detector <DETECTOR> \
  --manifest manifest.tsv \
  --outdir work/<run_name> \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8
```

> If you want “cells in parallel”, prefer `run-batch`, which exposes `--parallel` consistently.

### `run-multidetector` (multiple detectors)
- **Detectors are positional**
- **OUTDIR is positional** (the last required positional argument)
- Options can appear anywhere, but putting OUTDIR early avoids mistakes.

```bash
circyto run-multidetector \
  <DETECTOR_A> <DETECTOR_B> [<DETECTOR_C> ...] \
  work/<run_name> \
  --manifest manifest.tsv \
  --ref-fa ref/genome.fa \
  --gtf ref/genes.gtf \
  --threads 8 \
  --parallel 2
```

### `collect-matrix` (unified matrix collector)
- Uses a **detector flag** (`--detector`)
- Uses **explicit output file paths** (matrix + index files)

```bash
circyto collect-matrix \
  --detector <DETECTOR> \
  --indir work/<run_name> \
  --matrix work/<run_name>_matrix/circ_counts.mtx \
  --circ-index work/<run_name>_matrix/circ_index.txt \
  --cell-index work/<run_name>_matrix/cell_index.txt
```

---

## 🚀 Quick start (bundled chr21 mini example)

This repo includes a small chr21 reference + mini manifests (`manifest.tsv`, `manifest_2.tsv`) for smoke testing.

### 1) Run **find-circ3** on the bundled mini manifest (recommended: `run-batch`)

```bash
circyto run-batch \
  --detector find-circ3 \
  --manifest manifest_2.tsv \
  --outdir work/find_circ3_chr21 \
  --ref-fa ref/chr21.fa \
  --threads 4 \
  --parallel 2
```

### 2) Collect a circRNA × cell matrix

```bash
mkdir -p work/find_circ3_chr21_matrix

circyto collect-matrix \
  --detector find-circ3 \
  --indir work/find_circ3_chr21 \
  --matrix work/find_circ3_chr21_matrix/circ_counts.mtx \
  --circ-index work/find_circ3_chr21_matrix/circ_index.txt \
  --cell-index work/find_circ3_chr21_matrix/cell_index.txt
```

---

## 🧪 Golden path: 16-cell chr21 test (CIRI-full)

This is the recommended “sanity check” after installation. It runs a small Smart-seq2 subset on chr21 and produces a small circRNA × cell matrix.

### 1) Run CIRI-full (manifest → per-cell calls)

```bash
circyto run-batch \
  --detector ciri-full \
  --manifest manifest_2.tsv \
  --outdir work/ciri_full_chr21_16cells \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 4 \
  --parallel 2
```

### 2) Collect the MatrixMarket matrix

```bash
mkdir -p work/ciri_full_chr21_16cells_matrix

circyto collect-matrix \
  --detector ciri-full \
  --indir work/ciri_full_chr21_16cells \
  --matrix work/ciri_full_chr21_16cells_matrix/circ_counts.mtx \
  --circ-index work/ciri_full_chr21_16cells_matrix/circ_index.txt \
  --cell-index work/ciri_full_chr21_16cells_matrix/cell_index.txt
```

Expected outcome:
- A small sparse MatrixMarket matrix with non-zero entries
- `cell_index.txt` matches the cell IDs in `manifest_2.tsv`
- `circ_index.txt` contains circ loci in a consistent ID format

---

## 🔬 Multi-detector run + merging/comparison (chr21 example)

Run two detectors on the same manifest:

```bash
circyto run-multidetector \
  ciri-full find-circ3 \
  work/multidetector_chr21 \
  --manifest manifest_2.tsv \
  --ref-fa ref/chr21.fa \
  --gtf ref/chr21.gtf \
  --threads 4 \
  --parallel 2
```

Then:
- `circyto merge-detectors --help`
- `circyto compare-detectors --help`
- `circyto collect-multidetector --help`

(Those subcommands evolve fastest; the `--help` output is always the source of truth.)

---

## 📦 Installation

### System dependencies
- Python 3.10+
- `bowtie2`, `samtools`
- Detector binaries as needed:
  - **CIRI-full** JAR in `tools/` (see `tools/` layout)
  - **find-circ3** CLI installed (see below)

### From source (recommended)
```bash
git clone https://github.com/liuifrec/circyto
cd circyto
pip install -e .
```

### Install find-circ3 (recommended during development)
```bash
git clone https://github.com/liuifrec/find_circ3.git
cd find_circ3
pip install -e .
find-circ3 --help
```

---

## 🧫 Testing

### Unit tests
```bash
pytest -q
```

### Integration tests (detectors)
Run all integration tests:
```bash
pytest -m "integration"
```

Skip heavy tests (e.g. on laptops / CI):
```bash
export CIRCYTO_SKIP_INTEGRATION=1
```

---

## 🗺️ Roadmap

See `ROADMAP.md`.

Planned / under consideration for 2026 integration (STAR-based detectors):
- DCC
- circhunter
- CIRI3

---

## Citation

A methods manuscript is under preparation. In the meantime, please cite this repository:

> Liu, Y.-C. et al. “circyto: a unified CLI for single-cell circRNA detection and multimodal matrices.” GitHub repository: https://github.com/liuifrec/circyto
