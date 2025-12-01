<p align="center">
  <img src="assets/circCyto_logo.png" alt="circyto logo" width="440">
</p>


<p align="center">
  <a href="https://github.com/liuifrec/circyto/actions">
    <img src="https://img.shields.io/github/actions/workflow/status/liuifrec/circyto/ci.yml?label=CI&style=for-the-badge" alt="CI">
  </a>
  <img src="https://img.shields.io/badge/version-v0.8.0-blue?style=for-the-badge" alt="v0.8.0">
  <a href="https://github.com/liuifrec/circyto/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/liuifrec/circyto?style=for-the-badge" alt="License">
  </a>
</p>

---

# circyto

**circyto** is a unified Python CLI framework for **single-cell circRNA detection**,  
**detector orchestration**, **matrix generation**, and **multimodal export**.

It provides a clean, reproducible interface over multiple circRNA detectors and
produces standardized outputs suitable for downstream single-cell analysis
(**Scanpy, Seurat, scVI**, etc.).

---

## ✨ Features (v0.8.x)

### Core

- Unified CLI for single-cell circRNA workflows
- Detector drivers for:
  - **CIRI-full** — short-read, full-length circRNA assembly
  - **find_circ3** — modernized Python 3 rewrite of find_circ; junction detector
- Manifest-driven batch execution for many cells / samples
- Sparse **MatrixMarket circRNA × cell** export
- **Host-gene annotation** using GTF
- **Multimodal `.h5ad` export**  
  - `.X` = mRNA  
  - `obsm["X_circ"]` = circRNA counts  
  - `uns["circ"]` = circRNA metadata (including feature table & host-gene map)

### Pipeline overview

```text
FASTQ
  └─> circyto run-detector / run-manifest
         └─ per-cell circRNA calls
              └─> circyto collect
                    └─ circ_matrix.mtx  (+ circ_index / cell_index)
                          └─> circyto annotate-host-genes
                                  └─> circyto export-multimodal
                                          └─ AnnData with mRNA + circRNA
```

### Multi-detector support

| Detector           | Type                     | Status         | Notes                                           |
|--------------------|--------------------------|----------------|------------------------------------------------|
| **CIRI-full**      | Full-length detection    | ✔ Stable       | Best for Smart-seq2 / full-length short reads  |
| **CIRI2**          | Short-read               | ✔ Supported    | Adapter included                               |
| **find_circ3**     | Junction detector        | ⚗ Experimental | Validated on chr21 mini; Python 3 rewrite      |
| **CIRCexplorer2**  | Splice junction          | In progress    | Next integration target                        |

> Legacy **find_circ** has been removed in favor of **find_circ3**.

---

## 🚀 Quick Start

### 1. Run find_circ3 on the bundled chr21 mini example

```bash
circyto run-detector find-circ3 \
  --manifest manifest_2.tsv \
  --outdir work/find_circ3_chr21 \
  --ref-fa ref/chr21.fa \
  --threads 4 \
  --parallel 2
```

### 2. Collect circRNA counts

```bash
circyto collect-find-circ3 \
  --indir work/find_circ3_chr21 \
  --outdir work/find_circ3_chr21_matrix
```

### 3. Annotate host genes

```bash
circyto annotate-host-genes \
  --circ-tsv circ_feature_table.tsv \
  --gtf Homo_sapiens.gtf \
  --out circ_feature_table_annotated.tsv
```

### 4. Export multimodal AnnData

```bash
circyto export-multimodal \
  --rna-h5ad rna.h5ad \
  --circ-mtx circ_counts.mtx \
  --circ-index circ_index.txt \
  --cell-index cell_index.txt \
  --circ-features circ_feature_table.tsv \
  --out combined_multimodal.h5ad
```

---

## 🧪 Testing

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

## 🔬 Example: 16-cell chr21 Smart-seq2 test (CIRI-full)

This is a minimal integration test for a 16-cell Smart-seq2 subset mapped to
**chr21** using **CIRI-full** via the manifest interface.

1. **Prepare reference and manifest**

   - Place `chr21.fa` under `ref/` (e.g. `ref/chr21.fa`)
   - Use the provided `manifest_2.tsv` (16 Smart-seq2 cells)

2. **Run CIRI-full over the 16-cell manifest**

```bash
circyto run-manifest \
  --manifest manifest_2.tsv \
  --detector ciri-full \
  --ref-fa ref/chr21.fa \
  --outdir work/cirifull_chr21_manifest2 \
  --threads 8 \
  --parallel 4
```

3. **Collect a circRNA × cell matrix**

```bash
circyto collect \
  --manifest manifest_2.tsv \
  --detector ciri-full \
  --indir work/cirifull_chr21_manifest2 \
  --outdir work/cirifull_chr21_manifest2_matrix
```

Expected outcome (current v0.8.x behavior):

- A sparse MatrixMarket matrix of shape **≈12 × 16** with non-zero entries  
- Matching cell indices between `manifest_2.tsv` and `cell_index.txt`  
- A `circ_feature_table.tsv` containing coordinates and strand information for
  the detected circRNAs

This 16-cell chr21 run is the recommended **sanity check** after installation
to confirm that CIRI-full integration, collection, and indexing are working.

---

## 📦 Installation

### System dependencies

- Python **3.10+**
- `bowtie2`, `samtools`
- Detector binaries as needed (e.g. CIRI-full JAR, `find-circ3`)

### From source (recommended for now)

```bash
git clone https://github.com/liuifrec/circyto
cd circyto
pip install -e .
```

> A PyPI package is planned for the v1.0 release. Until then, please install
> from source as above.

---

## 🗺️ Roadmap

See [`ROADMAP.md`](ROADMAP.md) for details.

- **v0.8.x** – CIRI-full + find_circ3 integration, multimodal export
- **v0.9** – Multi-detector comparison + CIRCexplorer2 integration
- **v1.0** – PyPI release, docs site, deterministic test datasets

---

## 📚 Citation

A methods manuscript is under preparation.  
In the meantime, please cite this repository:

> Liu, Y.-C. *et al.* “circyto: a unified CLI for single-cell circRNA detection and multimodal matrices.”  
> GitHub repository: https://github.com/liuifrec/circyto
