<p align="center">
  <img src="https://raw.githubusercontent.com/liuifrec/circyto/main/assets/circCyto_logo.png"
       alt="circyto logo" width="440" />
</p>

<h1 align="center">circyto: Single-Cell circRNA Detection Suite</h1>

<p align="center">
  <a href="https://github.com/liuifrec/circyto/actions/workflows/ci.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/liuifrec/circyto/ci.yml?branch=main"
         alt="CI status" />
  </a>
  <img src="https://img.shields.io/badge/python-3.10%20|%203.11%20|%203.12-blue"
       alt="Supported Python versions" />
  <img src="https://img.shields.io/badge/status-experimental-orange"
       alt="Project status: experimental" />
</p>

---

`circyto` is a multi-detector circRNA discovery pipeline for single-cell datasets, supporting:

- **CIRI-full** (Smart-seq2 / full-length)
- **find_circ3**
- **CIRCexplorer2**
- Multi-detector comparison & merging
- Smart-seq2 / 10x support
- Matrix export (MatrixMarket + AnnData multimodal exporter)
- Host-gene annotation
- Basic demo + integration tests for each detector

------------------------------------------------------------------------

## 🚀 Features

### **Detector Support**

  ----------------------------------------------------------------------------
  Detector        Single‑Cell            Paired‑End           Notes
  --------------- ---------------------- -------------------- ----------------
  CIRI-full       ✅                     ✅                   Best for
                                                              full‑length &
                                                              Nanopore‑style
                                                              logic

  find_circ3      ✅                     ⚠️ (PE recommended)  Fast, sensitive

  CIRCexplorer2   ⚠️ Smart‑seq2 only     Requires STAR        Full
                                                              parse/annotate
                                                              workflow
  ----------------------------------------------------------------------------

### **Multidetector Pipeline**

``` bash
circyto run-multidetector ciri-full find-circ3 circexplorer2   --manifest manifest.tsv   --ref-fa genome.fa   multi_output/
```

### **Unified Matrix & AnnData Export**

-   circRNA matrix → `obsm["X_circ"]`
-   circRNA feature table → `uns["circ"]["feature_table"]`
-   Host‑gene map → `uns["circ_host_map"]`

------------------------------------------------------------------------

## 📦 Installation

``` bash
pip install circyto
```

Or development install:

``` bash
git clone https://github.com/your-org/circyto.git
cd circyto
pip install -e .
```

------------------------------------------------------------------------

## 🔧 Quickstart

### 1. Prepare a manifest

    cell_id    r1                         r2
    C1         fastq/C1_1.fastq.gz        fastq/C1_2.fastq.gz

### 2. Run a single detector

``` bash
circyto run-detector ciri-full   --manifest manifest.tsv   --ref-fa genome.fa   out_ciri/
```

### 3. Merge Detector Outputs

``` bash
circyto compare-detectors   --detector-dirs out_ciri out_findcirc out_circexp   merged_out/
```

### 4. Export Multimodal AnnData

``` bash
circyto export-multimodal   --rna input_gene.h5ad   --circ-dir merged_out/   --out combined_multimodal.h5ad
```

------------------------------------------------------------------------

## 📁 Repository Structure

    circyto/
      detectors/        # detector adapters
      pipeline/         # run, collect, export
      tools/            # bundled binaries
      tests/            # full CI-safe test suite

------------------------------------------------------------------------

## 🧪 Test Suite

Run all tests:

``` bash
pytest -vv
```

Run only integration tests:

``` bash
pytest -vv tests/test_*integration*.py
```

------------------------------------------------------------------------

## 📝 License

MIT License.

------------------------------------------------------------------------

## 👤 Author

James Liu (劉祐誠) --- 2025\
IFReC → Translational Data Science\
Single‑cell \| circRNA \| AI‑accelerated genomics
