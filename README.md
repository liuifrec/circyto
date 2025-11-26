# circyto

**circyto** is a unified Python CLI framework for **single-cell circRNA detection**,  
**detector orchestration**, **matrix generation**, and **multimodal export**.

It provides a clean, reproducible interface over multiple circRNA detectors  
and produces standardized outputs suitable for downstream single-cell analysis  
(Scanpy, Seurat, scVI, etc.).

---

## ✨ Features (v0.8.x)

### **Core**
- Unified CLI for single-cell circRNA workflows
- Detector drivers for:
  - **CIRI-full** (short-read, full-length circRNA assembly)
  - **find_circ3** (modernized Python rewrite of find_circ; junction detector)
- Manifest-driven batch execution
- Sparse MatrixMarket circRNA × cell export
- Host-gene annotation using GTF
- Multimodal `.h5ad` export  
  (`.X` = mRNA, `obsm["X_circ"]` = circRNA, `uns["circ"]` = metadata)

### **Pipeline Overview**

```
FASTQ → run-detector → per-cell circRNA calls → collect → circ_matrix.mtx
                                                  ↓
                                            annotate-host-genes
                                                  ↓
                                      export-multimodal → AnnData
```

### **Multi-detector Support**

| Detector | Type | Status | Notes |
|----------|------|--------|-------|
| **CIRI-full** | Full-length detection | ✔ Stable | Best for Smart-seq2 / long reads |
| **CIRI2**     | Short-read | ✔ Supported | Adapter included |
| **find_circ3** | Junction detector | ⚗️ Experimental | Validated on chr21 mini; Python3 rewrite |
| **CIRCexplorer2** | Splice junction detector | *In progress* | Next integration target |

> _Legacy find_circ has been removed in favor of **find_circ3**._

---

## 🚀 Quick Start

### 1. Run a detector
```bash
circyto run-detector find-circ3   --manifest manifest_2.tsv   --outdir work/find_circ3_chr21   --ref-fa ref/chr21.fa   --threads 4 --parallel 2
```

### 2. Collect circRNA counts
```bash
circyto collect-find-circ3   --indir work/find_circ3_chr21   --outdir work/find_circ3_chr21_matrix
```

### 3. Annotate host genes
```bash
circyto annotate-host-genes   --circ-tsv circ_feature_table.tsv   --gtf Homo_sapiens.gtf   --out circ_feature_table_annotated.tsv
```

### 4. Export multimodal AnnData
```bash
circyto export-multimodal   --rna-h5ad rna.h5ad   --circ-mtx circ_counts.mtx   --circ-index circ_index.txt   --cell-index cell_index.txt   --circ-features circ_feature_table.tsv   --out combined_multimodal.h5ad
```

---

## 📦 Installation

System dependencies:
- Python 3.10+
- `bowtie2`, `samtools`
- Detector binaries (e.g., CIRI-full JAR)

Install:
```bash
pip install circyto
```

Dev install:
```bash
git clone https://github.com/liuifrec/circyto
cd circyto
pip install -e .
```

---

## 🧪 Testing

### Unit tests
```bash
pytest -q
```

### Integration tests
(require bowtie2, samtools, find-circ3)

```bash
pytest -m "integration"
```

Skip heavy tests:
```bash
export CIRCYTO_SKIP_INTEGRATION=1
```

---

## 📅 Roadmap

See [`ROADMAP.md`](./ROADMAP.md).

- **v0.8.x:** CIRI-full + find_circ3 integration  
- **v0.9:** Multi-detector comparison + CIRCexplorer2  
- **v1.0:** PyPI release, docs site, deterministic datasets

---

## Citation

Please cite the repository; a manuscript is under preparation.
