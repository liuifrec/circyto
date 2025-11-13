# circCyto
<p align="center">
  <a href="https://github.com/liuifrec/circyto/actions/workflows/test.yml">
    <img src="https://github.com/liuifrec/circyto/actions/workflows/test.yml/badge.svg?branch=main" alt="CI status">
  </a>
  <a href="https://github.com/liuifrec/circyto/releases">
    <img src="https://img.shields.io/github/v/release/liuifrec/circyto?display_name=tag" alt="Latest release">
  </a>
  <a href="https://github.com/liuifrec/circyto/issues">
    <img src="https://img.shields.io/github/issues/liuifrec/circyto" alt="Issues">
  </a>
  <a href="https://github.com/liuifrec/circyto/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/liuifrec/circyto" alt="License">
  </a>
</p>

<p align="center">
  <img src="assets/circCyto_logo.png" alt="circyto logo" width="420">
</p>

<p align="center">
  <b>circyto</b> — a modular CLI toolkit for single-cell circRNA detection,<br/>
  integrating multiple circRNA detectors into a unified sparse matrix workflow.
</p>


---

## Overview

`circyto` is a Python CLI framework to:

- run external circRNA detectors (CIRI-full, CIRI-long, CIRCexplorer2, find_circ, …)
- normalize their outputs to a common circRNA schema
- assemble sparse circRNA × cell matrices
- export results for Scanpy / Seurat / ML downstream analysis.

Core ideas:

- Detector-agnostic adapters
- Manifest-based reproducible job execution
- Sparse matrix output (MatrixMarket)
- Planned `.h5ad` export and multi-detector benchmarking.
---
### Detector comparison

| Detector        | Input type              | Index / reference needed          | Output granularity         | Strengths                                                     | Limitations / Notes                                           | circyto support level |
|----------------|-------------------------|-----------------------------------|----------------------------|----------------------------------------------------------------|---------------------------------------------------------------|------------------------|
| **CIRI-full**  | Paired-end FASTQ        | Genome FASTA + BWA index + GTF    | Full-length circRNA + AS   | Full-length reconstruction, integrates RO/AS, detailed output | Heavy, multi-stage pipeline; requires good reference indices | ✅ Fully integrated (chr21 Smart-seq2, `run-manifest` + `collect`) |
| **CIRI-long**  | Long-read FASTQ (ONT)   | Genome FASTA + minimap2 index     | Full-length circRNA        | Designed for long-read; captures complex isoforms             | Long-read only; runtime and memory depend on read length     | 🔜 Planned integration |
| **CIRCexplorer2** | BAM (spliced alignments) | Genome FASTA + annotation        | Back-splice junctions      | Widely used; works from alignment BAMs                        | No full-length reconstruction; depends on upstream aligner   | 🔜 Planned integration |
| **find_circ**  | FASTQ (or BAM)          | Genome FASTA + BWA index          | Back-splice junctions      | Simple, fast; classic circRNA detector                        | Older; higher FP rate; no full-length info                   | 🔜 Planned integration |
| **circRNA_finder** (optional) | BAM    | Genome FASTA                      | Back-splice junctions      | Works directly on STAR outputs                                | STAR-specific; less widely maintained                        | ❓ Maybe (low priority) |

**Support levels**

- ✅ **Fully integrated**: wired through `run`/`run-manifest` + `collect`; tested end-to-end.
- 🔜 **Planned**: CLI + adapter design ready, implementation not finalized.
- ❓ **Maybe**: candidate for future integration based on demand.

## 🧬 Installation

```
git clone https://github.com/liuifrec/circyto.git
cd circyto
pip install -e .
```

Optional extras:

```
pip install -e .[detectors]
```

---

## ⚙️ Commands

| Command | Description |
|----------|--------------|
| `circyto prepare` | 10x BAM → FASTQs (`--chemistry tenx-3p|tenx-5p`) |
| `circyto run` | Run CIRI-full on batches produced by `prepare` |
| `circyto run-manifest` | Run CIRI-full using a manifest TSV (per-cell FASTQs) |
| `circyto collect` | Merge all CIRI outputs into a circ × cell matrix |
| `circyto convert` | Convert `.mtx` + index files to `.loom` or `.h5ad` |
| `circyto make` | All-in-one pipeline (detects manifest or BAM input) |

---


## 📖 Citation

Liu, Y.-C. *et al.* **circCyto** — a modular toolkit for single-cell circRNA profiling and integration with scRNA-seq workflows.  
*(Manuscript in preparation, 2025).*

---

## 📜 License

MIT License © 2025 **Yu-Chen (James) Liu**

