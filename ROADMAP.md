# circyto Roadmap

This roadmap reflects the current development state of **circyto** as of v0.8.x.

---

## ✅ Completed Milestones

### **v0.8.0 – Dual Detector Core**
- Full integration of **CIRI-full** for Smart-seq2 short-read data.
- New **find_circ3** detector:
  - Python3 rewrite of the original find_circ.
  - Two-pass Bowtie2 alignment.
  - Standardized per-cell circRNA output.
- New collector:
  - `collect-find-circ3` → MatrixMarket circ × cell matrix.
  - circ/cell index files.
- CI integration test:
  - Validates end-to-end mini-run on chr21.
- Multimodal exporter:
  - `export-multimodal` attaches circRNA to `.h5ad`.
- Host-gene mapping:
  - `annotate-host-genes` with GTF.

---

## 🎯 Upcoming Milestones

### **v0.9 – Multi-Detector Expansion**
- Integrate **CIRCexplorer2**:
  - FASTQ/BAM support
  - Splice-junction focused circRNA detection
  - New adapter + collector
- Unified multi-detector merge:
  - Combine outputs across:
    - CIRI-full
    - find_circ3
    - CIRCexplorer2
- Overlap metrics:
  - Jaccard (global and per-cell)
  - Venn-style summaries
- Deterministic test dataset:
  - Tiny FASTQs calling ≥1 circRNA for all detectors.
  - CI-safe with `CIRCYTO_ENABLE_HEAVY_TESTS=1`.

---

## 🚧 v0.10 – Multimodal & API Refinements
- Multiple circRNA layers in AnnData:
  - `obsm["X_circ_ciri_full"]`
  - `obsm["X_circ_find_circ3"]`
  - `obsm["X_circ_circ_explorer2"]`
- Detector API cleanup:
  - Standard run model
  - Output schema unification
- Improved logging + `circyto doctor`.

---

## 🧬 v1.0 – Stable Release
- PyPI release  
- Full documentation site  
- Example notebooks:
  - chr21 Smart-seq2 example
  - multi-detector comparison
- Methods manuscript preparation.

---

## 📌 Long-Term
- circyto-ST (spatial) extension  
- circRNA velocity prototype  
- GPU-accelerated junction caller  
- Multi-species detection

---

## Versioning
Semantic versioning:
- MAJOR = breaking changes  
- MINOR = new features  
- PATCH = bugfixes  
