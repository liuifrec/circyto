# circyto Roadmap

This document describes the planned development trajectory of **circyto**, a Python CLI toolkit
for single-cell circRNA detection, integration, and detector comparison.

The roadmap emphasizes **stability, reproducibility, and scientific credibility** over rapid
feature expansion.

---

## Current Status

**Latest stable version:** v0.8.2  
**Development phase:** Stabilized core, validated detector comparison

As of v0.8.x, circyto provides:

- Unified and locked CLI semantics
- Multi-detector execution (`run-detector`, `run-multidetector`)
- Per-detector and unified matrix collection
- Detector comparison with exact and fuzzy circRNA recovery
- Regression-tested recovery of CIRI-full calls by find-circ3
- AnnData multimodal export (`obsm["X_circ"]`)
- Fully passing test suite (unit + integration + regression)

No breaking CLI or data-format changes are planned within the v0.8 series.

---

## Versioning Policy

- **v0.8.x**: Stability and reproducibility
  - No breaking changes
  - No semantic drift in CLI behavior
  - Documentation, metadata, and ergonomics only
- **v0.9.0**: Scientific feature upgrade
  - New circRNA identity model
  - Enhanced detector comparison metrics
  - Publication-ready outputs

---

## Roadmap Overview

### v0.8.x — Stabilization & Trust (Now → Early 2026)

**Primary goal:** Make circyto boring, predictable, and reviewer-proof.

#### Scope
- Lock all CLI commands and flags
- Ensure README and examples exactly match tested behavior
- Improve discoverability and diagnostics without changing results

#### Planned Additions
- `circyto version`
- `circyto doctor`
- Embedded metadata in outputs
- Expanded documentation

#### Non-Goals
- No new detectors
- No changes to circRNA definition or matching logic
- No changes to matrix formats or AnnData layout

---

### v0.9.0 — Scientific Upgrade Release (Target: Spring 2026)

**Primary goal:** Elevate circyto from a pipeline to a scientific framework.

#### 1. CircRNA Identity Model
Formal detector-agnostic circRNA identity with exact and fuzzy equivalence modes.

#### 2. Advanced Detector Comparison
Recall, recovery rates, asymmetric matching, consensus circRNAs.

#### 3. Multimodal Integration Enhancements
Richer circRNA metadata and structured AnnData storage.

#### 4. Publication-Ready Reporting
Tables, plots, and supplementary-material–ready outputs.

---

## Tentative Timeline

| Period | Focus |
|------|------|
| Now → Dec | v0.8.2 stabilization |
| Jan | Diagnostics and documentation |
| Feb | Metadata & reproducibility |
| Mar–Apr | CircRNA identity model |
| Apr–May | Detector comparison upgrade |
| May–Jun | v0.9.0 release |

---

## Design Principles

1. Stability before features  
2. Explicit semantics  
3. Regression tests as scientific claims  
4. Detector-agnostic core  
5. Publication-oriented outputs  

---

## Long-Term Vision

circyto aims to become a reference framework for circRNA detector benchmarking
and multimodal circRNA analysis in single-cell genomics.
