# Circyto roadmap

This document is **not a tutorial**. It describes the product strategy, major milestones, and the “why” behind planned changes.

For step-by-step usage, see `docs/getting_started.md`.

---

## North star

**A new user can go from “I have FASTQs” → “I have a circRNA × cell matrix”** without reading developer logs, spelunking through source code, or guessing CLI conventions.

Current release line: **v0.10.0** (experimental).

---

## Guiding principles

1. **One golden path**  
   README stays short and points to the right doc pages. Detailed workflows live in `docs/`.

2. **Explicitness beats cleverness**  
   Prefer explicit subcommands (`circyto detectors`, `circyto doctor`) over help text that Typer can’t express well.

3. **Reproducibility is a feature**  
   Every run is manifest-driven, path-explicit, and produces standardized outputs that can be post-processed consistently.

4. **Detector integration is modular**  
   Detectors have different external dependencies. Circyto should report what is missing up front and keep the orchestration layer stable.

---

## Milestones

### v0.8.x — “Docs + usability to 100%” (completed)

**Goals**
- README becomes a single “golden path” entry point (install → minimal example → next steps).
- Introduce `docs/` pages that are the canonical reference for workflows and CLI semantics.
- Reduce user-facing friction caused by external dependencies and hidden detector lists.
- Make layout-specific detector behavior explicit instead of implying uniform support where upstream tools differ.

**Deliverables**
- `circyto detectors`: list detectors and their short descriptions + required external tools.
- `circyto doctor`: validate external dependencies (bwa/bowtie2/samtools/java; STAR optional) and print actionable messages.
- Document that `ciri-full` has two public execution modes:
  - paired-end rows: upstream bundled CIRI-full Pipeline
  - single-end rows: bundled CIRI2-based fallback normalized to the same TSV schema
- Move detailed workflows and edge cases into:
  - `docs/getting_started.md`
  - `docs/cli_policy.md`
  - `docs/detectors.md`
- CI hygiene: workflows run on **PRs** and **manual dispatch**, not on every push (avoid noisy notifications).

---

### v0.9.0 — “SMART-Seq3 workflow + QC + AnnData + annotation” (completed)

**Goals**
- Ship a usable end-to-end SMART-Seq3 workflow around CIRI3.
- Export analysis-ready QC tables and AnnData artifacts.
- Annotate circRNAs against external databases while preserving the existing matrix semantics.

**Completed deliverables**
- `workflow smartseq3-ciri3` with demux, alignment-cache preparation, detector execution, matrix collection, and resumable stage outputs.
- QC outputs including `cell_qc.tsv`, `circ_qc.tsv`, and workflow summaries.
- AnnData export via `circ_counts.h5ad`.
- `annotate-circs` for generic circRNA database annotation and known-versus-novel labeling.
- Documentation for real-data E-MTAB-8735 usage and annotation summaries.

---

### v1.0.0 — “Stable contract”

**Goals**
- Freeze a stable CLI contract and output schema that downstream tooling can rely on.
- Clearly version and document backward-incompatible changes.

**Possible deliverables**
- “Stable interface” guarantee for core commands:
  - `run-batch`, `run-detector`, `run-multidetector`
  - `collect-matrix`
  - multimodal export path (`annotate-host-genes`, `export-multimodal`)
- Versioned schemas for:
  - per-cell detector outputs
  - matrix artifacts (`circ_counts.mtx`, index files, metadata)

---

## Future / under consideration

STAR-based detectors are planned/under consideration for integration (e.g., DCC, circhunter, CIRI3). These will likely be gated behind optional dependencies and may arrive after the usability milestone.

---

## Release & maintenance notes (process)

- Keep git tags, `pyproject.toml` version, and GitHub releases aligned.
- Prefer short, actionable changelog entries aimed at users (not only developers).
