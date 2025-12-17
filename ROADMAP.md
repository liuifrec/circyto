# Circyto roadmap

This document is **not a tutorial**. It describes the product strategy, major milestones, and the “why” behind planned changes.

For step-by-step usage, see `docs/getting_started.md`.

---

## North star

**A new user can go from “I have FASTQs” → “I have a circRNA × cell matrix”** without reading developer logs, spelunking through source code, or guessing CLI conventions.

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

### v0.8.x — “Docs + usability to 100%” (current priority)

**Goals**
- README becomes a single “golden path” entry point (install → minimal example → next steps).
- Introduce `docs/` pages that are the canonical reference for workflows and CLI semantics.
- Reduce user-facing friction caused by external dependencies and hidden detector lists.

**Deliverables**
- `circyto detectors`: list detectors and their short descriptions + required external tools.
- `circyto doctor`: validate external dependencies (bwa/bowtie2/samtools/java; STAR optional) and print actionable messages.
- Move detailed workflows and edge cases into:
  - `docs/getting_started.md`
  - `docs/cli_policy.md`
  - `docs/detectors.md`
- CI hygiene: workflows run on **PRs** and **manual dispatch**, not on every push (avoid noisy notifications).

---

### v0.9.0 — “Detector plugin ergonomics”

**Goals**
- Make it easier to add detectors without touching unrelated code.
- Strengthen validation around input formats (manifest schema, reference files, etc.).

**Possible deliverables**
- A more explicit detector registration / metadata model (name, version, dependencies, capabilities).
- Structured output for `circyto detectors --json` to support tooling.
- Compatibility checks between detector outputs and downstream collectors.

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
