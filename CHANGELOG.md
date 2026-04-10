# Changelog

All notable changes to **circyto** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and semantic-ish versioning.

---
## [Unreleased]

### Changed

- Clarify public documentation for the current **v0.8.x** release line instead of older prototype/planning framing.
- Document that `circyto doctor` and `circyto detectors` are live commands and that bundled asset resolution is intended to be cwd-independent.
- Align README and workflow docs with the current CLI examples and flag names, including `--ref-fa`.
- Document alignment-first execution and real CIRI3 integration with explicit environment and setup guidance.
- Document CIRI3 readiness semantics for `circyto doctor` / `circyto detectors`, including direct `java -jar` readiness versus incomplete local contracts.
- Document validated local BWA + CIRI3 execution and distinguish it from STAR support that is implemented in code but not yet fully validated end-to-end.
- Document explicit CIRI3 external requirements: Java, a CIRI3 jar, `samtools`, and mode-specific `bwa` or `STAR`.
- Document validated local BWA + CIRI3 settings for single-cell alignment-first runs: unsorted SAM, `-S 0`, `-Ma 0`, and `bwa mem -k 15 -T 15`.
- Document STAR + CIRI3 as supported in code for alignment-first workflows, but not yet validated end-to-end in this release.

### Fixed

- Document the real `ciri-full` layout semantics:
  - paired-end rows use the upstream bundled **CIRI-full Java Pipeline**
  - single-end rows use a bundled **CIRI2-based fallback path**
  - both layouts normalize to the same TSV schema for downstream collection
- Record the recent `ciri-full` runtime diagnostics improvements in user-facing release notes.
- Fix CIRI3 alignment-first execution selection so an explicit `--command-template` takes precedence over direct `java -jar` mode.
- Clarify that explicit template execution does not require `--ref-fa` unless the template itself uses `{ref_fa}`, while direct execution still requires `--ref-fa`.

### Notes

- BWA + CIRI3 has been validated locally on a chr21 pilot with non-zero output.
- STAR + CIRI3 support is present in code for alignment-first workflows, but is not yet fully validated end-to-end in this release.

## [0.8.3] – 2025-12-17
- Add `circyto doctor` to validate external dependencies and report detector readiness.
- Add `circyto detectors` to list detectors with status and dependency requirements.

## [0.8.2] – 2025-12-16

### Added

* Stable **fuzzy circRNA recovery and comparison** utilities (`compare-ids`, fuzzy Jaccard).
* Robust circRNA ID normalization across detectors (CIRI-full, find-circ3, CIRCexplorer2).
* Unified multi-detector workflows: run → collect → merge → compare.

### Changed

* Locked CLI semantics for detector runners and collectors (`--outdir/-o`, `--indir`).
* Improved detector comparison logic to avoid strand and coordinate mismatches.
* README badges and metadata aligned with current release structure.
* CI workflow quieted and stabilized for local and PR-based development.

### Fixed

* Fuzzy matching regression where CIRI-full calls were not recovered by find-circ3.
* circRNA ID parsing inconsistencies across `:` / `-` / `|` separators.
* Missing LICENSE file and incorrect license badge rendering.

### Notes

* Heavy integration tests (STAR / CIRCexplorer2) remain optional and skippable via environment flags.
* This release represents a **stability + correctness milestone** rather than new detector features.

---

## [0.8.1]

* Previous incremental fixes and detector integration updates.

## [0.8.0]

* Initial unified detector API and multi-detector orchestration.
