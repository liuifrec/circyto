# Changelog

All notable changes to **circyto** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and semantic-ish versioning.

---
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
