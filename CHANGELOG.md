# 🧬 Changelog
# Changelog

## [0.4.0] - 2025-11-13 (unreleased)
### Added
- **CIRI-full integration**:
  - Support for running CIRI-full v2.x via `circyto run-manifest` with a manifest of paired-end Smart-seq2 FASTQs.
  - New `tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh` adapter that:
    - Handles gzipped FASTQs
    - Runs CIRI-full in Pipeline mode
    - Detects CIRI_output / CIRI-full_output files automatically
    - Normalizes outputs to a simple TSV schema: `circ_id, chr, start, end, strand, support`
  - End-to-end validation on a 16-cell chr21 Smart-seq2 subset (E-MTAB-6072) producing a 12 × 16 sparse matrix with 13 non-zero entries.

- **Integration test**:
  - `tests/test_cirifull_chr21_integration.py`: runs `run-manifest` on a 2-cell chr21 subset and `collect` to ensure:
    - Non-empty matrix (nnz > 0)
    - Number of matrix columns matches the manifest cell count.

### Changed
- Improved error messages and logging in the CIRI-full adapter, including:
  - Sanity checks for `bwa` and `samtools`
  - Reference index presence checks
  - Run directory and recursive file listings for debugging.
- MatrixMarket header uses `general` for `collect` outputs (compatible with Scanpy / scipy.io.mmread).

### TODO / Planned (0.4.x)

- Add `circyto export-h5ad` command to export circRNA × cell matrices as `.h5ad` for direct use with Scanpy.
- Normalize circRNA IDs to a consistent internal format: `chr:start|end|strand`.
- Add initial support and documentation for additional detectors (CIRI-long, CIRCexplorer2, find_circ).

## [v0.3.2] — 2025-11-07
**Fixes and Robustness**
- 🧩 Fixed `to_loom()` to safely handle all empty-matrix shapes (`0×0`, `0×N`, `N×0`) and produce valid `.empty.tsv` outputs instead of raising Pandas errors.
- 🧬 Improved H5AD export in `_to_h5ad()` — matrices are now correctly oriented as `(cells × circs)`, consistent with `AnnData` expectations.
- 🚦 All dummy and smoke tests (Smart-seq2 + PBMC v3) complete successfully with graceful empty fallbacks.

**Developer Improvements**
- CI workflow stabilized on `ubuntu-latest` runners.
- `black` auto-formatting integrated for consistent style.
- Added more robust optional-dependency handling (`loompy`, `anndata`).

---

## [v0.3.1] — 2025-11-06
**Additions**
- Introduced `circyto make` all-in-one command for both plate-based and 10x workflows.
- Added CI badge and version badge to `README.md`.
- Improved command help text and documentation examples.

---

## [v0.3.0] — 2025-11-05
**Initial Public Release**
- Core CLI commands (`prepare`, `run`, `run-manifest`, `collect`, `convert`) implemented.
- Modular detector support: CIRI-full, CIRI-long, find_circ, CIRCexplorer2.
- Output generation: circ×cell sparse matrix, `.loom` and `.h5ad` export.
- Ready for smoke testing on Smart-seq2 and 10x datasets.
