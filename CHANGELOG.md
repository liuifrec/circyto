# 🧬 Changelog
## [v0.8.0] - 2025-11-20

### Added
- End-to-end **CIRI-full integration** via `circyto run-manifest` + `circyto collect`, with a chr21 Smart-seq2 test plate.
- New **chr21 Smart-seq2 integration test** (`test_cirifull_chr21_integration.py`) validating a non-empty circ × cell matrix.
- Convenience scripts and docs for downloading ENA FASTQ subsets with `aria2c` (`make_fastq_urls.py` and usage notes).

### Changed
- Normalized and simplified `writers/convert.py` so matrix/AnnData conversion behaves consistently across commands.
- Cleaned and unified the collect pipeline (`pipeline/collect.py`), including consistent index writing and min-count handling.
- Updated command templates for CIRI-full (`run_cirifull_with_manifest`) to work cleanly under WSL/Ubuntu and local Linux.

### Removed
- Large meta-data files `E-MTAB-6072.idf.txt` and `E-MTAB-6072.sdrf.txt` from the repository.
- Kept only the **chr21 subset reference** (`ref/chr21.fa`, `ref/chr21.gtf`) and manifest files needed for tests.

### Documentation
- Updated `README.md` with WSL instructions and clarified how to run the chr21 integration tests locally.
- Extended `.gitignore` to exclude experiment metadata and runtime artifacts (e.g., `E-MTAB-6072*.txt`, `tests/out_*`, `tests/*.mtx`).

# Changelog — v0.7.0

## Added
- `run-multidetector` command
- `merge-detectors` TSV union + by-detector matrix builder
- `compare-detectors` Jaccard similarity and detector summary
- Proper summary.json generation
- Robust CLI argument handling

## Fixed
- Missing summary.json from run-multidetector
- Test-suite issues related to detector summary naming

## v0.6.0 — Multi-detector engine & CIRI2 integration (2025-11-18)

This release introduces the first major extension beyond single-detector workflows, adding a unified detector API, experimental CIRI2 support, and a fully functional multi-detector orchestration CLI.

### ✨ New Features

#### 1. Detector Engine API
- Added `DetectorBase`, `DetectorRunInputs`, and `DetectorResult`.
- Standardized normalized TSV output format:  
  `circ_id, chr, start, end, strand, support`.

#### 2. CIRI-full Engine
- Rewritten `CiriFullDetector` with:
  - Clean adapter scripting
  - Automatic output discovery
  - Robust logging + sanity checks
  - Normalization into the standard TSV format

#### 3. Experimental CIRI2 Support
- Added `circyto/detectors/ciri2.py`
- Wrapper for `CIRI2.pl`
- Supports low-stringency mode (`-w`)
- Automatically normalizes output

> Note: CIRI2 remains strict for single-cell RNA-seq and may produce sparse results.

#### 4. New CLI: `run-multidetector`
Run multiple detectors on the same manifest:

```bash
circyto run-multidetector ciri-full ciri2 \
  --manifest manifest.tsv \
  --outdir work/multi \
  --ref-fa ref.fa --gtf genes.gtf \
  --threads 8 --parallel 1
5. Utils

Added ensure_dir, read_tsv, write_tsv

Added build_default_engines() to enumerate detectors

🐛 Fixes & Improvements

Improved adapter environment handling

Path robustness in Codespaces

Deterministic output directory naming

Clearer detector availability messages

⚠ Known Issues

CIRI-full requires --parallel 1 for stability on Codespaces

CIRI2 strict filtering often eliminates marginal circRNAs

Multi-detector merging not yet implemented — planned for v0.7.0

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
