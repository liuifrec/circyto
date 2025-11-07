# 🧬 Changelog
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
