# Detector-free import validation

Date: 2026-10-13. Base: `17867f18f4ee806136ce465415ed075e732b9f7b`.
Branch: `feature/qc-report-usability`. Package version unchanged at 0.10.0.
This is one interoperability improvement, outside the frozen manuscript baseline.

## Empirical evidence

The executable audit is [`scripts/validate_import_circ_evidence.py`](../scripts/validate_import_circ_evidence.py).
It reads historical Git blobs into a new temporary directory, verifies the two
previously published object checksums before loading, exercises the actual CLI,
and compares sparse values, axis order, every original candidate annotation
column, source cell metadata, RNA values/metadata, and exact matched-cell overlap.
It never writes to the source files or regenerates manuscript analyses.

```bash
python scripts/validate_import_circ_evidence.py \
  --repo /path/to/circyto --out /tmp/circ-import-evidence
```

The repository must contain the historical commits. Run the copied script from
outside the checkout with the installed wheel on PATH to test distribution
behavior. Results are in `validation.json` and individual import bundles.

| Input | Cells | Candidates | Nonzero entries | Total support | Annotated host genes | RNA overlap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Smart-seq3 E-MTAB-8735 | 192 | 2,503 | 2,659 | 5,711 | 2,379 | 192 |
| IMR90 GSE278958 | 23 | 2,443 | 2,914 | 7,422 | 2,429 | 23 |

Both historical H5MU files come from
`c99cddaae2314d3af7d5232150c5b86bb6b52e25`:

- `load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu`, SHA-256
  `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.
- `load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu`, SHA-256
  `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.

Smart-seq3 median candidates/cell remains 12; median support remains 22.5.
`chr1:117402186|117420649` retains host gene MAN1A2, six detecting cells and 15
support. IMR90 medians are 119 candidates and 317 support. Source feature order,
chromosomes, starts/ends/strands, and all 14 annotation fields are preserved.
Stored integer counts become int64 without changing any value. Original RNA
`X`, observation/feature annotations, and cell order survive unchanged.

These are **processed-object round trips**, using explicitly named matrix
extractions, not reconstructed detector source files or independent CIRI3
acceptance. Genome paths identify hg38; the CIRI coordinate convention is
declared as 1-based closed and retained without arithmetic, not validated anew
against the reference. No per-cell detector logs were recovered for these
objects, so every imported status is `unknown` and detection fractions are
unavailable. `DIYHEK_192` retains its source RNA-only label and padded zero circ
row without becoming a completed zero. The 192-row overlap is physical object
membership, not 192 independently established circ-processed cells.

Additional historical collector acceptance:

| Source | Evidence identity | Result and scope |
| --- | --- | --- |
| CIRI2 ERR1041421 | `a2c934ec1a629528aae69dfb61df10af0f1407cd:tests/ciri2_err1041421_integration/ERR1041421.tsv`; SHA-256 `7fc1be53e7c141d47b99fce7dd4d034011b0515638af3ac7f66cf4755304977c` | Six candidates, 11 support; every coordinate and support value agrees with the committed normalized public-accession detector result. This tests the existing normalized collector, not a native CIRI2 parser. |
| find-circ ERR2139486 | `bfb52c5499f4c64ac0a360e180d7acd6bef29076:tests/find_circ3_manifest2_integration/ERR2139486/ERR2139486_splice_sites.bed`; SHA-256 `dc292574b65401ec06cdb7a10adad3f4716e3b8731fc560b82e726df8a77cd41` | Original 86,217-row archived output; 27 retained candidates and 27 unique-bridge support under the existing flag filter. An independent row scan checks filtering/support against collector → importer results. The feature table is taken from the original accepted rows. |

The current find-circ fixture was edited in a later commit and is **not** used
as original external evidence. Assemblies of these older chr21 detector examples
are explicitly labeled unverified, with their distinct coordinate/count units
retained. They are not combined. Per-cell completion remains unknown.

**Pending:** original native CIRI3 and CIRCexplorer call-file acceptance, and
independent external-user validation across additional detector versions.
No broad all-detector interoperability claim is supported. The tests demonstrate
specific archived collector compatibility and processed-value preservation.

## Distribution and automated checks

The complete default CI-style suite passes: **441 passed, 10 skipped, 7 warnings**
in 18.40 seconds, with `CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q .`.
This includes **39 importer tests**. Five warnings are pre-existing MuData
variable-name/AnnData index warnings; two are expected duplicate-cell warnings
from the deliberately invalid RNA input test. External detector execution is
not part of the default gate. No scientific thresholds or algorithm tests were
modified.

Tests cover all observation states, completed zeros, partial failed-cell support,
no-completed-cell denominators, literal IDs, first-appearance candidate order,
metadata order, square-matrix orientation, arbitrary source index order, missing
cells, duplicate IDs/pairs/matrix entries, invalid counts/intervals/strands,
conflicting declarations, integer precision above 2**53, overflow rejection,
source changes during import, provenance preservation, RNA mismatch/overlap,
RNA layers/embeddings, safe output staging, and escaped deterministic HTML.

Wheel and source-distribution builds pass with:

```bash
python -m build --no-isolation --outdir /tmp/circyto-import/packages
```

A fresh virtual environment without system-site packages installed the wheel
and its declared `[mudata]` dependencies. `pip check` passes. The executable
example, seven CLI/help checks, and the full archived-evidence audit pass from
`/tmp`, with imports confirmed under the environment's `site-packages`.
The installed Python stack is Python 3.10.20, AnnData 0.11.4, MuData 0.3.10,
NumPy 2.2.6, SciPy 1.15.3, pandas 2.3.3 and Typer 0.27.3. The repository test
stack uses Typer 0.27.0. Native MuData reads emit the documented upstream
FutureWarnings; neither a custom reader nor warning suppression is required.

Installed help checks: top-level help/version, `import-circ`, `collect-matrix`,
`export-mudata`, `export-multimodal`, and `report`. Ordinary `anndata.read_h5ad`
and `mudata.read_h5mu` load outputs. Repeating an import of the same inputs gives
identical JSON, TSV, index, Matrix Market and HTML bytes, with equal sparse
scientific values; HDF5 byte identity is not asserted. The demo's denominator
is 2 completed cells, of which 1 has support; the failed cell retains support 2
without entering the numerator or denominator.

The installed example's HTML was rendered and visually checked in headless
Chrome at 1440 × 1400 with hostname resolution disabled. It displays the
interpretation caveat, exact 1/2 fraction, six distinct states and original cell
order without external assets or scripts. HTML/JSON content is also checked
programmatically. Screenshot and local logs are under `/tmp/circyto-import/`.

Archive inspection confirms both importer modules in the wheel and the example,
audit script and documentation in the source archive; no historical H5MU or
manuscript recheck files are packaged. Focused Ruff checks pass (ignoring B008
for the repository's Typer Option default convention), and `git diff --check`
passes. The [committed audit JSON](examples/import_circ/validation.json) records
the installed empirical results and source identities.

## Preservation and submission decision

Existing detection algorithms, thresholds, source datasets, BioRender Figure
1A, main and manuscript branches are untouched. The untracked manuscript
recheck tree is excluded from staging and checked against its initial hashes,
sizes, modes and modification times: all 72 entries (including the root directory)
match. Main and manuscript refs also match the initial snapshot. No software version bump, merge or release.

The existing CIRI-full license and CIRI3 corresponding-source/embedded-notice
concerns in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) remain separate
release blockers; third-party distribution is unchanged. Detector-free use
requires no detector runtime, but does not resolve the wheel's existing bundled
third-party distribution concerns.

Recommendation: this supports a stronger **reuse and observation-aware
interoperability** argument, not a new calling method or biological result.
It does not by itself warrant changing the manuscript novelty claim. Original
external source acceptance and the limited count/coordinate contract must remain
visible if the feature is described in a later reviewed software snapshot.
