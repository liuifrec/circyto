# Real-object MuData interoperability

**PASS — remedy A, 2026-09-16.** Both checksum-matched originals load through
ordinary `mudata.read_h5mu` and pass native read/write/read equivalence in the
already documented stack. These originals remain the proposed public files;
no compatibility export is required. **The circyto helper separately FAILS on
Smart-seq3**, as detailed below; this is not a claim that all circyto read paths
work. [Machine-readable evidence](real_object_interoperability.json) includes
versions, warnings, the full failure traceback (with portable path prefixes),
input/output hashes, metadata changes, and equivalence results.

The worktree started clean on `manuscript/bioinformatics-application-note-v2`
at `fd52761dcf32436a53633ff7e0f9f30bf7b661cb`. Production `circyto/` and
`pyproject.toml` remain identical to baseline
`44697355bcab1c525ca7ef9b130e2ad0094d9e1b`, version 0.10.0. No biological
results, manuscript v2 assets, embeddings, or raw pipelines were regenerated.

## Inputs, environment and interfaces

Inputs were recovered only from the two [documented archive members](reproduce.md#obtain-the-two-processed-objects)
at `c99cddaae2314d3af7d5232150c5b86bb6b52e25`. Both SHA-256 values matched
before and after the audit:

```text
Smart-seq3  0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339
IMR90       bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8
```

The audit used a fresh Python process in an isolated conda environment:
Python **3.10.20**, MuData **0.3.10**, AnnData **0.11.4**, pandas **2.3.3**,
NumPy **2.2.6**, SciPy **1.15.3**, h5py **3.16.0**. All 13 pins in
[`requirements-results.txt`](requirements-results.txt) were verified; the
evidence also records every installed distribution. circyto was imported
from the checkout. No dependencies, global settings, or production code were
patched.

| Actual object | Ordinary MuData read | circyto read helper | Native read/write/read |
| --- | --- | --- | --- |
| Smart-seq3 | PASS | FAIL: `ValueError: cannot join with no overlapping index names` | PASS |
| IMR90 | PASS | PASS | PASS |

Ordinary reads use MuData 0.3.10's unchanged `pull_on_update=None` default.
Each native read, write, and reread emits two upstream `FutureWarning`s:
**12 warnings total**, all recorded and displayed. No warning filter declares
success. `read_elem` is used only as an independent stored-value reference;
it is not counted as a successful MuData read.

## Confirmed cause and remedy scope

MuData 0.3.10 routes its default through `_update_attr_legacy`; the circyto
helper's scoped `pull_on_update=True` selects the newer `_update_attr` path.
This distinction is visible in the installed source; upstream also documents
the retained legacy default in its [changelog](https://mudata.readthedocs.io/stable/changelog.html).

Smart-seq3 has unique string/object cell indices named `cell_id` in both
modalities, but an unnamed global observation index. The newer path builds
temporary MultiIndices with names `['cell_id', None]` and `[None, None]`, then
joins them. pandas rejects that join. IMR90's observation indices are all
unnamed. All original maps were checked against cell/feature identifiers and
are correct; an index dtype or stale-map diagnosis is not supported.

A controlled Smart-seq3 scratch copy changed **only the global observation
index name to `cell_id`**. All other stored semantics were checked unchanged;
the helper then read it successfully. Its hash and explicit diagnostic-only
purpose are recorded in the evidence. It is not a proposed public object.

The narrow remedy is **A: document and verify the native API on the originals**.
No production fix was made. A separate helper repair would need to reconcile
the temporary join-level names in the explicit synchronization path while
preserving stored axis names and metadata, with a named-index regression.
Changing the global synchronization policy would have a wider scope.

## Equivalence established

Stable identifiers and their ordering were checked before comparing matrices.
Modality inventory, CSR shapes/dtypes/indices/indptr/data (including stored
zeros), dense values, obs/var values/dtypes/categories/missingness, layers,
obsm/varm, obsp/varp, raw where present, and uns were compared exactly. IMR90's
CNV `mappabilitynorm` layer and modality/global provenance were included.
Global obs/var identities and existing values, obsmap/varmap values and their
identifier correspondence, and all global array/graph/uns slots passed.

- Smart-seq3: **192 cells**, RNA **63,187** features, circ **2,503** candidates,
  **2,659** nonzero circ entries, **2,379** annotated candidates, medians
  **12** detected candidates and **22.5** total support. Candidate
  `chr1:117402186|117420649`: **6 cells**, total support **15**.
- IMR90: RNA **23 × 63,187**, circ **23 × 2,443**, CNV **23 × 60,607**,
  overlap **23**, annotated candidates **2,429**.

Native reads pull extra modality columns into global metadata: Smart-seq3
obs **7 → 21** columns and var **4 → 17**; IMR90 obs **16 → 51** and var
**2 → 22**. No stored columns or values are lost. Smart-seq3 global
`var['chrom']` changes from categorical to object **in memory**, preserving
labels and missingness. The complete added-column lists and category metadata
are recorded. Before-write and after-reread global tables match exactly.

At the HDF5 level, **101 / 177 datasets** and **331 / 531 attributes** were
checked for Smart-seq3 / IMR90. Every dataset's shape, dtype and value is
unchanged; the only attribute changes are `encoder-version: 0.3.8 → 0.3.10`
at the root and each modality. Thus stored global tables (including the
categorical chromosome field), maps and all other attributes are preserved.
Scratch round-trip files have new SHA-256 hashes; byte identity is not claimed.

## Reproduce and finish

After the archive recovery in [reproduce.md](reproduce.md), run in the pinned
environment from the checkout root. Large outputs must be in a new directory
outside the checkout:

```bash
INTEROP_WORK=$(mktemp -d)
python scripts/manuscript/verify_real_object_interoperability.py \
  --smartseq3 "$MANUSCRIPT_INPUTS/load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu" \
  --imr90 "$MANUSCRIPT_INPUTS/load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu" \
  --output "$INTEROP_WORK/audit"
PYTHONPATH=. python -m pytest -q \
  tests/test_real_object_interoperability.py tests/test_manuscript_scripts.py
git diff --check
```

`PYTHONPATH=.` exposes the source package to the existing tests' subprocesses
without installing circyto into the analysis environment. The focused gate is
**25 passed, 4 warnings**; those four warnings are the existing intentional
overlapping-feature fixtures. The 13 new regressions exercise the audit,
named-index diagnosis, checksum/output guards, and detection of changed
support, ordering, missingness, layers, provenance, maps, and attributes.
The old **345 passed, 8 skipped, 5 warnings** full-suite result and subsequent
**12 passed, 4 warnings** manuscript result remain historical checkpoints,
not fresh full-suite runs.

**Remaining action:** none for the bounded native-MuData interoperability
gate. Before submission, confirm public access to the two checksum-identified
originals and the availability statement; no upload was authorized or made.
The disclosed Smart-seq3 circyto-helper defect remains a separate production
follow-up. No merge, push, release, or package-version change was made.
