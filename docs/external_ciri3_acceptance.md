# External per-cell CIRI3 + matched RNA acceptance

Status: **PENDING** (2026-10-13). This gate requires an independently generated
biological per-cell result and RNA from the same cells. The official native
single-sample demo establishes format compatibility only. A renamed RNA fixture,
a processed-H5MU-derived table, or a simulation cannot close this gate.

## Evidence searched

The integration review searched accessible project/data files, the recovered
historical objects in `/tmp/circyto-import/load_work`, official CIRI3 example
artifacts, and matching filenames across all local Git refs. It also checked
GitHub Actions artifacts (zero retained at search time) and release assets
(v0.8.3, v0.8.2 and v0.3.1: none). Earlier broader searches are recorded in
[reviewer readiness](reviewer_readiness.md). No authentic per-cell native CIRI3
result with matching RNA or completed public-data workflow QC directory was
found. This describes accessible evidence, not a claim that no such files exist.

Available evidence stays separate:

- Authentic upstream native demo: six candidates, 50 support; SHA-256
  `4885fc00bd30817ebdf81ba9f3344c48ddb18bcb56d2a077cee2297ff7f594ed`.
  Producing version unknown; format reviewed in official v3.0.1 tree
  `73108c4478a7e691cb93360921b65203e0701516`. Biological single-cell identity,
  assembly and matched RNA are not established.
- Official annotation example has simulated read IDs and is explicitly labeled.
- Smart-seq3 and IMR90 archived H5MU objects support checksum-verified processed
  matrix round trips and RNA/circ identity checks, not native-parser acceptance.
  See [exact checksums and values](import_circ_validation.md).

## Checklist for an external researcher and independent reviewer

- [ ] Record accession/laboratory provenance, consent for sharing where relevant,
  native file SHA-256, generating CIRI3 executable/version and reference assembly.
  Retain untouched files and run records; do not normalize them before testing.
- [ ] Establish that each native file represents exactly one biological cell.
  Review the original run organization: pooled files can share the same header.
  Provide unique `cell_id`, original `sample_id`, path, detector version and
  `result_scope=single-sample` in the [documented manifest](native_ciri3_import.md).
- [ ] Provide the full intended cell cohort in metadata order. Record completed,
  failed, skipped, unprocessed or unknown plus a status-evidence reference/hash.
  Completion needs run evidence bound to the cell and result, not file existence.
  Explicitly completed zero cells need their actual header-only native result.
  If unavailable, retain unknown rather than manufacturing an empty file.
- [ ] Supply actual genome build, `1-based-closed`, and `CIRI3_junction_reads`.
  Review that counts represent detector BSJ support, not UMIs, Score or FSJ.
  Known unsupported versions/layouts require separate review; never relabel them
  unknown merely to bypass validation. No silent shift, liftover or aggregation.
- [ ] Obtain matching RNA AnnData and its checksum. Verify original cell namespace,
  native ordering, genome build and expected shared/circ-only/RNA-only identities.
  Do not rename unrelated RNA cells to create overlap. Missing assembly metadata
  must be documented as unverified.
- [ ] Build/install the exact reviewed wheel with `[mudata]` in a new venv outside
  the checkout. Record Python/dependency versions, OS, `pip check`, import path
  in `site-packages`, CLI help, wheel checksum and installation wall time. Example
  files are in the source archive, not installed wheel commands.
- [ ] Run the command below in a new directory. Record native bytes, number of
  files/cells/candidates/nonzeros, wall time and peak RSS (GNU time on Linux).
  Record hardware and cache state; no detector execution is needed.
- [ ] Independently parse original tables (not circyto's adapter) and compare every
  candidate ID, chromosome, start/end, strand and integer support against CSR
  entries. Check metadata cell ordering and candidate first-appearance ordering.
  Check gene IDs/types and provenance without inventing symbols/GTF evidence.
- [ ] Load with standard `anndata.read_h5ad` and `mudata.read_h5mu`. Confirm sparse
  int64 circ values, modality axes, exact overlap, RNA counts/layers/annotations,
  failed-cell partial counts, and unknown status for RNA-only observations.
  Memory checks must not densify the full cell-by-candidate matrix.
- [ ] Compare copied native bytes and source hashes; retain the RNA source
  separately (hashed, not copied into evidence). Verify completion numerator and
  denominator from explicit statuses. No callability/true-negative/absence claim.
- [ ] Open HTML offline and confirm exact values, unavailable denominators and
  original/derived evidence separation. Repeat into a second new directory;
  compare JSON/TSV/MTX/HTML bytes and loaded HDF5 values, not HDF5 byte identity.
- [ ] Archive command, environment, checksums, comparison results, timing and
  reviewer sign-off with access instructions. Any mismatch keeps the gate open.

Run after preparing authentic inputs; these placeholders are not supplied data:

```bash
/usr/bin/time -v -o import.time.txt circyto import-circ \
  --ciri3-manifest native.tsv --cells cells.tsv --rna matched_rna.h5ad \
  --detector CIRI3 --genome-build GRCh38 \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads --outdir imported
```

Replace GRCh38 with the actual assembly. Importing requires Python >=3.10 plus
declared dependencies and optional MuData; Java, BWA, STAR and CIRI3 are not needed.
GNU time is only a measurement tool. For an immediately executable demonstration,
use the unchanged [official native example](native_ciri3_import.md#executable-example),
keeping its explicitly narrower evidence claim.
