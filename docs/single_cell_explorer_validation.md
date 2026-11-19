# Single-cell explorer validation

Validation date: 19 November 2026. Development branch:
`feature/single-cell-explorer`, based on the latest local checkpoint `8f729b3`
(`analysis/personal-rna-20261109`), in a separate `circyto-explorer` worktree.
The original worktree remained on `fix/external-detector-distribution`. No
applicable `AGENTS.md` was present. GitHub main was not used as the implementation
baseline. This branch is unmerged and has not been published.

## Artifacts

- [Working KYSE180 HTML](examples/explorer/kyse180/report.html),
  [selected data](examples/explorer/kyse180/explorer.json),
  [exact replay checks](examples/explorer/kyse180/validation.json), and
  [independently checked group summaries](examples/explorer/kyse180/verified_summary.tsv).
- [Synthetic HTML](examples/explorer/synthetic/report.html), including shuffled
  axes, partial overlap and distinct processing states.
- Screenshots of the actual interface:
  [paired plots](examples/explorer/validation/kyse180.png),
  [candidate detail](examples/explorer/validation/candidate.png),
  [narrow layout](examples/explorer/validation/narrow.png), and
  [synthetic demonstration](examples/explorer/validation/synthetic.png).
- [Guide and CLI examples](single_cell_explorer.md).

The KYSE180 replay uses the preserved `rna_analysis.h5ad` and
`circ_annotated.h5ad`, explicitly selecting the 181 SINGLE-layout cells: 92 at
12 Gy and 89 at 30 Gy. It reuses `X_umap_single_primary` and the accepted seven
curated junction identities. Six genes are embedded from 63,241 available
features; seven junctions are embedded from 3,513. All selected values and
coordinates were compared against native objects by exact cell ID. Both
completed-zero circRNA cells remain present. SHA-256 checks of objects, source
QC, candidate summaries and the curation TSV are unchanged.

The newer personal-RNA H5AD was recorded under
`/tmp/circyto-personal-rna-20261109/results/personal_RNA_181_annotated.h5ad`,
which is no longer present. The replay therefore verifies the preserved
accepted embedding; it does not claim to replay that later H5AD or recreate its
30-PC workflow. No alignment, gene counting, circRNA calling, normalization,
clustering or UMAP fitting was rerun. Expression units are raw gene-associated
single reads, not UMIs; support units retain the native CIRI3 declaration.

## Automated checks

The package suite passed with the existing Python 3.10.20 environment, AnnData
0.11.4, MuData 0.3.10, NumPy 2.2.6, pandas 2.3.3 and SciPy 1.15.3. The exact final
test counts and acceptance receipt are recorded in
[automated_checks.json](examples/explorer/validation/automated_checks.json).

```bash
CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q tests
python .github/scripts/check_acceptance_junit.py TEST_RESULTS.xml
```

Explorer tests cover shuffled cell axes; partial modality overlap; MuData and
external RNA precedence; explicit biological identity columns; missing RNA and
UMAP; duplicate gene symbols and stable IDs; exact BSJs and catalogue aliases;
explicit cohort lists; completed-zero, missing, failed, unprocessed and unknown
states; conflicting statuses; missing measurements; safe integer precision;
sparse CSR/CSC and backed files; a 100,000-gene matrix with a guard against full
densification; feature/value/byte limits; HTML injection and link schemes;
determinism; input preservation; output symlinks and source aliases. The
original 55 QC-report tests continue to pass.

Running `pytest -q .` across the entire inherited checkpoint fails collection
in `scripts/radiation_gse81812/test_recipe.py`: its unqualified `prepare_inputs`
import resolves to `scripts/cmpb_a02/prepare_inputs.py`. Those analysis scripts
and the existing root-level CI invocation are unchanged. The successful gate
above runs all package tests under `tests/`; it is not a claim that root-level
collection has been repaired. External detector integrations remain gated.

## Browser and exports

Headless Chrome 157 opened local `file://` reports with network name resolution
disabled and the report's strict content security policy. The test probe checks
real DOM interactions and the bytes/types of generated download Blobs; it does
not rely on a server or a mocked renderer. Receipts:

- [Synthetic controls and exports](examples/explorer/validation/synthetic_browser.json).
- [KYSE180 controls and exports](examples/explorer/validation/kyse180_browser.json).
- [Narrow layout](examples/explorer/validation/narrow_browser.json).
- [Missing UMAP fallback](examples/explorer/validation/missing_umap_browser.json).

Checks include stable-ID/alias/duplicate-symbol searches, exact paired point
positions, linked cell inspection, group filtering, exposure/cluster/QC
selectors, detection thresholds, denominator agreement, all-cell TSV exports,
SVG labels, PNG image generation, no JavaScript exceptions, and no HTTP/CDN
resource requests. The narrow browser has a 500-pixel effective viewport
(Chrome clamps the requested 390-pixel window). Screenshots were visually
reviewed for layout, labels and preserved point geometry. Browser timings in
receipts use virtual time and are **not wall-time performance benchmarks**;
Chrome's heap estimate is coarse.

## Resources and packaging

Measured in separate Python processes on this workstation. Peak RSS includes
imports, input readers and validation, not just serialized report arrays.

| Case | Embedded data | Generation | Peak RSS | HTML |
| --- | --- | ---: | ---: | ---: |
| KYSE180 | 181 cells × (6 genes + 7 junctions) | about 6.2 s | about 288 MiB | about 104 KiB |
| Large sparse axes | 1,000 cells; 200,000 genes and 20,000 junctions available; 8 + 8 selected | about 3.2 s | about 267 MiB | about 372 KiB |

Exact measurements are in the [KYSE180 receipt](examples/explorer/kyse180/validation.json)
and [sparse benchmark](examples/explorer/validation/scaling.json). These are
specific cases, not universal performance limits. Feature annotations and
AnnData layers can consume memory beyond the serialized feature budget.

Local wheel/sdist builds and an installed-wheel replay are documented in
[packaging.json](examples/explorer/validation/packaging.json). The installed
check uses a separate virtual environment with the existing environment's
dependencies, outside the repository. It does not claim a new resolver or
clean dependency installation. Assets are Python strings inside the wheel;
there are no external JavaScript/CSS package dependencies. Package version
remains the inherited `0.11.0rc3`; no release was uploaded.

## Preservation and review scope

Only reporting code, tests, examples, documentation and the acceptance inventory
were changed. No files under `manuscript/`, `review/`, existing analysis scripts
or analysis configurations were edited. Input hashes and branch/worktree checks
are retained in [preservation.json](examples/explorer/validation/preservation.json).
The SRRC manuscript/submission state remains untouched. This is a software
exploration interface, not new scientific evidence or a manuscript update.

Review the new reporting commits against `8f729b3`, rather than interpreting
all inherited workstation commits as part of this interface change. The CLI
requires explicit `--explorer`; default workflow reports retain their previous
behavior and schema. No merge, push, publication or submission was performed.
