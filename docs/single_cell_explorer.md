# Optional single-cell circRNA explorer

The explorer adds gene expression, native circRNA junction support and saved RNA
UMAPs to the existing offline report. It is an unmerged usability extension from
local checkpoint `8f729b3`, outside the frozen manuscript/submission baseline.
It never writes AnnData/MuData, normalizes expression, fits an embedding, or runs
alignment, gene counting or circRNA calling.

Open the [synthetic HTML](examples/explorer/synthetic/report.html) or
[181-cell KYSE180 example](examples/explorer/kyse180/report.html) in a browser.
Download the HTML itself when viewing through a repository web interface. All
data needed for interaction are inside that one file. No CDN, server, fetch,
remote font, telemetry or external data request is used. JavaScript must be
enabled. Supplied evidence URLs open externally **only when clicked**.

## CLI

The existing command remains valid and produces the same script-free QC report:

```bash
circyto report --workdir RESULTS
```

Add an explorer to a workflow report, with a separate report destination:

```bash
circyto report --workdir RESULTS \
  --explorer RESULTS/mudata/circyto_multimodal.h5mu \
  --rna processed_rna.h5ad --output-dir reports/combined \
  --umap-key X_umap --rna-layer log1p_cp10k \
  --expression-unit 'log1p CP10K' \
  --gene MAN1A2 --gene ENSG00000183023 \
  --circ 'chr1:117402186|117420649'
```

Use the matrix/layer and unit label appropriate to your RNA processing. **No
transformation is applied.** With no `--rna-layer`, the stored `X` is used; the
default expression label explicitly says units are undeclared. Native support
uses `support_definition` when present, or an undeclared-unit label. An explicit
`--support-unit` overrides that label. Native circRNA `X` must contain
nonnegative integer support, not normalized expression, scores or FSJ ratios.

A standalone explorer needs no workflow directory:

```bash
circyto report --explorer circ_counts.h5ad --rna processed_rna.h5ad \
  --output-dir reports/explorer --gene MAN1A2 \
  --circ hsa_circ_0000118 --candidate-metadata candidate_annotations.tsv \
  --expression-unit 'gene-associated single reads (not UMIs)'
```

Inputs are a circRNA AnnData `.h5ad` or MuData `.h5mu` with `circ` and/or `rna`
modalities. A separate `--rna` AnnData supplies **both** expression and its
coordinates/annotations, superseding MuData's RNA modality. Scanpy is not
needed to read an existing Scanpy output. MuData input requires the existing
`circyto[mudata]` extra; no other dependency was added.

Standalone output contains `report.html` and `explorer.json`. A combined
workflow report also contains the existing `metrics.json`, with an additional
explorer selection/provenance summary. Without `--output-dir`, a workflow
report still writes to `RESULTS/qc`. Import-circ bundles cannot be overwritten:
use a separate output directory to preserve their distinct QC schema. All
report content is prepared before staged replacement. Replacement of several
files is not a filesystem transaction; regenerate after an interruption.

## Exact identities and scope

Cell joins use unique, nonempty, exact `obs_names` strings. No row-position join,
barcode suffix removal, case folding, ID guessing, duplicate merging or
automatic ID renaming is performed. `--cell-id-key biological_cell_id` selects
an explicit identity column in every input `obs` instead. For MuData that
includes the global observation table. Prepare the biological cell mapping
before reporting if technical libraries use different namespaces. Disjoint
RNA/circRNA cell namespaces are rejected.

By default the modality union is retained. RNA-only and circ-only cells have
missing values in the absent modality. Missing coordinates omit a point, not a
cell from tables or detection denominators. Missing RNA or the requested UMAP
key omits the associated plots; no coordinates are invented. Non-finite
coordinates are treated as unavailable. The selected embedding must have two
columns. Use `--umap-key X_umap_single_primary` or another saved key explicitly;
the reporter does not guess between saved scientific representations.

For a declared cohort, `--cells cells.txt` takes one exact biological cell ID
per line without a header. Duplicates, blank lines and unknown IDs are rejected.
Only this explicit report scope affects the displayed population; source
matrices remain untouched. `--title` provides a visible dataset label.

## Feature selection and search

`--gene` and `--circ` are repeatable exact selection queries. Gene symbols,
stable gene IDs (including versionless Ensembl aliases), BSJ coordinates,
catalogue aliases and recorded host-gene fields are supported. BSJ aliases use
the recorded coordinates, e.g. `chr1:10|20`, `chr1:10-20` and a `:+`/`:-` suffix;
no genome-build or coordinate-convention conversion is performed. Searches in
the browser are case-insensitive substrings over the **embedded subset**.

Duplicate gene symbols retain every matched feature under separate stable IDs;
they are never summed. If duplicate `var_names` are symbols, a unique `gene_id`
column provides identities. Otherwise feature identities must be unique. A
query that exceeds the feature bound fails instead of choosing one duplicate.

The optional annotation TSV has one unique `circ_id` row per candidate:

```text
circ_id	circBase_id	uniform_id	annotation_source	evidence_url
chr1:10|20	hsa_circ_example	circGENE(2,3,4).1	curated fixture	https://example.org/evidence
```

This is an illustrative schema, not a catalogue assignment. Known scalar
annotation fields, source labels, coordinates, strand, build, conventions,
catalogue IDs and DOI/evidence URL fields are retained. Aliases can be separated
by semicolons; commas inside uniform circRNA names remain intact. Conflicting
metadata and rows absent from the circRNA object are rejected. No annotation
service is contacted. Catalogue identity, host assignment, read evidence and
functional confirmation remain separate claims.

Offline bounds are explicit:

| Option | Default | Behavior |
| --- | ---: | --- |
| `--max-genes` | 32 | Maximum embedded genes; zero disables gene values |
| `--max-circs` | 32 | Maximum embedded junctions; zero disables junction values |
| `--max-cells` | 50,000 | Oversized scopes fail; no silent cell sampling |
| `--max-values` | 2,000,000 | Maximum cells × total selected features |
| `--max-bytes` | 32,000,000 | Maximum UTF-8 explorer JSON size |

Explicit queries include only their matches. Without queries, circRNAs use
stable ID order, and genes prioritize recorded host genes of the selected
circRNAs, then stable ID order. This is a bounded preview, not a biological
ranking. Total versus embedded feature counts and the selection rule remain
visible. Rebuild with other queries to explore additional features.

Only selected matrix columns are materialized, one at a time. The complete
matrix is never densified. File loading uses backed read-only access where
supported; AnnData/MuData still read metadata and can load non-X layers into
memory. Bounds control serialized measurements, not all memory used by the
reader or full feature annotation metadata. Read-only source hashing scans
input bytes. Extremely large objects may need an explicitly scoped input.

Observation fields are prefixed by source (`rna:`, `circ:`, `global:`), so
conflicting annotations are not silently reconciled. Up to 48 fields are chosen
automatically, prioritizing clusters, annotations and condition/exposure. Use
repeatable `--obs-column rna:condition` to select up to 64 fields explicitly.
Categorical colour/group selectors accept at most 128 distinct labels; legends
show up to 12 labels with a repeating palette. Numeric QC fields are selectable
with stored-value labels. Metadata and text lengths also have bounds.

## Observation states and detection frequencies

| State for a selected junction | Plot/table behavior | Included in rate denominator |
| --- | --- | --- |
| Completed, support > 0 | Native value; `completed_detected` | Yes |
| Completed, support = 0 | Measured zero; `completed_zero` | Yes |
| Missing modality/value | Missing marker; absent values exported blank | No |
| Failed | Cross; supplied values retained in tables | No |
| Unprocessed/pending/planned/running/skipped | Triangle; supplied values retained | No |
| Unknown/unrecognized/no status | Diamond; supplied values retained | No |

Processing evidence comes from `circ_observation_status` and `detector_status`,
with `processing_status` as fallback. Both specific declarations must establish
completion when both exist; conflicting or missing statuses stay unknown. Any explicit failure or
unprocessed state overrides completion. An existing alignment status must also
establish completion. Recognized completion includes `completed`,
`completed_detected`, `completed_zero`, `success`, `empty`, `skipped_existing`.
Nonzero support alone never establishes completion. RNA status, when recorded,
is shown separately; absent RNA status means stored values are available, not
that a pipeline has been independently verified.

The selected candidate's group frequency is `cells with support ≥ threshold /
completed cells with a finite measurement`. Counts and denominators are
recomputed for the current report/group filter, not copied from cohort-level
feature annotations. Cells without UMAP coordinates remain eligible. A positive
support value below a raised threshold remains a measured positive value; it
does not become `completed_zero`. Completion does not establish locus-level
callability, and zero does not establish biological absence. Frequencies are
descriptive, with no inferred replicates, statistical test or causal claim.

## Export and reproducibility

Each plot exports SVG or PNG with the selected feature, units, linear colour
scale or category legend, and saved embedding key. Both plots use identical
coordinate limits. Hover/click identifies the same cell in both plots; entering
an exact cell ID also supports keyboard inspection.

The summary TSV exports candidate identity, grouping field, threshold, measured
denominator, detected and completed-zero counts, missing/failure state counts,
frequency, native support and its unit. The selected-cell TSV exports all
selected cells, coordinates, gene and circRNA IDs, untransformed values, units,
states and the current RNA plot's annotation/QC values. The on-screen cell table
is limited to 100 rows; TSV export is not. Missing values are blank, not zero.
Potential spreadsheet-formula text is apostrophe-prefixed in TSVs; JSON retains
the original strings. Values or support sums beyond JavaScript's exact integer
range are refused, avoiding silently rounded counts.

`explorer.json` uses `circyto.explorer.v1`. `cells` is the joined axis;
`coordinates`, observation `values`, and each feature's `values`/`states` follow
that exact axis. It also records the chosen layer, units, selection bounds,
source basenames/SHA-256/byte sizes and omitted-panel warnings. No arbitrary
`uns` objects or command lines are exported. Selected identifiers and
annotations are retained; this is not an anonymization tool.

For identical inputs, selection and package code, generation is deterministic.
The HTML embeds escaped JSON, uses DOM text insertion, permits only a hashed
local script and blocks network connections through its content security
policy. File/link guards protect sources and existing import bundles.

## Reproduce examples

```bash
PYTHONPATH=. python examples/single_cell_explorer_demo.py
PYTHONPATH=. python examples/replay_kyse180_explorer.py
PYTHONPATH=. python examples/benchmark_explorer.py
CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q tests
```

The KYSE180 script needs the preserved local data directory or `--source PATH`;
it only reads existing objects, writes reports and verifies exact values and
source hashes. It uses the saved 181-cell `X_umap_single_primary`, not the later
personal-RNA H5AD whose historical temporary path is unavailable. Its selected
genes use raw gene-associated single reads, explicitly not UMIs. The seven
candidate aliases reuse the existing curation TSV without rewriting it.

Browser checks use standard headless Chrome without npm packages:

```bash
python tests/explorer_headless.py prepare docs/examples/explorer/synthetic /tmp/explorer-test.html
google-chrome --headless --disable-gpu --disable-background-networking \
  --virtual-time-budget=30000 --run-all-compositor-stages-before-draw \
  --dump-dom file:///tmp/explorer-test.html > /tmp/explorer-dom.html
python tests/explorer_headless.py verify /tmp/explorer-dom.html /tmp/browser-results.json
```

The validation container additionally needs `--no-sandbox`; use the browser's
normal sandbox on a workstation. The probe checks generated Blob exports and
then replaces the document with its JSON receipt. It does not alter the actual
report. Screenshots are taken from the actual report:

```bash
google-chrome --headless --disable-gpu --disable-background-networking \
  --window-size=1440,1800 --hide-scrollbars --screenshot=/tmp/explorer.png \
  "file://$PWD/docs/examples/explorer/kyse180/report.html"
```

See [validation evidence](single_cell_explorer_validation.md) for measured
resources, preservation checks and the scope of test execution.
