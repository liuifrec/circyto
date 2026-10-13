# Import existing cell-resolved circRNA support

`circyto import-circ` converts an explicitly defined support table or an existing
circyto collector matrix into sparse AnnData, optional RNA+circ MuData, and an
offline HTML report. It runs no aligner or detector. Install `circyto[mudata]`
when using `--rna`; AnnData-only import uses the core installation.

For native 13-column single-sample CIRI3 results, use the additional
[`--ciri3-manifest` pathway](native_ciri3_import.md). It requires explicit mapping
and source semantics and reuses this import/export contract.

## Short executable example

With the installed wheel on PATH, run the distributed
[`examples/import_circ_demo.py`](../examples/import_circ_demo.py):

```bash
python examples/import_circ_demo.py /tmp/circ-import-demo
```

The destination must be new. This explicitly **synthetic** example creates six
cells, a supplied support table and a small RNA object. It checks count values
and order. Open `/tmp/circ-import-demo/result/qc/report.html`. No external
executable, reference download, network service or raw sequencing analysis is
needed. The script is in the source distribution/repository, not installed as a
wheel command; copy it anywhere and run it against the installed package.

For your own inputs:

```bash
circyto import-circ \
  --calls candidate_support.tsv --cells cells.tsv \
  --genome-build GRCh38 --coordinate-convention 1-based-closed \
  --support-definition BSJ_reads --detector CIRI3 \
  --rna matched_rna.h5ad --outdir imported
```

Omit `--rna` for AnnData only. Detector and support names are declarations,
not requests to execute or certify a detector. Supply exactly one assay/count
unit per import. The importer does not pool detectors, convert coordinates,
perform liftover, rename chromosomes, infer assemblies, or interpret normalized
expression as integer support.

## Input contract

All tables are tab-separated with unique headers. IDs are literal strings:
`001` and `NA` are retained. Empty IDs, duplicates and surrounding whitespace
are rejected. Extra cell metadata columns are retained as strings.

`cells.tsv`: required `cell_id`; optional `detector_status` and other metadata.
Row order defines the circ AnnData cell axis. Derived names `circ_observation_status`,
`circRNA_count`, `circRNA_total_support`, `membership`, `circ_present`, and
`rna_present` are reserved; rename historical versions (e.g. `source_membership`)
to retain them as input evidence. Include **all intended cells**,
including cells without candidates. Missing `detector_status` means `unknown`;
blank or unrecognized values are errors. Map legacy statuses explicitly using
the actual run evidence (e.g. confirmed `success` or `empty` to `completed`).
The importer does not infer successful processing from counts or file existence.

`candidate_support.tsv`: required columns:

```text
cell_id  circ_id  chrom  start  end  strand  support
```

Each row represents one cell/candidate pair; pairs must be unique. `support`
is a nonnegative integer (int64); zero is permitted. Duplicates are errors,
not silently summed. Every call cell must occur in metadata. `circ_id` is an
opaque stable identifier; start/end/strand are the authoritative coordinates.
Candidates follow first appearance order. One ID must have consistent coordinates
and host-gene annotations across cells. Strand is `+`, `-`, or `.`. Coordinate
validation follows the declared convention; no coordinate arithmetic is applied.
Integer support is never cast through floating point or narrowed to int32.
The total support must also fit int64 so summaries cannot overflow.

Optional annotation columns `host_gene`, `host_gene_source`,
`host_gene_from_gtf`, `host_gene_from_circatlas`, and
`host_gene_from_circatlas_id` are retained in `var`. Supplied values are preserved;
missing schema columns use the existing host-gene normalizer. A supplied gene
name alone is **not** labeled GTF evidence. Other row-specific evidence columns
remain in the exact source copy under `evidence/calls.tsv`.

Optional `genome_build`, `coordinate_convention`, `support_definition` and
`detector` columns in calls, features or cell metadata must exactly match the
command declarations in every row. Mixed or missing declarations in such
columns are rejected. When these columns are absent, the user is responsible
for verifying the declared source semantics; coordinates alone cannot identify
an assembly or a count unit. Use explicit labels such as
`historical_chr21_assembly_unverified` for uncertain historical assemblies;
this records uncertainty and does not make them suitable for cross-study joins.

## Reuse existing collectors

Existing `collect-matrix` commands already read saved outputs without detection:

```bash
circyto collect-matrix --detector ciri2 --indir saved_normalized_calls \
  --outdir collected --min-count-per-cell 0
circyto import-circ --matrix-dir collected --cells cells.tsv \
  --genome-build GRCh38 --coordinate-convention 1-based-closed \
  --support-definition BSJ_reads --detector CIRI2 --outdir imported
```

The matrix-directory contract is strictly **circ × cell**, including square
matrices. Required files are `circ_counts.mtx` (sparse integer Matrix Market),
`circ_index.txt`, `cell_index.txt`, and `circ_feature_table.tsv` with
`circ_id,chrom,start,end,strand`. Features must exactly cover the candidate index;
feature order is aligned to that index. Duplicate entries or index IDs are
rejected. All extra feature columns survive. Cell columns are reordered to the
metadata cohort. A metadata cell missing from the source matrix can be retained
as failed/skipped/unprocessed/unknown; marking it completed requires an explicit
source zero column. This prevents a filtered-out cell becoming a completed zero.

Collector-specific limits matter:

- CIRI-full/CIRI2/CIRI3 collectors read circyto's normalized per-cell TSV
  schema (`circ_id,chr,start,end,strand,support`), not every native vendor format.
- find-circ3 uses its existing flag filters and `uniq_bridges` support. It may
  drop cells without retained calls even with a zero threshold and does not
  emit the required feature table. Supply annotations from the original calls;
  do not infer zero-cell completion from a missing column.
- CIRCexplorer2 has legacy collection paths with different event/count handling;
  the standalone matrix parser treats rows as events. Inspect the particular
  collector output and declare the actual support unit. This import does not
  assert BSJ-read equivalence or repair collector behavior.
- `merge-detectors` aggregates over cells; its `circ_by_detector.tsv` is not a
  cell-resolved input to this command. `compare-detectors` compares those
  aggregates and does not recover the lost cell axis.

For unsupported native output layouts, create the documented long table using
a reviewed conversion and retain the original files. Do not claim that a
successful long-table import validates the original detector parser.

## Observation semantics and denominators

| Input detector_status | Derived circ_observation_status |
| --- | --- |
| completed, positive row | completed_detected |
| completed, all-zero row | completed_zero |
| failed | failed |
| skipped | skipped |
| unprocessed | unprocessed |
| unknown / status column absent | unknown |

The matrix is unchanged by status. Partial support in failed cells is retained.
`circRNA_count` counts strictly positive candidates;
`circRNA_total_support` sums supplied support. New metadata fields are additive;
legacy exporters and historical objects are unchanged. Older objects without
status cannot be interpreted as completed just because they contain zero rows.

Detection fraction = completed cells with positive support / all explicitly
completed metadata cells. Candidate `prevalence_completed` uses the same
completed-cell denominator and only positive calls in that subset. A denominator
of zero produces JSON `null` / AnnData `NaN`, not a zero rate. Raw totals and
`n_cells_detected` include all supplied values. Rates reflect the supplied table
and its filtering history, not unfiltered detector sensitivity. Cell-level
completion does not establish locus-level callability, true negatives or
biological absence.

## Outputs, RNA alignment and provenance

- `anndata/circ_counts.h5ad`: CSR int64 counts, metadata `obs`, annotated `var`,
  `uns['circyto']` declarations, source SHA-256 hashes and exact metrics JSON.
- `matrix/`: standard circ-by-cell Matrix Market, indexes and feature table.
- `qc/`: full `cell_qc.tsv`, `circ_qc.tsv`, `metrics.json`, `report.html`.
- `evidence/`: byte-for-byte input tables/matrix/indexes, including extra evidence.
- `mudata/circyto_multimodal.h5mu` with `--rna`: native `circ` and `rna`
  modalities and global union metadata. **Each modality retains its own cell
  order and coverage**; circ uses metadata order, RNA uses its original order.
  No RNA layers, annotations, embeddings or matrix values are reconstructed.
  RNA-only global rows have `unknown` circ status and missing circ metrics, not
  completed-zero status. `circ_present`, `rna_present`, and `membership` record
  modality membership. Match by ID, not row position. This uses MuData's native
  modality axes; unlike the older completed-workflow exporter it does not pad
  each modality to the union.

RNA IDs must be unique and have nonzero overlap. The report gives overlap and
unmatched counts. Any recorded RNA `genome_build` in `obs`, `var`, `uns`, or
`uns['circyto']` must match; otherwise import fails. If absent, the RNA build
remains explicitly unverified. RNA is hashed and preserved in MuData rather
than copied again into `evidence/`. Object-level provenance refers to the
supplied input files; archived upstream evidence remains the caller's
responsibility.

```python
import anndata, mudata
circ = anndata.read_h5ad('imported/anndata/circ_counts.h5ad')
paired = mudata.read_h5mu('imported/mudata/circyto_multimodal.h5mu')
```

The import-specific HTML reuses existing report styling, with a distinct import
schema and explicit denominator. It does not invent alignment or workflow
completion summaries to satisfy the older `circyto report` command. Reimport to
a new directory to regenerate it. JSON/TSV/HTML and scientific values are
reproducible for identical inputs; HDF5 byte identity is not promised. Existing
output directories are never overwritten. Reports include supplied cell metadata;
review their contents before sharing.
