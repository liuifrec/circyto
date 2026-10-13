# Native CIRI3 results to AnnData and MuData

`circyto import-circ --ciri3-manifest` accepts **existing native 13-column
single-sample CIRI3 result tables**. No CIRI3, Java, aligner, reference download,
or raw sequencing analysis is needed. Install the feature wheel with
`python -m pip install /path/to/circyto.whl`; use
`python -m pip install '/path/to/circyto.whl[mudata]'` for optional RNA integration.
Python >=3.10 and the existing declared dependencies are sufficient. The package
version remains 0.10.0; this feature is not part of the frozen release baseline.

## Executable example

The source archive/repository includes an unchanged official native output,
its original log, SHA-256 provenance and upstream license under
[`examples/native_ciri3`](../examples/native_ciri3). Copy that directory anywhere;
the CLI can run from outside the checkout against an installed wheel:

```bash
circyto import-circ \
  --ciri3-manifest examples/native_ciri3/manifest.tsv \
  --cells examples/native_ciri3/cells.tsv \
  --detector CIRI3 --genome-build upstream_demo_chr1_unverified \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads \
  --outdir /tmp/native-ciri3-import
```

Choose a new output directory. This imports one **demonstration sample**, six
candidates and 50 BSJ support. It is not presented as a biological single cell.
Open `/tmp/native-ciri3-import/qc/report.html` or load the sparse object:

```python
import anndata
circ = anndata.read_h5ad('/tmp/native-ciri3-import/anndata/circ_counts.h5ad')
```

The producing detector version and demo assembly are not established, so they
are explicitly labeled unknown/unverified. The artifact is an authentic
pre-existing upstream result, not a synthetic table constructed by circyto.
Its underlying read provenance is insufficient for biological benchmarking.

For your own files, provide the actual genome build and these two TSVs.
`native.tsv` has exactly these columns in order:

```text
cell_id  sample_id  path                  detector_version  result_scope
cell_A   sample_A   results/A/result.txt  3.0.1             single-sample
cell_B   sample_B   results/B/result.txt  3.0.1             single-sample
```

Use `unknown` only when the producing version is unavailable. **Different
versions must be imported separately**, including known versus unknown. Actual
files must use tabs. Relative paths resolve against the manifest directory. Each cell, sample
and physical source file can be mapped only once. IDs are literal strings;
filenames are never used to infer cell identities. `sample_id` is the original
sample identity supplied by the user; `cell_id` is the desired shared observation
namespace. A bulk sample must not be misrepresented as an individual cell.

`cells.tsv` defines the full cohort and ordering:

```text
cell_id  detector_status  batch
cell_A   unknown          batch1
cell_B   failed           batch1
cell_C   unprocessed      batch1
```

Only `cell_id` is required; absent status means `unknown`. Supported explicit
statuses are completed, failed, skipped, unprocessed and unknown. A completed
cell must have a mapped result, including a header-only result for a completed
zero. Extra metadata cells without native files stay failed/skipped/unprocessed/
unknown. Partial counts in a failed cell are retained, not normalized away.

To use a real matched RNA object, add `--rna /path/to/matched_rna.h5ad`. Matching
uses exact cell IDs and retains each modality's native axis, counts, metadata,
layers and embeddings. Zero overlap or conflicting recorded genome builds fail.
See [the general import contract](import_circ.md) for RNA-only cells and missing
build metadata. No biologically matched RNA accompanies the official native demo;
our optional demo MuData check is explicitly a technical ID-alignment test.

## Supported native schema and versions

Adapter schema: `ciri3.native-result-13.v1`. Header order and spelling are exact:

```text
circRNA_ID chr circRNA_start circRNA_end #junction_reads SM_MS_SMS #non_junction_reads junction_reads_ratio circRNA_type gene_id strand junction_reads_ID Score
```

The format is grounded in the official
[v3.0.1 tree, commit 73108c4478a7e691cb93360921b65203e0701516](https://github.com/gyjames/CIRI3/tree/73108c4478a7e691cb93360921b65203e0701516)
and its [output documentation](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/README.md#output).
The native result was first committed in
[27031610c234e7459bc44fb9bd9bded2050599b5](https://github.com/gyjames/CIRI3/commit/27031610c234e7459bc44fb9bd9bded2050599b5)
on 2023-08-16 and remains byte-identical in the tagged tree. **Presence in a
release is not evidence that that release generated it.** The file and log do
not record a producing version; the example manifest therefore says `unknown`.

Version declarations accepted by this narrow interface are `3.0.1` and `unknown`
(the latter only when producer metadata is unavailable). This is a schema
compatibility boundary, not certification of every CIRI3 mode or JAR. Other
versions require a separately reviewed long-table conversion; do not relabel a
known unsupported version as unknown.

- `circRNA_ID` must equal `chr:start|end` using the supplied coordinate columns.
  Chromosome labels are retained exactly. Strand must be `+` or `-` and remains
  a separate field; it is not appended to the native identifier.
- **Coordinates remain 1-based closed.** The official README calls the start
  1-based. Its [user-coordinate reader](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/src/com/zx/findcircrna/GetUserCircRNA.java)
  accesses flanking bases with Java `substring(start-3,start-1)` and
  `substring(end,end+2)`, consistent with an inclusive end. The
  [annotation writer](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/src/com/zx/hg38/Annotation.java)
  writes the same junction coordinates and compares them directly with GTF exon
  boundaries. The adapter does not perform shifts, liftover or chromosome mapping.
- `--support-definition CIRI3_junction_reads` means the integer
  `#junction_reads` field. The [summary source](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/src/com/zx/findcircrna/Summary.java)
  counts supporting read identifiers in sets. This is detector BSJ support,
  not UMIs, gene-expression counts, `Score`, FSJ counts or a junction ratio.
  The README describes `Score` separately as high-quality first-scan support.
  Raw fields are retained; the adapter neither recomputes counts from read IDs
  nor applies a detector-quality threshold.
- Duplicate candidate rows, malformed rows, invalid integer support, inconsistent
  candidate coordinates/strand or conflicting annotations across mapped samples
  fail with errors. No duplicate aggregation, strand coercion or count scaling.
- `gene_id` and `circRNA_type` are retained as `var['ciri3_gene_id']` and
  `var['ciri3_circRNA_type']`. Nonempty/non-NA native gene IDs also populate
  `host_gene_id` with `host_gene_id_source='CIRI3_gene_id'`. They are not converted
  to gene symbols or treated as evidence that circyto performed GTF annotation.
  `host_gene` stays empty unless supplied through another established workflow.

## Scope and observation limitations

CIRI3 produces single-sample files, pooled information tables, BSJ/FSJ matrices,
RNase-R outputs and user-specified-candidate outputs. The pooled information
file can have **the same header** as a single-sample result, and cannot be
identified from the header alone. Thus `result_scope=single-sample` is a required
user assertion backed by the original run organization. This adapter refuses
other declared scopes and does not accept multi-sample matrices, 12-column
variants or normalized circyto collector TSVs as native results. It cannot
detect a falsely labeled pooled file; inspect the upstream run before mapping.

The official log reports first/second scans and collation completion, but lacks
a detector version, stable sample identifier, artifact checksum or captured exit
code. It provides supporting execution evidence, not an independently bound
per-cell status certificate. The adapter does not promote status based on file
existence, log text, positive counts or a header-only table. A `completed` label
requires the user's independent run evidence. Unknown remains unknown.
Detection fractions exclude failed/skipped/unprocessed/unknown cells and never
establish locus-level callability, true negatives or biological absence.

## Evidence and outputs

The adapter uses the existing CIRI3 normalization function **after** strict
native validation, then delegates to the existing sparse long-table builder
and the same AnnData/MuData, provenance and HTML export pipeline. Existing
`--calls` and `--matrix-dir` behavior is unchanged; the three modes are mutually
exclusive.

Original bytes are copied to `evidence/ciri3_native_000001.tsv`, etc., with the
original manifest in `evidence/ciri3_manifest.tsv`. Converted outputs are under
`matrix/` and `anndata/`; no native file is rewritten. SHA-256 values are checked
before and after conversion. Native read IDs and all ancillary fields remain
in raw evidence. Logs are not required or imported; retain them with the
upstream run. The example's original log is distributed separately for auditing.

`obs` gains `ciri3_sample_id`, `ciri3_detector_version`, and `ciri3_source_file`.
`metrics.json['native_ciri3']` maps cell/sample/version to source evidence;
`uns['circyto']['native_ciri3_json']` preserves the same mapping. Count units,
genome and coordinate declarations remain in the existing schema. Review
read IDs, sample labels and manifest paths before sharing evidence bundles.

[Validation, measured performance and remaining gaps](native_ciri3_validation.md)
are reported separately. No new calling algorithm, statistical model, additional
detector parser or version bump is part of this milestone.
