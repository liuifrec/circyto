# SComatic Result Normalization

`circyto normalize-scomatic-results` is a read-only adapter for external
SComatic result files. It does not install, invoke, or wrap SComatic. Its only
job is to translate already-produced SComatic-like tables into circyto's
`scomatic_candidate_summary.tsv` schema for later DNA/RNA/circRNA integration.

Terminology guardrail:

- normalized SComatic rows are RNA-derived candidate variant signals
- they are not validated somatic mutations without orthogonal evidence

## Command

```bash
circyto normalize-scomatic-results \
  --scomatic-output external_scomatic.tsv \
  --cell-annotations cell_annotations.tsv \
  --provenance-metadata scomatic_run_metadata.json \
  --outdir work/scomatic_normalized
```

`--scomatic-output` can be repeated for multiple external files. The command
writes:

- `scomatic_candidate_summary.tsv`
- `scomatic_candidate_summary.tsv.provenance.json`
- `normalize_scomatic_results_summary.json`

Real SComatic Step 1/Step 2 calling files are accepted directly:

```bash
circyto normalize-scomatic-results \
  --scomatic-output hap1_batch10.HAP1.step2.calling.step2.tsv \
  --outdir work/scomatic_normalized
```

Step 1 files may contain callable non-candidate rows. The normalizer keeps only
rows with SComatic candidate allele/cell-type fields and streams the file rather
than loading the full table into memory.

Validation note: HAP1 batch10 SComatic BaseCellCounter, Step1, and Step2 execution has been validated as a technical smoke, and real Step1/Step2 output normalization is supported. Candidate rows remain RNA-derived candidate variant signals, not validated somatic mutations.

## Validation Status

Validated:

- real SComatic Step1 and Step2 output parsing
- normalized `scomatic_candidate_summary.tsv` schema
- HAP1 batch10 technical smoke output normalization

Partially validated:

- downstream `candidate_snv` MuData export through `merge-scomatic`
- cell-type-level semantics for real Step1/Step2 files
- provenance capture for external SComatic run metadata

Not yet validated:

- clinically validated variants
- validated somatic mutation discovery
- biological interpretation of genome-scale SComatic outputs
- robust calls from homogeneous datasets without meaningful cell-type/group contrasts

## Lessons Learned from HAP1 and IMR90

- HAP1 batch10 confirmed that real Step1/Step2 files can be normalized without treating rows as validated mutations.
- The HAP1 single-cell-type smoke produced no Step2 candidate calls, which reinforces the need to document grouping and cell-type limitations.
- IMR90 reinforces that processed DNA state, such as CNV, is the primary scRR DNA evidence; SComatic outputs remain RNA-derived candidate variant signals.

## Output Schema

`scomatic_candidate_summary.tsv` always uses these columns:

- `variant_id`
- `cell_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gene`
- `filter_status`
- `candidate_variant_class`
- `read_support`
- `vaf`
- `caller`

`variant_id` is preserved when present. If it is absent, circyto generates a
deterministic ID from `chrom`, `pos`, `ref`, `alt`, and `cell_id`.

`candidate_variant_class` defaults to
`RNA-derived candidate variant signal` when absent. `caller` defaults to
`SComatic` when absent.

## Accepted Minimal Input Fields

The external table must include columns that can be mapped to these fields:

- `cell_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gene`
- `filter_status`
- `read_support`
- `vaf`

Recognized aliases include SComatic-like names such as `CB`, `CHROM`, `POS`,
`REF`, `ALT`, `Gene`, `SComatic_filter`, `Read_count_ALT`, and `VAF`.

If required fields cannot be resolved, the command fails before writing a
candidate summary and reports the missing fields plus the available input
columns.

## Real SComatic Step 1/Step 2 Files

SComatic `BaseCellCalling.step1.py` and `BaseCellCalling.step2.py` outputs use a
VCF-like TSV structure:

- optional `##INFO` metadata lines
- a `#CHROM` header line
- calling fields such as `ALT`, `FILTER`, `Cell_types`, `Bc`, and `VAF`
- one or more trailing cell-type count columns

`circyto` detects this shape and normalizes one row per called
cell-type/allele observation. Multi-allelic values encoded with `|` are
flattened into separate normalized rows. For these real Step 1/Step 2 files,
`cell_id` is populated from SComatic `Cell_types`, because these files report
cell-type-level calls rather than unique-cell genotypes. The summary JSON records
this `cell_id` semantic mapping, the detected SComatic stage, source columns,
metadata line count, and column mappings.

Use the final Step 2 file for downstream integration when available. Step 1 is
supported for inspection and validation, but it predates RNA-editing, PoN, and
nearby-variant filters added by Step 2.

## Cell Annotation Table

When provided, the cell annotation table must contain:

- `Index`
- `Cell_type`

For generic per-cell candidate tables, normalized candidate cell IDs must be
present in `Index`. For real SComatic Step 1/Step 2 calling files, normalized
candidate labels come from `Cell_types`, so they are validated against
`Cell_type` instead. This catches accidental normalization against the wrong
SComatic annotation table while preserving the cell-type-level semantics of
real Step 1/Step 2 outputs.

## Provenance

`--provenance-metadata` can point to a JSON or text file. The command embeds the
metadata path and content in `normalize_scomatic_results_summary.json`, along
with input paths, column mappings, row counts, output paths, and standard circyto
provenance fields. The same provenance payload is also written to
`scomatic_candidate_summary.tsv.provenance.json`.

The normalized output can then be imported with:

```bash
circyto import-dna-snv-summary \
  --workdir WORKDIR \
  --dna-cell-summary dna_cell_summary.tsv \
  --scomatic-candidate-summary work/scomatic_normalized/scomatic_candidate_summary.tsv
```
