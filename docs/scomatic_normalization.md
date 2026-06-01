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
- `normalize_scomatic_results_summary.json`

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

## Cell Annotation Table

When provided, the cell annotation table must contain:

- `Index`
- `Cell_type`

The normalized candidate cell IDs must be present in `Index`. This catches
accidental normalization against the wrong SComatic annotation table.

## Provenance

`--provenance-metadata` can point to a JSON or text file. The command embeds the
metadata path and content in `normalize_scomatic_results_summary.json`, along
with input paths, column mappings, row counts, output paths, and standard circyto
provenance fields.

The normalized output can then be imported with:

```bash
circyto import-dna-snv-summary \
  --workdir WORKDIR \
  --dna-cell-summary dna_cell_summary.tsv \
  --scomatic-candidate-summary work/scomatic_normalized/scomatic_candidate_summary.tsv
```
