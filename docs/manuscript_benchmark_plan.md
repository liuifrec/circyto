# Manuscript Benchmark Plan

This page defines a lightweight manuscript-facing benchmark summary layer for completed `circyto` workflows.

Scope:

- workflow and output benchmarking only
- no rerunning workflows
- no biological claims

## Benchmark command

```bash
circyto summarize-benchmark \
  --workdirs WORKDIR1 \
  --workdirs WORKDIR2 \
  --output benchmark_summary.tsv \
  --json benchmark_summary.json
```

## Output columns

- `dataset_name`
- `workflow_type`
- `protocol`
- `read_layout`
- `n_cells`
- `n_rna_features`
- `n_circ_features`
- `median_rna_counts`
- `median_detected_genes`
- `median_circRNA_count`
- `median_circRNA_total_support`
- `h5mu_exists`
- `h5mu_size_bytes`
- `workdir_size_bytes`
- `cleanup_status`
- `workflow_succeeded`

## Intended manuscript use

This command is designed to support:

- Supplementary Tables S1-S3
- workflow reproducibility summaries
- disk and cleanup reporting across benchmark datasets

It is not intended to imply biological ranking or mechanistic interpretation by itself.
It does not replace the checksum-matched regeneration required by
`manuscript/application_note_evidence.md`.
