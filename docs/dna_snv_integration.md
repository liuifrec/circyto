# DNA/SNV Integration Scaffold

This page defines the current lightweight scaffold for future DNA/RNA integration in `circyto`.

Terminology guardrail:

- SComatic results in `circyto` are treated as `RNA-derived candidate variant signals`
- they are not called validated somatic mutations without orthogonal evidence

## Input contracts

### `dna_cell_summary.tsv`

Required columns:

- `cell_id`
- `dna_library_id`
- `cnv_burden`
- `replication_score`
- `cell_cycle_phase`
- `dna_variant_count`
- `notes`

### `dna_variant_summary.tsv`

Required columns:

- `variant_id`
- `cell_id`
- `chrom`
- `pos`
- `ref`
- `alt`
- `gene`
- `consequence`
- `evidence_type`
- `caller`
- `filter_status`

### `scomatic_candidate_summary.tsv`

Required columns:

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

## Import command

```bash
circyto import-dna-snv-summary \
  --workdir WORKDIR \
  --dna-cell-summary dna_cell_summary.tsv \
  --dna-variant-summary dna_variant_summary.tsv \
  --scomatic-candidate-summary scomatic_candidate_summary.tsv
```

Writes normalized copies under:

- `dna/dna_cell_summary.tsv`
- `dna/dna_variant_summary.tsv`
- `dna/scomatic_candidate_summary.tsv`
- `dna/dna_snv_import_summary.json`

## Joined summary command

```bash
circyto summarize-dna-rna-circ \
  --workdir WORKDIR \
  --write-summary
```

Writes:

- `qc/dna_rna_circ_cell_summary.tsv`
- `qc/dna_rna_circ_summary.json`

## Joined per-cell fields

- `total_rna_counts`
- `detected_genes`
- `circRNA_count`
- `circRNA_total_support`
- `cnv_burden`
- `replication_score`
- `cell_cycle_phase`
- `dna_variant_count`
- `scomatic_candidate_count`

This layer is integration bookkeeping only. It does not run DNA calling, CNV calling, or SComatic internally.
