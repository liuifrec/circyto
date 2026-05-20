# RNA Circ Integration

`circyto summarize-rna-circ` provides a lightweight integration summary across existing RNA profiling outputs and the circRNA matrix of a completed workflow directory.

## Command shape

```bash
circyto summarize-rna-circ \
  --workdir /path/to/completed_workflow
```

JSON mode:

```bash
circyto summarize-rna-circ \
  --workdir /path/to/completed_workflow \
  --json
```

Write QC summary files:

```bash
circyto summarize-rna-circ \
  --workdir /path/to/completed_workflow \
  --write-summary
```

## Inputs

The command reads:

- `rna/gene_counts.tsv`
- `qc/rna_qc.tsv` when present
- `rna/gene_feature_table.tsv` when `rna_qc.tsv` is missing and a fallback must be derived
- `matrix/circ_counts.mtx`
- `matrix/cell_index.txt`
- `matrix/circ_index.txt`
- `qc/circ_qc.tsv` is optional and not required

No alignments are required.

## What it reports

- `n_rna_cells`
- `n_circ_cells`
- `n_shared_cells`
- `n_rna_only_cells`
- `n_circ_only_cells`
- RNA-only cell IDs
- circ-only cell IDs
- per-cell:
  - `total_rna_counts`
  - `detected_genes`
  - `mitochondrial_fraction`
  - `ribosomal_fraction`
  - `circRNA_count`
  - `circRNA_total_support`
- a simple RNA-total-count versus circRNA-count relationship summary across shared cells

## Output files

With `--write-summary`, the command writes:

- `qc/rna_circ_cell_summary.tsv`
- `qc/rna_circ_summary.json`

Without `--write-summary`, the command is read-only.

## Downstream MuData Export

After RNA and circ summaries are in place, export a multimodal `.h5mu` bundle with:

```bash
circyto export-mudata \
  --workdir /path/to/completed_workflow
```

This preserves RNA-only cells in the shared observation table and zero-fills the circ modality for cells without circ matrix rows.

## RNA-only Cells

RNA-only cells are preserved in the joined summary. For cells such as `DIYHEK_192` that appear in RNA profiling outputs but not in the circRNA matrix:

- `membership = rna_only`
- `circRNA_count = 0`
- `circRNA_total_support = 0`

This is expected and does not require surviving alignment files.
