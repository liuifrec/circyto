# MuData Export

`circyto export-mudata` creates a scverse-style multimodal `.h5mu` bundle from an already completed RNA+circ workflow directory.

## Command shape

```bash
circyto export-mudata \
  --workdir /path/to/completed_workflow
```

Explicit output path:

```bash
circyto export-mudata \
  --workdir /path/to/completed_workflow \
  --output /path/to/full_length.h5mu
```

Overwrite an existing output:

```bash
circyto export-mudata \
  --workdir /path/to/completed_workflow \
  --overwrite
```

## Default output

When `--output` is omitted, the command writes:

- `WORKDIR/mudata/full_length.h5mu`

## Inputs

The exporter uses existing completed-workflow outputs only:

- `rna/gene_counts.tsv`
- `rna/gene_feature_table.tsv`
- `qc/rna_qc.tsv`
- `matrix/circ_counts.mtx`
- `matrix/cell_index.txt`
- `matrix/circ_index.txt`
- `matrix/circ_feature_table.tsv` when present
- `qc/rna_circ_cell_summary.tsv` when present
- `workflow_summary.json` when present
- `rna/rna_import_summary.json` when present
- `qc/rna_circ_summary.json` when present

No alignments are required.

## MuData structure

- `mdata["rna"]`
  - `X` = gene counts
  - `obs` = shared cell metadata aligned to the union of RNA and circ cells
  - `var` = `gene_feature_table.tsv`
- `mdata["circ"]`
  - `X` = circRNA counts
  - `obs` = same shared cell metadata
  - `var` = `circ_feature_table.tsv` when present, otherwise minimal `circ_id`
- `mdata.obs`
  - union of RNA and circ cells
  - includes `membership`, `total_rna_counts`, `detected_genes`,
    `mitochondrial_fraction`, `ribosomal_fraction`, `circRNA_count`,
    and `circRNA_total_support`
- `mdata.uns["circyto"]`
  - `workflow_summary`
  - `rna_import_summary`
  - `rna_circ_summary`
  - `command_name`
  - `circyto_version`
  - `source_workdir`

## RNA-only Cells

RNA-only cells such as `DIYHEK_192` are included in:

- `mdata.obs`
- `mdata["rna"]`

For the circ modality, these cells are retained with zero-filled circ rows. This preserves a stable shared cell axis across modalities without discarding RNA-only cells.

## Dependency behavior

`mudata` is optional. If it is not installed, the command fails clearly with:

```text
mudata is not installed; install circyto[mudata] or pip install mudata
```

## Read-only Downstream Inspection

After export, inspect the object structure and QC summaries with:

```bash
circyto inspect-mudata --input /path/to/full_length.h5mu
circyto summarize-mudata-qc --input /path/to/full_length.h5mu
```
