# Scanpy Downstream

`circyto` provides lightweight exploratory Scanpy helpers for completed `.h5mu` outputs.

## Dependency

Install the optional extra:

```bash
pip install -e .[scanpy]
```

If `scanpy` is missing, the commands fail clearly instead of modifying inputs.

## Commands

QC report:

```bash
circyto scanpy-qc-report \
  --input /path/to/full_length.h5mu \
  --output-dir work/scanpy_qc
```

PCA / UMAP / Leiden:

```bash
circyto scanpy-pca \
  --input /path/to/full_length.h5mu \
  --output-dir work/scanpy_pca
```

## `scanpy-qc-report`

Operates on the RNA modality and computes:

- `total_counts`
- `n_genes_by_counts`
- `pct_counts_mt`

Writes:

- `qc_violin.png`
- `qc_scatter_counts_vs_genes.png`
- `qc_mt_vs_counts.png`
- `scanpy_qc_summary.json`

## `scanpy-pca`

Operates on the RNA modality and runs:

- `normalize_total`
- `log1p`
- `highly_variable_genes`
- `PCA`
- `neighbors`
- `UMAP`
- `Leiden`

Writes:

- `rna_umap.png`
- `rna_leiden.tsv`
- `exploratory_rna_processed.h5ad`
- `scanpy_pca_summary.json`

## Scope

These outputs are exploratory only:

- the original `.h5mu` file is never overwritten
- no automatic filtering is applied
- no biological interpretation is attached
- no figures beyond basic QC and UMAP are generated
