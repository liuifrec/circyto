# MuData Downstream

`circyto` provides lightweight read-only inspection commands for completed `.h5mu` outputs.

## Commands

Inspect structure:

```bash
circyto inspect-mudata \
  --input /path/to/full_length.h5mu
```

JSON mode:

```bash
circyto inspect-mudata \
  --input /path/to/full_length.h5mu \
  --json
```

Summarize QC:

```bash
circyto summarize-mudata-qc \
  --input /path/to/full_length.h5mu
```

JSON mode:

```bash
circyto summarize-mudata-qc \
  --input /path/to/full_length.h5mu \
  --json
```

## `inspect-mudata`

Reports:

- modalities present
- `n_obs`
- RNA shape
- circ shape
- `obs` columns
- `rna.var` columns
- `circ.var` columns
- `uns["circyto"]` keys
- shared / RNA-only / circ-only cell counts when `membership` is present

## `summarize-mudata-qc`

Reports:

- `total_rna_counts` summary
- `detected_genes` summary
- `mitochondrial_fraction` summary
- `ribosomal_fraction` summary
- `circRNA_count` summary
- `circRNA_total_support` summary
- Pearson correlation between RNA total counts and circRNA count when both are available and non-constant

## Scope

These commands are descriptive only:

- no Scanpy normalization
- no dimensional reduction
- no plotting
- no biological claims
- no mutation of the `.h5mu` file
