# scverse Interoperability

`circyto` is designed to interoperate with the scverse ecosystem through explicit sparse matrices, `AnnData`, and `MuData`.

## Current compatibility

### AnnData

- `anndata/circ_counts.h5ad` remains the current default single-modality export
- the circ-focused `h5ad` is intended to be readable by standard `anndata` tooling

### MuData

- `circyto export-mudata` exports a multimodal `.h5mu`
- current modalities:
  - `rna`
  - `circ`
- future planned modalities:
  - `dna_cnv`
  - `dna_snv`

## h5mu expectations

Current `h5mu` outputs should provide:

- `mdata["rna"]`
- `mdata["circ"]`
- shared `mdata.obs`
- `mdata.uns["circyto"]` provenance

The `circ` modality should preserve circ feature coordinates and support downstream host-gene summarization and BED export.

## Scanpy interoperability

`circyto` provides optional Scanpy-facing downstream helpers for exploratory RNA analysis:

- `circyto scanpy-qc-report`
- `circyto scanpy-pca`

These are explicitly exploratory and do not overwrite the original `h5mu`.

## Genome-browser interoperability

`circyto export-circ-bed` exports a BED-like table with:

- `chrom`
- `start`
- `end`
- `circ_id`
- `support`

This supports lightweight genome-browser and interval-tool interoperability.

## Workdir validation

Use:

```bash
circyto validate-workdir --workdir WORKDIR --json
```

to validate:

- `workflow_summary.json`
- matrix outputs
- RNA outputs when present
- `h5ad` readability
- `h5mu` readability
- `obs` consistency across modalities
- expected circ `var` columns
- QC summary presence
