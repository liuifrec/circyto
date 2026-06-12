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
- validated scRR tri-modal extension:
  - `cnv`
- future exploratory modality:
  - `candidate_snv`

## h5mu expectations

Current `h5mu` outputs should provide:

- `mdata["rna"]`
- `mdata["circ"]`
- `mdata["cnv"]` when merged with `circyto merge-scrr-cnv`
- shared `mdata.obs`
- `mdata.uns["circyto"]` provenance

The `circ` modality should preserve circ feature coordinates and support downstream host-gene summarization and BED export.

CircRNA feature metadata should include:

- `host_gene`: final best host-gene annotation
- `host_gene_source`: `gtf`, `circatlas`, `circatlas_id`, or empty
- `host_gene_from_gtf`: GTF/GFF coordinate-overlap annotation only
- `host_gene_from_circatlas`: database-table host-gene annotation only
- `host_gene_from_circatlas_id`: fallback parsed from IDs such as `hsa-UBAP2_0052`

Priority is GTF/GFF overlap first, circAtlas host-gene table fields second, and circAtlas ID parsing third. Missing host genes are allowed for novel or unannotated circRNAs.

Existing files can be repaired without rerunning detection:

```bash
circyto repair-host-genes --input circ_counts.h5ad --output circ_counts.hostgene_fixed.h5ad
circyto repair-host-genes --input full_length.h5mu --output full_length.hostgene_fixed.h5mu --circ-mod circ
```

Add `--gtf annotation.gtf` to apply coordinate overlap before circAtlas fallbacks.

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
