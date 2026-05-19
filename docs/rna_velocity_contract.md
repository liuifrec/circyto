# RNA and Velocity Contract

This document defines the formal future input and output contracts for adding gene-expression and velocity-compatible layers to `circyto`.

The intent is to make multimodal full-length workflows feel closer to Cell Ranger:

- clear inputs
- predictable outputs
- strong validation
- explicit metadata
- reproducible summaries

## Scope

Current validated default output:

- `anndata/circ_counts.h5ad`

Future optional output:

- `mudata/full_length.h5mu`

## Future input contracts

### A. Gene count table import

Expected file: `gene_counts.tsv`

Required columns:

- `gene_id`
- `gene_name`
- one or more cell columns

Rules:

- `gene_id` must be unique
- cell columns must be non-empty
- counts should be integer-like
- cell column names must exactly match workflow cell IDs before MuData export

Current implemented behavior in `full-length-circrna`:

- `--gene-counts <TSV>` validates this contract
- writes a normalized snapshot under:
  - `rna/gene_counts.tsv`
  - `rna/gene_feature_table.tsv`
  - `rna/rna_import_summary.json`
- does not yet merge the RNA counts into `circ_counts.h5ad`

Example:

```tsv
gene_id	gene_name	cellA	cellB
ENSG000001	GENE1	10	0
ENSG000002	GENE2	4	7
```

### B. Spliced / unspliced / ambiguous layer import

Expected future directory layout:

```text
velocity_layers/
  barcodes.tsv
  features.tsv
  spliced.mtx
  unspliced.mtx
  ambiguous.mtx
```

Required files:

- `barcodes.tsv`
- `features.tsv`
- `spliced.mtx`
- `unspliced.mtx`

Optional file:

- `ambiguous.mtx`

Rules:

- `barcodes.tsv` cell IDs must be unique
- `features.tsv` gene IDs must be unique
- cell IDs must match the target workflow cell IDs exactly before final export

Example `features.tsv`:

```tsv
ENSG000001	GENE1
ENSG000002	GENE2
```

### C. Velocyto loom import

Future import contract:

- accept a loom file only through an explicit importer step
- normalize loom-derived layers to the same internal RNA contract:
  - total / spliced / unspliced / ambiguous
- do not expose raw loom assumptions downstream

Current status:

- not implemented
- no velocyto dependency is required

### D. featureCounts-style gene-count import

Future contract:

- support flat gene-count tables that can be normalized into:
  - `gene_id`
  - `gene_name`
  - cell columns

Preferred exported normalized form for `circyto` internals:

- `gene_counts.tsv` with explicit `gene_id` and `gene_name`

### E. Future internal counting from BAM/SAM

Future internal counting should remain optional.

Requirements before implementation:

- explicit method selection
- strict reference / annotation provenance
- reproducible summaries
- no hidden assumptions about strandedness or protocol

Current status:

- not implemented

## Future output contract

### Default output

Keep:

- `anndata/circ_counts.h5ad`

Contract:

- `X` = circRNA counts
- `obs` = cell-level metadata and QC
- `var` = circRNA feature table
- `uns["circyto"]` = provenance

### Optional future multimodal output

Planned file:

- `mudata/full_length.h5mu`

Recommended modalities:

- `rna`
- `circ`

### Required `obs` columns

Minimum shared `obs` columns:

- `sample_id`
- `protocol`
- `read_layout`
- `strandedness`
- `alignment_status`
- `detector_status`
- `circRNA_count`
- `total_circRNA_support`

Recommended additional columns when available:

- `dataset_id`
- `cell_type`
- `replication_state`
- `condition`
- `assigned_reads`

### Required `rna.var` columns

Minimum:

- `gene_id`
- `gene_name`

Optional:

- `gene_biotype`
- `chrom`
- `start`
- `end`
- `strand`

### Required `circ.var` columns

Minimum:

- `circ_id`
- `chrom`
- `start`
- `end`
- `strand`
- `host_gene`

### Required RNA layers

Minimum future `rna` modality expectations:

- `X` = total gene expression

Optional future velocity-compatible layers:

- `layers["spliced"]`
- `layers["unspliced"]`
- `layers["ambiguous"]`

### Required `uns["circyto"]` provenance fields

Minimum fields:

- `workflow_name`
- `protocol`
- `reference`
- `gtf`
- `manifest`
- `command_options`
- `export_h5ad`
- `export_mudata`
- `gene_expression_method`
- `velocity_layers`
- `cleanup_intermediates`

Recommended future fields:

- `dataset_id`
- `software_versions`
- `cell_id_join_policy`
- `counting_reference_contract`
- `velocity_source`

## Current validation helpers

Current lightweight preflight helpers:

- `validate_gene_expression_table_schema`
- `validate_velocity_layers_schema`
- `validate_cell_id_consistency`
- `validate_feature_id_uniqueness`

These helpers validate contract shape only. They do not count reads, run velocyto, or construct real RNA matrices.

## Future command shapes

Future multimodal export with normalized gene counts:

```bash
circyto workflow full-length-circrna \
  --manifest manifest.tsv \
  --outdir work/full_length_run \
  --protocol ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --export-h5ad \
  --export-mudata \
  --gene-expression-method featurecounts
```

Future velocity-compatible export:

```bash
circyto workflow full-length-circrna \
  --manifest manifest.tsv \
  --outdir work/full_length_run \
  --protocol ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --export-mudata \
  --gene-expression-method velocyto \
  --velocity-layers velocyto
```

Current status of both command shapes:

- documented only
- not implemented yet
