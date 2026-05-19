# Gene Expression and Velocity Integration

This document defines the future design for adding gene-expression and RNA-velocity-compatible layers to `circyto` full-length workflows.

It is a staged integration plan, not a production workflow.

## Current audit

Current state in `circyto`:

- `full-length-circrna` writes a circ-only `anndata/circ_counts.h5ad`
- that `h5ad` uses:
  - `X` = circRNA counts
  - `obs` = cell-level QC / metadata
  - `var` = circRNA feature metadata
  - `uns["circyto"]` = provenance
- circ-only matrix export is the default and validated behavior
- Smart-seq3 already has optional MuData export when a gene-count matrix is supplied
- current Smart-seq3 MuData support provides:
  - `mdata["rna"]` with total gene expression
  - `mdata["circ"]` with circRNA counts
  - shared `obs`
- current MuData support does not yet add velocity layers such as `spliced`, `unspliced`, or `ambiguous`

## Design options

### Option A: single AnnData

Possible designs:

- `X = circRNA`, with gene expression in `obsm`
- `X = gene expression`, with circRNA in `layers["circ_counts"]`
- mixed feature spaces in one object

Problems:

- mixes distinct biological feature spaces awkwardly
- makes RNA-velocity layer semantics less clear
- complicates downstream interpretation and tooling compatibility

### Option B: MuData

Recommended long-term design:

- `mdata["rna"]`
  - `X` = total gene expression
  - `layers["spliced"]`
  - `layers["unspliced"]`
  - `layers["ambiguous"]`
  - `var` = gene feature table
- `mdata["circ"]`
  - `X` = circRNA counts
  - `var` = circRNA coordinates, host gene, and circ-specific metadata
- shared `obs`
  - cell ID
  - protocol
  - dataset
  - QC
  - replication-state labels where available
- shared `uns["circyto"]`
  - provenance
  - reference metadata
  - workflow parameters

## Recommendation

Preferred long-term path: `MuData`.

Reason:

- RNA and circRNA are separate modalities
- velocity layers naturally belong under the RNA modality
- circRNA retains a clean modality-specific feature space
- shared cell metadata remains explicit

Backward compatibility rule:

- keep `anndata/circ_counts.h5ad` as the default circ-only output
- add multimodal outputs only as explicit opt-in

## Proposed future CLI

Future full-length workflow flags:

- `--export-h5ad`
- `--export-mudata`
- `--gene-expression-method none|simple-overlap|featurecounts|velocyto`
- `--velocity-layers none|velocyto`
- `--cleanup-intermediates alignments|demux|all`

Current scaffold status:

- `simple-overlap` is implemented for lightweight internal gene-expression sanity counting from workflow alignments plus GTF gene intervals
- `--gene-counts` external RNA snapshot import is implemented
- `featurecounts` and `velocyto` remain planned, not production
- cleanup execution is implemented as explicit opt-in post-success retention control
- velocity-compatible layers are still planned only

Formal contract reference:

- [RNA and velocity contract](rna_velocity_contract.md)

## Proposed workflow architecture

```text
FASTQ/BAM
-> circRNA detection with circyto
-> circRNA matrix + circ-only h5ad
-> optional gene-expression quantification
-> optional velocity-compatible layer generation
-> optional MuData export
-> optional cleanup planning
```

Proposed RNA-side quantification sources:

- `simple-overlap`
  - lightweight internal sanity profile
  - useful for chr21-scale testing and contract validation
  - not a replacement for production-grade RNA quantification
- `featurecounts`
  - practical total-gene matrix bridge
  - easier initial implementation
- `velocyto`
  - velocity-compatible layer source
  - must remain optional
  - should not become a hard dependency

## What is implemented now

Implemented safely now:

- design and schema documentation
- future flag validation for `full-length-circrna`
- external `gene_counts.tsv` RNA snapshot import
- internal `simple-overlap` RNA sanity counting
- explicit `NotImplementedError` for `featurecounts` / `velocyto` / velocity layers
- scoped cleanup planning and post-success cleanup execution

## What remains future work

### Phase 1: schema stabilization

- finalize RNA modality schema
- finalize circ modality schema
- finalize provenance block

### Phase 2: gene-expression import/export

- support `featurecounts`-style total gene matrices for full-length workflows
- align cell IDs robustly
- export MuData with `rna` + `circ`

### Phase 3: velocity-compatible layer support

- support `velocyto`-derived `spliced`, `unspliced`, and `ambiguous`
- preserve method metadata and references
- keep runtime optional

### Phase 4: workflow integration

- integrate multimodal export into validated public human workflows
- maintain circ-only `h5ad` as the default behavior

### Phase 5: study-level usage

- replication-state-aware scRR analyses
- circRNA plus gene-expression plus velocity interpretation
- no production claims until public datasets validate the contract

## Cleanup policy

Cleanup must obey all of the following:

- explicit opt-in only
- only after successful workflow completion
- preserve:
  - `workflow_summary.json`
  - detector summaries
  - `matrix/`
  - `anndata/`
  - `mudata/`
  - `qc/`
  - `manifests/`
  - `logs/`
  - provenance
- delete only large regenerable alignment intermediates
- dry-run must report planned cleanup without deleting

See also:

- [Intermediate cleanup policy](intermediate_cleanup_policy.md)

## Practical conclusion

The correct long-term architecture is:

- circ-only `h5ad` for backward compatibility
- MuData for future RNA + circ multimodal export
- optional velocity layers under the RNA modality only
- optional cleanup after successful export, never by default
