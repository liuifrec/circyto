# Intermediate Cleanup Policy

This document defines the intended cleanup policy for large regenerable workflow intermediates in `circyto`.

Current status:

- cleanup execution is implemented for explicit `full-length-circrna` opt-in
- dry-run still reports the cleanup plan without deleting files
- the current scaffold is limited to workflow-owned files under the output directory

## Goals

- reduce retained disk usage after successful export
- preserve user-facing outputs and provenance
- avoid deleting anything by default
- keep cleanup explicit, auditable, and Cell Ranger-like

## Retention Classes

### MUST_KEEP

Always preserve:

- `workflow_summary.json`
- detector summaries
- `matrix/`
- `anndata/`
- `mudata/` if present
- `qc/`
- `manifests/`
- `logs/`
- provenance JSON
- final detector TSVs

### SAFE_TO_DELETE_AFTER_SUCCESS

Delete only regenerable workflow-owned intermediates after a successful run, for example:

- `align/cache/`
- large `*.sam` / `*.bam` / `*.bai` intermediates
- STAR temporary outputs
- BWA rescue SAM/BAM
- demux per-cell FASTQs generated from pooled Smart-seq2/3 input
- temporary chunk files if reproducible

### USER_INPUTS_NEVER_DELETE

Never delete:

- original raw FASTQs
- original pooled Smart-seq FASTQs
- user-supplied manifests
- user-supplied `gene_counts.tsv`

`circyto` should treat these as immutable inputs even when they sit near workflow outputs. The current planner only inspects files inside the workflow output directory and does not mark external inputs for deletion.

## Planned Cleanup Scopes

The long-term CLI should support selective cleanup scopes:

- `--cleanup-intermediates alignments`
- `--cleanup-intermediates demux`
- `--cleanup-intermediates all`

Current implementation note:

- `full-length-circrna` accepts:
  - `--cleanup-intermediates alignments`
  - `--cleanup-intermediates demux`
  - `--cleanup-intermediates all`
- cleanup runs only after successful workflow completion
- failed runs skip cleanup entirely

## Guardrails

- cleanup is opt-in only
- cleanup should occur only after successful workflow completion
- dry-run should report planned deletions without removing files
- cleanup must never be applied to user inputs
- failure during cleanup must not invalidate already written final outputs
- cleanup should never touch files outside the workflow output directory

## Workflow-Specific Notes

### Smart-seq2/3 Pooled Workflows

Generated per-cell demultiplexed FASTQs may dominate disk usage after matrix and `h5ad` export. These demux FASTQs are workflow-owned intermediates and are reasonable `SAFE_TO_DELETE_AFTER_SUCCESS` candidates, but the original pooled FASTQs remain `USER_INPUTS_NEVER_DELETE`.

### Full-Length RamDA / scRR Workflows

Large STAR, BWA, SAM, BAM, and related rescue artifacts are regenerable and should be considered for cleanup only after the workflow summary, matrix outputs, QC tables, and `h5ad` export have been written successfully.

## Current Scaffold

The current implementation supports:

- dry-run cleanup planning for `full-length-circrna`
- summary of candidate regenerable files and bytes
- distinction between alignment-like intermediates and generated demux FASTQs
- explicit reporting that failed workflows should not plan deletion
- execution of scoped cleanup for workflow-owned regenerable intermediates
- workflow summary fields:
  - `cleanup_performed`
  - `cleanup_deleted_paths`
  - `cleanup_reclaimed_bytes`
  - `cleanup_scope`
