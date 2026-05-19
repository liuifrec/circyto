# Intermediate Cleanup Policy

This document defines the intended cleanup policy for large regenerable workflow intermediates.

Current status:

- cleanup execution is not implemented yet
- cleanup planning is dry-run only

## Goals

- reduce retained disk usage after successful export
- preserve user-facing outputs and provenance
- avoid deleting anything by default
- keep cleanup explicit and auditable

## Preserve

Always preserve:

- `workflow_summary.json`
- detector summaries
- `matrix/`
- `anndata/`
- `mudata/`
- `qc/`
- `manifests/`
- `logs/`
- provenance files

## Delete candidates

Only delete regenerable alignment intermediates, for example:

- unsorted SAM files
- temporary BAM files
- BAM indexes
- STAR chimeric junction intermediates
- STAR unmapped mate files
- BWA rescue SAM files

## Guardrails

- cleanup requires explicit `--cleanup-intermediates`
- cleanup should occur only after successful workflow completion
- dry-run reports planned deletions without removing files
- failure during cleanup must not invalidate already written final outputs
- cleanup should never touch files outside the workflow output directory

## Current scaffold

The current implementation supports:

- dry-run cleanup planning for `full-length-circrna`
- summary of candidate regenerable files and bytes

The current implementation does not yet perform deletion.
