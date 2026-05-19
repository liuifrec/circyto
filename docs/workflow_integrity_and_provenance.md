# Workflow Integrity and Provenance

`circyto` workflow outputs should be predictable, inspectable, restartable, and easy to validate after the fact.

This document defines the current workflow provenance contract and the read-only integrity checker.

## Standard provenance fields

`workflow_summary.json` should contain these top-level fields:

- `command_name`
- `circyto_version`
- `workflow_type`
- `workflow_uuid`
- `protocol`
- `read_layout`
- `genome_fasta`
- `gtf`
- `detector_backend`
- `started_at`
- `completed_at`
- `elapsed_seconds`
- `hostname`
- `python_version`
- `cleanup_summary`
- `completed_stages`
- `skipped_stages`
- `failed_stages`
- `partial_outputs_detected`

If RNA profiling is present, `workflow_summary.json` should also contain:

- `rna_import`

## Workflow UUID

Each workflow summary carries an immutable `workflow_uuid`.

- workflow commands generate a UUID when the summary is first written
- post-hoc commands such as `circyto add-rna-profile` preserve the existing UUID when updating the summary

## Restartability markers

`circyto` records:

- `completed_stages`
- `skipped_stages`
- `failed_stages`
- `partial_outputs_detected`

These fields are intended to make resumability and post-run inspection more explicit.

## Cleanup summary

`cleanup_summary` standardizes post-run retention information:

- `enabled`
- `scope`
- `performed`
- `deleted_paths`
- `reclaimed_bytes`

When cleanup planning or execution adds richer fields, they may still appear under the workflow-specific `cleanup` block, but `cleanup_summary` is the stable top-level contract.

## Integrity checker

Use:

```bash
circyto check-workflow --workdir /path/to/workdir
```

Behavior:

- read-only only
- validates `workflow_summary.json`
- validates required provenance fields
- validates required manifest and matrix outputs
- validates `h5ad` existence when recorded
- validates RNA outputs if `rna/` exists
- validates parseability of discovered `*.provenance.json` files
- reports missing or corrupted components

The command exits with:

- `0` when the workflow passes integrity checks
- `1` when required components are missing or corrupted

## Example

```bash
circyto check-workflow \
  --workdir /user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192
```
