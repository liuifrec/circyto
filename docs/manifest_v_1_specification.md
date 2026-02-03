# circyto Manifest v1 Specification (LOCKED)

## Status
**LOCKED – v1.0**

This document defines the **stable, backward-compatible manifest format** used by circyto for all detector execution, batching, and multimodal export. Once locked, this spec **must not be broken**; future changes require a new version (v2).

---

## Purpose

The manifest provides a **cell-resolved input contract** between upstream data preparation (demux, BAM splitting, user-provided FASTQs) and downstream detector execution (`run-detector`, `run-batch`, `run-multidetector`).

Design goals:
- Simple, tabular, human-readable
- Backward-compatible with legacy circyto manifests
- Flexible enough to support multiple platforms (Smart-seq2, 10x, custom)
- Explicit about *what is required* vs *what is metadata*

---

## File Format

- **Delimiter:** TAB (`\t`)
- **Header:** REQUIRED
- **Encoding:** UTF-8
- **Compression:** optional (`.tsv` or `.tsv.gz`)

---

## Required Columns (v1)

### `cell_id`
- **Type:** string
- **Description:** Unique identifier for a cell/sample
- **Constraints:**
  - Must be unique within a manifest
  - Used as output directory name and matrix column ID

### `r1`
- **Type:** path (string)
- **Description:** Path to read 1 FASTQ (gzipped or plain)
- **Required if:** `bam` is empty

### `r2`
- **Type:** path (string or empty)
- **Description:** Path to read 2 FASTQ
- **Notes:**
  - May be empty for single-end data
  - Required for paired-end detectors

### `bam`
- **Type:** path (string or empty)
- **Description:** Path to BAM file
- **Rules:**
  - If provided, `r1` and `r2` should be empty
  - Mutually exclusive with FASTQ inputs

**Exactly one of the following input modes must be satisfied per row:**
- FASTQ mode: `r1` (and optionally `r2`) provided
- BAM mode: `bam` provided

---

## Optional Metadata Columns

These columns are **ignored by detectors** but preserved for provenance and tooling.

### `platform`
- Example: `smartseq2`, `10x`, `plate`

### `library_id`
- Dataset or library identifier (e.g. GEO accession)

### `n_input_reads`
- Integer
- Informational only (e.g. from demux summary)

### Future-safe rule
Any **unknown columns**:
- MUST be ignored by circyto core
- MUST NOT cause validation failure

---

## Backward Compatibility Rules

circyto **MUST accept** the following legacy aliases:

| Legacy column | v1 canonical |
|--------------|-------------|
| `read1` | `r1` |
| `read2` | `r2` |

Internally, readers normalize to `r1` / `r2`.

This ensures all historical manifests remain valid.

---

## Example: FASTQ-based Manifest

```tsv
cell_id	r1	r2
sc01	fastq/sc01_R1.fastq.gz	fastq/sc01_R2.fastq.gz
sc02	fastq/sc02_R1.fastq.gz	fastq/sc02_R2.fastq.gz
```

## Example: BAM-based Manifest

```tsv
cell_id	bam
cellA	bams/cellA.bam
cellB	bams/cellB.bam
```

## Example: Rich Manifest (demux output)

```tsv
cell_id	platform	r1	r2	library_id	n_input_reads
sc01	smartseq2	fastq/sc01_R1.fastq.gz	fastq/sc01_R2.fastq.gz	SRR11140636	593
```

---

## Validation Semantics

`circyto manifest validate` checks:
- Header exists
- `cell_id` present and unique
- Exactly one input mode (FASTQ or BAM) per row
- Referenced files exist (unless `--no-check-files` is used)

Strict mode (`--strict`) enforces all of the above.

---

## Versioning Policy

- **Manifest v1**: this document
- Any breaking change (required columns, semantics) → **Manifest v2**
- v1 support MUST remain indefinitely

---

## Non-goals (Explicitly Out of Scope)

- Embedding detector-specific parameters
- Encoding barcodes or UMIs
- Representing multimodal matrices

Those belong in **detector configs** or **downstream AnnData**, not the manifest.

---

## Summary

Manifest v1 is:
- Minimal
- Stable
- Backward-compatible
- Platform-agnostic

This contract is now **locked** and forms the foundation of circyto’s execution model.

