# Circyto repo architecture (canonical)

This document is the canonical “map” of circyto’s internal structure and data contracts.
It is intentionally short, stable, and oriented around ownership of responsibilities.

---

## CLI entry points

### Main app
- `circyto/cli/circyto.py`
  - Root Typer app.
  - Registers subcommands (doctor/detectors/demux/manifest/smoke/…).
  - Defines CLI conventions (indir/outdir, default filenames).

### Smoke (one-command orchestration)
- `circyto/cli/smoke.py`
  - “One-command” flows intended for public data smoke tests and development.
  - Responsibilities:
    1) demux pooled FASTQ → per-cell FASTQ + manifest.tsv
    2) run detector over manifest
    3) collect matrix (mtx + indices)
    4) convert to h5ad (optionally loom)
  - Must be robust to: “0 assigned cells” (still write manifest + report).

---

## Data contract: Manifest v1 (locked)

Manifest is a TSV that describes per-cell inputs (fastq/bam) and minimal provenance.

**Required columns (v1):**
- `cell_id`
- `platform`
- `read1`
- `read2` (can be empty for single-end)
- `bam`   (can be empty if input_type=fastq)
- `library_id`
- `n_input_reads`

**Compatibility rule:**
- `pipeline/run_detector.py:read_manifest()` must accept BOTH:
  - legacy: `r1`, `r2`
  - v1: `read1`, `read2`

---

## Demux

### Smart-seq2 pooled demux
- `circyto/demux/smartseq2.py`
  - `SmartSeq2DemuxParams`
  - `demux_smartseq2_pooled()`
  - Writes:
    - `outdir/fastq/<cell>_R1.fastq.gz`
    - `outdir/fastq/<cell>_R2.fastq.gz`
    - `outdir/sink/unknown_barcode_R1.fastq.gz` (+ R2)
    - `outdir/demux_report.json`
    - `outdir/manifest.tsv` (Manifest v1; must be written even when 0 cells)

### CLI wrapper
- `circyto/cli/demux.py`
  - Parses CLI flags and calls `circyto.demux.smartseq2.demux_smartseq2_pooled()`
  - Optionally runs `circyto manifest validate`

---

## Detector running

### Manifest reader + runner
- `circyto/pipeline/run_detector.py`
  - `read_manifest(path)`
  - `run_detector_manifest(detector, manifest, outdir, ref_fa, gtf, threads, parallel)`
  - Constraints:
    - Respects `detector.max_parallel` if present.
    - Should treat manifest paths as authoritative (no inference).

### Engines registry
- `circyto/detectors/`
  - `build_default_engines()` returns name → engine
  - Engines implement `DetectorBase.run(inputs)`

---

## Collectors (detector outputs → matrix)

- `circyto/pipeline/collect.py`:
  - Generic/legacy collector for “per-cell TSV-like outputs”.
- `circyto/pipeline/collect_find_circ3.py`:
  - find-circ3 specific collector.
- `circyto/pipeline/collect_circexplorer2_matrix.py`:
  - circexplorer2 matrix collector.

**Collector output contract (default filenames):**
- `circ_counts.mtx`
- `circ_index.txt`
- `cell_index.txt`

---

## Convert / Multimodal export

- `circyto/writers/convert.py`
  - `convert_matrix_files(matrix_path, circ_index_path, cell_index_path, loom=None, h5ad=None)`
- `circyto/pipeline/export_multimodal.py`
  - Attach circRNA counts to AnnData (`obsm["X_circ"]`), plus metadata in `uns["circ"]`

---

## “Golden” end-to-end flow (smoke)

The canonical smoke flow is:

1) `circyto smoke smartseq2 ...`  
2) demux writes `demux/manifest.tsv`  
3) run-detector writes per-cell results under `run/<detector>/...`  
4) collect writes matrix under `matrix/<detector>/...`  
5) convert writes `h5ad/<detector>.h5ad`

Smoke is expected to work on public data (non-sensitive), locally or on permissive compute.
