# CIRI3 Alignment-First Workflow

## Purpose

Use the alignment-first path when you want to:

- align once and reuse those alignments across reruns
- run CIRI3 through the intended `prepare-alignment-cache -> run-detector-from-alignments -> collect-matrix` sequence
- validate workflow mechanics before scaling to a larger cluster run

This is the first-class `ciri3` execution path in `circyto`.

## Required inputs

Source manifest:

- one row per cell
- `read1` required
- `read2` optional for paired-end rows
- `cell_id` required
- source manifests can infer layout from `read2` presence or absence in simple cases
- `read_layout` is still strongly recommended in the source manifest and always required in the generated alignment manifest
- layout inference does not make single-end STAR + CIRI3 supported

Reference/runtime:

- `--ref-fa` matching the aligner index used for the run
- `bwa` and `samtools`
- Java plus a usable CIRI3 jar for direct execution, or an explicit CIRI3 command template

Useful preflight:

```bash
circyto doctor
circyto detectors --json
```

## Workflow

Prepare reusable alignments:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --detector ciri3 \
  --ref-fa /path/to/genome.fa \
  --outdir work/alignment_cache \
  --threads 8 \
  --parallel 4 \
  --chunk-size 48
```

Run CIRI3 from the emitted alignment manifest:

```bash
circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3 \
  --ref-fa /path/to/genome.fa \
  --threads 8 \
  --parallel 4 \
  --chunk-size 64
```

Collect the matrix:

```bash
circyto collect-matrix \
  --detector ciri3 \
  --indir work/ciri3 \
  --outdir work/ciri3_matrix
```

## Alignment manifest contract

The generated alignment manifest is the handoff between alignment preparation and CIRI3 execution.

Core columns:

```text
cell_id    bam    sam    group_id    read_layout    aligner    reference    cache_key    source_manifest
```

Current expectations:

- paths are normalized to absolute paths
- `read_layout` must be explicit
- BWA-mode CIRI3 rows point to unsorted `sam`
- `reference` must match the FASTA used later with `run-detector-from-alignments --ref-fa`

Validate before a large run:

```bash
circyto manifest validate-alignment --strict work/alignment_cache/alignment_manifest.tsv
```

## Local chr21 sentinel example

This is the exact bounded validation pattern used locally after the CIRI3 hardening work. It uses real local data subsets plus the local `chr21` reference.

Single-end sentinel:

```bash
circyto prepare-alignment-cache \
  --manifest work/local_validation_20260410/se_manifest.tsv \
  --aligner bwa-mem \
  --detector ciri3 \
  --ref-fa ref/chr21.fa \
  --outdir work/post_hardening_ciri3_se/align \
  --threads 2 \
  --parallel 1 \
  --chunk-size 1

circyto manifest validate-alignment --strict \
  work/post_hardening_ciri3_se/align/alignment_manifest.tsv

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/post_hardening_ciri3_se/align/alignment_manifest.tsv \
  --outdir work/post_hardening_ciri3_se/ciri3 \
  --ref-fa ref/chr21.fa \
  --threads 2 \
  --parallel 1 \
  --chunk-size 1

circyto collect-matrix \
  --detector ciri3 \
  --indir work/post_hardening_ciri3_se/ciri3 \
  --outdir work/post_hardening_ciri3_se/matrix
```

Paired-end sentinel:

```bash
circyto prepare-alignment-cache \
  --manifest work/local_validation_20260410/pe_manifest_10k.tsv \
  --aligner bwa-mem \
  --detector ciri3 \
  --ref-fa ref/chr21.fa \
  --outdir work/post_hardening_ciri3_pe/align \
  --threads 2 \
  --parallel 1 \
  --chunk-size 1

circyto manifest validate-alignment --strict \
  work/post_hardening_ciri3_pe/align/alignment_manifest.tsv

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/post_hardening_ciri3_pe/align/alignment_manifest.tsv \
  --outdir work/post_hardening_ciri3_pe/ciri3 \
  --ref-fa ref/chr21.fa \
  --threads 2 \
  --parallel 1 \
  --chunk-size 1

circyto collect-matrix \
  --detector ciri3 \
  --indir work/post_hardening_ciri3_pe/ciri3 \
  --outdir work/post_hardening_ciri3_pe/matrix
```

Structural success criteria:

- `alignment_manifest.tsv` exists and validates
- `alignment_prepare_summary.json` exists
- `detector_run_summary.json` exists
- per-cell normalized TSV exists
- matrix files exist

Biological emptiness is acceptable for this sentinel if the run completes cleanly.

## Cluster whole-genome recipe

Replace the paths below with your cluster locations and keep the command order unchanged.

Environment preflight:

```bash
circyto doctor
circyto smoke --detector ciri3 --aligner bwa-mem --outdir /scratch/$USER/circyto_smoke
```

Alignment preparation:

```bash
circyto prepare-alignment-cache \
  --manifest /project/manifests/project_manifest.tsv \
  --aligner bwa-mem \
  --detector ciri3 \
  --ref-fa /project/ref/genome.fa \
  --outdir /scratch/$USER/project_alignments \
  --threads 16 \
  --parallel 8 \
  --chunk-size 64 \
  --sentinel-cells 8
```

Alignment manifest validation:

```bash
circyto manifest validate-alignment --strict \
  /scratch/$USER/project_alignments/alignment_manifest.tsv
```

Detector execution:

```bash
circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest /scratch/$USER/project_alignments/alignment_manifest.tsv \
  --outdir /scratch/$USER/project_ciri3 \
  --ref-fa /project/ref/genome.fa \
  --threads 16 \
  --parallel 8 \
  --chunk-size 64 \
  --sentinel-cells 8
```

Matrix collection:

```bash
circyto collect-matrix \
  --detector ciri3 \
  --indir /scratch/$USER/project_ciri3 \
  --outdir /scratch/$USER/project_ciri3_matrix
```

After the sentinel succeeds, rerun the same `prepare-alignment-cache` and `run-detector-from-alignments` commands without `--sentinel-cells` for the full batch.

## STAR + CIRI3

Use STAR + CIRI3 only when you specifically want the official STAR hybrid path.

Current status:

- paired-end only
- real alignment-first execution is supported for paired-end rows
- single-end STAR + CIRI3 is not implemented in the current hybrid rescue path
- BWA + CIRI3 remains the baseline reliable path for sentinel and production use

Operator requirements:

- pass a usable STAR genome index explicitly
- use `--extra-flags "--genomeDir /path/to/star_index"`
- for gzipped FASTQs, also pass `--readFilesCommand zcat`
- keep `--ref-fa` consistent with the reference used for STAR and the downstream CIRI3 run

Paired-end STAR example:

```bash
circyto prepare-alignment-cache \
  --manifest /project/manifests/project_manifest.tsv \
  --aligner star \
  --detector ciri3 \
  --ref-fa /project/ref/genome.fa \
  --extra-flags "--genomeDir /project/ref/star_index" \
  --outdir /scratch/$USER/project_alignments_star \
  --threads 16 \
  --parallel 8 \
  --chunk-size 64 \
  --sentinel-cells 8

circyto manifest validate-alignment --strict \
  /scratch/$USER/project_alignments_star/alignment_manifest.tsv

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest /scratch/$USER/project_alignments_star/alignment_manifest.tsv \
  --outdir /scratch/$USER/project_ciri3_star \
  --ref-fa /project/ref/genome.fa \
  --threads 16 \
  --parallel 8 \
  --chunk-size 64 \
  --sentinel-cells 8

circyto collect-matrix \
  --detector ciri3 \
  --indir /scratch/$USER/project_ciri3_star \
  --outdir /scratch/$USER/project_ciri3_star_matrix
```

Paired-end STAR example for gzipped FASTQs:

```bash
circyto prepare-alignment-cache \
  --manifest /project/manifests/project_manifest.tsv \
  --aligner star \
  --detector ciri3 \
  --ref-fa /project/ref/genome.fa \
  --extra-flags "--genomeDir /project/ref/star_index --readFilesCommand zcat" \
  --outdir /scratch/$USER/project_alignments_star \
  --threads 16 \
  --parallel 8 \
  --chunk-size 64 \
  --sentinel-cells 8
```

Single-end STAR + CIRI3 is rejected deliberately. Use BWA + CIRI3 instead.

## Smoke command

Recommended install/runtime smoke:

```bash
circyto smoke --detector ciri3 --aligner bwa-mem --outdir work/smoke_ciri3_bwa
circyto smoke --detector ciri3 --aligner bwa-mem --read-layout single-end --outdir work/smoke_ciri3_bwa_se
circyto smoke --detector ciri3 --aligner star --outdir work/smoke_ciri3_star
```

What smoke does now:

- stages packaged demo FASTQs and a tiny packaged reference under the chosen `--outdir`
- builds lightweight local indexes as needed
- runs a packaged alignment-first demo path that validates CLI plumbing
- writes `smoke_summary.json` plus a final single-cell matrix

Limitations:

- smoke is not a real biological CIRI3 validation
- smoke validates packaged demo plumbing and CLI/runtime wiring
- the packaged demo is intentionally tiny
- empty outputs may still count as a PASS unless `--require-nonempty` is used
- the default smoke path uses a packaged CIRI3-compatible template for stability on fresh machines rather than the full real-data execution path

## Notes

- BWA mode is the baseline validated CIRI3 path.
- STAR + CIRI3 should currently be treated as a paired-end-only operator path.
- Keep the reference build consistent across alignment preparation, detector execution, and any downstream interpretation.
- For cluster runs, start with `--sentinel-cells` and inspect the emitted JSON summaries before scaling out.
