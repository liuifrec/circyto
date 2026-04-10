# Alignment-First Execution

## Why this exists

The original circyto detector orchestration path runs each detector per cell from FASTQ input. That means repeated alignment work across:

- every cell
- every rerun
- every detector that shells out to an aligner internally

For small smoke tests that is acceptable. For large demultiplexed datasets it becomes the dominant cost and makes reruns impractical.

## New flow

The alignment-first track splits execution into explicit phases:

1. source manifest -> prepare or reuse alignments
2. alignment manifest -> run alignment-capable detector backends
3. per-cell normalized TSV -> collect-matrix

This keeps the downstream TSV schema stable while making the expensive alignment stage reusable.

## Core commands

Prepare reusable alignments:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache \
  --chunk-size 48 \
  --sentinel-cells 8
```

Generate the alignment manifest via the alias:

```bash
circyto align-manifest \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --outdir work/alignment_cache
```

Run an alignment-capable detector:

```bash
circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_run \
  --ref-fa ref/genome.fa \
  --chunk-size 64
```

## Alignment manifest

The alignment manifest is a TSV with one row per cell. Core columns:

```text
cell_id    bam    sam    group_id    read_layout    aligner    reference    cache_key    source_manifest
```

Only one of `bam` or `sam` should be populated for a row.

## Cache behavior

Alignment outputs are cached by a hash over:

- source reads or input BAM
- read layout
- aligner strategy
- reference
- detector hint
- alignment flags

This is intended to let future detectors require different alignment caches without pretending one BAM always fits every backend.

## Resume and scheduling

The prepare phase writes:

- `alignment_manifest.tsv`
- `alignment_prepare_summary.json`
- `alignment_prepare_plan.json` when `--dry-run` is used
- per-alignment provenance JSON beside cached BAM/SAM files
- `chunks/chunk_*.json` per chunk

Useful controls:

- `--sentinel-cells`: run the first N cells before the rest
- `--chunk-size`: checkpoint the manifest progressively
- `--parallel`: bound concurrent work
- `--fail-fast`: stop after the first failed chunk
- `circyto summarize-chunks --indir ...`: inspect chunk completion state
- `circyto summarize-run-state --manifest ... --run-dir ... --mode prepare|detector`: audit missing and stale cells
- `circyto export-run-subset --manifest ... --run-dir ... --subset failed|missing|stale|incomplete --out subset.tsv`: export rerun manifests
- rerun the same command against the same `--outdir` to resume cached work; chunk history appends and detector summaries retain the whole-run cell set during subset retries

## CIRI3 execution modes

circyto supports two CIRI3 execution modes in the alignment-first path.

### BWA mode (validated locally)

- validated locally with BWA + CIRI3 on a chr21 pilot
- uses unsorted SAM input
- CIRI3 mapper mode: `-Ma 0`
- recommended single-cell stringency: `-S 0`
- recommended local BWA parameters: `bwa mem -k 15 -T 15`

Minimal example:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache/alignment_manifest.tsv \
  --outdir work/ciri3_run \
  --ref-fa ref/genome.fa
```

Assumptions:

- the BWA index matches `ref/genome.fa`
- the same reference build is used for alignment preparation and detector execution
- the manifest rows resolve to unsorted SAM for direct CIRI3 execution
- `samtools` is available for alignment handling and inspection

### STAR mode (supported in code)

- implemented in code for alignment-first workflows
- requires STAR alignment and a matching STAR index
- uses CIRI3 mapper mode: `-Ma 1`
- requires STAR chimeric outputs in the alignment manifest
- requires a paired BWA rescue SAM when using the official CIRI3 STAR hybrid input path
- requires `samtools` for alignment handling and inspection
- runs STAR in a local Linux temp workspace before copying artifacts back into the cache; use `CIRCYTO_STAR_TMPDIR` to force node-local scratch when needed

Minimal example:

```bash
circyto prepare-alignment-cache \
  --manifest manifest.tsv \
  --aligner star \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/alignment_cache_star

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/alignment_cache_star/alignment_manifest.tsv \
  --outdir work/ciri3_star_run \
  --ref-fa ref/genome.fa
```

Note: STAR mode now completes through the normal CLI on a small real chr21 subset with the official hybrid contract. BWA + CIRI3 remains the validated baseline, and larger STAR runs should still be rechecked on the target server environment.

## CIRI3 template validation

Template execution remains available when you need to pin an explicit local command contract.

Preflight before a large run:

```bash
circyto validate-ciri3-template \
  --template 'ciri3 --bam {alignment} --out {raw_output} --threads {threads} --cell {cell_id} --tmp {outdir} --fmt {alignment_format}'
```

Required placeholders:

- `alignment`
- `alignment_format`
- `cell_id`
- `threads`
- `raw_output`
- `outdir`

Optional placeholders:

- `ref_fa`
- `gtf`
- `extra_args`
- `read_layout`
- `group_id`
- `log_path`

Template mode does not require `ref_fa` unless your template explicitly uses `{ref_fa}`. Direct `java -jar` mode still requires `--ref-fa`.

Plan mode now records exact first-command previews for both alignment preparation and CIRI3 detector execution when a template is configured.

## Small local validation

To validate the full alignment-first plumbing with local repo assets:

```bash
circyto smoke --detector ciri3 --aligner bwa-mem --outdir work/smoke_bwa
circyto smoke --detector ciri3 --aligner star --outdir work/smoke_star
```

This uses a tiny local chr21 subset. It validates the end-to-end workflow shape:

- FASTQ manifest
- alignment preparation
- alignment manifest
- CIRI3 execution
- `collect-matrix`

This smoke path validates workflow mechanics, not biological CIRI3 correctness. Empty smoke outputs can still be a PASS unless `--require-nonempty` is requested.

## PRJNA607968-style workflow

For a large demultiplexed rerun, the goal is to pay alignment cost once, checkpoint aggressively, and resume only missing work.

1. Plan the run:

```bash
circyto plan-alignment-cache \
  --manifest prjna607968_manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/prjna607968_align
```

2. Run a sentinel-first alignment preparation:

```bash
circyto prepare-alignment-cache \
  --manifest prjna607968_manifest.tsv \
  --aligner bwa-mem \
  --ref-fa ref/genome.fa \
  --detector ciri3 \
  --outdir work/prjna607968_align \
  --sentinel-cells 8 \
  --chunk-size 48 \
  --parallel 8 \
  --index-bam
```

3. Validate the resulting alignment manifest:

```bash
circyto manifest validate-alignment work/prjna607968_align/alignment_manifest.tsv --strict
```

4. Validate the CIRI3 template once if you intend to use explicit template execution:

```bash
circyto validate-ciri3-template \
  --template 'ciri3 --bam {alignment} --out {raw_output} --threads {threads} --cell {cell_id} --tmp {outdir} --fmt {alignment_format}'
```

5. Run CIRI3 from alignments:

```bash
export CIRCYTO_CIRI3_CMD_TEMPLATE='ciri3 --bam {alignment} --out {raw_output} --threads {threads} --cell {cell_id} --tmp {outdir} --fmt {alignment_format}'

circyto run-detector-from-alignments \
  --detector ciri3 \
  --manifest work/prjna607968_align/alignment_manifest.tsv \
  --outdir work/prjna607968_ciri3 \
  --ref-fa ref/genome.fa \
  --chunk-size 64 \
  --parallel 8
```

For direct local execution instead of template mode, ensure Java and a CIRI3 jar are detected by `circyto doctor`, and keep `--ref-fa` consistent with the reference used during alignment preparation.

6. Collect the matrix:

```bash
circyto collect-matrix \
  --detector ciri3 \
  --indir work/prjna607968_ciri3 \
  --outdir work/prjna607968_ciri3_matrix
```

7. Resume only incomplete work:

```bash
circyto summarize-chunks --indir work/prjna607968_align
circyto summarize-chunks --indir work/prjna607968_ciri3

# rerun the exact same prepare or detector command
# completed rows are skipped via cache/provenance
```

This is the operational difference from the older week-scale per-cell FASTQ model: alignments are checkpointed and reused, and failed chunks do not force a full project restart.

## Resume after failure

Export only failed cells from an interrupted prepare run:

```bash
circyto export-run-subset \
  --manifest prjna607968_manifest.tsv \
  --run-dir work/prjna607968_align \
  --subset failed \
  --out work/prjna607968_failed_prepare.tsv
```

Export only failed detector rows:

```bash
circyto export-run-subset \
  --manifest work/prjna607968_align/alignment_manifest.tsv \
  --run-dir work/prjna607968_ciri3 \
  --subset failed \
  --out work/prjna607968_failed_detector.tsv
```

Export one failed chunk as a standalone rerun manifest:

```bash
circyto export-run-subset \
  --manifest work/prjna607968_align/alignment_manifest.tsv \
  --run-dir work/prjna607968_ciri3 \
  --chunk-index 3 \
  --out work/prjna607968_chunk3.tsv
```

## Sentinel-first overnight recipe

Recommended overnight pattern:

1. `circyto plan-alignment-cache ... --preview-rows 5`
2. `circyto validate-ciri3-template --template '...'`
3. `circyto prepare-alignment-cache ... --sentinel-cells 8 --chunk-size 48`
4. `circyto run-detector-from-alignments ... --sentinel-cells 8 --chunk-size 64`
5. next morning, inspect:
   `circyto summarize-chunks --indir work/prjna607968_align`
   `circyto summarize-run-state --manifest ... --run-dir ... --mode detector`

## Provenance

The detector summary and per-cell records now include:

- `read_layout`
- `execution_mode`
- `input_mode`
- `reused_alignment`
- `detector_backend`

This makes reused alignment behavior explicit rather than hidden.

## Current scope

Implemented now:

- alignment manifest format
- reusable alignment cache plumbing
- built-in `bwa mem` preparation with sorted BAM output
- resumable prepare phase with chunking, per-chunk summaries, and sentinel-first scheduling
- alignment-native detector execution path with chunking and skip-existing support
- first-class `ciri3` detector registration plus validated template execution
- manifest export for failed, missing, incomplete, and chunk-specific reruns
- manifest-aware stale-output auditing and chunk summaries

Still future work:

- richer built-in aligner recipes beyond `reuse-existing` and `bwa-mem`
- multisample slicing and detector-specific partitioners
- additional alignment-native detector backends
