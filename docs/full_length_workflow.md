# Full-Length CircRNA Workflow

`circyto workflow full-length-circrna` is the manifest-driven high-level workflow for full-length single-cell RNA-seq circRNA detection.

It is intended for per-cell FASTQ inputs such as RamDA and Shin-RamDA, where each manifest row already represents one cell or one library.

## Command

```bash
circyto workflow full-length-circrna \
  --manifest manifest.tsv \
  --outdir work/full_length_run \
  --protocol ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --detector ciri3 \
  --threads 8 \
  --export-h5ad
```

## Manifest

Supported columns include:

- `sample_id` or `cell_id`
- `fastq_1` or `read1`
- `fastq_2` or `read2`
- `protocol`
- `strandedness`
- `read_layout`

The current workflow expects one cell or one library per row.

## Stage Graph

For RamDA and Shin-RamDA:

1. manifest ingest
2. demux skipped
3. protocol-aware alignment prep
4. CIRI3 calling
5. matrix collection
6. circ-only AnnData export
7. workflow summary and QC tables

For SMART-Seq3:

- pooled SMART-Seq3 demux remains the existing validated path:
  `circyto workflow smartseq3-ciri3`
- `full-length-circrna` only supports SMART-Seq3 when the inputs are already demultiplexed per-cell FASTQs and `--skip-demux` is passed

## Protocol Behavior

### `ramda`

- demux skipped implicitly
- single-end rows use `bwa-mem` plus direct `CIRI3` SAM mode
- paired-end rows are only supported as `--dry-run` planning at present

### `shin-ramda`

- demux skipped implicitly
- single-end rows use `bwa-mem` plus direct `CIRI3` SAM mode
- paired-end rows are only supported as `--dry-run` planning at present

### `smartseq3`

- already-demultiplexed per-cell manifests can run with `--skip-demux`
- pooled SMART-Seq3 should continue using `circyto workflow smartseq3-ciri3`

## Outputs

The workflow writes:

- `align/`
- `ciri3/`
- `matrix/`
- `qc/cell_qc.tsv`
- `qc/circ_qc.tsv`
- `anndata/circ_counts.h5ad` when `--export-h5ad`
- `workflow_summary.json`

## Dry Run

Use `--dry-run` to plan alignment and detector stages without executing them:

```bash
circyto workflow full-length-circrna \
  --manifest manifest.tsv \
  --outdir work/full_length_plan \
  --protocol shin-ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --dry-run
```

The dry run still writes `workflow_summary.json` and preserves the stage graph, planned outputs, and underlying alignment and detector plans.
