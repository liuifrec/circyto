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
  --star-index /path/to/star_index \
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
- paired-end rows use the validated STAR+CIRI3 paired-end path
- real paired-end execution requires `--allow-paired-ramda`
- dry-run planning does not require `--allow-paired-ramda`
- paired-end rows also require `--star-index`
- the paired-end route is now locally validated on a real `GSE278952 / SRR30911454` chr21 subset

### `shin-ramda`

- demux skipped implicitly
- single-end rows use `bwa-mem` plus direct `CIRI3` SAM mode
- paired-end rows use the validated STAR+CIRI3 paired-end path
- real paired-end execution requires `--allow-paired-ramda`
- dry-run planning does not require `--allow-paired-ramda`
- paired-end rows also require `--star-index`
- the paired-end route is now locally validated on a real `GSE278952 / SRR30911454` chr21 subset

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

## Paired-End RamDA

The paired-end RamDA/Shin-RamDA execution path remains opt-in, but it is no longer only a dry-run scaffold:

- dry-run is always allowed so you can inspect the STAR+CIRI3 route
- real execution requires `--allow-paired-ramda`
- `--experimental-paired-ramda` still works as a deprecated alias
- the executable route was locally validated on a real `GSE278952 / SRR30911454` chr21 subset
- full hg38-scale biological validation against human paired-end scRR data is still in progress

Example paired-end dry-run:

```bash
circyto workflow full-length-circrna \
  --manifest paired_manifest.tsv \
  --outdir work/full_length_pe_plan \
  --protocol ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --star-index /path/to/star_index \
  --dry-run
```

Example paired-end real execution:

```bash
circyto workflow full-length-circrna \
  --manifest paired_manifest.tsv \
  --outdir work/full_length_pe_run \
  --protocol ramda \
  --genome-fasta ref.fa \
  --gtf genes.gtf \
  --star-index /path/to/star_index \
  --allow-paired-ramda \
  --export-h5ad
```

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
