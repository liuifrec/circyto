# Post-hoc RNA Profile

`circyto add-rna-profile` adds a lightweight RNA gene-count snapshot to an already completed workflow directory without rerunning alignment, detector, matrix, or `h5ad` export stages.

Current supported method:

- `simple-overlap`

## Command shape

```bash
circyto add-rna-profile \
  --workdir /path/to/completed_workflow \
  --gtf /path/to/annotation.gtf \
  --method simple-overlap
```

Dry-run:

```bash
circyto add-rna-profile \
  --workdir /path/to/completed_workflow \
  --gtf /path/to/annotation.gtf \
  --method simple-overlap \
  --dry-run
```

## What it does

1. discovers an existing alignment manifest under the completed workflow
2. reuses the existing aligned SAM or BAM paths from that manifest
3. runs lightweight `simple-overlap` gene counting with the supplied GTF
4. writes RNA snapshot outputs under `WORKDIR/rna/`
5. updates `workflow_summary.json` if present

It does not:

- rerun demux
- rerun STAR
- rerun BWA
- rerun CIRI3
- rerun matrix collection
- rerun `h5ad` export
- delete any files

`circyto add-rna-profile` still requires usable alignment inputs because it recomputes `gene_counts.tsv` from the completed workflow's alignment manifest.

## Manifest discovery

The command searches these paths in order:

- `WORKDIR/align/alignment_manifest.tsv`
- `WORKDIR/align/star/alignment_manifest.tsv`
- `WORKDIR/align/bwa_mem/alignment_manifest.tsv`

## Outputs

The command writes:

- `WORKDIR/rna/gene_counts.tsv`
- `WORKDIR/rna/gene_feature_table.tsv`
- `WORKDIR/rna/rna_import_summary.json`
- `WORKDIR/qc/rna_qc.tsv`
- `WORKDIR/qc/rna_gene_qc.tsv`

If `WORKDIR/workflow_summary.json` exists, it is updated in place:

- existing fields are preserved
- `command_name` becomes `circyto add-rna-profile`
- `rna_import` is added or replaced

## `simple-overlap` behavior

- parse `gene` features from the GTF when present
- if the GTF lacks explicit `gene` features, fall back to gene intervals aggregated from exon records sharing the same `gene_id`
- count one read or read-pair template if its primary alignment overlaps exactly one gene interval
- exclude ambiguous multi-gene overlaps
- group paired-end records by QNAME per cell to reduce mate double-counting when possible

This is a lightweight sanity profile, not a production replacement for `featureCounts` or velocity-aware quantification.

## Example server command

For the completed Smart-seq3 all192 workflow folder:

```bash
circyto add-rna-profile \
  --workdir /user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192 \
  --gtf /path/to/gencode_or_refseq_annotation.gtf \
  --method simple-overlap
```

Dry-run first:

```bash
circyto add-rna-profile \
  --workdir /user/ifrec/liuyuchen/circyto_redo/emtab8735/work/diySpike_workflow_all192 \
  --gtf /path/to/gencode_or_refseq_annotation.gtf \
  --method simple-overlap \
  --dry-run
```

## Inspecting RNA outputs

After profiling completes, inspect the RNA snapshot with:

```bash
python scripts/check_rna_profile_outputs.py \
  --workdir /path/to/completed_workflow
```

JSON mode:

```bash
python scripts/check_rna_profile_outputs.py \
  --workdir /path/to/completed_workflow \
  --json
```

The script reports:

- whether `rna/gene_counts.tsv` exists
- whether `rna/gene_feature_table.tsv` exists
- whether `rna/rna_import_summary.json` exists
- `n_genes`
- `n_cells`
- `total_counts_sum`
- `assigned_templates`
- `ambiguous_templates_excluded`
- `unassigned_templates`
- top expressed genes by total count
- lowest / highest total RNA count cells
- whether RNA cell IDs match `matrix/cell_index.txt` when present

## Refreshing RNA QC After Cleanup

If cleanup has already removed alignment SAM/BAM intermediates, regenerate RNA QC from the existing RNA outputs instead of rerunning `add-rna-profile`:

```bash
circyto refresh-rna-qc \
  --workdir /path/to/completed_workflow
```

This command:

- reads existing `rna/gene_counts.tsv`
- reads existing `rna/gene_feature_table.tsv`
- updates `rna/rna_import_summary.json`
- rewrites `qc/rna_qc.tsv`
- rewrites `qc/rna_gene_qc.tsv`
- updates `workflow_summary.json` when present

It does not require surviving alignment files and does not rewrite `gene_counts.tsv`.

## RNA Circ Integration Summary

After RNA profiling and circRNA matrix generation are both available, summarize overlap and per-cell joined metrics with:

```bash
circyto summarize-rna-circ \
  --workdir /path/to/completed_workflow \
  --write-summary
```

This writes:

- `qc/rna_circ_cell_summary.tsv`
- `qc/rna_circ_summary.json`

RNA-only cells remain in the summary with `circRNA_count = 0`.
