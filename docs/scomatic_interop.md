# SComatic Interoperability

`circyto` can now prepare lightweight interoperability tables for an external SComatic run and can join exploratory SNV candidates back to circRNA summaries.

This is intentionally an interoperability scaffold, not an in-process SComatic integration. `circyto` does not install, invoke, or validate SComatic here.

## Scientific framing

RNA-derived SComatic calls should be treated as exploratory candidate SNVs, not validated DNA somatic mutations unless orthogonal DNA validation exists.

RNA-based variant calls can reflect transcriptional noise, RNA editing, mapping artifacts, allele-specific expression, and protocol-specific biases. The exported and joined tables are therefore meant for hypothesis generation and prioritization.

## Commands

### `circyto export-scomatic-inputs`

Required inputs:

- `--bam-manifest <TSV>`
- `--cell-metadata <TSV>`
- `--outdir <DIR>`
- `--reference-fasta <FA>`
- `--protocol smartseq3|ramda|shin-ramda`

Outputs:

- `scomatic_bam_list.tsv`
- `scomatic_celltypes.tsv`
- `README_scomatic_next_steps.md`

Example:

```bash
circyto export-scomatic-inputs \
  --bam-manifest work/alignment_manifest.tsv \
  --cell-metadata work/cell_metadata.tsv \
  --reference-fasta ref/genome.fa \
  --protocol ramda \
  --outdir work/scomatic_export
```

### `circyto join-circ-snv-summary`

Required inputs:

- `--circ-matrix <TSV>`
- `--circ-feature-table <TSV>`
- `--scomatic-candidates <TSV>`
- `--cell-metadata <TSV>`
- `--outdir <DIR>`

Outputs:

- `circ_snv_cell_summary.tsv`
- `circ_snv_host_gene_summary.tsv`

Current join behavior:

- Per-cell circRNA count from a wide circRNA-by-cell TSV.
- Per-cell circRNA read-support sum when the matrix values are numeric.
- Per-cell candidate SNV count from the SComatic candidate table.
- Host-gene circRNA summary when `circ_feature_table.tsv` contains `host_gene`.
- Candidate SNV gene summary when the candidate table contains a gene-like column.
- Missing optional columns generate warnings instead of hard failures.

Example:

```bash
circyto join-circ-snv-summary \
  --circ-matrix work/circ_counts.tsv \
  --circ-feature-table work/circ_feature_table.tsv \
  --scomatic-candidates work/scomatic_candidates.tsv \
  --cell-metadata work/cell_metadata.tsv \
  --outdir work/circ_snv_join
```

## Scope limits

- No SComatic installation is required for these commands.
- No SComatic execution occurs inside `circyto`.
- The joined outputs are summary tables for exploratory interpretation, not validation artifacts.

For the longer-term integrated study direction, see:

- [SComatic-integrated circRNA study design](scomatic_circrna_study_design.md)
