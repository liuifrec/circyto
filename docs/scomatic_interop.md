# SComatic Interoperability

`circyto` can prepare full-length RNA inputs for external SComatic runs, import SComatic outputs, merge normalized RNA-derived candidate variant signals into MuData, and join exploratory summaries back to circRNA summaries.

This is intentionally an interoperability scaffold. `circyto run-scomatic` writes a command plan by default and only invokes external SComatic when `--execute` is explicitly supplied.

Recommended environment split:

1. Run `circyto prepare-scomatic-input` in the circyto environment.
2. Activate the SComatic environment, such as `conda activate scomatic-smoke`.
3. Run external SComatic BaseCellCounter, Step1, and Step2 from the generated plan or equivalent server script.
4. Reactivate the circyto environment.
5. Run `circyto import-scomatic` and `circyto merge-scomatic`.

Local WSL/native SComatic execution may be unstable. Prefer server, container, or HPC execution for real SComatic runs.

## Scientific framing

RNA-derived SComatic calls should be treated as exploratory RNA-derived candidate variant signals, not validated DNA somatic mutations unless orthogonal DNA validation exists.

RNA-based variant calls can reflect transcriptional noise, RNA editing, mapping artifacts, allele-specific expression, and protocol-specific biases. The exported and joined tables are therefore meant for hypothesis generation and prioritization.

## Validation Status

Validated:

- HAP1 batch10 technical smoke through BaseCellCounter, Step1, Step2, and normalization
- real Step1/Step2 output normalization into `scomatic_candidate_summary.tsv`
- BAM preparation mechanics for one-cell-per-BAM/SAM inputs

Partially validated:

- `candidate_snv` MuData export for RNA-derived candidate variant signals
- dry-run command planning for split circyto and SComatic environments
- circRNA summary joins using normalized candidate tables

Not yet validated:

- clinically validated variants
- validated somatic mutation discovery
- genome-scale biological interpretation of SComatic outputs
- robust candidate behavior in homogeneous datasets without meaningful cell-type/group contrasts
- local WSL SComatic execution as a supported path

## Lessons Learned from HAP1 and IMR90

- HAP1 batch10 proves the technical adapter path can run, including CB tag handling and Step1/Step2 normalization.
- HAP1 Step2 produced no candidate calls in the single-cell-type smoke test; this is a realistic limitation for homogeneous designs, not a circyto failure.
- IMR90 shows the main scRR DNA integration path should remain processed DNA state such as CNV, while SComatic remains an optional RNA-derived candidate variant signals sidecar.
- Both datasets argue for conservative manuscript language and explicit separation between DNA modalities and RNA-derived candidate variant signals.

## Commands

### `circyto prepare-scomatic-input`

Generic full-length RNA preparation command for Smart-seq2, Smart-seq3, RamDA-seq, ShinRamDA, and the scRR RNA arm.

It supports one-cell-per-BAM/SAM manifests and already merged multi-cell BAMs. See [Full-length RNA SComatic adapter](scomatic_full_length_adapter.md).

### `circyto run-scomatic`

Writes an external SComatic command plan for BaseCellCounter and optional Step1/Step2. Add `--execute` only when the SComatic environment and parameters are ready.

### `circyto import-scomatic`

Imports SComatic Step1, Step2, or generic candidate tables into `scomatic_candidate_summary.tsv`.

### `circyto merge-scomatic`

Adds normalized RNA-derived candidate variant signals to MuData as the exploratory `candidate_snv` modality.

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
- Per-cell RNA-derived candidate variant signals count from the SComatic candidate table.
- Host-gene circRNA summary when `circ_feature_table.tsv` contains `host_gene`.
- RNA-derived candidate variant signals gene summary when the candidate table contains a gene-like column.
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
- SComatic execution occurs only for `run-scomatic --execute`.
- The joined outputs are summary tables for exploratory interpretation, not validation artifacts.

For the longer-term integrated study direction, see:

- [SComatic-integrated circRNA study design](scomatic_circrna_study_design.md)
- [Full-length RNA SComatic adapter](scomatic_full_length_adapter.md)
