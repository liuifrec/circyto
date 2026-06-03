# scRR Replication Timing Modality

## Scope

This document defines circyto's lightweight integration path for processed HAP1 scRR-seq DNA replication timing or replication-state outputs from `GSE278952`.

The HAP1 files named for this integration are:

- `GSE278952_05_scRR-seq-DNA_HAP1_human_binarized_selectedsamples_all_geneintersect_hg38.txt.gz`
- `GSE278952_05_scRR-seq-DNA_HAP1_human_binarized_selectedsamples_all_sorted_hg38.txt.gz`
- `GSE278952_05_scRR-seq-DNA_HAP1_human_midS_Avg_RT_hg38.bedGraph.gz`

Local inspection status: these exact files were not present under the active workspace or `/tmp` during implementation, so circyto does not claim real-file validation for this importer yet. The parser is covered by synthetic tests that match the expected processed table structure: feature metadata columns followed by `DNA_*` or `RNA_*` sample columns.

## Biological Interpretation

The HAP1 scRR DNA branch should not be treated as CNV by default. The available HAP1 processed filenames indicate binarized replication-state profiles and an average replication-timing bedGraph. This is closer to a replication timing or replication-state modality than a copy-number modality.

Recommended modality name:

```text
rt
```

`rt` means processed scRR DNA replication timing/state. It is intentionally separate from:

- `cnv`: processed IMR90 copy-number state by genomic bin
- `candidate_snv`: optional RNA-derived SComatic candidate signals

## Difference From IMR90 CNV

IMR90 `GSE278958` processed DNA outputs include CNV state tables such as `summary_CNV_states_bin50kb` and mappability-normalized CNV signal tables. Those are imported as:

```text
MuData
|- rna
|- circ
`- cnv
```

HAP1 `GSE278952` processed DNA outputs listed for this task are replication timing/state outputs. They should be imported as:

```text
MuData
|- rna
|- circ
`- rt
```

Do not merge HAP1 replication-state tables into the `cnv` modality unless an upstream file is explicitly documented as CNV.

## AnnData Representation

`circyto import-scrr-rt` writes:

```text
rt.h5ad
rt_cells.tsv
rt_features.tsv
scrr_rt_import_summary.json
```

AnnData layout:

```text
obs = cells
var = genomic bins, gene-intersect features, or other processed RT features
X   = cell x feature replication-state or RT matrix
```

`X` semantics are detected from the processed table:

- binary values such as `0` and `1`: binary replication state
- integer numeric values: integer encoded replication state
- non-integer numeric values: numeric replication timing/state signal
- string categorical values: integer-encoded categories with the mapping stored in summary metadata

Feature metadata:

- coordinate-like tables keep `seqname`, `start`, `end`, and `bin_size`
- gene-intersect tables keep the source feature columns and set `feature_type`
- `feature_id` is unique and used as `var_names`

Average RT bedGraph:

- `--avg-rt-bedgraph` is optional
- if bedGraph coordinates exactly match `rt.h5ad.var[["seqname", "start", "end"]]`, values are stored in `var["avg_rt"]`
- otherwise the bedGraph path and row count are recorded in provenance only
- average RT is not stored as an AnnData layer because it is feature-level, not cell-by-feature

## Cell IDs

The importer reuses the scRR canonical cell convention:

```text
DNA_HAP1_A_001 -> HAP1_A_001
RNA_HAP1_A_001 -> HAP1_A_001
```

Supported obs ID strategies:

- `canonical`: default; strips leading `DNA_` or `RNA_`
- `dna`: uses DNA-style titles
- `rna`: uses RNA-style titles

When a `build-scrr-cell-map` TSV is supplied, exact `dna_cell_id`, `rna_cell_id`, or `sample_title` matches are used to populate paired RNA/DNA/canonical metadata.

## Commands

Import processed RT/state table:

```bash
circyto import-scrr-rt \
  --rt-table GSE278952_05_scRR-seq-DNA_HAP1_human_binarized_selectedsamples_all_sorted_hg38.txt.gz \
  --avg-rt-bedgraph GSE278952_05_scRR-seq-DNA_HAP1_human_midS_Avg_RT_hg38.bedGraph.gz \
  --outdir dna_rt
```

Metadata-only import:

```bash
circyto import-scrr-rt \
  --rt-table processed_rt_table.txt.gz \
  --outdir dna_rt \
  --no-h5ad
```

Merge with remapped RNA/circ MuData:

```bash
circyto merge-scrr-rt \
  --input rna_circ.remapped.h5mu \
  --rt dna_rt/rt.h5ad \
  --output rna_circ_rt.h5mu
```

## Limitations

- Real HAP1 processed RT files were not present locally during implementation, so real-file import validation is pending.
- The importer consumes processed GEO-style tables and does not rerun raw scRepli-seq DNA processing.
- Replication timing/state is not CNV unless the source file explicitly encodes CNV.
- `rt` is not an SNV modality.
- SComatic outputs remain RNA-derived candidate signals unless independently validated by DNA evidence.
