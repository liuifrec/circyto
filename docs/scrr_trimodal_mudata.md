# scRR Tri-Modal MuData

## Goal

Build a single scRR MuData object with:

```text
MuData
|- rna
|- circ
`- cnv
```

For HAP1 replication timing/state outputs, the parallel shape is:

```text
MuData
|- rna
|- circ
`- rt
```

The key requirement is that RNA/circ observations and DNA-modality observations must use the same biological-cell namespace.

## Step 1: Remap RNA/Circ MuData Obs IDs

```bash
circyto remap-scrr-mudata-obs \
  --input rna_circ.h5mu \
  --cell-map scrr_cell_map.tsv \
  --output rna_circ.remapped.h5mu
```

Default behavior:

- remaps `obs_names` from GSM IDs to `canonical_cell_id`
- preserves original GSM IDs in `obs["gsm_id"]`
- preserves original modality obs IDs in `obs["original_obs_id"]`
- adds modality-specific original ID columns such as `original_rna_obs_id`

Alternative target:

```bash
circyto remap-scrr-mudata-obs \
  --input rna_circ.h5mu \
  --cell-map scrr_cell_map.tsv \
  --target-id rna_cell_id \
  --output rna_circ.rna-title-ids.h5mu
```

Strictness:

- missing GSM IDs fail by default
- `--allow-partial` leaves unmapped obs IDs unchanged and reports them
- existing outputs require `--overwrite`

## Step 2: Merge CNV

```bash
circyto merge-scrr-cnv \
  --input rna_circ.remapped.h5mu \
  --cnv cnv.h5ad \
  --output rna_circ_cnv.h5mu
```

Default behavior:

- requires exact shared obs IDs across `rna`, `circ`, and `cnv`
- preserves `cnv.X` as integer CNV state
- preserves `cnv.layers["mappabilitynorm"]` when present
- writes a summary JSON beside the output h5mu

Use `--allow-partial` only when intentionally retaining modality-specific cell sets.

## Alternative Step 2: Merge HAP1 RT

Use the RT path for processed HAP1 replication timing/state outputs, not for IMR90 CNV summaries:

```bash
circyto import-scrr-rt \
  --rt-table GSE278952_05_scRR-seq-DNA_HAP1_human_binarized_selectedsamples_all_sorted_hg38.txt.gz \
  --avg-rt-bedgraph GSE278952_05_scRR-seq-DNA_HAP1_human_midS_Avg_RT_hg38.bedGraph.gz \
  --outdir dna_rt

circyto merge-scrr-rt \
  --input rna_circ.remapped.h5mu \
  --rt dna_rt/rt.h5ad \
  --output rna_circ_rt.h5mu
```

Default behavior:

- requires exact shared obs IDs across `rna`, `circ`, and `rt`
- preserves `rt.X` as the processed replication-state or RT matrix
- stores matching average RT bedGraph values as `rt.var["avg_rt"]`
- writes a summary JSON beside the output h5mu

## Merge Summary JSON

Default summary path:

```text
OUTPUT.summary.json
```

Fields:

```text
description
source_rna_circ_h5mu
source_cnv_h5ad
output_h5mu
output_summary_json
allow_partial
n_obs
modalities
modality_shapes
overlap_counts
```

`overlap_counts` includes:

```text
n_rna_obs
n_circ_obs
n_cnv_obs
n_rna_circ_overlap
n_rna_cnv_overlap
n_circ_cnv_overlap
n_trimodal_overlap
n_union_obs
```

For `merge-scrr-rt`, the CNV-specific keys are replaced with `rt` keys such as `n_rt_obs`, `n_rna_rt_overlap`, and `n_circ_rt_overlap`.

## Expected scRR Flow

```bash
circyto build-scrr-cell-map \
  --soft GSE278958_family.soft.gz \
  --out scrr_cell_map.tsv

circyto import-scrr-cnv \
  --cnv-states summary_CNV_states_bin50kb.txt.gz \
  --cnv-mappabilitynorm summary_CNV_mappabilitynorm_bin50kb.txt.gz \
  --outdir dna

circyto remap-scrr-mudata-obs \
  --input mudata/full_length.h5mu \
  --cell-map scrr_cell_map.tsv \
  --output mudata/full_length.scrr_remapped.h5mu

circyto merge-scrr-cnv \
  --input mudata/full_length.scrr_remapped.h5mu \
  --cnv dna/cnv.h5ad \
  --output mudata/full_length.scrr_trimodal.h5mu
```

## Interpretation

The CNV modality represents processed scRR/scRepli-seq DNA copy-number state by genomic bin for datasets whose GEO outputs are CNV summaries, such as the validated IMR90 path.

The `rt` modality represents processed replication timing/state for HAP1-style scRR DNA tables. Do not treat HAP1 binarized RT/state files as CNV by default.

RNA-derived candidate variant signals from SComatic remain a separate optional branch and should not be interpreted as orthogonally confirmed DNA variation.
