# scRR Tri-Modal MuData

## Goal

Build a single scRR MuData object with:

```text
MuData
|- rna
|- circ
`- cnv
```

The key requirement is that RNA/circ observations and CNV observations must use the same biological-cell namespace.

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

The CNV modality represents processed scRR/scRepli-seq DNA copy-number state by genomic bin. It should be treated as the primary scRR DNA branch. RNA-derived SComatic candidate signals remain a separate optional branch and should not be interpreted as validated DNA somatic mutations.
