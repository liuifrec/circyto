# scRR Cell Mapping

## Purpose

scRR-seq GEO metadata uses GSM accessions for sequencing samples, while the biological cell identity is encoded in `Sample_title`.

Example:

```text
GSM8558852 -> RNA_IMR90_A_100 -> IMR90_A_100
```

The matching DNA sample is:

```text
DNA_IMR90_A_100 -> IMR90_A_100
```

`circyto build-scrr-cell-map` creates the metadata bridge needed to align RNA/circ MuData objects with CNV AnnData objects.

## Command

```bash
circyto build-scrr-cell-map \
  --soft GSE278958_family.soft.gz \
  --out scrr_cell_map.tsv
```

Inputs:

- GEO `family.soft` or `family.soft.gz`

Output TSV columns:

```text
gsm_id
rna_cell_id
dna_cell_id
canonical_cell_id
sample_title
molecule
treatment
source_name
```

## Parsing Rules

- parse `^SAMPLE` blocks
- read `!Sample_geo_accession`
- read `!Sample_title`
- identify `RNA_*` and `DNA_*` titles
- `canonical_cell_id` is `Sample_title` with one leading `RNA_` or `DNA_` prefix removed
- pair RNA and DNA records by `canonical_cell_id`
- preserve RNA-only and DNA-only records

Validation:

- duplicate GSM IDs fail
- duplicate `Sample_title` values fail
- duplicate RNA canonical IDs fail
- duplicate DNA canonical IDs fail

RNA/DNA pairs intentionally share the same canonical ID and are not treated as duplicates.

## Example Row

```text
gsm_id       GSM8558852
rna_cell_id  RNA_IMR90_A_100
dna_cell_id  DNA_IMR90_A_100
canonical_cell_id  IMR90_A_100
sample_title RNA_IMR90_A_100
molecule     RNA
treatment    none
source_name  IMR90 cell
```

## Downstream Use

The map is consumed by:

- `circyto remap-scrr-mudata-obs`
- future scRR metadata/QC join commands

It is also accepted as a lightweight audit artifact because it makes the GSM-to-biological-cell mapping explicit.
