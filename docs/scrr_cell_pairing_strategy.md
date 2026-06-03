# scRR Cell Pairing Strategy

## Evidence

The scRR publication and GEO metadata describe paired DNA and RNA data from the same single cell. In `GSE278958`, sample names are listed in paired modality form, for example:

```text
RNA_IMR90_A_100
DNA_IMR90_A_100
```

The inspected processed `GSE278958` files follow the same pattern:

- CNV states header: `DNA_IMR90_A_100`, `DNA_IMR90_A_101`, ...
- RNA TPM header: `RNA_IMR90_A_100`, `RNA_IMR90_A_101`, ...

This supports direct one-to-one pairing by matching the suffix after the first modality prefix.

## Stable Pairing Rule

Define a canonical scRR physical-cell ID:

```text
canonical_cell_id = sample_id with one leading RNA_ or DNA_ prefix removed
```

Examples:

```text
DNA_IMR90_A_100 -> IMR90_A_100
RNA_IMR90_A_100 -> IMR90_A_100
DNA_RPE1_scRR1_MidS_010 -> RPE1_scRR1_MidS_010
RNA_RPE1_scRR1_MidS_010 -> RPE1_scRR1_MidS_010
```

Do not strip internal strings. Only strip one leading modality prefix.

## Validation Required

Before joining modalities, circyto should validate:

- all DNA sample IDs are unique
- all RNA sample IDs are unique
- all canonical IDs are unique within each modality
- each shared canonical ID maps to at most one DNA and one RNA record
- unmatched DNA-only and RNA-only cells are reported, not silently discarded

Column order must not be used as the primary pairing rule. Pair by exact canonical ID.

## Mapping Table

When sample IDs do not follow the `RNA_`/`DNA_` prefix convention, require a TSV mapping table. For GEO SOFT inputs, `circyto build-scrr-cell-map` creates this table from `^SAMPLE` blocks.

Current mapping columns:

```text
gsm_id    rna_cell_id    dna_cell_id    canonical_cell_id    sample_title    molecule    treatment    source_name
```

`canonical_cell_id` can be optional when it is inferable, but explicit values should take precedence.

## circyto Metadata Contract

Each modality should preserve both canonical and original IDs:

```text
obs_names                 canonical cell ID by default
obs["canonical_cell_id"]  modality-prefix-free physical-cell ID
obs["rna_cell_id"]        original RNA sample/cell ID when known
obs["dna_cell_id"]        original DNA sample/cell ID when known
```

For existing RNA/circ workflow outputs whose `obs_names` are not canonical, use `circyto remap-scrr-mudata-obs` before `circyto merge-scrr-cnv`.

## Failure Modes

Reject or warn clearly on:

- duplicate canonical IDs
- mapping entries for DNA IDs absent from the CNV table
- inferred RNA IDs that are not present in a provided RNA table
- mixed naming schemes in one import

## Current Implementation

`circyto import-scrr-cnv` implements the canonical-ID strategy for the CNV AnnData output. It also preserves original DNA IDs and inferred RNA IDs in `cnv.obs`.

`circyto build-scrr-cell-map`, `circyto remap-scrr-mudata-obs`, and `circyto merge-scrr-cnv` implement the GSM-to-biological-cell bridge needed for tri-modal RNA+circ+CNV MuData.
