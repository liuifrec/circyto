# scRR DNA Architecture Review

Reviewed: 2026-06-02

## Sources

- Nature Communications article: https://www.nature.com/articles/s41467-025-64688-1
- GEO SuperSeries `GSE278959`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278959
- GEO IMR90 G1 SubSeries `GSE278958`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278958
- scRR-seq code repository: https://github.com/mcbmieu/scRR-seq
- scRepliseq pipeline referenced by the scRR repository: https://github.com/kuzobuta/scRepliseq-Pipeline

Only paper text, GEO metadata, public repository metadata, and small processed GEO supplementary table headers were inspected. No raw sequencing data or hg38 workflows were run.

## What The DNA Modality Represents

The scRR-seq DNA component is a scRepli-seq-derived DNA readout from the same single cell as the RNA component. The publication describes DNA captured on beads, whole-genome amplification, and sequencing as in scRepli-seq. The DNA data are used to assess genome-wide copy number, karyotype G1 cells, and analyze replication state in S-phase cells.

The main DNA signal is therefore a genome-wide binned copy-number or replication-state profile, not a primary SNV-calling assay. The paper explicitly frames scRepli-seq as copy-number based and emphasizes replication timing, S-phase progression, CNV detection, haplotype-resolved analysis, and genome instability.

## GEO DNA Outputs

For `GSE278958` IMR90 G1, GEO exposes processed CNV summary files at multiple bin sizes:

- `summary_CNV_mappabilitynorm_bin50kb`
- `summary_CNV_mappabilitynorm_bin100kb`
- `summary_CNV_mappabilitynorm_bin200kb`
- `summary_CNV_mappabilitynorm_bin500kb`
- `summary_CNV_states_bin50kb`
- `summary_CNV_states_bin100kb`
- `summary_CNV_states_bin200kb`
- `summary_CNV_states_bin500kb`

Observed 50 kb file structure:

- first columns: `seqname`, `start`, `end`
- remaining columns: DNA sample IDs such as `DNA_IMR90_A_100`
- `states` values: categorical copy-number labels such as `0-somy`, `1-somy`, `2-somy`, `3-somy`
- `mappabilitynorm` values: numeric binned DNA signal values

The inspected `GSE278958` 50 kb states table has 60,607 genomic bins and 23 DNA sample columns. The matching RNA TPM table has RNA sample columns with the same suffixes, for example `RNA_IMR90_A_100`.

## Are CNV State Tables Directly Usable?

Yes, for a processed CNV modality. The GEO CNV tables are already binned, cell-indexed matrices with genomic coordinates. They are suitable for import into a `cnv` AnnData modality without rerunning the DNA pipeline.

The importer must preserve:

- source accession and filename
- bin size
- genomic coordinates
- original DNA sample IDs
- copy-number state labels or parsed integer copy numbers
- optional mappability-normalized numeric signal

The importer should not infer raw-read QC, ploidy model details, or new CNV calls unless those are provided by the upstream scRepliseq output.

## DNA/RNA Cell Linkage

`GSE278958` lists paired RNA and DNA samples with matching non-modality suffixes:

- `RNA_IMR90_A_100`
- `DNA_IMR90_A_100`

The processed CNV table headers and RNA TPM table headers follow the same pattern for the inspected IMR90 G1 samples. This supports deterministic cell pairing by stripping one leading modality prefix, then matching the remaining suffix.

This is evidence for direct cell-level linkage in the inspected GEO sample set. circyto should still validate one-to-one uniqueness rather than relying on column order.

## Resolution

For `GSE278958`, GEO provides processed CNV summaries at 50 kb, 100 kb, 200 kb, and 500 kb. The paper also reports scRR-seq DNA replication profiles at 40 kb or 80 kb resolution for S-phase analysis. These are related but not identical public-output resolutions.

Recommendation:

- import the bin size explicitly from the file or command metadata
- default examples to the processed GEO 50 kb table when a user wants maximum public-output resolution
- allow 100 kb, 200 kb, and 500 kb files as lower-memory alternatives

## Appropriate Data Structures

The native structure is a cell-by-genomic-bin matrix:

- observations: single cells
- variables: genomic bins with `seqname`, `start`, `end`, `bin_size`
- `X`: integer copy number parsed from `N-somy`
- optional layer: mappability-normalized numeric signal

This maps directly to AnnData and then MuData as `mdata.mod["cnv"]`.

## Limitations

- GEO processed files are outputs of an upstream DNA pipeline; they are not raw DNA reads.
- The CNV state labels should be treated as upstream calls, not recomputed circyto calls.
- Genome build should be taken from the GEO/pipeline metadata for the specific dataset; it should not be inferred from filename alone.
- S-phase replication-state interpretation requires phase-aware metadata. G1 CNV tables should not be reinterpreted as replication timing profiles.
- SNV calling is not the primary scRR DNA modality based on the publication and inspected public outputs.

## circyto Integration Implication

The evidence supports a long-term scRR branch centered on:

```text
MuData
|- rna
|- circ
|- cnv
`- candidate_snv  (future, optional, RNA-derived or externally validated)
```

The immediate DNA integration target should be CNV/copy-number state import, not SComatic as the core DNA branch.
