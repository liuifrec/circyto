# Human RamDA / Shin-RamDA Candidate Datasets

This note ranks human full-length total-RNA datasets that are plausible next-step validation targets for `circyto` on `hg38`.

The immediate purpose is pragmatic:

- pick a biologically correct human dataset after the mouse `GSE98664` engineering smoke test
- keep the first `hg38` validation small enough to debug
- prefer public sources with enough metadata to infer layout and reference compatibility

## Recommendation hierarchy

### BEST_FIRST_HG38_TARGET

`GSE278952`

- dataset_id: `GSE278952`
- protocol: `scRepli-RamDA-seq`
- species: `Homo sapiens`
- cell type/system: `HAP1 mid-S single cells`
- public source: GEO / SRA
- likely read layout: likely paired-end
- expected compatibility with CIRI3: moderate
- pros: clearly human; explicitly analyzed against `hg38`/human GENCODE in the associated GEO/sample metadata; public raw data are available in SRA; focused human cell-line system
- caveats: this is a DNA+RNA multi-omics dataset rather than plain single-cell RamDA; lightweight lookup surfaced a human sample/SRA experiment accession (`GSM8558695`, `SRX26315019`) but did not cleanly expose a matched RNA-side run accession in this session; likely paired-end rather than the single-end RamDA route already validated locally
- suitable for first hg38 validation: yes, if the first goal is true human full-length total-RNA compatibility rather than strict single-end continuity

### SECONDARY_TARGET

`GSE278958`

- dataset_id: `GSE278958`
- protocol: `scRepli-RamDA-seq`
- species: `Homo sapiens`
- cell type/system: `IMR-90 G1 single cells`
- public source: GEO / SRA
- likely read layout: likely paired-end
- expected compatibility with CIRI3: moderate
- pros: clearly human; associated paper and GEO metadata reference `hg38` and human GENCODE-based analysis; fibroblast system may be useful for later CNV/circRNA interpretation
- caveats: same scRR multi-omics complexity as `GSE278952`; accessions were easier to confirm at dataset/series level than at clean RNA run level in a lightweight lookup
- suitable for first hg38 validation: yes, but slightly less convenient than `GSE278952` because the HAP1 system is a tighter first-pass benchmark

### FUTURE_LARGE_SCALE_TARGET

`shin-ramda-riken` workflow resources

- dataset_id: `shin-ramda-riken-HEK293T_benchmark`
- protocol: `shin-ramda`
- species: `Homo sapiens`
- cell type/system: `HEK293T benchmark workflow`
- public source: GitHub workflow repository
- likely read layout: likely single-end for primary Shin-RamDA RNA use, but accession-level verification is still needed
- expected compatibility with CIRI3: potentially high once accession-level FASTQs are confirmed
- pros: directly relevant to Shin-RamDA; the repository explicitly contains `HEK293T_benchmark`, `Detection_of_genomic_breakpoints`, and `Lineage_estimation_using_etoposide_K562`; human cell systems match the next biological validation goal
- caveats: the GitHub repository says accession numbers are described in the manuscript rather than directly listing them in the repo README; lightweight lookup did not recover a clean public FASTQ run list in this session
- suitable for first hg38 validation: not yet; better treated as a second-phase target once accession-level FASTQs are pinned down

## Candidate table

### `GSE278952`

- dataset_id: `GSE278952`
- protocol: `scRepli-RamDA-seq`
- species: `Homo sapiens`
- cell type/system: `HAP1 mid-S single cells`
- public source: GEO / SRA
- likely read layout: paired-end is likely
- expected compatibility with CIRI3: workable, but probably through the paired-end route rather than the already-validated single-end BWA path
- pros: human; public; associated metadata explicitly references `hg38`; manageable focused cell-line system
- caveats: RNA and DNA are paired in one multi-omics workflow; not the simplest plain RamDA example
- whether suitable for first hg38 validation: yes

### `GSE278958`

- dataset_id: `GSE278958`
- protocol: `scRepli-RamDA-seq`
- species: `Homo sapiens`
- cell type/system: `IMR-90 G1 single cells`
- public source: GEO / SRA
- likely read layout: paired-end is likely
- expected compatibility with CIRI3: workable for human total-RNA validation, but not ideal for preserving the exact single-end route already debugged
- pros: human; public; biologically distinct fibroblast system
- caveats: likely larger and somewhat less convenient than `GSE278952` for an initial small validation
- whether suitable for first hg38 validation: yes, secondary

### `shin-ramda-riken-HEK293T_benchmark`

- dataset_id: `shin-ramda-riken-HEK293T_benchmark`
- protocol: `shin-ramda`
- species: `Homo sapiens`
- cell type/system: `HEK293T`
- public source: GitHub workflow repository plus manuscript
- likely read layout: likely single-end or protocol-specific; verify accession metadata before use
- expected compatibility with CIRI3: promising, especially if accession-level FASTQs prove to be ordinary stranded total-RNA libraries
- pros: directly aligned with the Shin-RamDA method family; human; benchmark framing is attractive for `circyto`
- caveats: public accession discovery remains incomplete from lightweight lookup
- whether suitable for first hg38 validation: not until accession-level runs are confirmed

### `shin-ramda-riken-etoposide-K562-lineage`

- dataset_id: `shin-ramda-riken-etoposide-K562-lineage`
- protocol: `shin-ramda`
- species: `Homo sapiens`
- cell type/system: `etoposide-treated K562`
- public source: GitHub workflow repository plus manuscript
- likely read layout: likely single-end or protocol-specific; verify accession metadata before use
- expected compatibility with CIRI3: potentially useful later, especially for integration with genome instability and breakpoint analysis
- pros: human K562 system is attractive for later SComatic-style follow-up
- caveats: more biologically complex than a first validation target; accession-level FASTQ confirmation still missing
- whether suitable for first hg38 validation: no, better as a future integration target

## Small curated accession note

Lightweight lookup in this session confirmed these public human accessions:

- `GSE278952`
- `GSE278958`
- `GSM8558695`
- `SRX26315019`

These are sufficient to justify human scRepli-RamDA follow-up planning, but they do not yet provide a clean, confirmed list of RNA-side `SRR` run accessions for a minimal first pull.

## Conclusion

For the first true human `hg38` RamDA-family validation, use `GSE278952` first.

Reasoning:

- human
- public
- explicit `hg38` compatibility in associated metadata
- focused cell-line system
- less ambiguous than the Shin-RamDA GitHub workflows, where the repo points back to the manuscript for accession numbers

For later SComatic integration exploration, prefer the human Shin-RamDA / K562 direction once accession-level FASTQs are pinned down, with `shin-ramda-riken-etoposide-K562-lineage` as the most interesting future target.

## Source notes

- The `GSE98664` GEO record is clearly mouse and should remain an engineering fixture rather than a human validation target.
- The `shin-ramda-seq-paper` repository explicitly documents `HEK293T_benchmark`, `Detection_of_genomic_breakpoints`, and `Lineage_estimation_using_etoposide_K562`, but says accession numbers are described in the manuscript.
- The scRepli-RamDA-seq paper and related GEO/SRA sample metadata indicate human `hg38` and human GENCODE-based analysis for the human subseries.
