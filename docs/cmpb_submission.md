# CMPB manuscript submission: scope and reproducibility

**Status:** Submitted on **20 November 2026** to *Computer Methods and Programs in Biomedicine* (CMPB) as a **Full Length Article**. This is a submission record, **not** a publication, acceptance notice, or peer-reviewed endorsement. A journal DOI is not available.

**Title:** circyto for integrated single-cell circRNA analysis in radiation research  
**Authors:** Yu-Chen Liu, Kengo Yoshida, Yoichiro Kusunoki  
**Research setting:** Radiation Effects Research Foundation (RERF), research protocol RP1-23-1  
**Software:** [circyto](https://github.com/liuifrec/circyto), a framework connecting junction-level circRNA evidence to cell-resolved gene expression and RNA-defined states.

## Submitted analysis at a glance

The manuscript assesses circyto's integration of CIRI3 back-splice junction (BSJ) read-support evidence and independently measured gene expression in AnnData/MuData objects. It retains cell identity, genomic coordinates, strand, support counts, and processing status, including completed libraries with zero detected junctions.

| Public dataset | Biological cells | Candidate BSJs in the manuscript analysis |
| --- | ---: | ---: |
| IMR90, scRepli-RamDA-seq (GSE278958) | 23 | 2,506 |
| HEK293T, Smart-seq3 (E-MTAB-8735) | 192 | 2,991 |
| KYSE180, Smart-seq2 (GSE81812; primary view) | 222 | 3,513 |
| **Distinct cells** | **437** | Dataset-specific candidate sets; **do not sum** |

In KYSE180, **181 irradiated cells** (92 cumulative 12-Gy and 89 cumulative 30-Gy) were analyzed for RNA-group and candidate associations. The three primary RNA groups contain **83, 60, and 38 cells**. Four junctions (three **MAN1A2** and one **SLC8A1**) were prioritized for future experimental investigation. These are **exploratory associations** and literature-linked hypotheses, not experimentally established circRNA functions in irradiated KYSE180 cells.

Software verification used **705 valid input configurations** and **10 invalid configurations** against an independent reference implementation. The largest **synthetic post-detection** workload contained 10,000 cells, 50,000 circRNA features and 1,000,000 nonzero entries (median approximately 47.9 s and 960 MiB peak RSS for the measured scope). These figures do **not** benchmark alignment or CIRI3 detection.

## What this study does *not* establish

- circyto is **not** a new BSJ detection algorithm; the study builds on CIRI3.
- Read support is not a Smart-seq3 UMI-deduplicated molecule count.
- Sequencing the 12-Gy libraries twice does not create independent biological replicates; the first-block comparison is a same-cell depth sensitivity analysis.
- Three cross-gene loci remain mapping/structure concerns rather than validated intergenic circRNAs.
- Association coefficients are exploratory adjusted estimates, not experimental proof of radiation-response mechanisms.
- Additional CNV, replication timing, SComatic, long-read, and explorer capabilities in this repository are **not** evidence that those modalities were evaluated in this CMPB manuscript.

## Software versus paper version

The current `main` source tree and its **0.10.0** package version have evolved beyond older tagged releases (**latest published GitHub release v0.8.3**). Some repository demonstrations and historical workflow benchmark counts differ from manuscript-specific runs. For example, an earlier Smart-seq3 workflow reports **2,503 circRNA features**, whereas the submitted analysis identifies **2,991 candidate BSJs**. They must remain separate provenance records.

For exact reproduction, record the commit, source data, read-depth/threshold choices, software environment, and manuscript analysis snapshot actually used. **Do not assume that a current checkout of `main` reproduces the submitted figures automatically.** Older documents in [`manuscript/`](../manuscript/) and [the original manuscript reproducibility plan](manuscript_reproducibility.md) reflect an earlier broader genome-state study concept, not the submitted CMPB experiment.

## Availability and citation

The three source sequencing datasets are public under the accessions above. The repository makes source code and documented examples available. No independent durable processed-data deposit or journal DOI is asserted here; access to manuscript-specific processed artifacts must be verified separately.

Until a journal article is published, cite **Yu-Chen Liu, circyto, GitHub repository**, with the precise software version and commit as indicated by [`CITATION.cff`](../CITATION.cff); cite upstream circRNA detectors separately. The submitted manuscript must not be described as accepted or published.

**Last documentation review:** 26 November 2026.
