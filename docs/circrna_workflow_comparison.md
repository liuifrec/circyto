# User tasks across circRNA workflows

Inspected on **2026-10-15**. This is a documentation and source comparison,
not an executed comparison of speed, sensitivity, accuracy, or user error rates.
The primary comparators are **nf-core/circrna and CirComPara2**.

## Inspected versions and primary sources

| Key | Exact snapshot / source | Material inspected |
| --- | --- | --- |
| C0 | circyto 0.10.0, `44697355bcab1c525ca7ef9b130e2ad0094d9e1b` | `docs/validated_workflows_summary.md`, `docs/full_length_workflow.md`, `docs/scrr_cell_mapping.md`, `docs/mudata_schema.md`, `docs/manuscript_software_baseline.md`; original-run limits in [the provenance audit](../manuscript/original_detection_provenance.md) |
| C1 | circyto development, version string still 0.10.0, `f841c296db8a6fda870eebc5b553eab10095fa0d` | `docs/import_circ.md`, `docs/native_ciri3_import.md`, `docs/pre_release_readiness.md`, import/native/report acceptance records; open draft [PR #19](https://github.com/liuifrec/circyto/pull/19), not a release |
| C2 | External-detector packaging, `6c280426ebf10abcbcaeab72f87c9d34060d6ee8` | `docs/external_detectors.md`, `docs/detector_distribution_validation.md`; separate unmerged distribution evidence |
| N | nf-core/circrna development commit `23fdf1d22698a34e06efa99bf5af68d1bb963b78` (2026-10-08) | Pinned [README](https://github.com/nf-core/circrna/blob/23fdf1d22698a34e06efa99bf5af68d1bb963b78/README.md), [usage](https://github.com/nf-core/circrna/blob/23fdf1d22698a34e06efa99bf5af68d1bb963b78/docs/usage.md), [outputs](https://github.com/nf-core/circrna/blob/23fdf1d22698a34e06efa99bf5af68d1bb963b78/docs/output.md); [Digby et al., 2023](https://doi.org/10.1186/s12859-022-05125-8) |
| P | CirComPara2 **v0.1.2.1**, `303c3621579c1a412493e40772c2ff2660148942` | Pinned [README and parameter/output descriptions](https://github.com/egaffo/circompara2/blob/303c3621579c1a412493e40772c2ff2660148942/README.md), [sample metadata](https://github.com/egaffo/circompara2/blob/303c3621579c1a412493e40772c2ff2660148942/test_circompara/analysis/meta.csv), [project settings](https://github.com/egaffo/circompara2/blob/303c3621579c1a412493e40772c2ff2660148942/test_circompara/analysis/vars.py), `src/sconstructs/main.py`, `src/utils/R/ccp_circrna_expression.R`; [Gaffo et al., 2022; online 2021](https://doi.org/10.1093/bib/bbab418) |

All links above were accessed on 2026-10-15. The GitHub tags API returned no
nf-core/circrna tags, so N is deliberately a commit pin, not an invented
release. Its current usage documents single- and paired-end samples; the 2023
paper's paired-end input description is not used to restrict the later snapshot.
CirComPara2's pinned README includes older Docker example tags; its source pin
here is v0.1.2.1. We did not install those containers. Published comparator
benchmarks are not transferred to circyto or treated as a matched benchmark.

## Task comparison

**D = documented support. ND = not documented in inspected sources.
NA = not assessed.** D does not mean executed here; circyto execution evidence
is separated in the [claim matrix](../manuscript/application_note_evidence.md#claim-to-evidence-matrix).
ND is a limit of this inspection, not proof of impossibility or absence.

| User task | circyto (C0; C1/C2 explicitly marked) | nf-core/circrna (N) | CirComPara2 (P) |
| --- | --- | --- | --- |
| Start from raw RNA reads | **D:** pooled Smart-seq3 R1/R2/I1/I2 demultiplexing; manifest-driven per-cell full-length routes. General adapter/quality trimming is not the central documented route. | **D:** samplesheet/FASTQ input, FastQC and Trim Galore; single/paired layout and strandedness in usage. | **D:** FASTQ + reference, optional Trimmomatic preprocessing, single/paired tests. |
| Detect, quantify and annotate candidates | **D:** established CIRI3, protocol-specific STAR/BWA routes, support matrices and host/known-circRNA annotation. | **D:** multiple BSJ callers, support/tool-count filters, GTF/database annotation and transcript quantification. | **D:** multiple callers, agreement/read thresholds, circRNA/host annotation and linear expression. |
| Manage sample metadata and execution | **D:** manifests, cell/library metadata, cached alignment/detector stages, resume and recorded summaries. | **D:** samplesheet, phenotype table, Nextflow profiles/resume; repeated sample IDs combine sequencing runs. | **D:** `meta.csv` associates files with sample IDs; `vars.py`, SCons task execution, dry run and parallelism. |
| Handle full-length single-cell protocols and demultiplex pooled reads | **D:** Smart-seq3 index-pair mapping; RamDA/scRR single-end route; paired RamDA requires explicit opt-in. Dataset execution scope is bounded. | **ND:** dedicated Smart-seq3 I1/I2 demultiplexing or scRR protocol policy. **D:** already separated samples can be supplied through its samplesheet. | **ND:** dedicated Smart-seq3 I1/I2 demultiplexing or scRR protocol policy. **D:** already separated samples can be supplied through `meta.csv`. |
| Map detector/sample IDs to biological cells | **D:** named cell axis and explicit GEO GSM-to-canonical scRR mapping. C1 native import requires an explicit unique sample/file/cell map. | **D:** custom sample IDs can carry externally established cell identities; spaces are normalized and repeated IDs merge runs. **ND:** scRR GSM-to-biological-cell mapping utility. | **D:** user sample IDs group read files and label results. **ND:** scRR GSM-to-biological-cell mapping utility. |
| Distinguish missing, failed and zero-support observations | **D, C1:** supplied completed/failed/skipped/unprocessed/unknown status; completed zero needs a mapped result; failed-cell support is retained. Historical C0 objects remain unknown without bound logs. | **D:** task failures, retries and execution reports. **ND:** a cell-level exported completion/zero-support contract equivalent to C1. | **D:** stop-on-error and explicit ignore-errors mode, with downstream caveat. **ND:** a cell-level exported completion/zero-support contract equivalent to C1. |
| Export analysis matrices/objects | **D:** Matrix Market, AnnData and optional RNA/circ MuData. | **D:** count matrices, HDF5/TSV quantification and merged `SummarizedExperiment` RDS. **ND:** native AnnData/MuData export. | **D:** reliable circRNA matrices, host linear expression and gene-expression tables. **ND:** native AnnData/MuData export. |
| Integrate RNA from exactly the same cells | **D:** explicit joins/remapping in C0; C1 exact-ID RNA alignment with modality axes retained and conflicts rejected. Genuine native-input + matched-RNA C1 acceptance is pending. | **D:** linear/circular quantification and sample phenotype integration. **ND:** exact-ID integration of an external RNA AnnData into MuData. | **D:** linear/circular and host-expression integration. **ND:** exact-ID integration of an external RNA AnnData into MuData. |
| Inspect QC and reports | **D:** C0 QC tables/summaries and inspection utilities. C1 offline workflow HTML and separate IMPORT HTML; original-workflow HTML acceptance pending. | **D:** FastQC/MultiQC HTML, Nextflow execution reports, timelines, traces and software versions. | **D:** FastQC/read-processing statistics and expression summary HTML. |
| Install and reproduce an analysis | **D:** Python wheel tests; external alignment/detector tools needed for raw workflows. C2 excludes Java binaries from distributions; detector-free import needs no detector runtime. C2 is unmerged. | **D:** Nextflow/configuration profiles, containers and version selection; execution provenance. | **D:** installer, Docker, project settings, reusable indices and test data; pinned README documents Python 2.7 dependencies. |
| Compare speed, accuracy or operator error on matched inputs | **NA:** no comparative execution or user study in this package. | **NA** | **NA** |

Ordinary sample identifiers can preserve cell identity when the input files and
mapping establish it. Neither comparator's use of the word “sample” proves an
inability to process per-cell files. The practical distinction described here
is which mapping, protocol and object-construction steps have documented
contracts, and which steps the user must supply or verify separately.

## Concrete work standardized by circyto

For pooled Smart-seq3, the user supplies the index-to-cell annotation and
reference inputs; circyto checks synchronized read identities, demultiplexes
index pairs, retains selected cell IDs through alignment/detection and collects
the outputs. For scRR, manifest layout and explicit paired-end opt-in choose the
documented route; a GEO-derived map connects library/GSM labels to biological
cells before multimodal merging. AnnData/MuData exports spare users a separate
matrix transpose, axis reconciliation and metadata assembly implementation.
These contracts target common manual mismatches; reduced user error has not
been measured in a comparative study.

C1 additionally requires explicit input semantics and status evidence. It
rejects duplicate identities and incompatible declarations and keeps unknown
processing status distinct from an asserted completed zero. These are tested
secondary extensions, not evidence that the original 192-row object had 192
independently verified successful detector runs.

CIRI3 remains an established detector with quantification and downstream
analysis, not merely a source of isolated calls. Its pinned
[v3.0.1 README](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/README.md)
and [primary paper](https://doi.org/10.1038/s41587-025-02835-1) document matrices
and differential analyses (accessed 2026-10-15). That format pin does not
identify the executable used for the manuscript runs. Orchestration, expression
matrices, host annotation and HTML QC are shared capabilities, not inventions
claimed for circyto.
