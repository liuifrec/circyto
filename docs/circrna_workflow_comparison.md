# Scope in relation to circRNA workflow software

Primary-source documentation reviewed 2026-10-09. This is a functionality
comparison, not a benchmark of sensitivity, specificity, speed, or biological
accuracy. Versions and datasets differ; no competing software was executed
during this readiness pass. An unlisted capability is **not** evidence that
another package lacks it.

| Software | Documented scope | Relevance to circyto's contribution |
| --- | --- | --- |
| **CIRI3** | BWA or STAR-derived alignment input; single/multiple-sample circRNA detection and quantification; BSJ and FSJ expression matrices; RNase R enrichment information; differential analysis of BSJ expression, junction ratios and relative expression within host genes. [Official documentation](https://github.com/gyjames/CIRI3). | circyto uses CIRI3 as its established detector. Its contribution is the protocol/cell workflow around the detector and interoperable cell-resolved outputs. Describing CIRI3 as only a detector, or claiming matrices as unique to circyto, would be inaccurate. |
| **circtools 2.0** | Modular circular RNA analysis: DCC/STAR detection and host-gene expression, optional CIRIquant integration, QC, CircTest, reconstruction, enrichment, alternative-exon analysis, primer/probe design, and new Nanopore/conservation modules. [Module overview](https://docs.circ.tools/en/latest/), [detection/metatool documentation](https://docs.circ.tools/en/latest/Detect.html). | A broad analysis toolkit with substantial downstream functionality. circyto's intended focus is protocol-aware full-length single-cell orchestration and AnnData/MuData integration. Do not frame QC, downstream analysis, or long-read processing as absent from circtools. |
| **nf-core/circrna** | Nextflow workflow with multiple BSJ callers, support-based call integration, annotation and quantification; the paper describes miRNA targeting and differential analysis. Current output documentation includes MultiQC and execution provenance. [Primary paper](https://doi.org/10.1186/s12859-022-05125-8), [official output documentation](https://nf-co.re/circrna/latest/docs/output). | Reproducible orchestration and reporting are shared concerns. circyto's comparison should emphasize its specific single-cell identity, protocol and multimodal object contracts, without claiming workflow automation or HTML QC is novel. The linked `latest` documentation is a moving development snapshot. |
| **CIRIquant** | CircRNA detection and quantification pipeline; its primary paper describes junction-based requantification, RNase R bias correction and differential analyses. [Official project](https://github.com/bioinfo-biols/CIRIquant), [primary paper](https://doi.org/10.1038/s41467-019-13840-9). | A relevant quantification/analysis comparator. circyto preserves its chosen detector's support semantics; its exported counts are not an implementation of CIRIquant bias correction or molecule counting. |

The circtools documentation currently requests its
[2019 paper](https://doi.org/10.1093/bioinformatics/bty948) and, for 2.0
functionality, the [2025 preprint](https://doi.org/10.1101/2025.02.16.638209).
Use the authors' requested citations rather than treating the original paper
as documentation of every later module.

For circyto, implemented capabilities and validation tiers must stay separate:

| Tier | Defensible scope |
| --- | --- |
| Established full-length workflows | Pooled SMART-Seq3 demultiplexing and CIRI3 orchestration; per-cell RamDA/scRR input; single-end BWA/direct-SAM and paired-end STAR/CIRI3 paths (paired RamDA requires opt-in). Repository evidence records public-data runs. This pass did not rerun them. |
| Cell-resolved downstream representation | Stable cell identities; candidate × cell Matrix Market outputs; cell × candidate AnnData; optional RNA/circ MuData; host-gene provenance and processed scRR CNV integration. See [validated workflow summary](validated_workflows_summary.md) and [object schema](mudata_schema.md). |
| Feature under review | Offline HTML assembled from existing workflow QC, with source checksums and conservative missing/failed-cell interpretation. Synthetic and workflow tests pass; original completed public-workflow report acceptance remains outstanding. |
| Experimental / limited | Heterogeneous detector adapters; processed RT import contract; candidate-SNV interoperability; generic Nanopore alignment/QC; chemistry-gated CIRI-long adapter. These do not establish HAP1 RT biology, biogenesis prediction, ordinary Nanopore circRNA validation, or general long-read support. |

Suggested positioning: **circyto provides protocol-aware, cell-resolved circRNA
analysis and interoperable AnnData/MuData representations across established
full-length single-cell workflows, while retaining the evidence and limitations
of the underlying detectors.** This is a scope statement, not a claim that
all comparator software lacks equivalent individual functions.
