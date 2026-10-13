# Detector-free interoperability assessment

Assessment base: `17867f18f4ee806136ce465415ed075e732b9f7b`,
`feature/qc-report-usability`. This change is outside the frozen manuscript
baseline; no detector algorithm, threshold or manuscript result is changed.

## Actual gap

Detector-free collection and export already existed. `collect-matrix` reads
saved detector outputs, `export-multimodal` attaches circ counts to an existing
RNA AnnData, and `export-mudata` reads completed-workflow artifacts without
alignments. Workflow helpers already export circ AnnData and MuData.
`merge-detectors` and `compare-detectors` supply aggregate comparisons, not a
cell-level import interface. Thus “detector-free export” alone would duplicate
existing functionality and is not a novelty claim.

The missing public contract connected those pieces to external cell-resolved
results while retaining processing status and declaring coordinate/count
semantics. Existing routes could filter zero cells, require workflow-specific
RNA tables, cast counts to int32, coerce string metadata to numbers, or omit
observation state. The previous HTML reporter appropriately withholds rates
when its aggregate QC cohort cannot establish completion, but cannot calculate
a completed-cell subset numerator from aggregate candidate totals.

The single improvement is `import-circ`: a strict input boundary plus
observation-aware AnnData/MuData export. It reuses collector outputs, host-gene
schema normalization, MuData construction/write policy and offline report
styling. It adds a public long-table/matrix contract, conservative validation,
explicit observation states, subset denominators, source retention/checksums,
and a compact import report. It does not add new detector parsing, calling,
ensemble consensus, normalization, reference annotation or statistical testing.
The new importer uses native AnnData construction to preserve int64 support and
literal metadata, avoiding legacy export helpers' int32/numeric coercions.

## Scientific utility and comparison

Separating failed/unprocessed/unknown cells from completed zero-support cells
prevents processing failures from inflating the denominator and depressing
observed detection rates. Retaining original support, coordinate declarations
and exact cell identities enables reviewable reuse alongside matched RNA
without conflating count units or silently recoding coordinates. Completion is
an observation-process label, not a model of locus-level callability.

- [CIRI3](https://github.com/gyjames/CIRI3) detects and quantifies circRNAs,
  provides BSJ/FSJ matrices and multilevel differential analysis. circyto's
  contribution here is the cell-resolved object/observation contract downstream
  of such results, not a replacement detector or a new quantification method.
- [circtools 2.0 detect](https://docs.circ.tools/en/latest/Detect.html) provides
  detection/metatool integration; the [2.0 paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11870457/)
  also describes probe design, conservation and long-read support. This change
  addresses faithful cell-state and scverse interchange, not those algorithms.
- [nf-core/circrna outputs](https://nf-co.re/circrna/latest/docs/output) include
  combined sample-specific BSJ calls, quantification and workflow reports.
  A portable detector workflow and a validated cell-observation object boundary
  serve different tasks; the importer complements saved outputs after a
  reviewed conversion, without claiming native nf-core compatibility.

These comparisons are scoped to the documented interfaces checked on
2026-10-13, not assertions that no comparable capability exists elsewhere.
The feature strengthens the software's reuse and reproducibility argument.
It is insufficient by itself to establish methodological novelty or justify
changing the frozen manuscript's claims. See the [validation record](import_circ_validation.md)
for empirical evidence and remaining external acceptance gates.
