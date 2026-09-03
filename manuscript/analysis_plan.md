# Application Note analysis plan

## Framing

The paper evaluates circyto as the missing workflow and data-model layer
between established circRNA detectors and scverse analysis. It does not
present a new detector or a circRNA biology survey.

The main paper must establish, in order:

1. full-length single-cell protocols require detector-aware preprocessing and
   cell-resolved collection;
2. circyto produces annotated cell-by-circRNA matrices with QC and provenance;
3. real Smart-seq3 data run through the workflow;
4. the resulting RNA and circRNA data behave as standard AnnData/MuData
   modalities;
5. processed DNA-state data can be added without changing modality semantics;
6. the frozen package is installable and reproducibly validated.

## Primary demonstration: Smart-seq3 / E-MTAB-8735

Use a checksum-matched copy of the audited Smart-seq3 RNA+circ H5MU object.
Regenerate all final numbers and plot data from that object in one controlled
run. The intended main-paper summaries are:

- 192 cells;
- 63,187 RNA features;
- 2,503 circRNA candidates;
- 2,379 host-gene annotations (95.0%);
- median 12 detected circRNAs per cell;
- median total circRNA support 22.5 per cell.

These values agree with the historical committed object audit, but the object
is not present in the frozen tree. They require pre-submission regeneration,
not a full raw-data workflow rerun, unless the archived object fails its
recorded SHA-256 check.

Figure analysis:

- compute the embedding from the RNA modality only with fixed preprocessing
  and random seed;
- overlay per-cell circRNA burden without using circRNA values to construct the
  embedding;
- overlay detection (or support, if clearly labelled) for the representative
  `chr1:117402186|117420649` MAN1A2-associated candidate;
- describe this as an illustrative software output, not a biological
  association or validation of MAN1A2 function.

## Secondary interoperability demonstration: IMR90 scRR

Use the validated 23-cell RNA+circ+processed-CNV MuData example to show that
modalities share an explicit cell axis and retain modality-specific features.
The committed evidence reports 23 cells in each modality and a trimodal
overlap of 23. CNV is imported from processed GEO summaries; circyto does not
infer CNV from raw DNA reads. Do not include CNV-burden correlations,
CNV-high host-gene programs, or locus-level claims in the main paper.

## Protocol breadth: HAP1 scRR

Use HAP1 only as evidence that the paired-end scRR RNA/circRNA route executes
on real data. The current authoritative workflow summary supports a validated
10-cell batch. Historical manuscript artifacts report a 63-cell RNA/circ
object and a 56-cell RT overlap, but current frozen-tree documentation says
real processed-file RT validation remains pending. Those historical counts and
all RT-circRNA association results therefore require reconciliation and are
not manuscript claims.

## Supplementary analyses

Keep the supplement limited to protocol/workflow validation, object schemas,
dataset inventory, software validation, versions, commands, and provenance.
Generic Nanopore and CIRI-long results may appear as a compact optional
interoperability note only if their chemistry and interpretation boundaries
remain explicit.

## Deferred biological follow-up

The following are not requirements for this Application Note:

- HAP1 RT-circRNA regression or discovery claims;
- IMR90 CNV-burden biology or CNV-high host-gene programs;
- local CNV-at-circRNA-locus inference;
- cross-dataset host-gene programs;
- radiation/exposure interpretation;
- predictive circRNA biogenesis modelling;
- CRR194209 or single-cell CIRI-long validation.
