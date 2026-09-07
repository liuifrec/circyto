# Application Note data inventory

This inventory distinguishes current committed evidence from values derived
from manuscript-scale objects that are no longer tracked. Exact provenance and
classification are maintained in `application_note_evidence.md`.

## Smart-seq3 / E-MTAB-8735

Role: primary real-data demonstration.

Audited object identity:

```text
load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu
SHA-256: 0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339
```

Checksum-verified regeneration (2026-09-07; `results_v2/smartseq3_summary.tsv`):

- RNA: 192 cells x 63,187 features;
- circRNA: 192 cells x 2,503 candidates;
- shared RNA/circRNA cells: 192;
- annotated circRNAs: 2,379 / 2,503 (95.0% to one decimal place);
- median detected circRNAs per cell: 12;
- median total circRNA support per cell: 22.5;
- representative candidate: `chr1:117402186|117420649`, annotated to
  `MAN1A2`, detected in 6/192 cells with total support 15 in the regenerated
  export.

Status: all manuscript quantities and Figure 1B/C are independently regenerated
and agree with the historical register. Source objects remain outside the
current tree; see `reproduce.md` for recovery and rerun commands.

## IMR90 scRR / GSE278958

Role: bounded multimodal interoperability demonstration.

Audited object identity:

```text
load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu
SHA-256: bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8
```

Committed validation evidence:

- RNA: 23 cells x 63,187 features;
- circRNA: 23 cells x 2,443 candidates;
- processed CNV: 23 cells x 60,607 bins;
- trimodal overlap: 23 cells;
- regenerated host-gene audit: 2,429 / 2,443 (99.4%).

Status: independently regenerated from checksum-matched bytes in
`results_v2/imr90_summary.tsv`; all historical comparisons agree. Do not
interpret processed CNV biology in the main paper.

## HAP1 scRR / GSE278952

Role: protocol validation in Supplementary Figure S1; optional RT contract in
Supplementary Figure S2.

Current authoritative status:

- the paired-end real-data RNA/circRNA route is validated on a 10-cell batch;
- RT import and RNA+circ+RT merge are implemented and covered by synthetic
  tests;
- current workflow documentation says real processed-file RT validation is
  pending.

Historical, unresolved manuscript-object values:

- RNA/circRNA: 63 cells;
- circRNA candidates: 3,209;
- RT: 56 cells x 56,881 bins;
- trimodal overlap: 56 cells;
- host-gene annotations: 3,117 / 3,209 (97.1%).

Status: historical/unverified for manuscript use and needs rerun/reconciliation.
Do not quote these counts or prior RT-circRNA correlations in the Application
Note until the source files, object SHA-256, and current workflow status are
reconciled.

## Long-read interoperability

Role: optional compact supplementary software-extensibility note.

- Generic ONT / SRR4048177: 52,696 input reads; 39,537 mapped primary
  queries (75.03%); 33,270 spliced primary queries (84.15% of mapped primary
  queries). This validates generic cDNA ingestion/alignment/QC only;
  `circRNA_call=false`.
- Official CIRI-long demo: 149 BSJs, 149 full-length isoforms, 149 expression
  rows, 149 usage rows, and 1,133 read assignments, with zero reconciliation
  discrepancies. This is bulk, chemistry-gated adapter interoperability, not
  single-cell biological performance.
