# circyto Bioinformatics Application Note workspace

Working title:

> circyto: a scverse-compatible framework for circular RNA analysis in full-length single-cell sequencing

Authors: Yu-Chen Liu, Kengo Yoshida, and Yoichiro Kusunoki.

This directory supports a compact Bioinformatics Application Note for
`circyto` 0.10.0. The manuscript's central claim is:

> circyto bridges circRNA detectors and the modern single-cell/scverse
> ecosystem by converting full-length single-cell sequencing into annotated
> cell-by-circRNA matrices and interoperable AnnData/MuData objects.

circyto is not a new back-splice detection algorithm. The primary real-data
demonstration is Smart-seq3 / E-MTAB-8735. IMR90 scRR supplies a bounded
multimodal interoperability example; HAP1, long-read interoperability,
candidate-variant signals, and biogenesis-ready exports are secondary or
future-facing evidence as specified in
[`application_note_evidence.md`](application_note_evidence.md).

## Main-paper scope

The planned paper has one main figure and this structure:

```text
Abstract: Summary; Availability and implementation; Contact; Supplementary information
1 Introduction
2 circyto workflow
3 Single-cell and multimodal outputs
4 Conclusions
Funding
Conflict of Interest
References
Figure 1
```

Figure 1A shows the workflow from full-length single-cell input through
protocol-aware alignment/detection, matrix and QC construction, annotation,
and AnnData/MuData export. Figure 1B overlays circRNA burden on an RNA-derived
Smart-seq3 UMAP. Figure 1C overlays detection of a representative
MAN1A2-associated circRNA. Optional CNV, RT, and candidate-variant integrations
must remain visually secondary.

## Files

- `submission_readiness.md`: verified numbers, assets, remaining submission gates.
- `reproduce.md`: two-command offline result/figure and supplement reproduction.
- `results_v2/`: checksum-verified numbers, Figure 1B/C, and exact provenance.
- `supplement_v2/`: Figure 1A, S1/S2, and Tables S1-S3 with provenance.
- `figure_legends.md`: conservative legends ready for manuscript v3 drafting.
- `application_note_evidence.md`: authoritative claim, number, source, and
  boundary register for manuscript v2/v3.
- `analysis_plan.md`: compact Application Note analysis strategy.
- `figures_plan.md`: Figure 1 and the minimum supplementary figures.
- `tables_plan.md`: minimum supplementary tables.
- `data_inventory.md`: dataset-level evidence and regeneration status.
- `results_checklist.md`: pre-submission evidence gates.
- `methods_commands.md`: command contracts; HAP1 RT association and IMR90 CNV
  biology sections are deferred biological follow-up, not main-paper tasks.
- `caveats.md`: required conservative language.
- `rerf_radiation_positioning.md`: deferred biological follow-up, outside the
  Application Note.

Reusable readers and analysis helpers live under `scripts/manuscript/`;
example command files live under `examples/manuscript/`.

## Reproducibility boundary

Large FASTQ files, references, GEO processed tables, circAtlas tables, and
manuscript-scale H5MU objects are intentionally not in the frozen tree. Earlier
commits contain checksum-backed audits and small derived tables, but the large
objects were removed for release hygiene. The two required objects have now
been recovered from historical Git and independently checksum-verified. Current
small outputs and figures are in `results_v2/`; follow `reproduce.md` to
regenerate them. Confirm public processed-data access and the final availability
statement before journal submission.
