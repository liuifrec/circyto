# Application Note results checklist

## Frozen software evidence

- [x] Baseline is circyto 0.10.0 at
  `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`.
- [x] Full suite recorded as 345 passed, 8 skipped, 5 warnings.
- [x] Wheel, sdist, clean wheel-only install, installed resource lookup, and
  H5MU round trip are recorded as passing.
- [x] MuData synchronization semantics and zero behavior-change warnings are
  documented.
- [ ] Replace any draft environment/version placeholders with the frozen values
  in Supplementary Table S3.

## Main Figure 1 and Smart-seq3 numbers

- [ ] Retrieve the Smart-seq3 object and verify SHA-256
  `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.
- [ ] Regenerate 192 cells, 63,187 RNA features, 2,503 circRNA candidates,
  2,379 annotated candidates, and 95.0% recovery in one report.
- [ ] Regenerate median detected circRNAs/cell = 12 and median total support =
  22.5.
- [ ] Regenerate the fixed-seed RNA-derived UMAP and record package versions.
- [ ] Confirm `chr1:117402186|117420649` is annotated to MAN1A2 and regenerate
  its detection/support overlay.
- [ ] Verify the embedding uses RNA only and both overlays reuse identical
  coordinates.

## Supplementary package

- [ ] Build Figure S1 from current workflow evidence for Smart-seq3, IMR90,
  and the validated HAP1 10-cell route.
- [ ] Build Figure S2 as an object/schema demonstration; label HAP1 RT as an
  optional contract pending reconciled real-file validation.
- [ ] Generate Table S1 dataset/workflow inventory.
- [ ] Generate Table S2 software/reproducibility validation.
- [ ] Generate Table S3 versions/commands/provenance.
- [ ] If long-read evidence is included, keep generic ONT and CIRI-long claims
  separate and preserve their scientific boundaries.

## Quantitative and wording audit

- [ ] Regenerate the IMR90 23-cell modality/overlap row from object SHA-256
  `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.
- [ ] Reconcile the historical HAP1 63-cell/56-cell-RT object with current
  documentation before quoting any HAP1 manuscript-scale count.
- [ ] Confirm every number in the abstract, main text, legends, and supplement
  appears in `application_note_evidence.md` with a source and status.
- [ ] Use “circRNA candidate” or “detected circRNA,” not validated biological
  circRNA, where detector calls have no orthogonal confirmation.
- [ ] State that circyto is not a new back-splice detection algorithm.
- [ ] Avoid CNV/RT biological interpretation, radiation claims, predictive
  biogenesis claims, and single-cell Nanopore circRNA validation claims.

## Deferred biological follow-up

- [ ] HAP1 RT-circRNA regression, IMR90 CNV programs, cross-dataset host-gene
  analyses, and biogenesis modelling remain outside the submission gate.
