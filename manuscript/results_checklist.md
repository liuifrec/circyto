# Application Note results checklist

## Frozen software evidence

- [x] Baseline is circyto 0.10.0 at
  `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`.
- [x] Full suite recorded as 345 passed, 8 skipped, 5 warnings.
- [x] Wheel, sdist, clean wheel-only install, installed resource lookup, and
  H5MU round trip are recorded as passing.
- [x] MuData synchronization semantics and zero behavior-change warnings are
  documented.
- [x] Replace draft environment/version placeholders with the frozen values
  in Supplementary Table S3.

## Main Figure 1 and Smart-seq3 numbers

- [x] Retrieve the Smart-seq3 object and verify SHA-256
  `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.
- [x] Regenerate 192 cells, 63,187 RNA features, 2,503 circRNA candidates,
  2,379 annotated candidates, and 95.0% recovery in one report.
- [x] Regenerate median detected circRNAs/cell = 12 and median total support =
  22.5.
- [x] Regenerate the fixed-seed RNA-derived UMAP and record package versions.
- [x] Confirm `chr1:117402186|117420649` is annotated to MAN1A2 and regenerate
  its detection/support overlay.
- [x] Verify the embedding uses RNA only and both overlays reuse identical
  coordinates.

## Supplementary package

- [x] Build Figure S1 from current workflow evidence for Smart-seq3, IMR90,
  and the validated HAP1 10-cell route.
- [x] Build Figure S2 as an object/schema demonstration; label HAP1 RT as an
  optional contract pending reconciled real-file validation.
- [x] Generate Table S1 dataset/workflow inventory.
- [x] Generate Table S2 software/reproducibility validation.
- [x] Generate Table S3 versions/commands/provenance.
- [x] Keep optional long-read evidence outside the minimal supplement.

## Quantitative and wording audit

- [x] Regenerate the IMR90 23-cell modality/overlap row from object SHA-256
  `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.
- [x] Formally exclude historical HAP1 full-object/RT counts from this paper;
  use only the validated 10-cell route and synthetic RT contract.
- [ ] Confirm every number in the abstract, main text, legends, and supplement
  appears in `application_note_evidence.md` with a source and status.
- [x] Use “circRNA candidate” or “detected circRNA,” not validated biological
  circRNA, where detector calls have no orthogonal confirmation.
- [x] State that circyto is not a new back-splice detection algorithm.
- [x] Avoid CNV/RT biological interpretation, radiation claims, predictive
  biogenesis claims, and single-cell Nanopore circRNA validation claims.

## Deferred biological follow-up

- [x] HAP1 RT-circRNA regression, IMR90 CNV programs, cross-dataset host-gene
  analyses, and biogenesis modelling remain outside the submission gate.

## Remaining editorial/submission gates

- [ ] Draft and coauthor-review the external manuscript v3; audit its final
  abstract/main-text numbers (the external draft was not supplied).
- [ ] Confirm public processed-object access and final availability statement.

Repository evidence, legends, and assets are ready for v3 drafting. See
`submission_readiness.md` for the verified result and file inventory.
