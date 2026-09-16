# Application Note submission readiness

Status: **READY FOR MANUSCRIPT V3 + COAUTHOR REVIEW**.
Finalization audit: 2026-09-08. The numerical evidence is complete; this is not
a claim that the external manuscript has already received coauthor approval.

2026-09-16 bounded interoperability update: **PASS for ordinary MuData APIs
on both actual archived objects**, including native read/write/read and full
semantic checks. See [real_object_interoperability.md](real_object_interoperability.md)
for the evidence and the separately confirmed Smart-seq3 circyto-helper failure.

## Frozen software

circyto **0.10.0**, baseline
`main@44697355bcab1c525ca7ef9b130e2ad0094d9e1b`.
Production `circyto/` and `pyproject.toml` are unchanged. The recorded frozen
gate remains **345 passed, 8 skipped, 5 warnings**, with zero MuData
behavior-change FutureWarnings, successful wheel/sdist builds, clean wheel
installation, packaged resource lookup, and H5MU round trip.

circyto bridges established circRNA detectors and the single-cell/scverse
ecosystem. It is not a new back-splice detection algorithm.

## Verified main-paper evidence

All values below were regenerated from checksum-matched processed objects.
The [comparison table](results_v2/evidence_comparison.tsv) contains **18/18
matching historical comparisons; no numerical discrepancies**.

| Smart-seq3 / E-MTAB-8735 quantity | Final value |
| --- | ---: |
| Cells | 192 |
| RNA features | 63,187 |
| circRNA candidates | 2,503 |
| Nonzero circRNA matrix entries | 2,659 |
| Host-gene annotated candidates | 2,379 |
| Host-gene annotation fraction | 0.950459448661606 (95.0%) |
| Median detected circRNA candidates/cell | 12 |
| Median total circRNA support/cell | 22.5 |
| MAN1A2 candidate detecting cells | 6/192 |
| MAN1A2 candidate total support | 15 |

The illustrative candidate is `chr1:117402186|117420649`, annotated to MAN1A2.
No regulation, function, independent validation, or biological absence claim
follows from this overlay.

| IMR90 scRR / GSE278958 quantity | Final value |
| --- | ---: |
| RNA cells × features | 23 × 63,187 |
| circRNA cells × candidates | 23 × 2,443 |
| CNV cells × bins | 23 × 60,607 |
| Shared RNA/circ/CNV cells | 23 |
| Host-gene annotated candidates | 2,429/2,443 |
| Host-gene annotation fraction | 0.994269340974212 (99.4%) |

IMR90 establishes processed RNA+circ+CNV interoperability only. HAP1
supplementary evidence uses the committed **10-cell RNA/circ route**;
historical full-object/RT numbers are excluded.

## Figure assets and supplement

All figures have SVG originals and 300-dpi PNG previews; conservative text is
provided in [figure_legends.md](figure_legends.md).

| Asset | Source and rendered output |
| --- | --- |
| Figure 1A | Existing `docs/manuscript_figure_skeleton.md` scaffold; rendered by `scripts/manuscript/build_application_note_supplement.py` to [figure1a_architecture.svg](supplement_v2/figure1a_architecture.svg). No prior rendered panel was overwritten. |
| Figure 1B | [smartseq3_umap_cells.tsv](results_v2/smartseq3_umap_cells.tsv) → [circRNA burden SVG](results_v2/figure1b_smartseq3_circrna_burden.svg) |
| Figure 1C | [smartseq3_selected_candidate.tsv](results_v2/smartseq3_selected_candidate.tsv) → [candidate detection SVG](results_v2/figure1c_smartseq3_man1a2_detection.svg); identical coordinates to 1B |
| Figure S1 | [Protocol/workflow validation](supplement_v2/figure_s1_protocol_workflows.svg): Smart-seq3, IMR90, HAP1 batch10 |
| Figure S2 | [Multimodal schema](supplement_v2/figure_s2_multimodal_schema.svg): real RNA+circ/CNV evidence; optional RT explicitly synthetic-tested |
| Table S1 | [Dataset/workflow inventory](supplement_v2/table_s1_dataset_workflow_inventory.tsv) |
| Table S2 | [Software/reproducibility validation](supplement_v2/table_s2_software_reproducibility.tsv) |
| Table S3 | [Versions, commands, checksums, provenance](supplement_v2/table_s3_versions_commands_provenance.tsv) |

No RNA UMAP was stored in the source objects. The historical RNA-only method
was regenerated with Scanpy **1.11.5**, normalization to 10,000, log1p,
2,000 Seurat-flavor variable genes, 50 ARPACK PCs, 15 Euclidean neighbors,
UMAP min_dist 0.5/spread 1.0, spectral initialization, seed **17**, and one
numerical thread. CircRNA values are overlays only. Complete parameters and
the executed script digest are in [results provenance](results_v2/provenance.json).

## Reproducibility and verification

Source archive: Git commit `c99cddaae2314d3af7d5232150c5b86bb6b52e25`.
Both inputs were checksum-verified before/after regeneration and again during
finalization:

```text
Smart-seq3: 0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339
IMR90:      bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8
```

Use the pinned [analysis requirements](requirements-results.txt) with Python
3.10.20. From the repository root, after obtaining the two objects:

```bash
python scripts/manuscript/regenerate_application_note_results.py \
  --smartseq3 <SMARTSEQ3_H5MU> \
  --imr90 <IMR90_H5MU> \
  --output manuscript/results_v2_reproduced
python scripts/manuscript/build_application_note_supplement.py \
  --results manuscript/results_v2_reproduced \
  --output manuscript/supplement_v2_reproduced
```

Replace the bracketed placeholders with local paths; output directories must
be new. [reproduce.md](reproduce.md) gives exact archive-member paths and
recovery commands. Both scripts are offline; source objects are never written.
Source-tree circyto imports are reported truthfully as `0.10.0 (source tree)`.
The regeneration script uses read-only `anndata.io.read_elem` for numerical
extraction from stored modality groups. The subsequent
[real-object audit](real_object_interoperability.md) corrects the earlier
blanket reader limitation: ordinary MuData 0.3.10 reads both originals in this
stack; circyto's explicit-pull helper fails only on Smart-seq3's mismatched
global/modality observation index names. Production code is unchanged. The
historical Scanpy PCA-argument deprecation warning is separate.

The independent `results_v2_recheck` run already establishes determinism:
**all 10 result/table/figure artifacts are byte-identical**. Both runs' script,
input, and output digests verify. The supplement's **nine** output digests and
source digests verify. See [regeneration_verification.json](regeneration_verification.json)
and [recheck provenance](results_v2_recheck/provenance.json). Duplicate recheck
figures and runtime caches remain preserved locally and are excluded from Git;
recheck tables and provenance are retained. No analysis was rerun at finalization.

2026-09-08 focused checks: `python -m compileall -q scripts/manuscript`,
`python -m pytest -q tests/test_manuscript_scripts.py` (**12 passed, 4 warnings**),
and `git diff --check`. The four test warnings concern intentional overlapping
feature names; none is a MuData behavior-change FutureWarning. No raw pipeline
was executed. Large inputs and caches are excluded from the manuscript commits.

2026-09-16 focused checks: `PYTHONPATH=. python -m pytest -q
tests/test_real_object_interoperability.py tests/test_manuscript_scripts.py`
(**25 passed, 4 warnings**) and `git diff --check`. The native real-object
audit separately records **12 upstream FutureWarnings**, without suppression.
The frozen full production suite was not rerun.

## Deferred / nonblocking

HAP1 real RT reconciliation; CRR194209; predictive biogenesis modelling;
broad AnnData 0.13 migration; new detector work. These do not block v3 drafting.

## Remaining submission blockers

- The external manuscript v3 must be drafted/updated, checked against these
  results, and reviewed by coauthors; its document was not supplied in this task.
- Public processed-object access and the final availability statement must be
  confirmed before journal submission. Historical Git recovery is documented;
  a dedicated data deposit or DOI is not yet confirmed.

There are no remaining repository-evidence blockers to manuscript v3 drafting.
