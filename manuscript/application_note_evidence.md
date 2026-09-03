# Authoritative evidence for the Bioinformatics Application Note

This file is the bridge between the frozen circyto repository and manuscript
v2/v3. It governs manuscript claims, numbers, figures, and supplementary
material. When another manuscript planning file disagrees with this file or
with the frozen software reports cited here, use this file and resolve the
disagreement before drafting.

## Authority and evidence classes

- Software baseline: circyto 0.10.0 at
  `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`, the merged PR #13/#15 tree.
- **committed/reproducible**: supported by reports, code, or machine-readable
  summaries committed in the baseline history.
- **local processed-data derived**: calculated from a named, checksum-audited
  local object; the large object is not present in the frozen tree.
- **historical/unverified**: appears in prior manuscript artifacts but
  conflicts with, or is not established by, the current authoritative
  workflow evidence.
- **needs rerun**: must be independently regenerated before the final
  manuscript is submitted. This normally means a small read-only pass over a
  checksum-matched processed object, not an expensive raw-data workflow.

Historical committed result files are identified as
`COMMIT:path/to/file`. They remain traceable Git evidence even though generated
objects and result tables were removed at `c3971ec` for release hygiene.

## Central manuscript claim

> circyto bridges circRNA detectors and the modern single-cell/scverse
> ecosystem by converting full-length single-cell sequencing into annotated
> cell-by-circRNA matrices and interoperable AnnData/MuData objects.

Required adjacent sentence:

> circyto is not a new back-splice detection algorithm; it orchestrates
> protocol-aware detector workflows and standardizes cell-resolved outputs for
> downstream analysis.

## Claim-to-evidence matrix

| Claim | Evidence source | Validated? | Placement | Wording boundary |
| --- | --- | --- | --- | --- |
| circyto provides protocol-aware full-length single-cell routes, cell-by-circRNA matrices, QC, annotation, and AnnData/MuData export | `docs/manuscript_software_baseline.md`; `docs/validated_workflows_summary.md`; workflow and export tests | Yes, committed/reproducible | Main paper | Framework/orchestration claim only; not a novel detector |
| The pooled Smart-seq3 / E-MTAB-8735 route ran end to end on real data | `README.md` “Validated benchmark result”; `docs/validated_workflows_summary.md` | Yes, committed/reproducible | Main paper | Call outputs circRNA candidates/detections unless orthogonally validated |
| The audited Smart-seq3 object contains 192 cells, 63,187 RNA features, and 2,503 circRNA candidates | `c99cdda:manuscript/results/manuscript_object_audit.md`; object SHA-256 below; current `README.md` independently records 192 cells and 2,503 candidates | Historically checksum-audited; needs rerun for final text | Main paper | RNA feature count and final object shape come from the processed object, not a current tracked artifact |
| Smart-seq3 host-gene annotation recovered 2,379/2,503 candidates (95.0%) | `c99cdda:manuscript/results/host_gene_annotation_recovery.tsv`; historical object audit | Local processed-data derived; needs rerun | Main paper or legend | Host-gene annotation, not validation of circRNA biology; 95.0% is one-decimal rounding of 0.9504594487 |
| Smart-seq3 has median 12 detected circRNAs and median total support 22.5 per cell | `c99cdda:manuscript/results/manuscript_object_audit.md` | Local processed-data derived; needs rerun | Main paper or Figure 1 legend | Define “detected” as nonzero candidate support; do not call support normalized expression |
| RNA-derived UMAP coordinates can be overlaid with cell circRNA burden | `c99cdda:manuscript/results/smartseq3_umap_cells.tsv`; `c99cdda:scripts/manuscript/export_smartseq3_figure1_data.py`; method record | Historically generated; needs rerun with fixed environment | Main Figure 1B | UMAP is derived only from RNA; circRNA burden is an overlay, not an embedding input |
| A representative MAN1A2-associated candidate can be overlaid on the same UMAP | Historical top-candidate and selected-feature tables at `c99cdda:manuscript/results/`; candidate `chr1:117402186|117420649` | Historically generated; needs rerun | Main Figure 1C | Illustrative detection, not MAN1A2 functional biology or independent circRNA validation |
| circyto outputs operate as scverse AnnData/MuData objects | `docs/manuscript_software_baseline.md`; `docs/mudata_compatibility.md`; H5MU round-trip and multimodal tests | Yes, committed/reproducible | Main paper | Claim interoperability and preserved semantics; do not claim the deferred broad AnnData 0.13 dtype migration |
| IMR90 demonstrates RNA+circ+processed-CNV integration with 23 shared cells | `docs/validated_workflows_summary.md`; `docs/current_project_status.md`; `c3971ec^:load_work/scrr_imr90/full_length_rna_circ_cnv.summary.json` | Yes, committed/reproducible; regenerate supplement row | Brief main text and Supplement | Processed GEO CNV import; no raw-DNA CNV inference and no CNV biological conclusion |
| The HAP1 paired-end scRR RNA/circ route executes on real data | `docs/validated_workflows_summary.md`; `docs/current_project_status.md` | Yes, on a 10-cell batch | Supplement | Protocol validation only; do not substitute unresolved 63-cell/RT claims |
| HAP1 RT import/merge contracts exist | `docs/validated_workflows_summary.md`; RT import/merge tests and documentation | Yes, implemented with synthetic tests | Supplementary schema only | Do not claim completed real-file RT biology or RT-circRNA discovery |
| Generic Nanopore cDNA interoperability is validated on SRR4048177 | `docs/validation/srr4048177_nanopore_interoperability.md`; expectation profile and smoke script | Yes, committed/reproducible | Optional Supplement | Alignment/QC/provenance only; `circRNA_call=false`; no single-cell Nanopore circRNA validation |
| The CIRI-long adapter preserves official-demo outputs | `docs/validation/ciri_long_official_demo.md`; adapter tests and demo runner | Yes, committed/reproducible | Optional Supplement | Bulk RCRT/circRNA-enriched official demo; not single-cell performance or biological accuracy |
| circyto provides biogenesis-ready export schemas | `docs/biogenesis_schema.md`; `circyto/biogenesis/`; tests | Yes, schema/export only | Defer or one future-work sentence | No predictive biogenesis model; positive/unlabelled does not define true negatives |
| MuData behavior-transition compatibility preserves historical circyto synchronization semantics | `docs/mudata_compatibility.md`; synchronization regression tests | Yes, committed/reproducible | Supplementary reproducibility | Scoped `pull_on_update=True`; broader dtype migration remains deferred |
| The frozen package is reproducible/installable | `docs/manuscript_software_baseline.md`; `docs/release_readiness_audit.md` | Yes, committed/reproducible | Main availability sentence and Supplement | Report the frozen gate, not the earlier pre-compatibility warning count |

## Main-paper evidence

The main paper should use only these evidence blocks:

1. **Need and workflow.** Existing detectors do not themselves provide the
   complete cell-resolved matrix, annotation, QC, provenance, and scverse
   object layer supplied by circyto.
2. **Real-data execution.** E-MTAB-8735 Smart-seq3 is the primary demonstration
   of real full-length single-cell data passing through the workflow.
3. **Standard objects.** AnnData/MuData outputs retain modality-specific axes,
   metadata, sparse values, and cell mappings under the explicit
   synchronization policy.
4. **Bounded multimodality.** The validated IMR90 23-cell
   RNA+circ+processed-CNV object demonstrates extensibility without inviting a
   CNV biology story.
5. **Reproducibility.** The final source/package gate, clean wheel-only install,
   installed resource lookup, and H5MU round trip all passed.

## Proposed main Figure 1

- **A — architecture:** full-length scRNA-seq → protocol-aware alignment /
  circRNA detection → cell-by-circRNA matrix + QC → host-gene / known-circRNA
  annotation → AnnData / MuData. Show RNA and circRNA prominently. Show CNV,
  RT, and candidate-variant integrations as optional secondary branches.
- **B — Smart-seq3 burden:** RNA-derived UMAP of 192 E-MTAB-8735 cells,
  coloured by detected circRNA candidates per cell.
- **C — representative output:** the same coordinates overlaid with binary
  detection of `chr1:117402186|117420649`, a MAN1A2-associated candidate.

The figure must not imply that circyto is a detector, that circRNA values were
used to construct the RNA UMAP, or that the MAN1A2 example establishes a
biological association.

## Smallest useful Supplementary package

### Supplementary Figure S1 — protocol/workflow validation

Compare Smart-seq3, IMR90 scRR, and HAP1 scRR by input layout,
alignment/detector route, output type, and validation status. The HAP1 row is
the validated 10-cell paired-end route, not the unresolved historical RT
object.

### Supplementary Figure S2 — multimodal object/schema

Show RNA+circ objects, the validated IMR90 RNA+circ+CNV object, and the
optional RT integration contract. Mark RT as implemented with synthetic tests
and pending reconciled real processed-file validation.

### Supplementary Table S1 — dataset/workflow inventory

Record accessions, protocol/layout, route, reference, modalities, validation
scale, manuscript role, and limitations.

### Supplementary Table S2 — software/reproducibility validation

Record workflow, dataset/fixture, cells, circRNAs, annotation recovery,
validation status, evidence source, and limitation. Leave unresolved cells
blank or label them “not regenerated.”

### Supplementary Table S3 — versions, commands, and provenance

Record circyto/baseline identity, Python/dependency/external-tool versions,
commands or scripts, input accessions/checksums, output provenance, and status.

Long-read evidence is optional. If included, append a compact subsection or
rows to S1-S3 rather than creating another figure or table.

## Quantitative claim register

### Frozen software and packaging

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| Test suite | 345 passed, 8 skipped, 5 warnings | committed/reproducible | `docs/manuscript_software_baseline.md`; four intentional overlapping-feature warnings and one AnnData index-conversion warning |
| MuData behavior-change warnings | 0 | committed/reproducible | `docs/mudata_compatibility.md` |
| Version | 0.10.0 | committed/reproducible | `pyproject.toml`, package version, baseline report |
| Baseline commit | `44697355bcab1c525ca7ef9b130e2ad0094d9e1b` | committed/reproducible | merged remote `main` verified before this manuscript branch |

### Smart-seq3 / E-MTAB-8735

Audited object SHA-256:
`0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| Cells | 192 | committed/reproducible and local processed-data derived | Current `README.md` plus historical checksum-backed audit; independently regenerate for final table/legend |
| RNA features | 63,187 | local processed-data derived; needs rerun | Historical object audit; regenerate from checksum-matched H5MU |
| circRNA candidates | 2,503 | committed/reproducible and local processed-data derived | Current `README.md` plus historical audit; regenerate for final table/legend |
| Nonzero circ matrix entries | 2,659 | committed/reproducible | Current `README.md`; supplementary-only if useful |
| Host-gene annotated candidates | 2,379 | local processed-data derived; needs rerun | Historical recovery table/object audit |
| Host-gene recovery | 95.0% | local processed-data derived; needs rerun | 2,379/2,503 = 95.0459%, rounded to one decimal |
| Median detected circRNAs/cell | 12 | local processed-data derived; needs rerun | Historical object audit |
| Median total circRNA support/cell | 22.5 | local processed-data derived; needs rerun | Historical object audit |
| MAN1A2-associated candidates in historical host aggregate | 3 | local processed-data derived; needs rerun | Historical top-host table; not required in main text |
| Cells detecting representative `chr1:117402186\|117420649` | 6/192 | local processed-data derived; needs rerun | Historical top-candidate/selected-feature exports |
| Total support for representative candidate | 15 | local processed-data derived; needs rerun | Historical export; optional legend value only |

### IMR90 scRR / GSE278958

Audited object SHA-256:
`bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| RNA cells/features | 23 / 63,187 | committed/reproducible | Current workflow docs plus historical machine-readable merge/object audit; regenerate S1/S2 row |
| circRNA cells/candidates | 23 / 2,443 | committed/reproducible | Same sources; 2,443 is the authoritative IMR90 circRNA count |
| CNV cells/bins | 23 / 60,607 | committed/reproducible | Same sources; processed GEO CNV summaries |
| Trimodal overlap | 23 cells | committed/reproducible | `docs/validated_workflows_summary.md` and merge summary |
| Host-gene annotated candidates | 2,429/2,443 (99.4%) | local processed-data derived; needs rerun | Historical object audit; supplementary-only |

### HAP1 scRR / GSE278952

Historical object SHA-256:
`6d3f460371000dde2f307c67887cfb5939ea591527c4e813d8daf2ef2b1bece5`.

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| Currently validated RNA/circ route scale | 10-cell batch | committed/reproducible | `docs/current_project_status.md`; suitable for S1 protocol validation |
| Historical RNA/circ cells | 63 | historical/unverified; needs rerun | Historical object audit conflicts with current “HAP1 full pending” status; do not quote |
| Historical circRNA candidates | 3,209 | historical/unverified; needs rerun | Same conflict; do not quote |
| Historical RT cells/bins | 56 / 56,881 | historical/unverified; needs rerun | Current workflow summary says real processed-file RT validation is pending |
| Historical trimodal overlap | 56 | historical/unverified; needs rerun | Do not use in manuscript until reconciled |
| Historical host-gene annotations | 3,117/3,209 (97.1%) | historical/unverified; needs rerun | Do not use in manuscript until reconciled |
| Prior RT-circRNA correlations/regression | various | historical/unverified; deferred biological follow-up | Exclude from Application Note |

### Optional long-read evidence

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| SRR4048177 input reads | 52,696 | committed/reproducible | Generic ONT report; not circRNA validation |
| Mapped primary queries | 39,537 (75.03%) | committed/reproducible | Generic ONT report |
| Spliced primary queries | 33,270 (84.15% of mapped primaries) | committed/reproducible | Generic ONT report |
| Exploratory circRNA calls | 0 (`circRNA_call=false`) | committed/reproducible | Diagnostic patterns are explicitly not circRNA calls |
| Official-demo BSJs | 149 | committed/reproducible | CIRI-long validation report |
| Full-length isoforms | 149 | committed/reproducible | CIRI-long validation report |
| Expression / isoform-usage rows | 149 / 149 | committed/reproducible | CIRI-long validation report |
| Read assignments | 1,133 | committed/reproducible | CIRI-long validation report |
| Reconciliation discrepancies | 0 | committed/reproducible | All five normalized artifacts matched official inputs |

## Source files and regeneration routes

Current-tree sources:

- `README.md`
- `docs/manuscript_software_baseline.md`
- `docs/validated_workflows_summary.md`
- `docs/release_readiness_audit.md`
- `docs/mudata_compatibility.md`
- `docs/validation/ciri_long_official_demo.md`
- `docs/validation/srr4048177_nanopore_interoperability.md`
- `scripts/manuscript/summarize_mudata_inventory.py`
- `scripts/manuscript/known_novel_circ_summary.py`

Historical checksum/object and Figure 1 sources:

- `c99cdda:manuscript/results/manuscript_object_audit.md`
- `c99cdda:manuscript/results/manuscript_object_manifest.tsv`
- `c99cdda:manuscript/results/host_gene_annotation_recovery.tsv`
- `c99cdda:manuscript/results/smartseq3_umap_cells.tsv`
- `c99cdda:manuscript/results/smartseq3_top_circRNA_candidates.tsv`
- `c99cdda:manuscript/results/smartseq3_top_hostgene_features.tsv`
- `c99cdda:manuscript/results/smartseq3_selected_feature_candidates.tsv`
- `c99cdda:scripts/manuscript/audit_manuscript_objects.py`
- `c99cdda:scripts/manuscript/export_smartseq3_figure1_data.py`
- `c3971ec^:load_work/scrr_imr90/full_length_rna_circ_cnv.summary.json`
- `c3971ec^:load_work/scrr_hap1/full_length_rna_circ_rt.summary.json`

For final regeneration, first verify the archived object's SHA-256, then run
the current inventory script. Reuse or restore the historical fixed-seed
Smart-seq3 export logic in a controlled manuscript environment and record its
full package versions. Do not regenerate from an object with a different hash
without documenting why it changed.

## Scientific boundaries and claims that must not be made

- Do not claim circyto is a new circRNA/back-splice detector.
- Do not describe detector candidates as biologically validated circRNAs
  without independent evidence.
- Do not imply the MAN1A2 overlay establishes a biological mechanism,
  association, or benchmark truth.
- Do not claim a finished HAP1 RT-circRNA discovery analysis.
- Do not interpret the IMR90 processed-CNV example as a demonstrated
  CNV-circRNA mechanism or overinterpret 23 cells.
- Do not describe public HAP1/IMR90 analyses as radiation-exposure results.
- Do not call RNA-derived candidate variant signals validated DNA variants.
- Do not claim generic ONT cDNA validates circRNA detection or absence.
- Do not claim the official bulk CIRI-long demo establishes single-cell
  performance or biological accuracy.
- Do not claim a predictive biogenesis model or treat unlabelled candidates as
  true negatives.
- Do not claim broad MuData 0.4 / AnnData 0.13 dtype migration.

## Remaining pre-submission tasks

1. Retrieve the Smart-seq3 and IMR90 processed objects from the manuscript data
   archive and verify their recorded SHA-256 values.
2. Independently regenerate the Smart-seq3 quantitative row, RNA-only UMAP,
   burden overlay, and MAN1A2-candidate overlay; capture package versions and
   commands in Supplementary Table S3.
3. Regenerate the IMR90 23-cell modality/overlap row from its object.
4. Reconcile or formally retire the historical HAP1 63-cell/RT object. Until
   then, use only the validated 10-cell RNA/circ route and the synthetic RT
   contract.
5. Generate Supplementary Figures S1-S2 and Tables S1-S3; add long-read rows
   only if they fit without creating a second story.
6. Audit every manuscript and legend number against this register. Record
   “not regenerated” rather than filling a gap from old Markdown.
7. Confirm the data/code archive and public accession statements used in
   “Availability and implementation” before submission.

The repository is ready for v3 text drafting from a claim-control perspective.
Final quantitative prose, Figure 1, and the supplementary numeric tables are
not submission-ready until tasks 1-4 are completed.
