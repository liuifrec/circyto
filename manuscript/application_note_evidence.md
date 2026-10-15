# Authoritative evidence for the Bioinformatics Application Note

This file is the bridge between the frozen circyto repository and manuscript
v2/v3. It governs manuscript claims, numbers, figures, and supplementary
material. When another manuscript planning file disagrees with this file or
with the frozen software reports cited here, use this file and resolve the
disagreement before drafting.

2026-10-15 Nikaido response: [response register and Introduction draft](nikaido_response.md),
[versioned task comparison](../docs/circrna_workflow_comparison.md), and
[original-detection provenance](original_detection_provenance.md). The detailed
claim matrix below supersedes the former broad feature table. Historical raw-run
versions/settings are not established by the 0.10.0 software baseline or the
processed-object regeneration. The current RERF submission, BioRender Figure 1A,
numerical evidence and approved Conclusions remain unchanged.

Regeneration completed on 2026-09-07: both expected H5MU hashes match and all
18 quantitative historical comparisons agree. Authoritative current outputs
are [`results_v2/regeneration_summary.json`](results_v2/regeneration_summary.json),
[`results_v2/evidence_comparison.tsv`](results_v2/evidence_comparison.tsv), and
[`results_v2/provenance.json`](results_v2/provenance.json). See
[`reproduce.md`](reproduce.md) for the offline entry point and the documented
legacy-object reader/embedding decisions. S1/S2 and Tables S1-S3 are in
[`supplement_v2/`](supplement_v2/).

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
- **checksum-verified regeneration**: independently recalculated from the
  expected processed-object bytes; small tables and figures are committed,
  with input/output digests and the execution environment recorded.
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

Updated 2026-10-15 for Nikaido comments 1–3. This is the single detailed matrix;
the editable [TSV](claim_evidence.tsv) contains full source identities. “Pass”
describes the cited historical test, not an execution repeated during this
documentation revision. The numerical registers below remain unchanged.

Reference keys in the readable table:

- C0 = `44697355bcab1c525ca7ef9b130e2ad0094d9e1b`, frozen 0.10.0 implementation.
- M0 = `c29831c3f6e7a31b708cce493c0be8036176077a`, manuscript evidence snapshot.
- O0 = `c99cddaae2314d3af7d5232150c5b86bb6b52e25`, historical object archive.
- C1 = `f841c296db8a6fda870eebc5b553eab10095fa0d`, development/PR #19.
- C2 = `6c280426ebf10abcbcaeab72f87c9d34060d6ee8`, external-detector packaging.
- SS_SHA256 = `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339`.
- IM_SHA256 = `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.

Files on C1/C2 are evidence read with `git show <full-ref>:<path>`; their software
was not merged into this manuscript branch. Original objects record 0.9.0 at
export/refresh; C0 is the later validated implementation baseline, not a claim
that 0.10.0 originally generated the alignments. See
[original detection provenance](original_detection_provenance.md).

| Feature | User-facing claim | Implementation/version | Input/provenance | Test class | Operation tested | Result | Limitation | Placement | Evidence record |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C01 Smart-seq3 workflow | Standardizes pooled index-to-cell processing and CIRI3 output collection | Available in 0.10.0@C0; original executing SHA unknown; archived export/refresh says 0.9.0 | E-MTAB-8735 diySpike; O0:load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu SHA256=SS_SHA256 | original workflow execution (reported); embedded original summary | R1/R2/I1/I2 -> demux -> STAR/BWA rescue/CIRI3 -> matrix/QC | Recorded 192 selected/aligned; 191 success + 1 empty; final 192 x 2503 object | Exact raw run IDs, tool versions and per-cell completion binding missing; DIYHEK_192 source membership RNA-only | Main workflow; detection details supplement | C0:docs/validated_workflows_summary.md; O0:load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu /uns/circyto/workflow_summary_json |
| C02 IMR90 workflow and cell mapping | Connects per-cell RNA/circ outputs to canonical scRR cells and processed CNV | Available in 0.10.0@C0; original executing SHA unknown; export/refresh says 0.9.0 | GSE278958; O0:load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu SHA256=IM_SHA256 | original workflow execution (reported); stored mapping | Single-end BWA/direct-SAM CIRI3; GSM remap; processed CNV merge | 23 RNA/circ/CNV shared cells; stored remap and workflow UUID recovered | Only 2/23 run mappings recovered; no raw-DNA CNV inference; thresholds/executable hashes unknown | Brief main text; supplement | C0:docs/validated_workflows_summary.md; O0:load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu /uns/circyto_scrr_remap and /uns/circyto_scrr_trimodal |
| C03 HAP1 protocol route | Executes paired-end scRR RNA/circ workflow with explicit paired-RamDA opt-in | 0.10.0 implementation@C0; original execution version/SHA unknown | GSE278952; reported batch10, exact ten cells not recovered | original workflow execution (project report) | Paired FASTQs -> STAR/BWA rescue/CIRI3 -> matrix/QC/h5ad | 10-cell RNA/circ route reported executable | Not the pilot or unresolved 63-cell/RT object; native batch10 files and selection absent | Supplement only | C0:docs/current_project_status.md; C0:docs/validated_workflows_summary.md |
| C04 Processed statistics and Figure 1B/C | Reproduces the fixed numerical evidence and RNA-only embedding overlays | Manuscript record@M0; run HEAD 7e9c1cb9a09a4f0bfd114bb66c0c34df4f31fdbc plus recorded working script SHA256=4f97474384bb975f3428d49d4effe62a4dad839af0c46c8027672a61147aec32 | O0:load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu SHA256=SS_SHA256; O0:load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu SHA256=IM_SHA256 | checksum-verified processed-object regeneration | Read stored modalities; calculate statistics; RNA normalize/log/HVG/PCA/UMAP seed 17; overlay circ | 18/18 historical values match; 10 recheck artifacts byte-identical; Smart-seq3 192/2503, IMR90 23/2443 | No raw detector rerun or recovery of detector settings; read_elem extraction is not native MuData acceptance | Main results/1B/C; supplement provenance | M0:manuscript/results_v2/provenance.json; M0:manuscript/regeneration_verification.json |
| C05 Native MuData round trips | Archived objects load and retain scientific semantics through ordinary MuData APIs | Audit record@M0; run HEAD fd52761dcf32436a53633ff7e0f9f30bf7b661cb; MuData 0.3.10/AnnData 0.11.4 | O0:load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu SHA256=SS_SHA256; O0:load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu SHA256=IM_SHA256 | native object interoperability | mudata.read_h5mu -> write_h5mu -> read_h5mu; values/axes/maps/metadata/layers checks | Both pass; 12 upstream FutureWarnings recorded | Separate circyto explicit-pull helper FAILS on Smart-seq3 index names; no all-reader or cross-version dtype claim | Main interoperability sentence; supplement details | M0:manuscript/real_object_interoperability.md; M0:manuscript/real_object_interoperability.json |
| C06 Detector-free import of archived matrices | Preserves explicitly extracted circ values, annotations and matched RNA axes | 0.10.0 development@5bc38f293e66d4051c1badac254f4501f6983acf; retained@C1 | O0:load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu SHA256=SS_SHA256; O0:load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu SHA256=IM_SHA256 | processed-object extraction/import round trip | Installed CLI matrix import + exact-ID RNA -> h5ad/h5mu/IMPORT report | Smart-seq3 2659 nnz/5711 support/192 physical overlap; IMR90 2914/7422/23; original values preserved | Not genuine native CIRI3 input; all completion statuses unknown; padded RNA-only row retained | Secondary supplement, not central novelty | C1:docs/import_circ_validation.md; C1:docs/examples/import_circ/validation.json |
| C07 Authentic native-format import | Reads the inspected 13-column CIRI3 single-sample output without value changes | 0.10.0 development@9a80ec537aed2a81495c2e14b58f9f9c20e51375; retained@C1 | gyjames/CIRI3@73108c4478a7e691cb93360921b65203e0701516:data/circRNA/Single/result.txt; SHA256=4885fc00bd30817ebdf81ba9f3344c48ddb18bcb56d2a077cee2297ff7f594ed | authentic upstream native-format acceptance | Installed importer; independent CSV/value/annotation checks; ordinary AnnData; IMPORT HTML | 1 sample x 6 candidates; support [9,9,5,6,14,7], sum 50; repeated text outputs agree | Producer version unknown; v3.0.1 pins format only; assembly and biological single-cell identity unestablished; status unknown | Secondary supplement | C1:docs/native_ciri3_validation.md; C1:docs/examples/native_ciri3/acceptance.json |
| C08 Native annotation fixture | Retains native gene/type fields where supplied | 0.10.0 development@C1 | Official data/DE/BSJ/With_RE/case1_CIRI3.txt, copied as annotated_simulation.txt; provenance.json pins source/digest | simulation-derived annotation contract | Import native gene/type/strand fields and test conflicts | Available annotation preserved; conflicting declarations rejected in regressions | Read IDs explicitly simulate:; not biological acceptance or new host annotation validation | Supplement technical validation only | C1:examples/native_ciri3/provenance.json; C1:tests/test_import_ciri3_native.py |
| C09 Synthetic RNA and observation states | Preserves exact IDs and distinguishes recorded completion from missingness | 0.10.0 development@C1 | Synthetic import demo/status cases; optional official demo uses relabeled synthetic RNA with 3 cells x 1 feature | synthetic mapping/status and technical RNA alignment | Completed-zero/failed/skipped/unprocessed/unknown; exact-ID joins; RNA layers/order/value checks | Demo denominator 2, supported completed cells 1; failed support 2 retained; technical native/RNA overlap 1 | Synthetic RNA relabeling is a test only; not real matched RNA; status does not establish locus callability | Secondary supplement | C1:docs/import_circ_validation.md; C1:docs/native_ciri3_validation.md; C1:docs/examples/native_ciri3/rna_acceptance.json |
| C10 Workflow HTML | Summarizes supplied workflow QC with conservative denominators | 0.10.0 development@C1 | Synthetic 12-cell/6-candidate report, dry-run and failed/missing fixtures | synthetic workflow/report integration | Installed report generation, offline rendering, workflow hooks, import-bundle overwrite guard | Fixture support 77; report values/repeat text outputs checked; failed/incomplete fractions unavailable | Reporter interprets supplied QC; not a fresh detector audit or executed public-workflow HTML acceptance | Secondary supplement | C1:docs/qc_reporting_validation.md; C1:docs/reviewer_readiness.md; C1:docs/pre_release_readiness.md |
| C11 Original-workflow HTML acceptance | Would test standalone report on original public-workflow QC | Candidate implementation@C1; no completed original-data test | Original Smart-seq3/IMR90/batch10 workflow_summary.json + qc/cell_qc.tsv + qc/circ_qc.tsv not locally recovered | pending original-workflow acceptance | Planned read-only source copy -> standalone workflow report | PENDING; embedded H5MU summary alone is insufficient | Do not reconstruct missing QC TSVs from H5MU and call it original workflow acceptance | Gap/plan only | C1:docs/external_ciri3_acceptance.md; 2026-10-15 bounded search in original_detection_provenance.md |
| C12 Frozen package gate | The documented baseline passed source and wheel checks | 0.10.0@C0 | Hermetic tests and clean wheel fixtures; Python 3.10.20 | installed-package/source regression | Full default tests, wheel/sdist, fresh install, resource lookup, H5MU round trip | 345 passed, 8 skipped, 5 warnings; package checks pass | Historical validation; not a new raw-data run or current public-release clearance; later packaging reviewed separately | Availability and supplement | C0:docs/manuscript_software_baseline.md |
| C13 Integration CI | Importer/MuData/report success paths execute in reviewed CI | 0.10.0 development@C1 | Hermetic fixtures plus authentic upstream demo; Python 3.10/3.12 | CI/source and installed package | Critical-module JUnit gate, builds, fresh wheel, native demo outside checkout | Local 475 passed/10 skipped/9 warnings; exact-head remote run 38351002530 success (verified 2026-10-15) | Six external detector tests remain gated; CI is not real matched-cell or detection accuracy evidence | Supplement technical validation | C1:docs/pre_release_readiness.md; https://github.com/liuifrec/circyto/actions/runs/38351002530 |
| C14 External-detector distribution | Separates detector-independent Python use from external Java installations | Implementation c9c7ac68c5b2956c38da7c1bbfd83dad7520678b; record@C2; version still 0.10.0 | Clean wheels/sdist, injected stale assets, native demo and synthetic RNA | installed-package/distribution and CI | Check omitted JAR/class/manual; keep adapters; resolve external paths; import outputs | Local 480 passed/10 skipped/9 warnings; CI run 38355386066 success on c9c7ac6 | Separate unmerged packaging branch; no release/legal blanket clearance; no historical execution version inferred | Supplement/installation decision | C2:docs/detector_distribution_validation.md; https://github.com/liuifrec/circyto/actions/runs/38355386066 |
| C15 External CIRI3 runtime | Installed adapter can execute the reviewed external JAR on official input | Packaging implementation c9c7ac68c5b2956c38da7c1bbfd83dad7520678b; external JAR SHA256=569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0 | Official demo SAM SHA256=4872d2cc5c73f95f18cf404cebc6f9cc4941106a5eb64e6c68c2a6f4a9e986c9; reference fragment SHA256=0223fb4957af56602f37de9d8c09c3a3cc21879c9becc366ad36bb164927f97b | actual external adapter execution on upstream demo | Java 17.0.18, 1 thread, -S 0 -Ma 0, no supplied GTF; independent output comparison | Six candidates, support 50; IDs/coordinates/strand/support match original upstream output | No new FASTQ alignment; demo is not established single-cell biology; these settings do not recover manuscript-run settings | Supplement technical validation | C2:docs/detector_distribution_validation.md |
| C16 Genuine per-cell native input plus matched RNA | Would establish a small real-cell import-to-RNA-to-output/report demonstration | Proposed reviewed installed importer@C2; not executed | No verified per-cell native file + RNA pair recovered | pending genuine matched-cell acceptance | Plan in nikaido_response.md: availability/identity-selected cells, independent value and ID audit | PENDING; route A recommended | Route B requires completion of this example and review; unit tests, native demo or renamed synthetic RNA cannot substitute | Gap/plan only | C2:docs/external_ciri3_acceptance.md; manuscript/nikaido_response.md |

Route **A** is recommended: the validated workflows and existing object results
remain central; C1/C2 import/report capabilities are limited secondary extensions.
The real-cell gate and original-workflow report gate remain open. Optional
long-read/RT/biogenesis material retains the restrictions below and is not
promoted by these new technical checks. No new biological claim follows.

## Main-paper evidence

The main paper should use only these evidence blocks:

1. **Need and workflow.** Describe the protocol, biological-cell mapping and
   scverse object-construction steps circyto standardizes. The task comparison
   recognizes existing orchestration, expression matrices, annotation, linear
   RNA integration and QC in other workflows; no uniqueness or superiority
   inference follows from a documentation gap.
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
| Cells | 192 | checksum-verified regeneration | `results_v2/smartseq3_summary.tsv` |
| RNA features | 63,187 | checksum-verified regeneration | Same summary |
| circRNA candidates | 2,503 | checksum-verified regeneration | Same summary |
| Nonzero circ matrix entries | 2,659 | checksum-verified regeneration | Same summary |
| Host-gene annotated candidates | 2,379 | checksum-verified regeneration | Same summary |
| Host-gene recovery | 95.0% | checksum-verified regeneration | Exact fraction 0.950459448661606; one-decimal percentage |
| Median detected circRNAs/cell | 12 | checksum-verified regeneration | Same summary |
| Median total circRNA support/cell | 22.5 | checksum-verified regeneration | Same summary |
| MAN1A2-associated candidates in historical host aggregate | 3 | local processed-data derived; needs rerun | Historical top-host table; not required in main text |
| Cells detecting representative `chr1:117402186\|117420649` | 6/192 | checksum-verified regeneration | `results_v2/smartseq3_selected_candidate.tsv` |
| Total support for representative candidate | 15 | checksum-verified regeneration | Same per-cell export; illustrative only |

### IMR90 scRR / GSE278958

Audited object SHA-256:
`bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`.

| Quantity | Value | Classification | Source/action |
| --- | ---: | --- | --- |
| RNA cells/features | 23 / 63,187 | checksum-verified regeneration | `results_v2/imr90_summary.tsv` |
| circRNA cells/candidates | 23 / 2,443 | checksum-verified regeneration | Same summary |
| CNV cells/bins | 23 / 60,607 | checksum-verified regeneration | Same summary; processed GEO CNV |
| Trimodal overlap | 23 cells | checksum-verified regeneration | Same summary; modality cell-ID intersection |
| Host-gene annotated candidates | 2,429/2,443 (99.4%) | checksum-verified regeneration | Exact fraction 0.994269340974212; supplementary-only |

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

The authoritative entry point is now
`scripts/manuscript/regenerate_application_note_results.py`; follow
[`reproduce.md`](reproduce.md). It verifies both hashes before reading stored
modality groups, regenerates the historical RNA-only embedding method, and
records full provenance. Legacy global metadata synchronization fails for these
objects in the analysis stack; the read-only modality reader avoids that step.
No input repair, production change, raw pipeline, or detector run was needed.

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

1. Draft manuscript v3 and audit its abstract/text/legends against the committed
   regenerated tables; the external manuscript document was not supplied here.
2. Confirm a public processed-data deposit and final code/data availability
   statement before journal submission. Historical Git recovery is documented;
   no new release or archive DOI is claimed.
3. Obtain coauthor review of text, figures, author/citation details, funding,
   conflicts, and contact information.

The Smart-seq3 and IMR90 regeneration and minimal supplement are complete.
Historical HAP1 full-object/RT numbers are formally excluded from this
Application Note and remain deferred, nonblocking work. The historical count
of three MAN1A2-associated candidates was not needed or regenerated and is not
a manuscript claim. This register supports v3 drafting and coauthor review;
final journal submission still requires the availability and editorial steps.
