# CODEX 10 — CIRCYTO CMPB integrated manuscript draft v1

**Task status:** Staged for owner-triggered Codex execution only. This is a **writing / evidence-integration task**, not an experiment, journal submission, public release, repository merge, or authorization to change scientific results. Target journal: *Computer Methods and Programs in Biomedicine* (CMPB), traditional/subscription route as discussed with RERF; TCBB remains a fallback. Venue choice and the institutional non-OA permission do not authorize manuscript submission, software publication, coauthor approval or data deposition.

## 1. Start from the accepted local figure lineage

Base the new, isolated local worktree on the **locally accepted Figures 2–4 figure-review commit**
`6e5a1a0fbc5d2a6a9361814057ec31bdc9fcb882`.
Its parent/accepted Figure 4: `81ffd7509fcaa4cbb44b7f9e82064edfd40e9837`; A06 `f858edff019b295c8e92e1c8645884ab9d2b721e`; A05 `90381340cc0365bc66c815d06ff01b2978dc3d34`; A03/A04 `02ee3d4acfed63c5b2f02c4baef6d22c822d4dac`; A02 `307758f57a39c009393d03fc0f69e4638b2498c6`. The GitHub task-planning branch is a *separate branch based on A05 for task discovery*, not a substitute for the newest local figure lineage. Verify actual local branch tips and ancestry before writing. If the accepted latest local commit cannot be found, **stop** rather than silently using older figures. Never reset, force-push or replace prior branches. Create `manuscript/cmpb-manuscript-v1-review` if free (otherwise a uniquely named alternative).

Freeze and reference the following accepted evidence, each at its original scope:
- A02 real four-cell technical acceptance: `manuscript/cmpb/a02_complete/cohort_acceptance.md`, cell trace, native checks, report checks. Preserve A100 original failed alignment-status result and separate repaired status replay.
- A03/A04: `manuscript/cmpb/a03_a04/results_for_review.md`, reference-comparison run table, complete robustness outcomes, denominators/error aggregates and source manifest.
- A05: `manuscript/cmpb/a05_descriptive/a05_results_for_review.md` and `cnv_eligibility_decision.md`.
- A06: `manuscript/cmpb/a06_scaling/results_for_review.md` and scaling run/summary tables.
- Figures 2–3: `manuscript/cmpb/figures2_3_review/` including captions, figure-source tables, validation and revised Figure 4 copies.
- Frozen earlier original manuscript: `manuscript/bioinformatics-application-note-v2@c29831c3f6e7a31b708cce493c0be8036176077a` with `manuscript/application_note_evidence.md`, `manuscript/README.md`, `manuscript/figure_legends.md`, `manuscript/reproduce.md`, `manuscript/submission_readiness.md`, existing original Figure 1B/C sources and archival provenance.
- `manuscript/nikaido_response.md`, `docs/circrna_workflow_comparison.md` at pinned Nikaido audit `42459975b055d771197fbe181c6c56ffec7647cb`, `manuscript/cmpb/study_design.md`, `benchmark_plan.yaml`, `claim_evidence.tsv`, `figure_plan.md`.

The external circulated manuscript and approved BioRender Figure 1A may not be present in the repo: **locate and inventory only supplied/current authorized copies; do not make up, silently reconstruct, edit or replace them**. The original authors/title and coauthor-approved Conclusions content are authoritative only where actually documented; report external-source gaps explicitly. The original published/project reports are not permission to invent original execution hashes.

## 2. Scientific story and contribution claims

Working direction: **"circyto: protocol-aware and provenance-preserving integration of cell-resolved circular RNA for multimodal single-cell analysis"** (draft title, pending author approval).

The central contribution is a reusable, protocol-aware **scientific data and observation-state contract**: exact biological ID matching, native caller-evidence integrity, explicit completed/failed/unknown/zero states, separate modality axes, annotated sparse scverse-compatible scientific objects, QC/report provenance and reproducible exports. **It is not a new back-splice caller, superior quantification, a missingness estimator, a speed breakthrough, nor a claim that existing competitors cannot do the job.** Existing CIRI3 is an external caller. Competent bespoke code was as correct as circyto in all tested valid/invalid A04 configurations and was faster on the tiny A03 endpoint. The novelty argument is availability, explicitness and integration of the validated contract into a reusable workflow, not mathematical uniqueness.

Preserve the Nikaido competitor response. Compare actual **user tasks** fairly with nf-core/circrna, CirComPara2, circtools and the CIRI3 native caller based on their inspected docs and accepted version pins. These packages do substantial workflow, counting, QC and annotation work; do not invent their absence of functionality from a lack of AnnData output documentation. **Named-tool full pipelines were not executed.** No nf-core/CirComPara2 speed, accuracy or biological benchmark claim. The nf-core input-only task was deferred; record it as unexecuted.

## 3. Manuscript outline, content and exact evidence boundaries

Prepare a coherent, full-length CMPB **original research manuscript draft**, **not** the compact Bioinformatics Application Note with its old one-figure structure. A body of roughly 3,500–4,500 words is a flexible internal writing target only: **do not claim it is an official CMPB requirement**. Confirm available current official author guidance if accessible; if not, mark journal-format checks as pending and avoid inventing article limits. Use scientific reader-facing headings throughout: no "A02", "A03" through "A07" work-package codes in title, abstract, Results headings, figure headings, captions or concluding text. Keep internal IDs in source/evidence provenance appendix.

Suggested narrative:
1. Structured abstract with Background and Objectives / Methods / Results / Conclusions using only observed accepted outcomes. A clear software/data-system contribution, not biological accuracy.
2. Introduction: biological importance of cell-resolved circRNA and the specific data-integration challenge; acknowledge existing callers and complete pipelines fairly. Distinguish problem definition from previously solved tasks.
3. Methods:
   - Architecture, supported input protocols, external detector dependencies, exact biological cell keys, 13-column native CIRI3 records, annotation provenance, sparse count semantics, independent RNA axis and MuData export, observation status/denominator contract.
   - Four genuine IMR90 new full-library single-end runs fixed by accession rank; **mixed execution lineage** A100 original scientific run + separate status-only repair, A101–A103 repaired workflow. New reference GENCODE release 45 GRCh38 primary assembly for the NEW rerun only. Four-cell historical RNA overlap slice with **unverified historical reference digest** and different row order. Separate three-cell fresh-workflow versus four-cell import report denominators.
   - A03 genuinely matched four-cell post-detection task versus competent independently implemented Python/AnnData/MuData reference: same scientific endpoint but different output footprints, three fresh-process timings and warm cache, not a human-usability or full-pipeline comparison.
   - A04 all 715 configurations/arm, 705 valid and 10 invalid; independent reference and deliberately unsafe illustration. Twenty fixed seeds are simulated repeats, not independent biological cells. Same-completed-cohort candidate fractions, undefined when no completed cells, no dropout correction.
   - A05 192 Smart-seq3 stored rows (191 recorded circ members plus RNA-only padded row), historical IMR90 23 stored matched RNA/circ rows, separate support thresholds 1,2,3,5, archived source facts, unknown original cell-bound completion and reference/tool hashes. Do not compare different assays' yields as protocol efficiency.
   - A06 five synthetic post-detection workloads S01–S05, one host, warm-cache three measured repeats, explicit input size/read-ID/RNA controlled contrasts, no detector time. Report process RSS separately from cgroup cache peak.
4. Results, mapping explicitly to accepted Figures 1–4 (Figure 1 source and BioRender preservation gate below):
   - Figure 1: existing original BioRender architecture (without replacement), original RNA-derived Smart-seq3 embeddings/overlays and qualified descriptive statistics.
   - Figure 2: real four-cell native evidence, exact-ID modality matching, the **distinct** historical 23-cell RNA/circ view, complete technical reconciliation. Four-cell union 435 unique candidates; 472 positive cell-feature entries; 1,699 CIRI3 junction-read support; per-cell native 142/474, 95/317, 112/419, 123/489. These are not UMIs or independent biological truth.
   - Figure 3: both competent arms exact on 705 valid and reject 10 invalid, unsafe illustration 172 valid exact / 533 valid failed checks, unsafe invalid rejections inherited. Complete table 5,159,619 candidate/configuration/arm rows, 781,440 nonzero unsafe candidate-fraction errors, 2,443 undefined fractions retained by both competent arms in no-completed control. Neither the unsafe illustration nor 5M rows represents independent biological replicates or a typical competing system.
   - Figure 4: A03 median wall 1.199 s circyto and 0.654 s custom reference, peak process RSS 187,572 vs 156,000 KiB, **reference faster and smaller on this input**. A06 synthetic 10k cells/50k circ candidates/1M nonzero entries median 47.879 s, max measured process RSS 960.46 MiB, plus S04 long-ID and S05 synthetic RNA controls. No generalized linear scaling or cross-tool superiority.
5. Discussion: contribution and fair alternatives, strengths/reuse, independent evidence, original historical provenance limits, small genuine-cell pilot, simulated missingness limits, RNA simple-overlap interpretation, CNV gate, non-availability of independent lab adoption or transcript accuracy, feature/compatibility boundaries.
6. Conclusion: preserve the previously approved sentence *verbatim*, if confirmed in the accepted original record:
   "By representing circRNAs as an explicitly annotated single-cell modality, circyto enables future investigations into how circRNA profiles relate to cellular transcriptional states, host-gene expression, and genomic alterations within individual cells."
   Do not amplify future potential into achieved biological association.
7. Availability/declarations, bibliography and supplementary inventory. Distinguish historical frozen v0.10.0, evaluated repaired but unpublished implementation commit `e6843a38469d49f96d44758a4e37893ebbe2734b`, locally tested wheel hash `28c7cbde3c31cf90dd7c247487ba5f410c5b82f6f97d7a2281a0ce52912f7383`, and eventual future release (not performed). External CIRI3 JAR excluded from the distributable artifact. Never invent an active PyPI/Bioconda release or data/archival DOI.

**CNV gate:** A05 verified 23 × 60,607 stored GEO state values but not source callability/missingness semantics; no Figure 2 CNV heatmap, calculated CNV effect, or invented zero interpretation. The 23 historical profiles are *descriptive stored measurements* with unknown per-cell completion; the four-cell new pilot does not recover their original detection provenance.

**Technical test counts** are assertions/fixtures, not experimental sample sizes or accuracy confidence intervals. No inferential p-values or CIs from three process repeats. Do not fabricate biological differences, treatment/radiation effects, or independent lab replication.

## 4. Two tables, captions, bibliography and exact figure handling

Propose **Table 1: supported and verified computational tasks across circyto and documented external workflows**, with labels like documented/tested/user-code needed/not assessed (definitions explicit, no "impossible" claims). **Table 2: dataset/input lineage and acceptance evidence** (separate historical stored and genuinely newly processed cohorts, accession/source hashes, biological denominator, limitations), or justify another clearer pair. Every cell/number has an evidence path. Supplement contains accepted full metrics/threshold/resource/robustness tables linked not recomputed. Do not copy tens of thousands of test rows into manuscript.

Reuse **Figure 1A only as the currently approved BioRender original**, if actually available and authorized; if absent create a clearly labeled placeholder/source request, not a replacement scaffold or AI-drawn figure. Reuse accepted RNA-derived Figure 1B/C with their original scientific data. Link **latest reader-facing Figures 2/3 and title-revised Figure 4**, SVG as authoritative, without rescaling scientific data or relabeling quantities. Include full numbered captions source-linked, but do not embed machine-generated graphics in place of missing Figure 1A.

The bibliography must be real: check available exact DOI/author/year/title sources in the original manuscript, Nikaido evidence, and reputable literature. Mark missing citation identifiers rather than hallucinating. Use citations to CIRI3, nf-core/circrna, CirComPara2, circtools, original Smart-seq3 and GSE278958 source, AnnData/MuData and any justified relevant prior framework literature. Do not copy abstract paragraphs from sources. Provide a separate bibliography audit table with `verified`, `incomplete`, `unverified` and exactly where checked.

## 5. Deliverables, safety, owner stop gate

Create a **new** `manuscript/cmpb/manuscript_v1_review/` directory containing at least:
- `manuscript.md` complete draft with clearly marked author/affiliation placeholders only when original approval cannot be verified.
- `abstract.md`, `figures_legends.md`, `table1_*.tsv`, `table2_*.tsv`, `supplementary_inventory.md`, `references.bib` (only verified entries plus separately tracked missing sources).
- `claim_to_sentence_audit.tsv` (every quantitative manuscript sentence → accepted exact source, version, cohort, check); `reference_audit.tsv`; `manuscript_readiness.md` with unresolved scientific/reviewer decisions and source/figure availability.
- `editorial_decisions.md` explaining differences from the original Application Note and no internal A-code reader-facing headings; `diff_preservation.json` confirming no original text/figures/production/data altered.

Check the draft for errors of scope/denominator, false biological/competitor claims, reference integrity, figure labels, internal A-code leaks, software lineage and consistent units; use a separate lightweight linter, but manuscript logic needs human review. No new analysis, no reloading large read files or full 715-case grid, no science alteration, no correcting historic sources to fit story.

Use existing local tools; target a bounded **writing-focused** session (suggest <= 90 minutes with no large computations). Commit **new manuscript review directory only** on isolated local branch, clean diff; **do not push, merge, tag, publish, submit, contact coauthors, send emails or change release artifacts** without a separate owner instruction.

Export a compact, checksum-verified ZIP to `~/Downloads/CIRCYTO_CMPB_MANUSCRIPT_V1_<timestamp>/` containing the complete new draft, tables, citations and audit packet; do not include H5MU/FASTQ/JAR/environments/bulk objects/private details. Report exact path, size, SHA256, local commit, outstanding source gates, draft word count and whether all four figures have accepted/located source files. Stop for owner review. Coauthor approval and final CMPB journal formatting remain subsequent distinct steps.
