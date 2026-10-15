# CMPB first-session study registration

2026-10-15. **Planning and read-only inventory completed; new experiments pending review.**
This is the execution registration for the supplied CMPB expansion blueprint,
not a second manuscript or a revision of the circulated Application Note.
CMPB is primary through the subscription/non-open-access route; TCBB is the
fallback. Venue selection does not authorize scientific claims, release,
submission, or simultaneous submission. The RERF form and private correspondence
are outside this session.

## Checkpoints and evidence authority

| ID | Exact checkpoint | Role |
| --- | --- | --- |
| B0 | `44697355bcab1c525ca7ef9b130e2ad0094d9e1b` | Historical software baseline, version 0.10.0 |
| M0 | `c29831c3f6e7a31b708cce493c0be8036176077a` | Circulated manuscript evidence, statistics and figure regeneration |
| I0 | `f841c296db8a6fda870eebc5b553eab10095fa0d` | PR #19, verified open/draft/unmerged |
| C2 | `6c280426ebf10abcbcaeab72f87c9d34060d6ee8` | External-detector packaging; base for this evaluation branch |
| N0 | `42459975b055d771197fbe181c6c56ffec7647cb` | Completed Nikaido response, original provenance and task comparison |
| P0 | `0f13ad4f21966f65f809ffba3bacaaebd09ec1c3` | Completed local pip/conda preparation and installation evidence |
| O0 | `c99cddaae2314d3af7d5232150c5b86bb6b52e25` | Archive containing the two checksum-identified processed objects |

Local refs and fetched remote refs agree with these checkpoints. The new branch
is `manuscript/cmpb-evaluation`, based on C2; no other branch was merged or moved.
The original checkout had no tracked changes. Its three supplied planning files
and entire `manuscript/results_v2_recheck` tree are preserved: **61 files,
71 descendant entries, 2,004,895 bytes**, using per-file SHA-256 comparisons.
An old count of 72 included a different counting convention; it is not enforced.

Reuse the following immutable records rather than treating this registration as
a replacement audit:

- [Nikaido response and contribution boundaries](https://github.com/liuifrec/circyto/blob/42459975b055d771197fbe181c6c56ffec7647cb/manuscript/nikaido_response.md),
  [91-field original provenance](https://github.com/liuifrec/circyto/blob/42459975b055d771197fbe181c6c56ffec7647cb/manuscript/original_detection_provenance.tsv),
  [215 stored cell identities](https://github.com/liuifrec/circyto/blob/42459975b055d771197fbe181c6c56ffec7647cb/manuscript/original_cell_identity.tsv).
- [Manuscript evidence register](https://github.com/liuifrec/circyto/blob/c29831c3f6e7a31b708cce493c0be8036176077a/manuscript/application_note_evidence.md),
  [reproduction guide](https://github.com/liuifrec/circyto/blob/c29831c3f6e7a31b708cce493c0be8036176077a/manuscript/reproduce.md),
  submission-readiness, figures/tables plans, regeneration and native MuData audits
  at M0. The guide's generated architecture scaffold does **not** replace the
  current BioRender Figure 1A.
- C2's [native/import/report records](../../docs/native_ciri3_validation.md),
  [integration review](../../docs/pre_release_readiness.md),
  [external acceptance gap](../../docs/external_ciri3_acceptance.md),
  [distribution record](../../docs/detector_distribution_validation.md) and
  [notices](../../THIRD_PARTY_NOTICES.md).
- [Pip/conda preparation](https://github.com/liuifrec/circyto/blob/0f13ad4f21966f65f809ffba3bacaaebd09ec1c3/docs/release_preparation.md)
  and its machine-readable record. Seven prior clean installations and the
  Python 3.10/3.12 test results are reused, not rerun. No registry publication
  is inferred. The 480-test packaging CI belongs to implementation
  `c9c7ac68c5b2956c38da7c1bbfd83dad7520678b`, not documentation head C2.

## Questions and success criteria

The contribution under evaluation is preservation and integration of cell-resolved
circRNA evidence through protocol handling, explicit observation states, biological
identifiers and annotated scientific objects. Circyto is **not a new back-splice
detector**, missing-data estimator or circRNA accuracy improvement.

| Analysis | Question | Prespecified endpoint | State |
| --- | --- | --- | --- |
| A01 | What is already established, and what original run evidence is missing? | Classified source fields, exact object identity and version-linked claims | Inventory completed; historical gaps remain |
| A02 | Do genuine native calls and same-cell RNA survive import/export and reporting? | Every value, axis, annotation, state, source hash and report value checked independently | Pilot proposed; not executed |
| A03 | What work is required for the same scientific output? | Numerical equivalence with a competent ID-aware reference; preparation operations, errors, wall time/RSS | Registered; not executed |
| A04 | What happens when identities or observation states are perturbed? | Assignment/value errors, invalid-input detection, false rejection, state preservation and same-cohort prevalence error | Registered; not executed |
| A05 | What can the stored sparse measurements show descriptively? | Within-dataset distributions, separate support thresholds and same-cell CNV/RNA/circ display | Registered; missing CNV mask information noted |
| A06 | How does import/report overhead scale? | Bounded synthetic configurations, actual file/string sizes, three process repeats, all failures retained | Registered; not executed |
| A07 | Can users install exact artifacts and read the outputs? | Prior artifact-specific pip/conda/native object checks | Existing evidence reused; final release still pending |

Acceptance never requires a positive association, a preferred circRNA count,
circyto beating the competent reference, or independent laboratory adoption.
Author-generated reruns are technical acceptance. Detection truth, clinical
utility, radiation effects and protocol superiority are outside the design.

## Available inputs and historical gaps

The two source objects were rehashed read-only in this session:

| Source | Bytes | SHA-256 | Established scope |
| --- | ---: | --- | --- |
| Smart-seq3 O0:`load_work/emtab8735_smartseq3/full_length.hostgene_fixed.h5mu` | 30,911,142 | `0ecd36bb0a74455db7f0affb9ade5023c1934c1dd234aca975365c0b69d8b339` | 192 stored identities; 2,503 candidates; 2,659 nonzeros; 5,711 support; 2,379 annotated; medians 12 and 22.5 |
| IMR90 O0:`load_work/scrr_imr90/full_length_rna_circ_cnv.hostgene_fixed.h5mu` | 41,036,354 | `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8` | 23 matched identities; 2,443 candidates; 2,914 nonzeros; 7,422 support; 2,429 annotated; imported CNV 23 × 60,607 |

Numbers above are reused audited evidence, not newly recalculated results.
Smart-seq3's `DIYHEK_192` remains RNA-only with a padded circ zero row; neither
the 192 stored rows nor the aggregate 191-success/one-empty summary establishes
cell-bound completion. Historical statuses remain unknown. HAP1 evidence remains
the reported paired-end batch10 route; its exact ten cells/native files/object
have not been recovered. Unverified 63-cell/56-cell RT numbers are excluded.

[data_provenance.tsv](data_provenance.tsv) reuses all 91 historical fields from
N0 and adds this session's verified bytes and public pilot mappings. Its four
classes are `actual run record`, `contemporaneous script with execution unconfirmed`,
`inference`, and `missing`. **Read `evidence_kind` and `chain` with the class**:
directly recorded object metadata, public metadata and a project execution report
are records at those scopes, not original per-cell detector logs. Recorded export
version 0.9.0 is not the executing workflow SHA. Current C2 defaults are not
historical settings. Original executable/ref/GTF/index hashes, executed argv,
full selection/exclusion rules and per-cell status binding remain missing.

The previous N0 search is reused. This session checked the named original root,
the two exact local object directories and the existing local work/legacy subset
directories once. The server root is unmounted; the local legacy subset contains
only a manifest. No server was contacted. No original per-cell native calls,
alignments or full workflow QC directory were recovered. This is a bounded local
finding, not an assertion about all storage. [missing_inputs.md](missing_inputs.md)
gives the finite recovery request and two alternatives.

## Fixed matched-cell pilot

Choose the first four numeric GEO accessions among the 23 stored IMR90 identities,
before viewing circ support. All four remain in the planned cohort, including
failed or unavailable outcomes. The other 19 are outside this technical pilot
solely by accession rank; this is not their historical exclusion rule.

| Canonical cell | GEO / experiment / run | ENA single FASTQ | Compressed bytes | Reads |
| --- | --- | --- | ---: | ---: |
| IMR90_A_100 | GSM8558852 / SRX26321174 / SRR30918126 | SRR30918126.fastq.gz | 225,876,952 | 3,824,354 |
| IMR90_A_101 | GSM8558853 / SRX26321183 / SRR30918117 | SRR30918117.fastq.gz | 230,818,973 | 3,935,595 |
| IMR90_A_102 | GSM8558854 / SRX26321184 / SRR30918116 | SRR30918116.fastq.gz | 236,808,421 | 4,002,138 |
| IMR90_A_103 | GSM8558855 / SRX26321185 / SRR30918115 | SRR30918115.fastq.gz | 233,399,236 | 3,954,925 |

Public GEO titles and experiment links agree with stored canonical/RNA/DNA IDs;
ENA lists one SINGLE-layout FASTQ for each. Exact URLs, MD5s, bytes, source-response
SHA-256s and acquisition state are in [benchmark_plan.yaml](benchmark_plan.yaml).
These are metadata-verified public inputs, **not downloaded local FASTQs**.
GEO's generic extraction paragraph says paired-end/150 bp, while these run records
say SINGLE and report 151 bases/read. Both are retained; inspect acquired files
before execution and use the run-specific layout. Sources: [GEO A100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8558852),
[A101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8558853),
[A102](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8558854),
[A103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8558855),
and the exact ENA API query in the YAML; accessed 2026-10-15.

Genuine same-cell RNA **is locally available** in the verified IMR90 object's
`rna` modality: 63,187 features, source-relative pilot order A103/A102/A101/A100.
Keep that order when extracting a new four-cell H5AD; circ metadata will use
A100/A101/A102/A103 to exercise real ID matching. Preserve original values,
features, layers and metadata. It is the historical lightweight `simple-overlap`
profile, not an independently generated full transcript quantification. Original
reference hashes remain uncertain. Do not add a false verified-build declaration.

Route A retrieves original per-cell native calls with status and source mapping.
Route B is a **NEW** full-library SE rerun: 15,717,012 reads, 926,903,582 compressed
FASTQ bytes, no trimming/demux/subsampling, BWA direct SAM and the pinned external
CIRI3 artifact. Proposed reference is the complete GRCh38 primary assembly and
matching comprehensive primary-assembly GTF from [GENCODE release 45](https://www.gencodegenes.org/human/release_45.html),
with new source/index hashes before execution. This choice does not recover the
historical `hg38.fa` identity. A separate new `simple-overlap` RNA profile from
the rerun is labeled separately from the existing matched RNA.

Budget for B, **estimated and unapproved**: 8 threads, one cell at a time,
24 GiB process-group RAM, 16 GiB Java heap, 80 GiB total disk and 16 hours maximum
wall time; anticipated 4–12 hours including reference/index preparation. Review
after the first cell if it exceeds 4 hours or 20 GiB RSS. Stop at enforceable
caps, retaining failed/timeout/OOM outcomes. No jobs are left running. Paired-end
extension is deferred until its identities and inputs are established.

## Proposed commands and independent checks

`benchmark_plan.yaml` is the authoritative argv/configuration list. Commands
below are **not executed**. `PILOT_ROOT` means a new approved absolute directory
named `imr90-four-new-v1`, outside historical workdirs. A reviewed helper must
first create the four-row manifest, the genuine RNA subset, and later native/status
tables from actual run evidence. Those helpers are not implemented in this session.

```bash
circyto manifest validate "$PILOT_ROOT/inputs/manifest.tsv" --strict
circyto doctor
circyto workflow full-length-circrna \
  --manifest "$PILOT_ROOT/inputs/manifest.tsv" --protocol ramda --skip-demux \
  --genome-fasta "$PILOT_ROOT/ref/GRCh38.primary_assembly.genome.fa" \
  --gtf "$PILOT_ROOT/ref/gencode.v45.primary_assembly.annotation.gtf" \
  --detector ciri3 --threads 8 --parallel 1 --chunk-size 1 \
  --export-h5ad --gene-expression-method simple-overlap --no-report \
  --outdir "$PILOT_ROOT/workflow"
circyto report --workdir "$PILOT_ROOT/workflow"
circyto import-circ \
  --ciri3-manifest "$PILOT_ROOT/inputs/native.tsv" \
  --cells "$PILOT_ROOT/inputs/cells.tsv" \
  --rna "$PILOT_ROOT/inputs/rna/imr90_original_overlap_profile.four.h5ad" \
  --detector CIRI3 --genome-build GRCh38 \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads \
  --outdir "$PILOT_ROOT/imported"
```

Reuse the existing audited C2 wheel, SHA-256
`eca5f45f35b351ec44df4e5b0f73eb4e5e0564e2691631f8eb417b0e8a7ca03b`:
all 87 packaged Python members match C2 in this session. Its 0.10.0 metadata
does not identify the frozen manuscript release. No new distribution is built
here. The prior distinct 0.11.0rc1 proposal remains an owner release decision.
Use an isolated installed environment outside source, record origins and dependency
versions, and retain external-detector exclusions. The external JAR SHA-256 is
`569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0`;
the prior help banner contains literal `$version`, so the pinned upstream release
artifact and hash establish the new executable identity. Do not infer versions
from `Java_18.0.1` or label a known unsupported original version `unknown`.

The independent oracle uses ordinary CSV parsing, exact integer arithmetic and
SciPy CSR construction, without circyto normalization/import/report functions.
It compares every native support value and coordinate, strand, gene/type field,
source byte hash, cell/feature axis, RNA field/layer, overlap and native MuData
read/write/read result. Keep ancillary native fields in raw evidence; do not
demand that every read-ID string become a feature annotation. Workflow collector
filtering, including a possible zero-cell omission, is audited separately; no
native record is dropped to make the two outputs look equal.

Audit all HTML-displayed values against independent records, TSV and JSON;
record the import report's first-100-candidate display limit. Workflow reporting
requires authentic workflow QC and conservative whole-cohort completion;
import reporting can use an explicitly completed subset. Neither report is a
status certificate. Genuine completed-zero files, if present, remain header-only;
partial failed-cell support is retained. If no such cell occurs, say not observed.

## Fair comparison and controlled analysis

The primary reference converts the **same native files and RNA** using competent
ID-aware Python/AnnData/MuData code and equivalent evidence copies, QC tables and
offline HTML. Validate both independently. Compare preparation operations and
bespoke transformations with scripts retained, and time the same complete
conversion/export/report endpoint. This is a scripted task comparison, not a
human usability study. Agreement is a successful result.

The minimum proposed named-pipeline empirical task is nf-core/circrna's supported
sample-sheet/input-channel construction on the same four real libraries. A small
reviewed harness will invoke its unmodified `PIPELINE_INITIALISATION`, including
`nf-schema@2.5.1`, and capture emitted IDs, layouts and lane grouping. The harness
is not written or run. Preserve ordinary sample IDs as biological cell IDs through
an explicit map. Repeated samples can validly mean additional lanes in nf-core;
do not score that documented behavior as an ID-loss defect. This limited test
does not establish complete workflow performance. Full nf-core/CirComPara2
execution remains not assessed and requires its own frozen task and budget.

The task assessment below extends the existing [N0 comparison](https://github.com/liuifrec/circyto/blob/42459975b055d771197fbe181c6c56ffec7647cb/docs/circrna_workflow_comparison.md).
Labels mean: **D** supported/documented, **T** observed in an existing test,
**U** requires additional user code for the specified target, **ND** not documented
in inspected sources, **NA** not assessed, **N/A** not applicable. An ND entry is
not a claim of impossibility or loss of cell identity. No speed/accuracy benchmark
is implied.

| User task | circyto C2/B0 | CIRI3 v3.0.1 | nf-core/circrna pinned commit | CirComPara2 v0.1.2.1 | circtools 2.0.4 |
| --- | --- | --- | --- | --- | --- |
| Reads/alignment → detection, counts, annotation | D/T: protocol workflows invoke existing callers; host annotation | D/T: alignment-native calling; upstream demo tested | D: preprocessing, multi-caller BSJ/linear counts, GTF annotation | D: preprocessing, multiple callers, expression/annotation | D: DCC detect/metatool, CIRIquant integration; modules have their own inputs |
| Sample metadata/workflow management | D/T: manifests, protocol/layout routes | D: sample/SAM lists; caller-specific outputs | D: sample sheet, repeated-library grouping, Nextflow profiles/resume | D: sample metadata/vars.py, SCons execution | D: module inputs and sample grouping; generic scheduler behavior NA |
| Pooled full-length single-cell protocol handling | D/T: Smart-seq3 I1/I2 mapping; scRR per-cell SE/PE routes | ND: protocol-specific demux; U for this task | D: SE/PE samples; ND: Smart-seq3 demux/scRR mapping | D: SE/PE samples; ND: Smart-seq3 demux/scRR mapping | ND: these two protocol-specific routes |
| Biological cell identity mapping | D/T: explicit cell IDs and scRR canonical map; original logs incomplete | U: caller sample → chosen cell namespace | D: sample IDs can be cell IDs; U: supply biological map | D: sample IDs can be cell IDs; U: supply biological map | D: sample grouping; U: supply biological map |
| Failed/unknown/unprocessed versus completed-zero observations | D/T: explicit import state contract; original cell status missing | ND: equivalent observation-state object contract | D: execution records; ND: equivalent multimodal state contract | D: workflow records; ND: equivalent multimodal state contract | ND: equivalent multimodal state contract |
| Native AnnData/MuData plus exact-ID RNA | T: object/matrix round trips; genuine native+RNA pilot PENDING | ND: native export; U: reference conversion | D: SummarizedExperiment; ND: native AnnData/MuData; U: conversion | ND: native AnnData/MuData; U: conversion | ND: native AnnData/MuData; U: conversion |
| QC/reporting | T: synthetic workflow/IMPORT HTML; real workflow report PENDING | D: logs and native output metrics | D: FastQC/MultiQC and execution reports | D: FastQC/HTML summary | D: quickcheck STAR/DCC QC figures |
| Installation/reproducibility | T: exact local pip/conda artifacts on P0; external tools separate | T: pinned external JAR on upstream example | D: version/config/container routes; local execution NA | D: installer/container/config; local execution NA | D: pip/Docker and R dependencies; local execution NA |

Pins inspected/reused on 2026-10-15:
[CIRI3](https://github.com/gyjames/CIRI3/tree/73108c4478a7e691cb93360921b65203e0701516),
[nf-core](https://github.com/nf-core/circrna/tree/23fdf1d22698a34e06efa99bf5af68d1bb963b78),
[CirComPara2](https://github.com/egaffo/circompara2/tree/303c3621579c1a412493e40772c2ff2660148942),
[circtools](https://github.com/jakobilab/circtools/blob/3f941e3d5ccc311caa76d8c674736da063c40d34/README.md).
Primary papers: [nf-core/circrna](https://doi.org/10.1186/s12859-022-05125-8),
[CirComPara2](https://doi.org/10.1093/bib/bbab418),
[circtools 2.0](https://doi.org/10.1186/s12859-026-06517-w).
Automation, matrices, RNA quantification, host annotation and QC are not unique
to circyto. The circtools paper also describes nanopore, conservation and padlock
probe modules; these are outside this matched short-read integration task.

Robustness uses distinct labeled perturbations and 20 fixed seeds, with fractions
0/5/10/25/50% where meaningful; deterministic invalid controls are not needlessly
crossed with every fraction. Historical statuses are never rewritten. For the
controlled simulation only, initial completion is an explicit invented condition.
Compare circyto, the competent reference and an unsafe positional/zero-fill
illustration. The heuristic does not represent nf-core or CirComPara2.

For candidate j, completed-cohort prevalence is
`sum_i[e_i * I(x_ij > 0)] / sum_i[e_i]`, undefined with no completed cells.
Compare against the **same** completed cohort, separately for uniform and
RNA-profile-depth-associated missingness. Difference from the full cohort is a
selection effect, not bias corrected by metadata. Source tables retain every
configuration, failure and null result; seeds are not biological replicates.

Descriptive analyses use the unchanged originals: per-cell support thresholds
1/2/3/5 separately from pooled-candidate thresholds, annotation and sparsity,
and RNA complexity within each dataset. CNV display requires a verified missing/
uncallable mask before its fixed 1-Mb aggregation. Numeric zero is not automatically
missing. Existing CNV importer code rejects missing numeric input; a successful
old import therefore cannot prove callability. No host-gene association is selected
in this first pilot. Synthetic scaling uses the five bounded configurations in
YAML, not a full Cartesian grid; native string/file sizes and RNA dimensions are
measured alongside nnz. All later analyses still require review.

## Evidence map and proposed manuscript outline

[claim_evidence.tsv](claim_evidence.tsv) retains separate test classes and exact
versions rather than presenting test counts as experiments.

| Claims | Evidence | Present limitation / proposed location |
| --- | --- | --- |
| C01–C03 | Historical workflow reports / stored summaries | Original provenance incomplete; Methods and dataset table |
| C04–C05 | Figure regeneration and native MuData round trips | Separate provenance chains; preserve Figure 1 and disclose helper-reader limitation |
| C06–C10 | Processed import, authentic 6/50 upstream format, synthetic RNA/state/report tests | No real matched-native cell acceptance; technical supplement |
| C11 / C16 / C18 | Original-workflow report / genuine matched-native pilot / real reports | PENDING; Figure 2 requires new evidence |
| C12–C15 / C17 | Baseline, CI, packaging, external-demo runtime and local installations | Not biological validation or public release; Availability/Supplement |
| C19–C22 | Comparison, robustness, descriptive utility and scaling | PENDING; Figures 2–4 and full source tables |

The proposed CMPB outline is a drafting choice. A body target of approximately
3,200–3,500 words, four figures and two tables is **not a verified journal limit**.
The official detailed [author guide](https://www.sciencedirect.com/journal/computer-methods-and-programs-in-biomedicine/publish/guide-for-authors)
returned HTTP 403 in this session. Verify it before final formatting; do not use
CMPB Update rules. [Official scope](https://www.us.elsevierhealth.com/computer-methods-and-programs-in-biomedicine-0169-2607.html)
includes biomedical computing/software. Acceptance and novelty remain editorial
judgments, and charges require owner confirmation before submission.

1. Abstract: Background and Objective / Methods / Results / Conclusions;
   **new Results content PENDING**.
2. Introduction: biomedical researcher task and related work; no new detector claim.
3. Methods: software/data contract; datasets and original-versus-new provenance;
   independent checks, comparators, perturbations, descriptive and scaling methods.
4. Results: established workflow/objects (existing evidence); **PENDING matched-cell
   acceptance** (`pilot_cell_trace`, `native_value_checks`, `report_value_checks`);
   **PENDING task comparison** (`reference_comparison`, `nf_input_mapping`);
   **PENDING robustness** (`robustness_outcomes`, `completed_cohort_prevalence`);
   **PENDING descriptive utility** (`descriptive_cell_metrics`, `threshold_sensitivity`,
   `cnv_display_bins`); **PENDING scaling** (`scaling_runs`). No Results prose drafted.
5. Discussion: sparse/conditional observation, historical provenance and input-format
   limits, implementation effort, comparison scope and outcomes including failures.
6. Conclusions: review after evidence; retain the approved sentence exactly below.
7. Availability/declarations: exact eventual reviewed release, data access and
   component notices; no unperformed publication or institutional clearance claimed.

> By representing circRNAs as an explicitly annotated single-cell modality, circyto enables future investigations into how circRNA profiles relate to cellular transcriptional states, host-gene expression, and genomic alterations within individual cells.

Original title/authors, BioRender 1A, audited 1B/C and circulated manuscript remain
untouched. [figure_plan.md](figure_plan.md) maps each proposed panel to source
tables and gates; no figures or previews were generated. Session output hashes,
validation argv and private path aliases are retained in the local audit packet;
the Git commit identifies these six text files. One next decision: **retrieve A's
four original native/status records (recommended), or approve B's bounded new
four-cell rerun**. No further execution follows automatically from this registration.
