# Response to Nikaido's three manuscript comments

Prepared **2026-10-15** for review. This is a revision package, not an edited or
submitted RERF manuscript. The current RERF submission, BioRender Figure 1A,
numerical evidence and Yoshida's approved Conclusions wording are preserved.
The external submission/BioRender/approved Conclusions source files were not
supplied here; no replacement wording for the Conclusions or Figure 1A was
created. Final supplementary numbering is deferred.

## Verified starting state

Local heads, origin-tracking heads and GitHub branch heads all matched the
four supplied snapshots on 2026-10-15:

| Purpose | Verified ref |
| --- | --- |
| Production implementation | `main@44697355bcab1c525ca7ef9b130e2ad0094d9e1b` |
| Manuscript base | `manuscript/bioinformatics-application-note-v2@c29831c3f6e7a31b708cce493c0be8036176077a` |
| Import/report integration | `review/pre-release-integration@f841c296db8a6fda870eebc5b553eab10095fa0d` |
| External-detector distribution | `fix/external-detector-distribution@6c280426ebf10abcbcaeab72f87c9d34060d6ee8` |

[PR #19](https://github.com/liuifrec/circyto/pull/19) remained **open, draft,
unmerged**, targeting main. Its exact-head
[CI run 38351002530](https://github.com/liuifrec/circyto/actions/runs/38351002530)
completed successfully. That result does not close the biological acceptance
or original-provenance gaps.

The original checkout remained on the distribution branch, with no tracked
changes and four untracked duplicate figures under `manuscript/results_v2_recheck/`.
The entire recheck tree, including ignored files, was inventoried (61 files,
2,004,895 bytes) and retained. Work was done in the separate
`/tmp/circyto-nikaido-response` worktree on `manuscript/nikaido-response`, based
on the verified manuscript commit. Other refs were read with `git show`; none
was merged or changed. The response uses the existing evidence register,
submission record, reproduction guide and figure/table plans rather than a
second readiness system.

## Point-by-point response register

“Addressed” means addressed in this reviewable repository package. It does not
mean the external manuscript has been revised or approved by collaborators.

| Comment / requested change | Manuscript or supplement location | Evidence and exact ref | Action completed | Unresolved information | Status |
| --- | --- | --- | --- | --- | --- |
| **1. Compare user tasks with nf-core/circrna and CirComPara2; explain protocol, identity and scverse work.** | Introduction draft below; descriptive comparison for supplement, number unassigned. | nf-core/circrna `23fdf1d22698a34e06efa99bf5af68d1bb963b78`; CirComPara2 v0.1.2.1 `303c3621579c1a412493e40772c2ff2660148942`; primary papers; circyto C0/C1/C2 above. | Added pinned task comparison and bounded Introduction text; recognizes comparator matrices, linear-RNA integration, workflow management and reporting. | External manuscript integration/coauthor review; comparative speed/accuracy/user-error rates are not assessed. | **addressed** |
| **2. Connect claims to implementation, dataset and tested behavior; genuine matched-cell demonstration if import/report becomes central.** | Existing `application_note_evidence.md` claim matrix; `claim_evidence.tsv`; route decision below. | Manuscript `c29831c3f6e7a31b708cce493c0be8036176077a`; archive `c99cddaae2314d3af7d5232150c5b86bb6b52e25`; C0/C1/C2 acceptance reports, with individual milestone refs in matrix. | Replaced broad feature matrix with one 16-row feature/version/input/test/operation/result/limit/placement matrix. Separated original workflows, regeneration, native MuData, processed import, authentic demo, synthetic cases, installation/CI and report gates. Recommend route A. | Genuine per-cell native CIRI3 + same-cell RNA and original-workflow HTML acceptance remain pending. A plan is prepared; no such demonstration was executed. | **partially addressed** |
| **3. Document ORIGINAL detection conditions separately from figure regeneration.** | Methods provenance note; original-detection table and exact stored-cell list; `methods_commands.md` and `reproduce.md` warnings. | O-SS/O-IM at `c99cddaae2314d3af7d5232150c5b86bb6b52e25`, verified hashes; S-May `9e028d8dd6b054bbf34fb77debb4ec095916c359`; S-June `777dc3406e9793703aba79bac7dfa2473911d739`; C0 server/workflow records. | Read embedded original summaries, recovered original source paths, options, versions and 215 stored cell identities; classified fields D/R/U; retained conflicts; prepared narrow read-only collection list for James. | Original launch/executable records, full accession/selection lists, FASTA/GTF identities/digests and native-result→object checksum links; batch10 exact identities. | **partially addressed** |

Proposed response to comment 1: We now compare circyto with both requested
workflows by the work a user performs. The revised framing recognizes existing
orchestration, quantification, annotation, expression matrices and QC, and
specifies the protocol, cell-mapping and AnnData/MuData steps circyto standardizes.
Documentation gaps in comparators are labelled as such, not claimed absences.

Proposed response to comment 2: Each claim now identifies the tested software,
input source, operation and limitation. We recommend retaining the established
workflows as the main contribution. Import/report extensions have narrower
technical evidence; a genuine native per-cell input plus matched-RNA example
is still required before making them central.

Proposed response to comment 3: The original and reproduction chains are now
separate. We recovered original metadata from the verified objects but cannot
yet state all original detector versions/settings or complete input/reference
checksums. The table marks those fields unknown and lists the precise server
records needed. Script defaults and pilot settings are not reported as executed
manuscript-scale conditions.

## Revised Introduction passage — draft for review

Existing workflows, including nf-core/circrna and CirComPara2, automate circRNA
detection, quantification and annotation and provide expression summaries and
quality-control outputs. [Digby et al., 2023](https://doi.org/10.1186/s12859-022-05125-8);
[Gaffo et al., 2022](https://doi.org/10.1093/bib/bbab418).
For full-length single-cell analyses, users must also connect protocol-specific
read processing and detector outputs to biological cell identities and align
those identities with RNA measurements. circyto standardizes these steps through
pooled Smart-seq3 index-to-cell demultiplexing, manifest-driven single- and
paired-end scRR/RamDA routes, explicit cell mapping, and annotated cell-by-circRNA
AnnData and RNA/circRNA MuData outputs. These interfaces support reuse of detector
support values and cell metadata in single-cell analysis. **circyto is not a new
back-splice detection algorithm; it orchestrates established detectors and
standardizes their cell-resolved outputs for downstream analysis.**

The software assertions in this draft are tied to C01–C05 and C12 in the
[claim matrix](application_note_evidence.md#claim-to-evidence-matrix). The
[task comparison](../docs/circrna_workflow_comparison.md) supplies versioned
documentation and the access date. The passage makes no unique/first/superior
claim and does not imply a measured reduction in operator errors.

## Manuscript route recommendation

**Route A: retain the validated workflows as the main contribution.** Use the
existing Smart-seq3 and IMR90 numerical evidence and the bounded HAP1 batch10
execution report. Describe detector-free/native import and offline HTML as
limited secondary software extensions, with their distinct tests and pending
real-cell/report acceptance visible. The authentic upstream demo is one sample
of unestablished biological single-cell identity; its synthetic RNA alignment
test is not real matched RNA. Historical completion status remains unknown.

**Route B remains conditional**, requiring the genuine matched-cell example below,
review of the updated implementation and its report semantics. Passing unit
tests or CI is insufficient. A final reviewed release can eventually be the
installation target while Methods retains the original analysis versions and
checksums. The response does not set a new version, release a wheel or make
internal branch chronology the scientific narrative. Original metadata gaps
still need attention under either route; route A does not turn unknown Methods
fields into known ones.

Suggested Methods provenance note for review:

> Original workflow execution and later processed-object regeneration were
> audited separately. Archived Smart-seq3 and IMR90 objects retain circyto 0.9.0
> export/refresh metadata; the subsequently validated implementation baseline
> is 0.10.0. Original detector executable versions, full input manifests and
> reference checksums remain under source-record verification. Historical
> script settings are therefore reported as unconfirmed reconstructions, not
> as executed conditions. Figures and summary statistics were regenerated from
> the specified checksum-matched objects using the separately recorded analysis
> environment.

This is transparent interim wording, not a substitute for retrieving the
missing Methods records before final external manuscript revision.

## Small real-cell acceptance plan

**Prepared only; inputs incomplete; no execution authorized or performed.**
Existing genuine RNA is available inside the checksum-matched IMR90 object.
The missing inputs are original native per-cell CIRI3 files and run-bound cell
identities/status/reference records. A normalized table extracted from H5MU
cannot serve as native acceptance.

Prefer two IMR90 cells because public GSM/SRX/SRR mappings and original RNA
identities are already recorded for them. Fix this proposed subset by identity
and availability, without inspecting circRNA support:

| Original library/sample | Run | Shared biological cell ID | Candidate native-file locator (unverified) | RNA input |
| --- | --- | --- | --- | --- |
| GSM8558852 / SRX26321174, RNA_IMR90_A_100 | SRR30918126 | IMR90_A_100 | `work_hg38_full23_v2/ciri3/GSM8558852.ciri3_run/ciri3_raw.tsv` | O-IM `/mod/rna`, row IMR90_A_100 |
| GSM8558853 / SRX26321183, RNA_IMR90_A_101 | SRR30918117 | IMR90_A_101 | `work_hg38_full23_v2/ciri3/GSM8558853.ciri3_run/ciri3_raw.tsv` | O-IM `/mod/rna`, row IMR90_A_101 |

Candidate locators come from the historical adapter naming contract, not a
file-existence claim. Confirm them against original summaries and manifests;
keep original bytes unchanged. If either identity/file cannot be established,
record the gap and stop this fixed example. Do not replace a cell based on
its signal. RNA extraction would retain actual existing cell IDs, values,
features and annotations and be described as processed same-cell RNA from the
archived workflow, not independent expression quantification.

Before an execution request, fill the two source paths, sizes, SHA-256 values,
producer versions, single-cell scope, reference identity and status evidence.
Install a wheel built from the reviewed C2 implementation in a new environment
outside the checkout; record wheel digest, dependencies and `site-packages`
import path. Version string 0.10.0 alone cannot identify this unreleased code.
Supported native format is `ciri3.native-result-13.v1`, reviewed at upstream
v3.0.1 `73108c4478a7e691cb93360921b65203e0701516`. Producer version must be
truthful: `unknown` only if genuinely unavailable, not to bypass a known
unsupported format/version.

Prepare `native.tsv` with original sample IDs and the two canonical cell IDs;
`cells.tsv` fixes that order and uses only evidence-backed status. Header-only
native files may establish a completed zero only with completion records.
Extract the two genuine RNA rows to `matched_rna.h5ad`, hash it and document
the extraction. Proposed local storage is a fresh
`/tmp/circyto-imr90-real-cell-acceptance-<date>/`, never an existing result tree.
For two files and the two-cell RNA slice, propose a **1 GiB output allowance,
4 GiB RAM and 10-minute import/report limit**; these are conservative planning
ceilings, not measured runtime/scaling results. Environment installation and
source acquisition are separate.

Run the **installed** importer after those inputs are reviewed and
`VERIFIED_BUILD` has been set to the documented assembly:

```bash
/usr/bin/time -v -o import.time.txt circyto import-circ \
  --ciri3-manifest native.tsv --cells cells.tsv --rna matched_rna.h5ad \
  --detector CIRI3 --genome-build "$VERIFIED_BUILD" \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads --outdir imported
```

Independently parse native inputs and check every candidate ID, coordinate,
strand, support and annotation; compare sparse values/axis order, retained RNA
data/layers/annotations and exact overlap with ordinary AnnData/MuData readers.
Verify statuses and unavailable denominators, source hashes and offline IMPORT
HTML. Repeat into a second new directory to compare JSON/TSV/MTX/HTML and HDF5
semantics. This is an **IMPORT report**; never run `circyto report` on that
bundle. If the original workflow summary and both QC TSVs are separately
recovered, copy them to a distinct directory and test standalone workflow HTML
there, preserving source checksums and checking the displayed original QC.

### Bounded fallback if original native files are unavailable

Prepare a **new validation**, subject to explicit approval, for the same two
public runs, single-end, without signal-based selection or read subsampling.
Fetch only their RNA-side FASTQs and accession/checksum metadata, use a new
output root, retain native CIRI3 results and run-bound logs, and pair them to
the verified same-cell RNA IDs above. No raw DNA or other samples are needed.
The proposed adapter is C2 with the reviewed external CIRI3 JAR SHA-256
`569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0`; the documented
single-end route uses BWA/direct SAM, so STAR is unnecessary for this example.

The approval proposal must first pin BWA/Java executables and digests,
FASTA/GTF release/bytes, source FASTQ sizes/checksums and the complete commands.
They are deliberately not filled from supposed historical defaults. Proposed
resource ceilings are **two runs, 16 threads, 64 GiB RAM, 120 GiB scratch and
8 hours wall time**, with input-size/reference-index preflight and a stop if
these ceilings cannot accommodate full inputs. These are budgets, not a claim
of measured feasibility. Record all new settings as new validation metadata;
even successful new detection cannot recover the original manuscript settings.
No download, reference indexing, detector execution or acceptance run was done.

## Review artifacts and minimum next action

- [Task comparison](../docs/circrna_workflow_comparison.md), with primary papers
  and immutable software pins.
- [Single readable claim matrix](application_note_evidence.md#claim-to-evidence-matrix)
  and [editable TSV](claim_evidence.tsv).
- [Original-detection table, conflicts and source-record request](original_detection_provenance.md),
  [field-level TSV](original_detection_provenance.tsv) and
  [exact stored-cell list](original_cell_identity.tsv).
- Existing [submission readiness](submission_readiness.md),
  [reproduction guide](reproduce.md), [figure plan](figures_plan.md) and
  [table plan](tables_plan.md) now point to these bounds. No S1/S2/S3 asset is
  overwritten or renumbered.

**Minimum next action:** James retrieves the named original manifests,
command/version records and reference/object checksums using the scoped
read-only list. That request is prepared here, not sent. Review the Introduction
and route A while those records are collected. Stop broader recovery until
the records arrive; do not start a raw rerun to fill historical blanks.

## Validation of this documentation package

Completed 2026-10-15: TSV structure and unique keys checked; all 16 readable
claim rows agree with the TSV; all 215 cell rows/order/mapping fields checked
against the unchanged original objects; both object SHA-256 values match;
31 cited Git ref/file pairs resolve; local document links and fenced command
syntax checked; `git diff --check` passes. The 91 provenance fields comprise
38 directly recorded, nine script reconstructions and 44 unknowns.

The quantitative register and all result/supplement assets are unchanged.
All 61 inventoried original recheck files retain their sizes and SHA-256 values;
the original checkout status and four protected branch heads match the initial
state. No production code, package version, test or executable evidence script
changed, so no production test suite or raw-data workflow was rerun. Syntax
checks of proposed commands are not acceptance execution.
