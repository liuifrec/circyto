# Finite input recovery request and pilot alternatives

2026-10-15. **Original per-cell native calls and bound status records are not
locally recovered.** Both original processed H5MU byte identities and the four
proposed cells' stored RNA identities are verified. Public FASTQ metadata is
available; FASTQs, full reference and new index have not been acquired for this
pilot. No SSH, server download, analysis or report generation occurred.

This request is for James to review, not a message sent to a collaborator.
Server/mount paths use aliases; their private local mapping is retained outside
Git. No hostname, credentials, private correspondence or full server mirror is
requested. Only records needed for the stated cells/cohorts should be returned.

## A. Recover original records first (recommended)

The existing object-linked leads are:

| Alias and relative workdir | Binding key | What it establishes |
| --- | --- | --- |
| `HISTORICAL_SMARTSEQ3_ROOT/work/diySpike_workflow_all192` | UUID `03d6099d-7800-4ac2-be15-6bb93e9debb3`; Smart-seq3 source SHA-256 in data_provenance.tsv | Directory recorded in archived object, not confirmation that remaining files are intact |
| `HISTORICAL_IMR90_ROOT/work_hg38_full23_v2` | UUID `6c62bf6c-8fa5-49fe-a30e-bc6db7e19232`; IMR90 source SHA-256 in data_provenance.tsv | Use the original object-bound workflow, not the two-cell pilot |
| `HISTORICAL_HAP1_ROOT` | Original validated batch10 manifest/log required | Exact workdir and ten cells still unknown; do not substitute `work_hg38`, `work_hg38_fullstack`, a pilot or the 63-cell/RT object |

The UUIDs above are the exact values in data_provenance.tsv (D004 and D037).
Match them to the returned summaries; no fuzzy UUID match may join runs.

For each of the two object-linked workdirs request these existing relative files,
with size and SHA-256. Missing entries should be marked missing, not recreated:

```text
workflow_summary.json
manifests/all_manifest.tsv
manifests/single_end_manifest.tsv
manifests/paired_end_manifest.tsv
demux/manifest.tsv
demux/demux_summary.json
align/alignment_manifest.tsv
align/alignment_prepare_summary.json
ciri3/detector_run_summary.json
logs/run_ciri3_summary.json
matrix/cell_index.txt
matrix/circ_index.txt
matrix/circ_feature_table.tsv
qc/cell_qc.tsv
qc/circ_qc.tsv
rna/rna_import_summary.json
```

For the four fixed IMR90 pilot cells, use the detector manifest/summary to resolve
the actual single-sample native outputs and their per-cell logs. Current adapter
path forms give these **candidate locations, not recovered files**:

```text
ciri3/GSM8558852.ciri3_run/ciri3_raw.tsv
ciri3/GSM8558853.ciri3_run/ciri3_raw.tsv
ciri3/GSM8558854.ciri3_run/ciri3_raw.tsv
ciri3/GSM8558855.ciri3_run/ciri3_raw.tsv
ciri3/GSM8558852.ciri3.log
ciri3/GSM8558853.ciri3.log
ciri3/GSM8558854.ciri3.log
ciri3/GSM8558855.ciri3.log
```

If original filenames use canonical IDs, return the explicit mapping instead of
renaming evidence. The needed cells are A100/A101/A102/A103 and runs
SRR30918126/117/116/115 respectively. Do not send unrelated raw read files.
Retain every native ancillary field and any genuine header-only result. Request
the command record, exit code and log that bind each cell and input to that
native file; file existence or a positive count does not establish completion.

Also request the following narrowly scoped original records:

| Needed record | Exact purpose / scope |
| --- | --- |
| IMR90 `manifest_all.tsv`, original selected manifest and `scrr_cell_map.tsv` | All 23 intended/selected cells, original accession/library/cell mapping, exclusions and reasons; retain original file hashes |
| Smart-seq3 original normalized plate annotation, run-to-R1/R2/I1/I2 list and selection record | Connect pooled source runs, index pairs and 192 stored identities; verify top-N/min-read/subsampling parameters and exclusion reasons |
| HAP1 batch10 selected/input manifest and its exact workflow summary | Identify the ten validated cells and layout; bounded immediate candidate-workdir listing only if its name is unknown |
| Original job script/launch log and checked-out revision plus patch, if any | Exact argv, threads, chimeric/rescue/strand/layout settings, index generation, filters, count units and duplicate treatment; no entire shell history |
| Run-bound STAR/BWA/CIRI3 version output and executable/JAR checksums | Original tool identity only. A checksum of today's executable is a present-day record unless linked to the preserved original binary |
| Original reference/GTF acquisition manifest, FASTA/GTF hashes and index build record | Resolve assembly/provider/release/contigs and bytes; `hg38.fa` or a v45 filename alone is insufficient |
| Original annotation repair, remap, CNV merge and object export commands plus input/output hashes | Link native calls and intermediates to the exact archived H5MU, rather than merely a similar object |
| `GEO_lists/GSE278958_06_scRR-seq-DNA_IMR90_G1_summary_CNV_states_bin50kb.txt.gz` and `..._mappabilitynorm_bin50kb.txt.gz` provenance/headers and missing/callability definition | Determine original CNV missing-state semantics before new heatmap aggregation; request full public tables only if needed in the separately approved descriptive analysis |

Start with metadata/manifests/status logs and the four small native results.
Inspect sizes first; if a selected log/result is unusually large, return its size
and relevant command/status excerpt for review before copying. Reference/FASTQ/
alignment **checksums and manifests** suffice for this first recovery request;
do not copy raw data or entire reference trees. Retain missing records explicitly.

Optional read-only checks for the person who already has authorized server access
(not run here; set the aliases privately):

```bash
stat -- "$HISTORICAL_IMR90_ROOT/work_hg38_full23_v2/workflow_summary.json"
sha256sum -- "$HISTORICAL_IMR90_ROOT/work_hg38_full23_v2/workflow_summary.json"
sha256sum -- "$HISTORICAL_IMR90_ROOT/manifest_all.tsv" \
  "$HISTORICAL_IMR90_ROOT/scrr_cell_map.tsv"
sha256sum -- "$HISTORICAL_REFERENCE_ROOT/hg38.fa" \
  "$HISTORICAL_REFERENCE_ROOT/gtf/gencode.v45.annotation.gtf"
```

These commands are only inventory/hash collection. The last pair describes
current available bytes unless immutable run records establish their historical
identity. Do not start or resume a workflow. Original QC files must be copied
to a **new** review directory before any later standalone report writes HTML/JSON.
Never construct missing workflow QC tables from H5MU and call that recovery.

Conflicting evidence stays visible: the two-cell server pilot is not the final
23-cell run; public ramdaq/GENCODE v38 metadata is not proof of circyto's reference;
stored v45 paths still lack digests; current defaults are only code inspection;
GEO generic paired-end prose differs from these run-specific single-end files.

## B. New technical rerun if A cannot supply the records

The exact fixed cohort and HTTPS FASTQ URLs/MD5/bytes are in
[benchmark_plan.yaml](benchmark_plan.yaml), with public mappings in
[data_provenance.tsv](data_provenance.tsv). All four library identities are verified;
their genuine same-cell RNA profiles are present in the original H5MU. The
926,903,582-byte public FASTQ total is an acquisition estimate from ENA metadata,
not a claim of local download or authenticated local SHA-256.

Before running, require approval of this new validation and its budget, acquire
complete libraries, check ENA MD5/size and FASTQ headers/counts, record SHA-256,
and prepare the specified full GENCODE v45 primary assembly/GTF with provider
checksums and local decompressed/index SHA-256s. The source links are fixed; the
bytes have not been downloaded. BWA 0.7.17-r1188 and Java 17 are proposed dependency
versions, not historical claims. Verify exact executables and the already local,
hash-matched external CIRI3 JAR before execution. Capture the real emitted argv
and status; retain the native 13-column schema and single-sample scope.

Use the installed C2 artifact, a new workdir and the commands/config in the YAML.
The four-cell RNA slice retains original cell names and relative order. The
workflow may also generate a separately labeled new lightweight RNA profile.
Never manufacture native input by extracting the old circ matrix. No chromosome
restriction or read subsampling is proposed. Any changed cohort/reference/layout
returns for review before execution.

Estimated maximum: **8 threads, 24 GiB RAM, 80 GiB disk, 16 hours**, with first-cell
review triggers and limits in the YAML. Enforce caps before running. Actual
runtime and memory are unknown; preserve OOM/timeout/setup failures. Report all
four cells without replacing failed or low-signal cells. This budget covers the
pilot only; later comparator, robustness, descriptive and scaling jobs remain
separately pending. No obligation to leave jobs running during APBC.

What B can establish: an author-generated, documented native→matched-RNA→object
and report technical acceptance, including a genuine NEW workflow report.
What B cannot establish: original detection conditions, independent laboratory
validation, circRNA accuracy, historical per-cell completion or new biology.

## Remaining gates outside the immediate pilot decision

- Historical Methods still need exact original selection/versions/settings and
  native-to-object linkage; a successful B does not fill those fields.
- Figure 2C needs source CNV missing/callability semantics. If unavailable, return
  a narrower display proposal; do not invent a mask or replace source bins.
- The independent conversion and nf-core input-stage harness are planned, not
  delivered executable scripts. Review their source/contract before timing.
- The genuine paired-end extension is deferred; HAP1 batch10 and Smart-seq3 raw
  provenance are insufficient for a fixed additional cohort today.
- Current CMPB detailed guide returned HTTP 403. Check official formatting/charges
  before submission; the institutional form is not edited in this task.
- Package release, public data deposit, DOI, tags, journal submission and coauthor
  scientific clearance remain separate decisions. Prior local pip/conda acceptance
  does not authorize publication.

**Single next decision:** choose A (recommended finite original-record retrieval)
or approve B's fixed four-cell NEW rerun and resource cap. Stop at this checkpoint.
