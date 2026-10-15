# A02 first-cell technical acceptance

**FAILED: workflow alignment status is lost before QC/AnnData/report generation.**
The new A100 alignment and CIRI3 run completed successfully. Its native import
preserved all checked counts, annotations and same-cell RNA fields. Nevertheless,
the workflow QC and AnnData label alignment `missing`; the standalone workflow
report consequently says `incomplete` and withholds completed-cohort prevalence.
This fails the complete first-cell acceptance. No production repair, second
biological cell, A03–A06 analysis or publication was performed.

Execution: 2026-10-15. Handoff finalized 2026-10-16 Japan time after a server
restart interrupted audit writing. No managed job remained unfinished and no
scientific stage was resumed. Registration base:
`17b733c8f21bb128a035e7b100e704099e2bdb0e`, branch
`manuscript/cmpb-evaluation`. Evaluated package/source:
`6c280426ebf10abcbcaeab72f87c9d34060d6ee8`. The commit containing this packet is
the review head; its exact hash is reported in the final handoff.

## Actual cohort status

| Cell | GSM / SRR | New execution | New detector evidence | Acceptance |
| --- | --- | --- | --- | --- |
| IMR90_A_100 | GSM8558852 / SRR30918126 | Finished; workflow exit 0 | Completed; BWA and CIRI3 exit 0, input/output hashes bound | Failed: workflow status propagation |
| IMR90_A_101 | GSM8558853 / SRR30918117 | Not started | Unprocessed | Not started |
| IMR90_A_102 | GSM8558854 / SRR30918116 | Not started | Unprocessed | Not started |
| IMR90_A_103 | GSM8558855 / SRR30918115 | Not started | Unprocessed | Not started |

Original detector completion remains **unknown** for these historical objects.
The four-cell pilot is incomplete. See [pilot_cell_trace.tsv](pilot_cell_trace.tsv)
and [pilot_progress.json](pilot_progress.json). No cell was replaced based on
signal, and no combined four-cell workflow summary was synthesized.

## Exact failure and retained evidence

All paths in this section are below the durable evaluation root
`/home/liuyuchen/data/circyto_cmpb_a02_20261015`.

| Evidence | Actual observation |
| --- | --- |
| `workflow_firstcell/align/bwa_mem/alignment_prepare_summary.json` | A100 `status=aligned`; `status_counts={"aligned":1}`; no failed cells |
| `workflow_firstcell/align/alignment_prepare_summary.json` | Root routing summary has manifests/layout counts but no per-cell alignment records or status counts |
| `workflow_firstcell/workflow_summary.json` | One planned/completed cell; completed alignment/detector stages; empty alignment status counts |
| `workflow_firstcell/qc/cell_qc.tsv` and `anndata/circ_counts.h5ad` | A100 alignment `missing`, detector `success`, 142 candidates, 474 support |
| `workflow_firstcell/qc/metrics.json` and `qc/report.html` | Workflow `incomplete`; reason: alignment status missing; prevalence denominator unavailable |

Source inspection of the installed, wheel-matched C2
`circyto/pipeline/workflow_full_length_circrna.py` (`_write_outputs`, around
lines 469–486) shows that QC reads the root alignment summary's per-cell map.
The successful per-cell records are instead present in the `bwa_mem` summary.
This identifies the status propagation mismatch. The report correctly preserves
the incomplete evidence it receives; editing its HTML or filling in QC status
by hand would conceal the defect and was not done.

Four explicit failures appear in
[report_value_checks.tsv](report_value_checks.tsv): recorded alignment versus
report status, recorded alignment versus workflow AnnData status, overall report
completion, and completed-cohort denominator. The independent status audit,
source hashes and unchanged output bytes remain in
`records/acceptance/status_acceptance.json`. No raw alignment/detection rerun is
required merely to inspect this evidence.

## Completed technical checks

| Check | Actual result and boundary |
| --- | --- |
| Full public input | 3,824,354 reads; 577,477,454 bases; gzip CRC/FASTQ structure/unique read IDs/read counts passed; ENA accession, layout, URL, bytes and MD5 matched registration |
| Deposited layout | One SINGLE file; every preserved original flow-cell header has `/1`, none `/2`. This establishes the deposited single-end input; generic GEO paired-end prose remains a documented metadata conflict |
| Full reference | Complete GENCODE v45 GRCh38 primary-assembly FASTA and matching primary-assembly GTF; provider compressed MD5s verified, compressed/uncompressed SHA-256 recorded, contig/annotation compatibility checked; all five BWA index files hashed |
| Execution | BWA 0.7.17-r1188; external CIRI3 3.0.1 JAR; Java 17.0.18-internal; registered BWA `-t 8 -k 15 -T 15`, CIRI3 `-T 8 -S 0 -Ma 0`; workflow `simple-overlap` and `--no-report` |
| Native evidence | 142 native candidates, 142 nonzeros, 474 `#junction_reads` support units. Native source SHA-256 `4c807adc516e5367a315b4ac3d6d41c9d6d573a3bf222904c05865c24c5964f6`; ancillary fields/read IDs retained byte-for-byte outside Git |
| Historical same-cell RNA | Original H5MU SHA-256 `bb2e12f7c3b36f9fa72d98cd71e8bea905a67f50e22af1d6b713550ee92b60c8`; A100 slice 1 × 63,187; all available values, features, layers, annotations, embeddings and original metadata preserved |
| Native import / objects | 149 independent value checks passed, including every native coordinate/ID/strand/count/annotation, source bytes, axes, MatrixMarket and ordinary AnnData/MuData read/write/read. Genuine RNA identity matched exactly |
| Numeric/rendered reports | 478 numeric/HTML consistency checks passed before the four semantic status failures were appended. Both separate reports opened locally with DNS disabled; browser tables/cards/SVG values match source HTML |
| New overlap RNA | Separate same-run lightweight simple-overlap profile, 1 × 63,241; 2,880,944 assigned templates. It does not replace or retrospectively verify the original RNA profile |

The original RNA reference identity remains incompletely established. Import
metadata correctly keeps `rna_genome_build_verified=false`. Neither RNA profile
is an independent transcript quantification. The importer preserves native gene
IDs/type fields; workflow host annotations are separately GTF-derived. Native
and workflow count values agree here, with workflow AnnData using int32 and the
native import int64. No claim is made about legacy int32 counts beyond its range.
No zero cell column was omitted in this positive-support cell; the registered
legacy collector's cell-total filter is documented in the reproduction record.

## Resources, failures and audit limitations

Managed shared preparation took 4,344.31 s (72.41 min), including full acquisition
and reference checking (1,203.09 s) and BWA index preparation (3,050.57 s).
First-cell managed stages took 691.77 s; the workflow itself took 669.29 s
(11 min 9 s), native import 1.63 s, standalone report 0.88 s. The external-tool
wrapper measured BWA at 493.38 s and CIRI3 at 39.47 s; these times include
post-execution hashing and are not isolated algorithm benchmarks.

All managed work, including control tests, totaled 5,046.32 s (84.11 min).
The calendar execution envelope through the last stage was 2 h 4 min 32 s;
the first-cell envelope was about 28 min 28 s. Peak sampled cgroup memory was
18.14 GiB, summed process-tree RSS 11.31 GiB, and retained evaluation storage
19.23 GiB. Sampling and accounting limits are explicit in
[resource_usage.json](resource_usage.json). The later server-restart/handoff gap
is separately recorded; it is not detector time. The conservative absolute
16-hour launch deadline has now expired: no budget was reset, and future work
requires review of remaining resources.

User-service memory/swap and process-tree wall limits were enforced. CPU quota
was unavailable on this host; inherited affinity to CPUs 0–7 and a watchdog
checking every descendant thread were tested. The disk watchdog stops at 79 GiB,
below the 80 GiB cap; it is not a filesystem quota. Tiny control tests verified
descendant termination on memory, wall, disk and CPU-affinity violations. Their
intentional nonzero exits are retained and are not biological-cell failures.

Retained issues and corrections:

- The combined environment solve failed on incompatible BWA/Python zlib
  requirements. Separate new `env` and `env-bwa` prefixes retained the registered
  versions; existing research environments were untouched. Python 3.10.22,
  AnnData 0.11.4, MuData 0.3.10, NumPy 2.2.6, SciPy 1.15.3 and pandas 2.3.3 were
  used. Module origins, full versions, `pip check` and all 87 installed circyto
  Python files were checked against the reviewed wheel.
- The first FASTQ helper did not count `/1` suffixes in the second header field.
  Its record was retained; the corrected helper rescanned the complete file.
- BWA stdout uses a Python copy pipe. Its wrapper did not directly hash the
  stdout stream, and the original record retains that null hash. The tested
  binding uses the fresh one-cell manifest, BWA exit 0, CIRI3's pre-execution SAM
  hash and unchanged SAM bytes. This limitation is explicit in the binding record.
- The first workflow value check rejected int32 versus int64 storage widths.
  That failed stage/log and all its audit outputs remain in
  `records/acceptance_attempt1`. The corrected read-only comparison requires
  exact integer equality after lossless widening; floats and altered values
  fail. Native-import/RNA round-trip checks remain strict.
- MuData emitted future-default warnings, retained in logs. No critical check
  was skipped. The final analysis-helper suite passed 17 tests. Counts of tests
  are engineering checks, not biological evidence or comparative performance.
- The workflow status defect remains **unfixed**. Source HTML, JSON, TSV and
  AnnData still contain the observed discrepancy, with four explicit audit FAILs.

## Provenance, files and next decision

The named historical `work_hg38_full23_v2` directory was unavailable and no
authorized remote location was supplied. No new broad search or SSH connection
was attempted. All new execution records are **Route B validation**, not recovered
historical Methods or completion states. The six-file registration and protected
recheck paths/bytes were compared with their preserved manifests; the two original
H5MU hashes remained unchanged. Circulated manuscripts, figures, Conclusions,
original branches and data were not edited.

[commands.json](commands.json), [input_manifest.json](input_manifest.json),
[environment_manifest.json](environment_manifest.json),
[output_artifact_manifest.json](output_artifact_manifest.json),
[native_value_checks.tsv](native_value_checks.tsv) and
[reproduce_first_cell.md](reproduce_first_cell.md) provide exact commands,
versions, evidence classes and hashes. Large inputs, JAR, references, indices,
SAM and scientific objects remain at the durable root. Both offline previews
remain under `records/acceptance/browser_*`; no manuscript figure was generated.

The single next decision is whether to authorize a narrowly scoped repair and
revalidation of workflow alignment-status propagation. Review that defect before
authorizing the remaining three cells. This session establishes neither complete
four-cell acceptance, biological detection accuracy, independent-laboratory
validation, all-version compatibility nor historical provenance recovery.

**FIRST-CELL TECHNICAL ACCEPTANCE FAILED — workflow alignment-status propagation;
all evidence retained, remaining three not started.**
