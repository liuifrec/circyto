# Offline QC reports and external-laboratory quick start

This is the unreleased `feature/qc-report-usability` feature, based on
`main@44697355bcab1c525ca7ef9b130e2ad0094d9e1b` (v0.10.0). It does not change
the manuscript branch or its scientific software baseline. Established
high-level workflows retain their experimental CLI status. The legacy top-level
`run` command is preserved; a new manifest-based `run` interface is deferred
because overloading its existing required arguments would create a compatibility
and maintenance burden.

## Try a report without running an analysis

From this branch's repository root, after installing the package as described
in the [README](../README.md):

```bash
python examples/qc_report_demo.py --outdir work/qc_report_demo
circyto report --workdir work/qc_report_demo
```

The demo generator refuses to overwrite an existing directory. To regenerate
an existing demo report, repeat only the `circyto report` command. The generator
uses the Python standard library and writes deterministic synthetic QC tables,
a small Matrix Market matrix, and an illustrative workflow summary. It also
creates separate toy FASTQ/reference inputs for the planning example below.
No detector, alignment, annotation, or manuscript analysis is executed.

Open `work/qc_report_demo/qc/report.html` by double-clicking it or using the
browser's Open File action. The HTML contains all styling and SVG charts; it
does not fetch external resources. On GitHub, download the
[checked-in HTML](examples/qc_report/qc/report.html) instead of viewing its source.
The [sample metrics](examples/qc_report/qc/metrics.json) are machine readable.
The report is also reproducible directly from the checked-in sample directory:

```bash
circyto report --workdir docs/examples/qc_report
```

Expected synthetic values: **12 cells; 6 candidates; 77 total support; median 2
candidates/cell; 2/12 zero-count cells; 4/6 recorded host-gene coverage**. Candidate
prevalence numerators are 8, 5, 3, 1, 2, 2 cells, respectively, with denominator
12. The “Completed” example status is illustrative and is prominently identified
as synthetic; it does not document a real analysis run.

## Plan and run an established workflow

First run `circyto doctor` and inspect both core readiness and its CIRI3 entries.
The doctor lists
missing tools and asset locations. A single-end RamDA/Shin-RamDA run uses BWA,
SAMtools, Java >=12, and the configured CIRI3 jar/wrapper. Paired-end routes
add STAR and a compatible STAR genome index. Install these external tools
separately; `pip install circyto` does not install them. See
[full-length workflows](full_length_workflow.md) and the README external-tool
section for the supported setup and reference requirements.

Once BWA and SAMtools are on `PATH`, this command is an executable planning-only
example using the toy inputs created above:

```bash
circyto workflow full-length-circrna \
  --manifest work/qc_report_demo/inputs/manifest.tsv \
  --protocol ramda \
  --genome-fasta work/qc_report_demo/inputs/genome.fa \
  --gtf work/qc_report_demo/inputs/genes.gtf \
  --outdir work/qc_workflow_plan \
  --dry-run
circyto report --workdir work/qc_workflow_plan
```

This plan's report is **Dry run** with unavailable detection metrics. Automatic
reporting is skipped for dry runs. Planning still checks alignment-tool
availability; without configured CIRI3 runtime it may show that the detector
command preview is unavailable. Do not remove `--dry-run` for these toy inputs:
they are not suitable for a biological analysis.

For your real data, create a tab-separated manifest with one row per cell, such
as the following single-end RamDA example (replace the paths with real files):

```text
sample_id	fastq_1	fastq_2	protocol	strandedness	read_layout
cell_A	/absolute/path/cell_A.fastq.gz		ramda	unstranded	single
cell_B	/absolute/path/cell_B.fastq.gz		ramda	unstranded	single
```

Keep the empty `fastq_2` column for single-end reads. Cell IDs must be unique.
Use absolute FASTQ paths to avoid working-directory ambiguity; supply the
genome FASTA and matching annotation GTF for your organism/build. Validate the
manifest with `circyto manifest validate manifest.tsv --strict`, then run:

```bash
circyto workflow full-length-circrna \
  --manifest manifest.tsv --protocol ramda \
  --genome-fasta ref/genome.fa --gtf ref/genes.gtf \
  --outdir results/my_run --threads 8 --export-h5ad
```

These real-data paths are user inputs, unlike the self-contained demo above.
For paired-end RamDA/Shin-RamDA, supply both FASTQs, use `paired` read layout,
and supply `--star-index` plus the existing `--allow-paired-ramda` opt-in.
For already demultiplexed SMART-Seq3, use `--protocol smartseq3 --skip-demux`
and the required STAR index. Pooled SMART-Seq3 continues to use
`circyto workflow smartseq3-ciri3 --help`, with transcript R1/R2, I1/I2, plate
annotation, index/cell-column names, reference FASTA, and STAR index. Neither
reporting nor the help improvements extend protocol support.

## Outputs and metric definitions

All paths are relative to the workflow's `--outdir`. The original collected
matrix is `matrix/circ_counts.mtx`, in **candidate × cell** orientation. Its row
and column identities are `matrix/circ_index.txt` and `matrix/cell_index.txt`.
AnnData is **cell × candidate** at `anndata/circ_counts.h5ad` when exported.
Optional SMART-Seq3 MuData remains at `mudata/circyto_multimodal.h5mu`.
No matrix, annotation, count, threshold, normalization, or object semantics
are changed by reporting.

The report population is the rows of `qc/cell_qc.tsv`. These can include
selected zero-count cells absent from the collected matrix. Do not infer the
QC cell denominator from the matrix's nonempty column count. The report does
not load AnnData/MuData or recount the matrix.

| Metric | Definition and exact source |
| --- | --- |
| QC cells | Unique `cell_id` rows in `qc/cell_qc.tsv` |
| circRNA candidates | Unique `circ_id` rows in `qc/circ_qc.tsv` |
| Per-cell candidate counts | `cell_qc.tsv:circRNA_count`, unchanged |
| Per-cell support | `cell_qc.tsv:total_circRNA_support`, unchanged |
| Total support | Sum of `total_circRNA_support`; checked independently against sum of `circ_qc.tsv:total_support` |
| Median candidates/cell | Median of every QC cell's `circRNA_count`, including zeroes |
| Candidate prevalence | `circ_qc.tsv:n_cells_detected` divided by the QC cell count |
| Zero-count frequency | Number of QC rows with `circRNA_count == 0` divided by all QC cells |
| Recorded host-gene coverage | Nonblank, non-placeholder `host_gene` rows divided by all candidate QC rows |

Support retains the detector's count semantics; it is not automatically a UMI,
molecule, or expression measure. An `empty` detector result records no detected
candidates. A failed/missing detector result cannot establish biological
absence, even if its QC row contains zeroes. Host-gene annotation is not
orthogonal confirmation. A present but entirely blank `host_gene` column yields
0 recorded annotations; an absent column yields **Not available**. Placeholders
`NA`, `N/A`, `NaN`, `None`, `null`, `unknown`, `unassigned`, `-`, and `.` are
treated as unannotated (case-insensitive).

Missing or invalid numeric values invalidate the affected aggregate; they are
never imputed or silently excluded from a sum. An empty, valid table records
zero rows, while a missing/malformed table records `null`. Zero denominators
produce unavailable fractions. Duplicate identities or malformed rows invalidate
that table. Counts are parsed as exact nonnegative integers, including integral
decimal spellings such as `2.0`. Source-value inconsistencies are surfaced and
never repaired; prevalence is suppressed if cross-table consistency checks fail.

`qc/metrics.json` uses schema `circyto.qc_report.v1`. Scalar metric entries have
`value`, `unit`, `source`, `definition`, and `reason`; fractions also retain
integer `numerator` and `denominator`. `per_cell` and `per_candidate` include
all rows in source order. Candidate prevalence uses the shared
`candidate_prevalence_denominator`. Percentages in HTML are rounded; original
integer counts and denominators remain exact. HTML tables show at most 100 rows
to keep the document manageable; JSON and source TSVs retain all rows.

## Completion, failures, and optional data

`workflow_summary.json`, `qc/cell_qc.tsv`, and `qc/circ_qc.tsv` are the primary
inputs. When present, `ciri3/detector_run_summary.json`,
`align/alignment_prepare_summary.json`, and `rna/rna_import_summary.json` supply
additional recorded evidence. The reporter uses these fixed paths under
`--workdir`; it does not follow obsolete absolute paths in copied summaries.
Original matrix/index file presence is reported without opening or validating
their contents. A completion timestamp, usable QC, consistent populations,
completed alignment/detector evidence, and matrix/index file presence are
required for the report's **Completed** label.

Explicit failures override completed stage labels. Planned, pending, missing,
or conflicting statuses prevent a successful label. Missing summaries produce
**Unknown** unless other sources explicitly record failure. A dry-run summary
is labelled **Dry run**. Recorded status reflects available files, not a live
process monitor or a full scientific integrity audit. A historically completed
run copied without its matrices is conservatively labelled **Incomplete** in
that copy. For a deeper existing artifact check, use `circyto check-workflow`.

`circyto report` exits 0 when it successfully writes the two report artifacts,
even when the workflow itself failed or is incomplete; inspect
`metrics.json:workflow.status` or the prominent HTML status. A reporting I/O or
output-path error exits nonzero. Missing input values ordinarily produce a
diagnostic report instead of preventing it from being written.

Both high-level workflows generate the report after persisting their existing
scientific outputs and summary. `--no-report` disables automatic generation.
If reporting raises an error, the workflow command surfaces the error and exits
nonzero with an explicit statement that scientific results were retained.
Fix the report write problem (for example, directory permissions), then rerun
only `circyto report --workdir RESULTS_DIRECTORY`. Do not rerun detection merely
to regenerate HTML. If a workflow stops before its final summary, this standalone
command can still produce a diagnostic report from whatever evidence exists.

## Reproducibility, privacy, and screenshot

Reporting writes only `qc/report.html` and `qc/metrics.json`, replacing previous
report artifacts. It reads source files without modifying them. Both report
contents are prepared before writing and each destination is replaced atomically;
the pair is not a transactional filesystem update. An interruption between
replacements can leave artifacts from different attempts; regenerate both with
the same report command. Symlinked output files/directories are refused, and
source symlinks outside the workflow directory are not read.

Identical inputs and file availability produce byte-identical reports with
the same reporter code/version. Source byte lengths and SHA-256 checksums,
available run version/timestamps, protocol, detector, workflow ID, and reference
basenames are recorded. Hostnames, full paths, command lines, and environment
variables are not exported; paths and network locations in free-text warnings
are redacted. Cell/candidate IDs and annotations are retained and HTML-escaped,
so review those identifiers before sharing. This is not a patient anonymization
tool. No full reference checksum is inferred from a filename.

The [documentation screenshot](images/qc_report.png) is a browser rendering of
the checked-in synthetic HTML, not a manuscript scientific figure. To reproduce
it with installed Chrome (optional documentation tooling only):

```bash
circyto report --workdir docs/examples/qc_report
google-chrome --headless --disable-gpu --no-first-run \
  --disable-background-networking --hide-scrollbars \
  --window-size=1440,1400 --force-device-scale-factor=1 \
  --user-data-dir=/tmp/circyto-report-screenshot-profile \
  --screenshot="$PWD/docs/images/qc_report.png" \
  "file://$PWD/docs/examples/qc_report/qc/report.html"
```

The validation environment required Chrome's `--no-sandbox` switch because it
runs inside an isolated execution environment. A normal desktop install should
use its default Chrome sandbox. Font rendering may vary by OS/browser version.

## Submission positioning

The feature can be described as an optional usability extension once its own
code and artifacts are reviewed. With the manuscript software baseline frozen
at v0.10.0, the recommended scope is a **post-review usability enhancement**.
Do not attribute this report to the existing manuscript baseline, use synthetic
demonstration counts as scientific evidence, or replace any scientific figures.
See [validation results](qc_reporting_validation.md) for the tested scope.
