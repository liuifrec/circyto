# Native CIRI3 acceptance record

Base: `5bc38f293e66d4051c1badac254f4501f6983acf`.
Branch: `feature/native-ciri3-import`. Version stays 0.10.0.
This record describes one native CIRI3 import milestone, not a release.

## Authentic input and version

The unchanged official `data/circRNA/Single/result.txt` is copied under
[`examples/native_ciri3`](../examples/native_ciri3). It is 1,913 bytes, six
candidates, 50 BSJ support, SHA-256
`4885fc00bd30817ebdf81ba9f3344c48ddb18bcb56d2a077cee2297ff7f594ed`.
The original 299-byte run log is SHA-256
`f2aabca9870f621c3e1c520b6ade555677ec57d8edbb0c33cfd31ef224464b6f`.
Both are byte-identical to the official v3.0.1 tree
`73108c4478a7e691cb93360921b65203e0701516`. The result first appeared in
`27031610c234e7459bc44fb9bd9bded2050599b5` (2023-08-16), before that release.
Its producing version is therefore recorded as **unknown**, not 3.0.1.
The tag establishes a reviewed format snapshot, not the producing executable.
The `Java_18.0.1` JAR suffix describes a Java target, not a CIRI3 version.

The official repository retains the corresponding sample SAM, reference fragment,
result and log. Independently, all 50 listed junction-read IDs were found in the
original SAM, and each native count equals its unique listed read-ID count.
The SAM SHA-256 is
`4872d2cc5c73f95f18cf404cebc6f9cc4941106a5eb64e6c68c2a6f4a9e986c9`.
This checks correspondence, not detector correctness. No CIRI3 or alignment was rerun for this milestone. No H5MU-derived
native table or normalized circyto fixture was substituted. The project history's
`tests/fixtures/detectors/ciri3_raw.tsv` at `0645b68` has a different seven-column
mock layout and was explicitly excluded from authentic acceptance. The earlier
Smart-seq3/IMR90 H5MU round trips remain distinct evidence.

The supplied sample reference is a 10-Mb chr1 fragment and lacks an established
assembly/accession in the inspected evidence. It is labeled
`upstream_demo_chr1_unverified`, and the observation is a demonstration sample,
not asserted to be a biological single cell. This is native-output compatibility
validation, not new biological validation or a sequencing performance benchmark.

A second unchanged official output, `data/DE/BSJ/With_RE/case1_CIRI3.txt`, is
included as `annotated_simulation.txt` to test available native annotations.
Its read IDs explicitly say `simulate:`; it is labeled as simulation-derived
annotation-contract evidence. It is not counted as biological acceptance.
Hashes, sizes, original member paths and upstream license are provided in the
[fixture provenance](../examples/native_ciri3/provenance.json).

## Reproduce without a detector

Install the wheel (with `[mudata]` only if using RNA). Copy the example directory
and [`examples/validate_native_ciri3.py`](../examples/validate_native_ciri3.py)
outside the checkout. Then run:

```bash
python validate_native_ciri3.py --fixtures native_ciri3 --out /tmp/native-acceptance
```

The Linux audit uses GNU `/usr/bin/time`; importing itself has no GNU-time
requirement. It checks fixture hashes before ordinary CSV reads, runs the
installed CLI twice, uses standard AnnData/MuData readers and checks every native
ID, coordinate, support, chromosome, strand, native annotation and status.
It also checks copied evidence, object provenance and repeatability. Optional
`--rna existing.h5ad` exercises exact-ID alignment to an existing RNA object
containing `upstream_single`; it does not imply biologically matched RNA.

## Installed acceptance and performance

The full default CI-style suite passes: **474 passed, 10 skipped, 9 warnings**
in 19.71 seconds (`CIRCYTO_SKIP_INTEGRATION=1 python -m pytest -q .`).
Focused native/general importer regressions pass: **72 passed**, including 33
native-interface cases. External detector executions remain gated. Warnings are
existing AnnData/MuData warnings and the expected ordinary MuData reader
FutureWarnings; none were suppressed.

New checks cover authentic native values, explicit mapping/version/scope,
normalized-versus-native rejection, unsupported headers and malformed rows,
duplicates, conflicting coordinates/strands/annotations/versions, exact support,
header-only zero results, missing/failed/skipped/unprocessed/unknown cells,
metadata ordering and optional RNA integration. A guard makes sparse
`toarray`/`todense` raise during import; native conversion/export still passes.
Code inspection confirms coordinate-format CSR construction with O(stored calls)
data and one-dimensional per-cell/per-candidate summaries, not a dense
cell-by-candidate array. Long tables/read-ID strings are still memory-resident.

A fresh venv without system-site packages installed the built wheel and its
`[mudata]` dependencies. Imports resolve to that venv's `site-packages` from
`/tmp/circyto-native/independent`, outside the repository. `pip check` passes.
Top-level help/version, `import-circ --help` (including the new option), and
existing collector/export help pass. No detector binary or raw analysis is
needed for the example. The installed stack is Python 3.10.20, AnnData 0.11.4,
MuData 0.3.10, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3 and Typer 0.27.3.

Standard `anndata.read_h5ad` preserves all six candidate IDs and coordinates,
strands, native gene/type fields and `[9, 9, 5, 6, 14, 7]` support values, with a
CSR int64 matrix. Status remains unknown and detection fractions remain
unavailable. The raw native file is retained byte-for-byte separately from the
converted matrix, with its checksum and explicit sample mapping in provenance.
Repeating identical input gives identical JSON, TSV, HTML and Matrix Market
bytes; HDF5 byte identity is not promised.

Optional `mudata.read_h5mu` acceptance used an existing **synthetic RNA fixture**
from the preceding import milestone. For this technical test only, its `cell_A`
label was mapped to `upstream_single` and its synthetic genome label removed;
RNA values remained unchanged. The saved RNA object had three cells and one
feature, one shared observation and two RNA-only observations. RNA values,
metadata, ordering and circ support survive; RNA-only rows retain unknown circ
status. The original fixture was not overwritten. Unit tests separately cover
RNA layers. This is not cross-study biological pairing or matched-RNA evidence.

Measured with GNU `/usr/bin/time`, two fresh CLI processes per configuration,
warm filesystem cache, on Linux x86_64 (kernel 6.8.0-138-generic, glibc 2.35):

| Input and output | Native input | Run 1 wall / peak RSS | Run 2 wall / peak RSS |
| --- | --- | --- | --- |
| Official native demo → AnnData + QC | 1,913 bytes; 1 × 6; 6 nonzeros; support 50 | 0.51 s / 120,888 KiB | 0.48 s / 120,996 KiB |
| Same native demo + technical RNA → MuData + AnnData + QC | same native data; RNA fixture 3 × 1 | 0.62 s / 124,272 KiB | 0.59 s / 123,696 KiB |

These include interpreter/library startup and file I/O. They are measurements
of the actual supplied native input, not biological benchmarking, detector
performance or evidence of large-cohort throughput. No synthetic scaling result
is presented as biological evidence. Machine-readable installed acceptance
records are [AnnData](examples/native_ciri3/acceptance.json) and
[technical MuData](examples/native_ciri3/rna_acceptance.json).

The final installed HTML was rendered and visually checked at 1440 × 1300 in
headless Chrome with hostname resolution disabled. It clearly shows six
candidates, support 50, unknown completion and unavailable detection fraction.
A small import-report style correction keeps coordinates/identifiers on one
line with horizontal scrolling; it changes no scientific value. All report
assets are local; no script/CDN dependency is introduced.

Wheel/sdist builds pass with `python -m build --no-isolation`. Archive checks
confirm the adapter in the wheel and complete examples, unchanged native files,
original log, upstream license and guides in the sdist. Neither large SAM/H5MU
inputs nor protected manuscript files are packaged. Focused Ruff checks
(ignoring B008 for Typer's Option convention) and `git diff --check` pass.
Local logs and browser evidence are under `/tmp/circyto-native/`.

The current milestone commit is the commit introducing this record; obtain its
full SHA with `git log -1 --format=%H -- docs/native_ciri3_validation.md`.
The exact pushed SHA is also reported in the delivery response.

## Limits and remaining gaps

- Scope is the exact 13-column single-sample result table, reviewed at v3.0.1;
  native producer version can be explicitly unknown. No claim that all outputs
  or versions named CIRI3 share this format.
- The pooled information table can share the same header. The manifest's
  single-sample assertion must be established by the user from upstream run
  organization; a falsely labeled pooled table cannot be identified by syntax.
- Multi-sample BSJ/FSJ matrices, 12-column variants, user-candidate quantification,
  enrichment/ratio outputs and other known versions remain unsupported here.
- File/log existence is insufficient for independently bound per-cell completion.
  The official example remains unknown; detection fractions are unavailable.
  Status does not change counts or establish locus-level callability/absence.
- Real matched RNA and independent biological single-cell CIRI3 output from
  another laboratory are not supplied by this demo. Optional RNA acceptance is
  a technical object-alignment test; no RNA–circ biological conclusion follows.
- Gene IDs are retained as IDs and do not become gene symbols or invented GTF
  provenance. Conflicting annotation across files is rejected rather than merged.
- Sparse matrices are used throughout, but native/normalized long tables and
  read-ID fields are read in memory. This is not a streaming or out-of-core
  importer. Tiny-input measurements do not predict large-cohort scaling.

Only this adapter, associated integration, tests and documentation are changed.
No detector execution code, algorithm, threshold, original dataset, manuscript
claim or package version changes. No additional detector parser or model work.
Original protected branches and manuscript recheck files are verified separately
before committing; the feature branch alone is pushed, with no merge or release.
