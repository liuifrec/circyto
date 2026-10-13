# External-detector packaging validation

Base: `review/pre-release-integration@f841c296db8a6fda870eebc5b553eab10095fa0d`.
Branch: `fix/external-detector-distribution`. Version remains 0.10.0.
Obtain this record's commit SHA with
`git log -1 --format=%H -- docs/detector_distribution_validation.md`.
PR #19, its head, all original branch tips and all 72 protected manuscript
recheck entries remain unchanged. No merge, tag, publication or manuscript edit.

## Findings and changes

The upstream permission/source evidence and precise artifact boundary are in
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md). No new permission grant was
obtained. CIRI3's GPL source availability is established, but a complete
corresponding build/embedded notice review is not. CIRI-full v2.0 JAR/manual
permission remains unresolved. The original two binary hashes are retained in
the audit; files in the source checkout are not deleted or modified.

Both Python distribution formats exclude Java JAR/class files and the CIRI-full
manual. Retained source tools: two shell adapters, CIRI2 and CIRI-AS Perl scripts,
their two manuals and two GPL texts, plus CIRI3 README/license (10 tool assets).
A setuptools build guard refuses a contaminated build tree instead of emitting
an apparently clean wheel from stale files. It does not delete user files.
Build from a clean tree after a refusal. Tests also rebuild a wheel from the
sdist and check that adapter/source/license assets remain present.

External JAR paths and HOME discovery are unchanged. Doctor and runtime errors
now explain that Java detectors are external and identify configuration
variables. The reviewed Java requirement is relabeled accurately. CIRI-full's
existing single-end fallback still preflights the JAR; use the standalone CIRI2
adapter if no permitted CIRI-full JAR is available. No detection algorithm,
threshold or runtime command construction was changed.

## Local tests and installed acceptance

- Full default suite: **480 passed, 10 skipped, 9 warnings**, 21.16 s on Python
  3.10.20, with `CIRCYTO_SKIP_INTEGRATION=1`. Focused packaging/path/doctor tests:
  **24 passed**. Existing six detector integration skips, two Scanpy skips and
  two missing-MuData negative-test skips remain explicit. All import/MuData
  success paths execute.
- New artifact tests build a direct wheel, sdist and sdist-derived wheel and
  inspect members. Stale `.jar`, `.class` and CIRI-full manual injections all
  produce build errors. This is a package-content test, not detector evidence.
- Fresh venv without system-site packages installed the actual wheel from `/tmp`,
  outside the repository. `pip check`, top-level/import/report/workflow help,
  doctor and detector listing pass. No JAR/class/manual exists in the installed
  package. Packaged shell adapters and licensed Perl source resolve. Missing
  detectors produce actionable errors while core doctor readiness succeeds.
- Explicit external CIRI3 JAR, CIRI3 HOME, Java override and CIRI-full JAR paths
  resolve correctly alongside the installed packaged adapter. Real Java 17.0.18
  and CIRI3 `-H` execute. CIRI-full discovery is tested with the existing local
  historical binary; **no CIRI-full execution or permission clearance is claimed**.
- The manifest-driven full-length workflow dry-run uses the configured external
  Java/JAR paths in its command preview. Inputs are labeled synthetic planning
  fixtures; no BWA or STAR execution occurred in that check.
- Detector-independent native import still reads the authentic official result
  into a 1 × 6 CSR object with support 50 and unknown completion. Two imports
  preserve original evidence and repeat text/MTX outputs; standard AnnData
  loading passes. The separate synthetic matched-RNA demo passes MuData loading
  and value checks. Neither establishes real per-cell + matched RNA acceptance.

## Actual CIRI3 adapter execution

The installed wheel's existing `Ciri3Detector` executed the external JAR with
one thread, existing default `-S 0 -Ma 0`, no supplied GTF, and copies of the
original upstream demo SAM/reference. This is real process execution, not a
mock, and writes exclusively to a temporary output directory.

| Input | SHA-256 |
| --- | --- |
| CIRI3 JAR | `569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0` |
| Official demo SAM | `4872d2cc5c73f95f18cf404cebc6f9cc4941106a5eb64e6c68c2a6f4a9e986c9` |
| Official reference fragment | `0223fb4957af56602f37de9d8c09c3a3cc21879c9becc366ad36bb164927f97b` |

The detector returned **six candidates, support 50**, in approximately 0.65 s
(includes adapter/process startup; warm filesystem, not a scaling benchmark).
An independent CSV comparison matched all candidate IDs, chromosomes, starts,
ends, strands and junction support to the authentic original upstream result
SHA-256 `4885fc00bd30817ebdf81ba9f3344c48ddb18bcb56d2a077cee2297ff7f594ed`.
Normalization produced six rows. No biological single-cell provenance or genome
assembly is inferred for this demo, and this does not validate full FASTQ
alignment workflows or biological sensitivity/specificity. Established manuscript
runs were not repeated. Original input data and historical outputs are untouched.

To reproduce the adapter runtime check after external installation, using copies
of the pinned official `data/circRNA/Single/sample.sam` and `data/circRNA/ref.fa`:

```python
from pathlib import Path
from circyto.detectors.base import DetectorRunInputs
from circyto.detectors.ciri3 import Ciri3Detector
result = Ciri3Detector().run(DetectorRunInputs(
    cell_id="official_demo", outdir=Path("new_runtime_output"),
    sam=Path("sample.sam").resolve(), input_mode="alignment",
    ref_fa=Path("ref.fa").resolve(), threads=1,
))
print(result.tsv_path)
```

## Remote CI evidence

[Workflow run 38355386066](https://github.com/liuifrec/circyto/actions/runs/38355386066)
passed on implementation commit `c9c7ac68c5b2956c38da7c1bbfd83dad7520678b`:
Python 3.10 and 3.12 each ran **480 passing tests, 10 explicit skips**. Both
passed the critical-module execution audit, source/wheel builds, fresh isolated
wheel installation, missing-detector/package-content checks and authentic native
import. Logs and retained XML/acceptance artifacts were retrieved and inspected.
This is separate from the local actual external CIRI3 process execution above;
remote CI did not download or execute CIRI3. Later documentation commits do not
change detector behavior; subsequent run results are identified separately.

## Decision and stopping boundary

**Engineering GO:** newly built wheel/sdist exclude unresolved Java distribution,
preserve adapters and detector-independent import, and run the reviewed external
CIRI3 on authentic demo input. The CI workflow includes artifact regressions,
installed-package absence/discovery checks and the importer execution audit.
Local logs and machine records are under `/tmp/circyto-distribution`; remote CI,
if run, is reported separately with its exact tested head.

**NO-GO for blanket public-release/legal clearance.** Redistribution of the
excluded JARs/manual, historical Git/release assets, and any complete upstream
license review require rights-holder/maintainer/institutional decisions. No
additional engineering work is proposed to infer those permissions. Existing
real-cell/RNA acceptance, original workflow report acceptance and public-data
archive gates from PR #19 also remain open. Stop this milestone here; RERF
internal review remains independent.
