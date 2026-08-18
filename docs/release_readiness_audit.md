# Release-readiness audit

Audit date: 2026-08-18

Scope: GitHub issue #16, circyto 0.10.0

Branch: `feature/nanopore-feasibility`

## Verdict

**Ready to freeze after the stacked PRs are reviewed and the release environment is provisioned with the external tools required by its intended workflows.** The Python package, source distribution, wheel, resource lookup, and clean-install CLI smoke checks pass. Missing workflow executables and this host's Java 8 are reported as optional readiness gaps, not package failures.

## Checks performed

- Confirmed the initial worktree was clean. The current branch is the head of draft PR #15, whose base is `feature/biogenesis-foundation` (draft PR #13); no unrelated merge or rebase was performed.
- Exercised the installed editable CLI, `doctor`, `detectors`, and help for CIRI3, alignment-first, collection/conversion, MuData, SComatic, scRR/RT, Nanopore, CIRI-long, and biogenesis commands.
- Generated [the machine-readable CLI inventory](cli_inventory.json) from the installed command tree: 90 command paths, including 83 current commands, two compatibility aliases, and five legacy commands. No duplicate purposes or documented-but-missing commands were confirmed.
- Audited imports, dependency groups, package data, entry points, runtime resources, documentation wording, version declarations, warnings, tracked artifacts, absolute paths, and release-impact TODO/FIXME items.
- Ran the full test, build, compile, diff, artifact-content, and clean-wheel checks listed below. CRR194209 was not run.

## Defects fixed

- Wheel installations could not locate the repository-bundled CIRI3/CIRI-full/CIRI2 runtime assets. Required jars, adapters, scripts, manuals, and licenses are now package data, and path resolution falls back to `circyto/resources/tools` outside a source checkout.
- `doctor` conflated missing optional tools with package failure and omitted core-package, minimap2, CIRI-long, protocol-boundary, and useful Java checks. Its output now separates core readiness, packaged assets, optional executables, detector readiness, and protocol limitations; remediation names each missing tool and its purpose. It remains hermetic and performs no network checks.
- `detectors` omitted the supported CIRI2 adapter. CIRI2 now has an isolated readiness entry.
- The build extra was undeclared, `pysam` was not exposed as an optional BAM dependency, unused `tqdm` burdened the base environment, and the setuptools license form was deprecated. These declarations were corrected without a broad dependency refactor.
- README first impressions did not clearly separate short-read detection, generic Nanopore alignment/QC/provenance, and CIRI-long RCRT support. The capability boundaries, minimal quick start, readiness check, and output object are now visible early.
- Workstation-specific SComatic and validation paths were replaced with environment variables or portable placeholders. The stale current-release label in `ROADMAP.md` was corrected from 0.9 to 0.10.
- Removed 140 tracked generated/snapshot files (231.0 MiB), including H5MU outputs, detector logs/outputs, BWA indexes, PDFs, and repository/CLI snapshots; ignore rules now cover them. Intended compressed FASTQ fixtures and licensed runtime detector assets remain.

Focused hermetic regressions cover optional-versus-core doctor status, detector isolation, and packaged-resource fallback/source labeling.

## Packaging and clean installation

`python -m build` produced:

- `circyto-0.10.0-py3-none-any.whl` (2.5 MiB, 103 entries)
- `circyto-0.10.0.tar.gz` (2.5 MiB, 204 entries)

Both contain the Nanopore expectation profile, smoke resources, CIRI3 jar/license/readme, and the CIRI-full/CIRI2 adapters and runtime files. Neither contains generated H5AD/H5MU/BAM/SAM outputs, BWA indexes, prior detector-output directories, or large test data.

A new environment under `/tmp`, outside the repository, installed only the final wheel and its declared dependencies. `pip check` passed. `circyto --version`, `--help`, `doctor`, `detectors`, `nanopore --help`, and `ciri-long --help` all exited zero. Imports resolved from the environment's `site-packages`; repository discovery returned `None`; the packaged CIRI3/CIRI-full/CIRI2 assets and the SRR4048177 expectation profile loaded successfully.

## Dependencies, warnings, and versions

- Core imports are covered by the base dependencies. Guarded integrations remain optional; `pysam` is available through the `bam` extra. `build` is in the development extra. No clear undeclared Python runtime dependency remains.
- Test result: **337 passed, 8 skipped, 131 warnings**. Of the warnings, 126 are MuData 0.4 behavior-change `FutureWarning`s, four are expected overlapping-variable-name warnings from manuscript fixtures, and one is an AnnData string-index conversion warning. No actionable circyto warning was found.
- MuData's future `pull_on_update=False` default was not adopted or globally suppressed because doing so could change AnnData/MuData synchronization semantics. It needs a dedicated compatibility task and explicit semantic regression tests.
- Version declarations are consistent at 0.10.0 across `pyproject.toml`, `circyto.__version__`, CLI output, README/current docs, and generated provenance. No version bump was made.

## Intentionally deferred

- Fifteen advanced, compatibility, or legacy command paths are discoverable in CLI help but lack dedicated prose documentation. They remain recorded in the inventory; no public command was renamed.
- Installing Java 17 and workflow-specific BWA, Bowtie2, STAR, minimap2, find-circ3, or CIRI-long is an environment concern. This audit did not install or execute those external workflows.
- Consolidating source-tree and packaged copies of licensed detector assets may reduce maintenance duplication later, but is not required for release correctness.
- MuData 0.4 behavior adoption is deferred as described above.
- `prepare-public-dataset` retains explicitly labeled `TODO_REAL_ACCESSION` fallback rows for planning-only targets. The command and its documentation already prevent those placeholders from being presented as validated accessions.

## Final validation

- `python -m pytest -q .`: 337 passed, 8 skipped
- `python -m build`: passed
- `git diff --check`: passed
- `python -m compileall -q circyto scripts`: passed
- Final isolated wheel smoke and resource load: passed

Recommended final small task before the weekly freeze: review and merge the PR #13/#15 stack in order, then run `circyto doctor` in the intended release execution environment with Java 17 and only the external detector tools that environment is expected to support.
