# Pip and conda release preparation

Status: **PREPARED FOR OWNER REVIEW**; no publication authorized or performed.
Proposed version: **0.11.0rc1**, pending owner approval. Only disposable local
build copies change version; `pyproject.toml` and `circyto.__version__` in Git
retain 0.10.0. This does not change the frozen manuscript release or move a tag.

Source lineage: `fix/external-detector-distribution` at
`6c280426ebf10abcbcaeab72f87c9d34060d6ee8`; PR #19 remains draft/open at
`f841c296db8a6fda870eebc5b553eab10095fa0d`. Preparation branch:
`release/pip-conda-preparation`. The tested source commit is
`48dd8d2bc45af4978f238a6fa5f31031c44549d3`, plus the two disposable version
literal substitutions below. Later commits complete this evidence record and
installation instructions; they do not change candidate runtime code. The
matching sdist therefore contains the earlier preparation instructions, not this
post-validation record. Use this record with the exact artifacts below.

## Owner publication steps (not performed)

1. Review the installation evidence, retained-source license expression/notices,
   version proposal and open scientific acceptance gaps. Approve a distinct
   release version, then update both version declarations on a reviewed release
   commit. Do not reuse 0.10.0 or upload artifacts from automatic Git archives.
2. Check account ownership/permission on PyPI and TestPyPI separately. Public
   API 404 responses do not establish name availability. Configure each account
   with its own pending/existing Trusted Publisher: owner `liuifrec`, repository
   `circyto`, workflow `release.yml`, environment `pypi` or `testpypi` respectively.
   No pasted API tokens are needed or requested.
3. Configure both GitHub environments with required reviewers, prevent self-review
   and restrict deployment to main. Protect main/release workflow edits. These
   account/environment actions are **not supplied by YAML**. Leave repository
   variable `CIRCYTO_PUBLISH_ENABLED` unset until that configuration is reviewed.
4. The prepared workflow runs only by manual dispatch on main, checks the exact
   approved SHA/version, runs full tests and fresh installs, and retains artifacts
   for review before its environment-gated upload job. Only that job receives
   `id-token: write`; it has no checkout or build. It cannot publish from this
   preparation branch, arbitrary pushes or PR events. All actions are SHA-pinned.
5. Obtain separate approval even for TestPyPI, because it is a public upload.
   After an approved upload, create a fresh environment, preinstall reviewed
   dependencies from the intended index (normally PyPI), then install the exact
   downloaded TestPyPI wheel with `--no-deps --no-index`, checking the JSON API's
   filename/hash. Do not resolve dependencies using mixed public/test indices.
   Run `pip check` and the same base/multimodal installed acceptance. A TestPyPI
   result is not PyPI publication; inspect PyPI's exact files/version separately.
6. For Bioconda, replace the local recipe source variables with the final stable
   **clean PyPI sdist** URL and literal SHA-256, update version, rerun recipe lint,
   build and mulled/container acceptance under Bioconda's then-current pins.
   No GitHub auto-generated archive: it can include excluded binaries/manuals.
   Only after owner approval submit the recipe PR. Local conda-build success is
   not Bioconda acceptance or publication. No Anaconda upload/channel is created.

The setuptools backend, CLI entry point and optional extras remain in place.
The aggregate metadata is `MIT AND GPL-2.0-only`, with individual license texts
included; it does not relicense circyto's MIT code or grant permission for excluded
assets. JAR/class files and the unresolved CIRI-full manual stay excluded in all
artifacts. Historical Git assets and upstream licensing questions retain the
boundaries in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md). This engineering
preparation is not legal clearance.

## Official guidance inspected on 2026-10-15

- [PyPA packaging tutorial](https://packaging.python.org/en/latest/tutorials/packaging-projects/):
  isolated builds, wheel/sdist distinction and installation verification.
- [PyPA license expression](https://packaging.python.org/en/latest/specifications/license-expression/):
  distribution metadata must account for shipped components.
- [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/):
  job-scoped OIDC and separate TestPyPI configuration.
- [PyPI name retention](https://docs.pypi.org/project-management/name-retention/):
  missing project pages are not proof of claimable names.
- [Bioconda recipe guidelines](https://bioconda.github.io/contributor/guidelines.html):
  stable checksummed source, noarch Python and separation of tests requiring
  optional dependencies from mulled runtime-only commands.
- [Bioconda local testing](https://bioconda.github.io/contributor/building-locally.html):
  local builds and container/mulled validation are distinct steps.

## Executed validation

Machine-readable results: [release_preparation_evidence.json](release_preparation_evidence.json).
All artifacts, complete archive inventories, dependency records, JUnit XML,
acceptance outputs and logs remain local under `/tmp/circyto-install-lab`.
Nothing in that directory has been uploaded. Retain that directory privately
for review; `/tmp` is not a durable release archive.

### Source, version and artifacts

The disposable `candidate-source` tree was exported from the exact commit above.
Only `version = "0.10.0"` in `pyproject.toml` and `__version__ = "0.10.0"` in
`circyto/__init__.py` were changed to `0.11.0rc1` before building. Their modified
file hashes are in the machine record. No old tag or historical version changed.
Builds used setuptools' existing PEP 517 backend in isolation; `python -m twine
check --strict` passed for the wheel, sdist and rebuilt wheel.

Paths below are relative to `/tmp/circyto-install-lab`:

| Artifact | SHA-256 | Inspection |
| --- | --- | --- |
| `dist/circyto-0.11.0rc1-py3-none-any.whl` | `42a9f6a63f0d9a9761cc1872744701a2a7d639ddd705c93a000b1c89aa9d04c2` | 113 members; installed directly |
| `dist/circyto-0.11.0rc1.tar.gz` | `99c46f476c402ad6894cad00f5b6af49b48428c6edee00104fb5476d20aa813f` | 273 members; clean conda source; examples included |
| `rebuilt/circyto-0.11.0rc1-py3-none-any.whl` | `dd7ecbb1b9636872402bf05a4780e8d86a9d2856afba8e347bfdfda54b292471` | 113 members; all payload bytes identical to direct wheel |
| `conda-artifacts/noarch/circyto-0.11.0rc1-py_0.conda` | `fa0741c4f1e56de1cc3d80c3911d72c75d92b7506ebe61c698989c710064311a` | 145 members across package/info components; installed independently |

All archives and the extracted sdist have **zero JARs, classes or excluded
CIRI-full manuals**. Both adapters, CIRI2/CIRI-AS source/manuals/licenses and
CIRI3 README/license remain present. Python metadata and conda license records
include the aggregate SPDX expression and individual notices. No compiled
extension or platform-specific executable is installed by the recipe. Noarch
describes the Python distribution plus text resources; shell/Perl adapters and
external Java/aligners still need their own runtimes/platform checks.

Build commands (fresh directories only):

```bash
python -m build --sdist --wheel --outdir /tmp/circyto-install-lab/dist
# Extract the audited sdist into a new directory before this second build:
python -m build --wheel /tmp/circyto-install-lab/from-sdist/circyto-0.11.0rc1 \
  --outdir /tmp/circyto-install-lab/rebuilt
python -m twine check --strict /tmp/circyto-install-lab/dist/* \
  /tmp/circyto-install-lab/rebuilt/*.whl
```

The first command's working directory was the disposable versioned source tree.
The installed example files were copied from the extracted audited sdist into
`/tmp/circyto-install-lab/examples`, then run with isolated Python from outside
every source checkout. No editable install or source-tree PYTHONPATH was used.

### Installed acceptance and dependency observations

| Route | Python | AnnData / MuData | NumPy / SciPy / pandas | Result |
| --- | --- | --- | --- | --- |
| Direct wheel, base | 3.10.22 | 0.11.4 / absent | 2.2.6 / 1.15.3 / 2.3.3 | PASS, including actionable missing-MuData error |
| Direct wheel, `[mudata]` | 3.10.22 | 0.11.4 / 0.3.10 | 2.2.6 / 1.15.3 / 2.3.3 | PASS |
| Direct wheel, base | 3.12.15 | 0.13.4 / absent | 2.5.4 / 1.18.1 / 3.0.6 | PASS, including actionable missing-MuData error |
| Direct wheel, `[mudata]` | 3.12.15 | 0.13.4 / 0.4.1 | 2.5.4 / 1.18.1 / 3.0.6 | PASS |
| Conda environment + exact wheel | 3.12.15 | 0.13.4 / 0.4.1 | 2.5.4 / 1.18.1 / 3.0.6 | PASS; 96 conda dependencies, wheel installed `--no-deps --no-index` |
| Native conda artifact, base | 3.12.15 | 0.13.4 / absent | 2.5.4 / 1.18.1 / 3.0.6 | PASS; package + 95 conda dependencies |
| Native conda artifact + MuData | 3.12.15 | 0.13.4 / 0.4.1 | 2.5.4 / 1.18.1 / 3.0.6 | PASS; package + 96 conda dependencies |

Typer 0.27.3 and Rich 15.0.0 were used throughout. All seven installs passed
`pip check`; module origins resolve under each fresh environment's `site-packages`.
User-site/PYTHONPATH were isolated. PATH contained only the selected environment's
bin directory for installed acceptance: Java, BWA and STAR were absent, core
doctor succeeded, and optional detector readiness was reported missing.

The authentic upstream 13-column example loaded as CSR int64, one observation by
six candidates, support **50**, unchanged axes/coordinates/annotations, and
**unknown** processing status. Its source SHA-256 remains
`4885fc00bd30817ebdf81ba9f3344c48ddb18bcb56d2a077cee2297ff7f594ed`.
With MuData, an explicitly synthetic 3-by-1 RNA matrix (values 3, 4, 9) supplied
one shared ID and two RNA-only IDs. Standard MuData loading preserved both
modalities' values and axes and RNA annotations. Base installs report
`Optional RNA export requires circyto[mudata]`; conda users add `mudata>=0.3`
to their solve as documented, rather than treating pip extras as conda syntax.

Offline HTML rendered in headless Chrome from a local file with DNS blocked.
The report displays six candidates/support 50, unknown completion and the warning
that import success is not detector completion. Standalone workflow-report
regeneration correctly refuses this IMPORT bundle. Authentic format acceptance
does not establish a biological single cell, synthetic RNA is not real matched
RNA, and original completed-workflow reporting remains untested here.

The actual conda build used **conda-build 26.9.1**. Its host and recipe-test solve
selected **Python 3.15.0** under the recipe's `python >=3.10`; the separate fresh
artifact installs above explicitly selected 3.12.15. A requested 3.12 variant
file did not pin the explicit host constraint. Recipe-test MuData loading passed
on 3.15, but the full suite was run only on 3.10/3.12. Record the observed solve,
not an inferred version. The recipe-test process could see system Java; the
seven separate installed checks above explicitly excluded it from PATH.

All conda dependency URLs were verified from `conda-meta`: conda-forge supplied
the resolved dependencies, with bioconda second and strict channel priority.
circyto itself came from the exact local artifact. The recipe ran
`pip install . --no-deps --no-build-isolation --no-index`; its declared host
dependencies supplied pip/setuptools/wheel. It did not fetch Python runtime
dependencies with pip. The conda environment YAML is a separate installation
route and is not described as a package.

The recipe source URL was the local `file:///tmp/circyto-install-lab/dist/`
sdist with its exact checksum above. An environment-specific condarc, cache and
build root were under `/tmp/circyto-install-lab`; the user's global condarc and
research environments were not modified. The build command used:

```bash
export CIRCYTO_SDIST_URL=file:///tmp/circyto-install-lab/dist/circyto-0.11.0rc1.tar.gz
export CIRCYTO_SDIST_SHA256=99c46f476c402ad6894cad00f5b6af49b48428c6edee00104fb5476d20aa813f
CONDARC=/tmp/circyto-install-lab/condarc \
CONDA_PKGS_DIRS=/tmp/circyto-install-lab/mamba/pkgs \
CONDA_ENVS_PATH=/tmp/circyto-install-lab/build-envs CONDA_CHANNEL_PRIORITY=strict \
  conda build packaging/conda --no-anaconda-upload --package-format 2 \
  --override-channels -c conda-forge -c bioconda \
  --croot /tmp/circyto-install-lab/conda-build \
  --output-folder /tmp/circyto-install-lab/conda-artifacts \
  -e /tmp/circyto-install-lab/variants.yaml
```

No Bioconda lint/container/mulled run or recipe submission is claimed. The
environment-dependent local source URL must become a stable checksummed release
sdist URL before submission. The builder emitted a parsed-output-block warning;
the inspected rendered metadata had one intended `circyto` output, correct
runtime requirements, optional `mudata >=0.3` constraint and noarch entry point.

Full resolved dependency records are local `pip*-dependencies.txt` and
`conda-*-dependencies.json`, including conda package URLs/builds/hashes. These
are observed environments, not promises about every combination admitted by
the existing lower bounds. No runtime dependency bounds or extras were changed.
Only Linux x86_64 was tested. macOS/Windows and external detector portability
remain outside this execution evidence; [Bioconda's supported platforms](https://bioconda.github.io/)
do not include Windows.

### Full suite, focused checks and preservation

- Python 3.10.22: **480 passed, 10 skipped, 9 warnings**, 24.44 s.
- Python 3.12.15: **480 passed, 10 skipped, 55 warnings**, 26.19 s.
- Critical-module JUnit audit passed on both: all import/MuData/report success
  paths executed. The 10 explicit skips are six external integrations, two
  optional Scanpy cases, and two missing-MuData negative tests (MuData installed).
  The separate base-installed acceptance exercised missing-MuData errors.
- Four packaging tests passed, including aggregate license metadata, retained
  source notices, sdist examples and refusal of stale JAR/class/manual assets.
  Direct wheel, sdist-derived wheel and conda archive inventories passed too.
- Release workflow passed actionlint 1.7.12 and explicit trigger/permission/gate
  checks. Shell syntax, Python compilation and `git diff --check` passed. The
  workflow was **not dispatched**, and no GitHub publishing environment, Trusted
  Publisher or enable variable was configured by this preparation.
- The first full-suite invocation lacked `python` on child-shell PATH, causing
  seven harness failures on each interpreter; retained `*-initial` logs document
  them. Prepending the fresh environment's bin directory fixed the invocation.
  No production code changed. The newer SciPy stack emits `mmread` future sparse
  array default warnings; current assertions pass, and future stacks need testing.
- Optional external CIRI3 execution was not repeated. Its prior reviewed JAR hash
  and actual runtime evidence remain in [detector distribution validation](detector_distribution_validation.md).
  This preparation neither bundles that JAR nor expands its permission status.
- The original checkout, all 61 files in `manuscript/results_v2_recheck` including
  ignored files, and existing branch tips are checked against the initial snapshot.
  No production algorithms, filters, figures, manuscript files or original branches
  were changed. Only the new preparation branch is committed/pushed.

Minimum owner action: review the proposed version, exact artifacts and notices,
then decide whether to authorize a separate TestPyPI trial and account/environment
configuration. PyPI publication and Bioconda acceptance remain later decisions.
