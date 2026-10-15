# Pip and conda release preparation

Status: local validation in progress; no publication authorized or performed.
Proposed version: **0.11.0rc1**, pending owner approval. Only disposable local
build copies change version; `pyproject.toml` and `circyto.__version__` in Git
retain 0.10.0. This does not change the frozen manuscript release or move a tag.

Source lineage: `fix/external-detector-distribution` at
`6c280426ebf10abcbcaeab72f87c9d34060d6ee8`; PR #19 remains draft/open at
`f841c296db8a6fda870eebc5b553eab10095fa0d`. Preparation branch:
`release/pip-conda-preparation`. Exact candidate source, artifact inventory,
dependency versions and results will be recorded after execution below.

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

To be completed with exact local results before handoff.
