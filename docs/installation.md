# Installing the prepared candidate

This branch prepares an installation; it does not publish a release. On
2026-10-15 the PyPI, TestPyPI and Bioconda package APIs returned HTTP 404 for
`circyto`. No public version or owner was listed. A missing listing does not
establish that an account can claim the name. The owner must check account-side
rights separately. See [release preparation](release_preparation.md) for exact
artifact hashes, executed checks, dependency versions and remaining owner work.

## NOW-TESTED: exact local Python artifacts

Use Python 3.10 or 3.12 on Linux x86_64 for the tested configurations. The local
candidate version is **0.11.0rc1**, proposed for owner approval. The version was
changed only in a disposable build tree; tracked historical 0.10.0 is unchanged.
Obtain the audited wheel **and matching clean sdist** from the preparation output
directory documented in the evidence record. They have not been uploaded.
For this local handoff that directory is `/tmp/circyto-install-lab/dist/`.
Verify their recorded SHA-256 before installing. For example, from a new working
directory outside any checkout:

```bash
python3 -m venv .venv
source .venv/bin/activate
unset PYTHONPATH
export PYTHONNOUSERSITE=1
python -m pip install /absolute/path/circyto-0.11.0rc1-py3-none-any.whl
python -m pip check
circyto --help
circyto doctor
```

For optional RNA+circ MuData instead, install that same wheel with its existing
extra (quote the complete path):

```bash
python -m pip install '/absolute/path/circyto-0.11.0rc1-py3-none-any.whl[mudata]'
```

The base installation exports AnnData and HTML. MuData is optional; the existing
`loom`, `scanpy`, `bam` and `dev` extras are retained. Java/CIRI3, BWA and STAR are
not required for import/report. Doctor can report core readiness while detector
readiness is missing. Install the external runtimes only for detector workflows
using [external detector instructions](external_detectors.md). No JAR is bundled
or downloaded during build, installation or import. circyto code is MIT; the
retained upstream source assets have separate GPLv2 notices.

## Supplied native example, without a source checkout

Examples live in the **sdist**, not in the installed wheel. Extract only the
example directory from the exact matching audited sdist into the current working
directory. This is distribution data, not an editable/source installation:

```bash
tar -xzf /absolute/path/circyto-0.11.0rc1.tar.gz \
  --strip-components=1 circyto-0.11.0rc1/examples
python -I examples/installed_release_check.py \
  --fixtures examples/native_ciri3 --out native-check
```

The script runs help/doctor and the installed `import-circ` command, reads the
result with AnnData and verifies values, identities, coordinates and annotations.
Open `native-check/native/qc/report.html` directly in a browser. It is offline
HTML with no server or JavaScript. Detailed output: `native-check/acceptance.json`.
The equivalent direct command is:

```bash
circyto import-circ --ciri3-manifest examples/native_ciri3/manifest.tsv \
  --cells examples/native_ciri3/cells.tsv --detector CIRI3 \
  --genome-build upstream_demo_chr1_unverified \
  --coordinate-convention 1-based-closed \
  --support-definition CIRI3_junction_reads --outdir native-direct
```

For a MuData installation, use a new output directory:

```bash
python -I examples/installed_release_check.py \
  --fixtures examples/native_ciri3 --out multimodal-check --multimodal
```

This creates explicitly **synthetic RNA** and tests one exact matching ID, two
RNA-only IDs, matrix preservation and real `mudata.read_h5mu` loading. Output:
`multimodal-check/matched/mudata/circyto_multimodal.h5mu` and its IMPORT HTML.
The authentic upstream format sample has six candidates and support 50; its
producer version, assembly and processing completion are unknown. It is not an
established biological single cell. The synthetic RNA is not biological matched
RNA. These checks do not close real matched-cell or original workflow-report
acceptance gaps. Standalone `circyto report` correctly refuses the import bundle.

## NOW-TESTED: conda environment around the local wheel

This is an **environment installation**, not a native conda package. Extract
`packaging/environment.yml` from the matching sdist, or use the reviewed file in
this branch. The file includes MuData. Use a new environment prefix:

```bash
tar -xzf /absolute/path/circyto-0.11.0rc1.tar.gz \
  --strip-components=1 circyto-0.11.0rc1/packaging/environment.yml
micromamba create --no-rc -y -p /absolute/path/circyto-env \
  --override-channels -c conda-forge -c bioconda --strict-channel-priority \
  -f packaging/environment.yml
micromamba run -p /absolute/path/circyto-env python -m pip install \
  --no-deps --no-index /absolute/path/circyto-0.11.0rc1-py3-none-any.whl
micromamba run -p /absolute/path/circyto-env python -m pip check
micromamba run -p /absolute/path/circyto-env python -I \
  examples/installed_release_check.py --fixtures examples/native_ciri3 \
  --out conda-wheel-check --multimodal
```

The scientific dependencies come from conda-forge/bioconda, with strict priority;
pip installs only the exact local wheel. No global `.condarc` changes are needed.

## NOW-TESTED: native local conda artifact

The separate recipe at `packaging/conda/meta.yaml` builds from the audited sdist,
not a GitHub source archive. See the evidence record for its exact filename/hash.
The local artifact is
`/tmp/circyto-install-lab/conda-artifacts/noarch/circyto-0.11.0rc1-py_0.conda`.
Use the extracted environment file to solve scientific dependencies, then install
that exact local package offline into a new prefix:

```bash
micromamba create --no-rc -y -p /absolute/path/circyto-conda \
  --override-channels -c conda-forge -c bioconda --strict-channel-priority \
  -f packaging/environment.yml
micromamba install --no-rc -y --offline -p /absolute/path/circyto-conda \
  --override-channels -c conda-forge -c bioconda --strict-channel-priority \
  /absolute/path/circyto-0.11.0rc1-py_0.conda
micromamba run -p /absolute/path/circyto-conda python -m pip check
micromamba run -p /absolute/path/circyto-conda python -I \
  examples/installed_release_check.py --fixtures examples/native_ciri3 \
  --out conda-package-check --multimodal
```

MuData is an optional conda dependency: include `mudata>=0.3` in the environment
solve for multimodal work, as this file does. For a base environment, omit its
MuData line and the check's `--multimodal` flag. Both variants were tested.
`circyto[mudata]` is pip syntax, not a conda package
name. The recipe's MuData acceptance runs as a test dependency; core recipe tests
use only runtime dependencies. Noarch eligibility covers the pure Python package
and text resources, not portability of external shell/Perl/Java/aligner tools.
Only Linux x86_64 was executed in this preparation.

## AFTER-PUBLICATION: commands to verify later

These are future commands, **not verified registry routes today**. Once the
owner publishes and checks the exact version and files on PyPI:

```bash
python -m pip install circyto
python -m pip install 'circyto[mudata]'
```

For an approved prerelease use an exact version (`circyto==0.11.0rc1`), since pip
does not normally select prereleases automatically. After Bioconda accepts the
recipe and publishes the package, verify the version before using:

```bash
micromamba create --no-rc -n circyto --override-channels \
  -c conda-forge -c bioconda --strict-channel-priority circyto
# Add mudata to that solve for optional multimodal support.
```

The same channel order and strict priority apply with conda. External detector
installation remains an additional, explicitly selected step.
