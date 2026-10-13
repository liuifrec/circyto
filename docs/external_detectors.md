# External Java detectors

This packaging branch does **not** ship CIRI3 or CIRI-full JARs in wheels or
sdists. It also excludes the CIRI-full manual whose redistribution permission
has not been established. Adapters, licensed CIRI2/CIRI-AS Perl source and their
notices remain packaged. No automatic detector download occurs. The version
remains 0.10.0; use the exact packaging-branch wheel, not an older same-version
wheel. Historical wheels and Git source archives are not made compliant by
these changes.

`import-circ` (including native CIRI3 import), AnnData/MuData export and HTML
reporting require **no Java or detector executable**. Install the reviewed wheel
in a fresh environment; add `[mudata]` if you need RNA integration:

```bash
python -m venv circ-env
. circ-env/bin/activate
python -m pip install '/absolute/path/circyto-0.10.0-py3-none-any.whl[mudata]'
python -m pip check
circyto doctor
circyto import-circ --help
```

The example assumes POSIX shell and Python >=3.10. Use a new environment or
`--force-reinstall --no-deps` for an existing installation of the same circyto
version. Record the wheel SHA-256. Ordinary pip upgrades cannot distinguish
same-version feature wheels. Do not copy a new wheel over an old installed tree.

## CIRI3 installation for established workflows

Obtain CIRI3 directly from its official repository, retaining its license.
The following pins the same upstream binary reviewed in the existing workflows;
it does not substitute a newer Java-8 build or alter biological parameters:

```bash
git clone https://github.com/gyjames/CIRI3.git CIRI3-external
git -C CIRI3-external checkout --detach 73108c4478a7e691cb93360921b65203e0701516
export CIRCYTO_CIRI3_JAR="$PWD/CIRI3-external/CIRI3_Java_18.0.1.jar"
printf '%s  %s\n' \
  '569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0' \
  "$CIRCYTO_CIRI3_JAR" | sha256sum --check
java -version
java -jar "$CIRCYTO_CIRI3_JAR" -H
circyto doctor
circyto detectors --json
```

Use Java >=12 for this reviewed binary, preferably Java 17 (the pinned source
classpath targets JavaSE-12). The filename's `18.0.1` is not the detector version
or a requirement for Java 18. Other binaries may have other requirements and
are not certified by this check. `CIRCYTO_CIRI3_JAVA=/absolute/path/to/java`
selects a non-default Java executable.

The JAR's actual upstream location must match the command above; see the pinned
repository tree. Alternatively set `CIRCYTO_CIRI3_HOME` to the directory directly
containing the desired `CIRI3*.jar`. For deterministic selection prefer the exact
JAR path, especially if several JARs exist. Existing `CIRCYTO_CIRI3_BIN`,
`CIRCYTO_CIRI3_SEARCH_ROOTS`, `CIRCYTO_CIRI3_CMD_TEMPLATE` and extra-argument
settings remain supported; a wrapper alone requires its existing execution
contract/template. Inspect resolved paths in `doctor`; existing resolver fallback
behavior is unchanged, so a misspelled override may fall back to another local
installation. Do not assume that setting a variable proves which binary ran.

For full workflows also install BWA, SAMtools, and STAR for STAR mode. The
[Linux environment recipe](../environment.ramda-se.yml) provides BWA/SAMtools/
Java 17; it does not download CIRI3. Supply real reference/index/annotation files
and use the existing [workflow instructions](../README.md). `--dry-run` verifies
planning, not biological execution. `doctor` reports core Python readiness
separately from optional detector readiness; exit zero alone is not a detector
readiness assertion.

## CIRI-full

Obtain an appropriately permitted copy from the rights holder or your
institution. This project does not supply a download command for the unresolved
v2.0 binary or infer a permission grant from publicly visible source:

```bash
export CIRCYTO_CIRI_FULL_JAR=/absolute/path/to/your/permitted/CIRI-full.jar
circyto doctor
circyto detectors --json
```

The packaged shell adapter is still discovered automatically and receives the
resolved JAR path. BWA and Java are needed; CIRI2/CIRI-AS sources/licenses remain
beside the adapter. Do not replace `CIRCYTO_TOOLS_DIR` with a directory containing
only JARs: that setting overrides discovery of the packaged adapters too.
Use the individual JAR variables above. Existing source checkouts may still
resolve their historical `tools/` copies; an installed wheel outside the checkout
is the authoritative packaging test.

CIRI-full's existing single-end fallback also performs its JAR preflight before
using CIRI2. That behavior is unchanged; the standalone CIRI2 adapter remains
available without a CIRI-full JAR. No claiming that paired-end CIRI-full biological
execution was validated by path discovery, mocks or `doctor`.

## Attribution and distribution limits

Cite circyto ([CITATION.cff](../CITATION.cff)) and each detector actually used:

- CIRI3: [Detecting and quantifying circular RNAs in terabyte-scale RNA-seq
  datasets with CIRI3](https://doi.org/10.1038/s41587-025-02835-1), and record the
  exact upstream revision/executable hash.
- CIRI-full: [Reconstruction of full-length circular RNAs enables isoform-level
  quantification](https://doi.org/10.1186/s13073-019-0614-1).
- CIRI2: [Circular RNA identification based on multiple seed matching](https://doi.org/10.1093/bib/bbx014).
- CIRI-AS when used: [Comprehensive identification of internal structure and
  alternative splicing events in circular RNAs](https://doi.org/10.1038/ncomms12060).

See [third-party notices](../THIRD_PARTY_NOTICES.md) for exact permission findings.
A direct upstream installation does not clear permission to mirror or redistribute
that installation. No institutional permission or new upstream grant was obtained.
The unresolved binary/source/notice questions require maintainer or legal review;
the package exclusion is a technical boundary, not certification of upstream
compliance or retroactive clearance of old releases.
