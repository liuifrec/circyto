# Third-party distribution policy: external Java detectors

Packaging branch audit: 2026-10-13, based on PR #19 head
`f841c296db8a6fda870eebc5b553eab10095fa0d`. No package version change.
The MIT license covers circyto-authored code, not third-party assets.

| Asset | Current wheel/sdist policy | Confirmed evidence / unresolved questions |
| --- | --- | --- |
| CIRI3 JARs, embedded HTSJDK classes | **Excluded**; external user installation | Official pinned CIRI3 tree supplies GPLv2 and Java source. Exact corresponding build/dependency notice route is not certified. |
| CIRI-full JAR and CIRI-full manual | **Excluded** | No applicable permission found for bundled v2.0 JAR or manual; source visibility/newer binaries do not resolve it. |
| CIRI2 / CIRI-AS Perl scripts, their manuals and adjacent GPLv2 texts | Retained unchanged | Script copyright headers refer to accompanying license; complete supplied Perl source is shipped. Preserve notices and GPL obligations when redistributing. |
| CIRI3 README / license and authentic example / provenance / license | Retained unchanged | Supplied upstream GPLv2 notice retained. These are documentation/input evidence, not executable JARs or new biological validation. |
| circyto shell adapters and Python code | Retained | External JAR resolution remains available; no calling algorithms or thresholds changed. |

`MANIFEST.in` excludes Java binaries/classes and the unresolved manual from the
sdist. Setuptools package-data exclusions apply to wheels too. A build guard
rejects those assets if stale build trees bypass the exclusion lists. Regression
tests inspect direct wheels, sdists and wheels rebuilt from sdists, and seed
stale files to verify refusal. This policy applies to **newly built Python
artifacts only**. Historical tracked binaries, previous releases, and GitHub
source archives still require separate maintainer/legal review; preserving
existing branches/history is not retroactive redistribution clearance. Do not
publish old or automatic Git source archives as if this policy covered them.

## Upstream evidence checked in this pass

- CIRI3 pinned tree `73108c4478a7e691cb93360921b65203e0701516` (v3.0.1)
  contains GPLv2, Java source, `.classpath` with JavaSE-12 and
  `lib/htsjdk-3.0.4.jar`; current upstream main queried as
  `a19560452a94a63ee23757fb4a7032582d052b45`. Its updated Java-8 binary is not
  substituted for the existing reviewed runtime. Binary hashes are below.
- All 1,181 embedded HTSJDK class entries match the 3.0.4 library already
  identified below. Official HTSJDK tag 3.0.4 resolves to
  `e2943b78ff595e144281b8a3c126be5a4d70085d`.
  Its [SAMRecord source header](https://github.com/samtools/htsjdk/blob/3.0.4/src/main/java/htsjdk/samtools/SAMRecord.java)
  supplies the Broad Institute MIT permission and notice requirement. This
  confirms permission for that source, not a completed inventory of every
  embedded component. Neither inspected JAR has LICENSE/NOTICE archive entries.
- CIRI-full current upstream tree remains
  `03a5719ba0b28ae31c17a4029310ffbfac5b7241`: source/Makefile and newer v2.1.2
  binary present, no applicable root license found. The locally bundled v2.0
  manual and JAR contain no permission text. Nested Perl licenses are not
  assumed to license either asset.

**Decision:** engineering GO for exclusion of the unresolved Java distribution
and preservation of external adapters, subject to the recorded artifact tests.
**NO-GO for redistributing those JARs/manual or claiming blanket legal clearance.**
No maintainer permission was requested/received and no institutional legal
opinion was obtained. Further clearance of these assets or past releases is a
maintainer/legal decision, not another code change in this task.
See [installation and detector citations](docs/external_detectors.md) and
[validation record](docs/detector_distribution_validation.md).

# Historical detector distribution audit (before packaging exclusion)

Audit date: 2026-10-09. The root MIT license covers circyto-authored code.
It does **not** relicense the detector assets below. They are already present
in the wheel and source distribution under `circyto/resources/tools/`;
installing circyto therefore distributes these assets even when a detector
is not executed. Original notices, manuals, and detector code are unchanged.

| Packaged asset (relative to `circyto/resources/tools/`) | Attribution and supplied license | Distribution finding |
| --- | --- | --- |
| `CIRI3/CIRI3_Java_18.0.1.jar` | [CIRI3 upstream](https://github.com/gyjames/CIRI3); GPL version 2 text in adjacent `LICENSE.txt` | Binary and README are packaged. Exact upstream binary identity is verified below; corresponding source/build material is not included in the Python distribution. |
| `CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl` | BIOLS, Chinese Academy of Sciences; adjacent `LICENSE.txt` supplies GPL version 2 | Perl source, manual, and license are packaged together. |
| `CIRI-full_v2.0/bin/CIRI_AS_v1.2/CIRI_AS_v1.2.pl` | BIOLS, Chinese Academy of Sciences; adjacent `LICENSE.txt` supplies GPL version 2 | Perl source, manual, and license are packaged together. |
| `CIRI-full_v2.0/CIRI-full.jar` | [CIRI-full upstream](https://github.com/bioinfo-biols/CIRI-full); bundled manual identifies Yi Zheng, BIOLS, CAS as contact | No license specifically covering this JAR was found in the packaged directory, JAR, manual, or inspected upstream tree. The subdirectories' CIRI2/CIRI-AS licenses must not be assumed to cover it. Redistribution permission remains unresolved. |

The bundled CIRI3 JAR has SHA-256
`569428a4fe0d6573dbeeda6d3d4a0d457ac01d01bbf0b49ab1a42b84e24f1be0`.
Its Git blob `c2d3c23322a734a14bebe20a7cdfc58f996a9dd1` matches
[upstream commit 73108c4](https://github.com/gyjames/CIRI3/tree/73108c4478a7e691cb93360921b65203e0701516).
That tree contains source and library directories. A source URL alone is not
a completed audit of the corresponding-source distribution obligations in
[the supplied GPLv2 terms](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/LICENSE.txt).
The JAR also embeds `htsjdk/*` classes and contains no LICENSE/NOTICE entries;
the exact embedded library version and its required notices remain to be
verified and packaged before a release. This pass does not certify the bundle
as license-complete.

The existing CIRI-full v2.0 JAR has SHA-256
`5cdba739ca9ec5cb3dea145872916899598747fe6cb1af4566a93843a577b3ac`.
It includes Java source files but no license text. The inspected
[upstream tree 03a5719](https://github.com/bioinfo-biols/CIRI-full/tree/03a5719ba0b28ae31c17a4029310ffbfac5b7241)
contains a newer v2.1.2 JAR; it is not evidence of permission for the bundled
v2.0 binary. Maintainers should confirm permission with upstream or separately
review external-tool acquisition before a new public software release.
No upstream code or binary was replaced during this audit.

BWA, SAMtools, Java, STAR, Bowtie2, find-circ3, minimap2, and CIRI-long are
external runtimes, not executables shipped in this wheel. Their own
distributions provide their licenses. The optional environment recipe installs
BWA, SAMtools and OpenJDK through conda-forge/bioconda; it does not change their
licenses. `CIRI-vis.jar` exists in the development checkout but is not included
in the inspected Python wheel/sdist.

Cite the actual detector and protocol used, in addition to circyto's
[`CITATION.cff`](CITATION.cff). Tool attribution and a software citation do not
establish biological validation of a circRNA candidate.

## Pre-release packaging recommendation (2026-10-13)

The integration audit reconfirmed both JAR hashes above. Neither archive has
LICENSE/NOTICE entries. All **1,181** embedded `htsjdk/*.class` entries in CIRI3
are byte-identical to the corresponding entries in upstream
`lib/htsjdk-3.0.4.jar` (SHA-256
`41c98336c12f5df81d6a3bef024d98dbb8344396cc140c3df955edf79ef726d6`).
This identifies matching class bytes; it does not complete a transitive
license audit or establish the reproducible source/build correspondence of
CIRI3 as a whole. The library archive also lacks LICENSE/NOTICE entries.

**Recommendation, not implemented:** make the narrow packaging change that
excludes the two unresolved JARs from both wheel and sdist, while retaining the
Python import/report commands, detector adapters, and separately reviewed
source scripts/notices. Preserve detector execution through the existing
`CIRCYTO_CIRI3_JAR` / `CIRCYTO_CIRI3_HOME` and `CIRCYTO_CIRI_FULL_JAR`
external-install resolution. Document user acquisition from upstream, supported
versions/checksums and Java requirements; test installed adapters with externally
supplied tools before shipping this migration. Do not add automatic downloads,
mirror the unresolved binaries, or call a Python extra a licensing remedy.
Do not remove an adapter or silently change its default availability.

If maintainers instead retain CIRI3 redistribution, first provide a reviewed
GPLv2 section 3 compliance route with exact corresponding source, build
instructions and dependency notices; a repository link alone is not the
completed review. See the supplied [GPLv2 text](https://github.com/gyjames/CIRI3/blob/73108c4478a7e691cb93360921b65203e0701516/LICENSE.txt)
and the [FSF GPLv2 distribution FAQ](https://www.gnu.org/licenses/old-licenses/gpl-2.0-faq.html).
CIRI-full v2.0 additionally needs an applicable permission grant from the rights
holder before redistribution; newer source visibility is not that grant.
Existing GPL Perl sources and the native example/license still need their
notices retained and a final artifact inventory review. This is a release
engineering recommendation pending maintainer/legal review, not legal clearance.
No third-party files or packaging inclusion rules are changed by this audit.
CI retains test records only, not built detector-containing distributions.
