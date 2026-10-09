# Third-party detector distribution and attribution

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
