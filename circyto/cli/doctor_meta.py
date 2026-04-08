from __future__ import annotations

from dataclasses import dataclass
from typing import List

from circyto.paths import PathResolution, find_ciri_full_jar


@dataclass(frozen=True)
class DetectorSpec:
    name: str
    det_type: str  # "CLI" | "JAR" | "STAR" | etc.
    required_cmds: List[str]
    required_assets: List[str]
    hint_lines: List[str]


DETECTOR_SPECS: List[DetectorSpec] = [
    DetectorSpec(
        name="find-circ3",
        det_type="CLI",
        required_cmds=["bowtie2", "samtools"],
        required_assets=[],
        hint_lines=[
            "conda install -c bioconda bowtie2 samtools",
            "or use apt/brew to install bowtie2 and samtools",
        ],
    ),
    DetectorSpec(
        name="ciri-full",
        det_type="JAR",
        required_cmds=["bwa", "java"],
        required_assets=["CIRI-full-jar"],
        hint_lines=[
            "Install deps: conda install -c bioconda bwa openjdk",
            "Set CIRCYTO_CIRI_FULL_JAR=/abs/path/CIRI-full.jar or place it under tools/CIRI-full_v2.0/",
        ],
    ),
]


def resolve_asset(asset: str) -> PathResolution:
    if asset == "CIRI-full-jar":
        return find_ciri_full_jar()
    return PathResolution(label=asset, resolved_path=None, checked_paths=(), source=None)
