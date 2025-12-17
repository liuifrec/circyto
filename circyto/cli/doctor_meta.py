from __future__ import annotations

from dataclasses import dataclass
from typing import List


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
            "Place CIRI-full*.jar under: tools/",
        ],
    ),
]
