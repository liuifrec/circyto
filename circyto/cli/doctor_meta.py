from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from typing import List

from circyto.paths import PathResolution, find_ciri2_adapter, find_ciri_full_adapter, find_ciri_full_jar
from circyto.detectors.ciri3 import inspect_ciri3_runtime


@dataclass(frozen=True)
class DetectorSpec:
    name: str
    det_type: str  # "CLI" | "JAR" | "STAR" | etc.
    required_cmds: List[str]
    required_assets: List[str]
    hint_lines: List[str]
    optional_cmds: List[str] = field(default_factory=list)
    validation_status: str = "experimental adapter"
    notes: List[str] = field(default_factory=list)


DETECTOR_SPECS: List[DetectorSpec] = [
    DetectorSpec(
        name="find-circ3",
        det_type="CLI",
        required_cmds=["find-circ3", "bowtie2", "samtools"],
        required_assets=[],
        hint_lines=[
            "conda install -c bioconda bowtie2 samtools and install find-circ3 from upstream",
            "or use apt/brew for bowtie2 + samtools and make sure `find-circ3` is on PATH",
        ],
        validation_status="experimental adapter",
        notes=[
            "Adapter and parser are available for synthetic and small-fixture use.",
            "Not part of the manuscript-scale validated detector set.",
        ],
    ),
    DetectorSpec(
        name="ciri2",
        det_type="CLI",
        required_cmds=["perl", "bwa"],
        required_assets=["CIRI2-adapter"],
        hint_lines=[
            "Install deps: conda install -c bioconda bwa perl",
            "Keep tools/CIRI-full_v2.0/bin/ciri2_adapter.sh available, or set CIRCYTO_CIRI2_ADAPTER.",
        ],
        validation_status="experimental adapter",
        notes=[
            "Implemented with synthetic/unit coverage and a gated chr21 single-end known-positive regression.",
            "No manuscript-scale validation claim is made for CIRI2.",
        ],
    ),
    DetectorSpec(
        name="ciri-full",
        det_type="JAR",
        required_cmds=["bwa", "java"],
        required_assets=["CIRI-full-jar", "CIRI-full-adapter"],
        hint_lines=[
            "Install deps: conda install -c bioconda bwa openjdk",
            "Set CIRCYTO_CIRI_FULL_JAR=/abs/path/CIRI-full.jar or place it under tools/CIRI-full_v2.0/",
            "Keep tools/CIRI-full_v2.0/bin/ciri_full_adapter.sh available; circyto uses it at runtime.",
        ],
        validation_status="experimental adapter",
        notes=[
            "Implemented with synthetic/unit coverage and a gated chr21 Smart-seq2 integration example.",
            "Single-end rows use a CIRI2 fallback; upstream CIRI-full Pipeline behavior requires paired-end reads.",
        ],
    ),
    DetectorSpec(
        name="ciri3",
        det_type="JAR",
        required_cmds=[],
        required_assets=[],
        hint_lines=[
            "CIRI3 requires Java plus a local CIRI3 jar; set CIRCYTO_CIRI3_JAR or CIRCYTO_CIRI3_HOME.",
            "Install samtools for alignment handling; install bwa for BWA mode or STAR for STAR mode.",
            "If direct java -jar execution is unsuitable locally, set CIRCYTO_CIRI3_CMD_TEMPLATE or CIRCYTO_CIRI3_BIN.",
            "Set CIRCYTO_CIRI3_JAVA if circyto should use a non-default java executable.",
            "Set CIRCYTO_CIRI3_EXTRA_ARGS for local tuning such as -S 0 during single-end validation.",
            "STAR is only required for STAR-based CIRI3 workflows.",
        ],
        optional_cmds=["bwa", "STAR", "samtools"],
        validation_status="primary validated detector",
        notes=[
            "Primary detector for the validated Smart-seq3 and RamDA/scRR manuscript workflows.",
            "Direct jar mode requires Java plus a local CIRI3 jar; template mode may be configured separately.",
        ],
    ),
    DetectorSpec(
        name="circexplorer2",
        det_type="CLI",
        required_cmds=["CIRCexplorer2", "STAR"],
        required_assets=[],
        hint_lines=[
            "Install CIRCexplorer2 and STAR only when using this optional adapter.",
            "Set STAR index and refFlat inputs through the detector workflow options.",
        ],
        validation_status="optional experimental adapter",
        notes=[
            "Optional adapter with parser/collector coverage; not part of manuscript-scale validation.",
        ],
    ),
]


def resolve_asset(asset: str) -> PathResolution:
    if asset == "CIRI-full-jar":
        return find_ciri_full_jar()
    if asset == "CIRI-full-adapter":
        return find_ciri_full_adapter()
    if asset == "CIRI2-adapter":
        return find_ciri2_adapter()
    return PathResolution(label=asset, resolved_path=None, checked_paths=(), source=None)


def detector_runtime_status(name: str) -> dict:
    if name == "ciri3":
        details = inspect_ciri3_runtime()
        errors = details.get("template_errors", []) if details.get("template_configured") else []
        if details["direct_ready"]:
            status = "READY"
            reason = (
                f"direct java -jar contract available (Java {details.get('java_version')}, "
                f"need >={details.get('required_java_major')})"
            )
        elif details["template_configured"] and not errors:
            status = "READY"
            reason = "template contract configured"
        elif details.get("jar") and details.get("java") and not details.get("java_version_ok"):
            status = "NOT READY"
            reason = (
                f"Java version too old (found {details.get('java_version')}, "
                f"need >={details.get('required_java_major')} for bundled CIRI3 jar)"
            )
        elif details.get("jar") or details.get("bin"):
            status = "PARTIAL"
            reason = "found local CIRI3 assets but runtime contract is incomplete"
        else:
            status = "NOT READY"
            reason = "no local CIRI3 jar or wrapper detected"
        return {
            "status": status,
            "reason": reason,
            "details": details,
        }

    spec = next(spec for spec in DETECTOR_SPECS if spec.name == name)
    missing_cmds = [cmd for cmd in spec.required_cmds if shutil.which(cmd) is None]
    missing_assets = [asset for asset in spec.required_assets if not resolve_asset(asset).found]
    if not missing_cmds and not missing_assets:
        status = "READY"
        reason = "all required commands and assets found"
    elif missing_cmds:
        status = "NOT READY"
        reason = f"missing commands: {', '.join(missing_cmds)}"
    else:
        status = "PARTIAL"
        reason = f"missing assets: {', '.join(missing_assets)}"
    return {
        "status": status,
        "reason": reason,
        "details": {
            "missing_cmds": missing_cmds,
            "missing_assets": missing_assets,
        },
    }
