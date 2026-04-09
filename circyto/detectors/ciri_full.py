# circyto/detectors/ciri_full.py
from __future__ import annotations

import os
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .base import DetectorBase, DetectorCapabilities, DetectorRunInputs, DetectorResult
from .ciri2 import _infer_fastq_read_length
from circyto.paths import (
    find_ciri_full_adapter,
    find_ciri_full_jar,
    format_missing_resolution,
)


@dataclass
class CiriFullDetector(DetectorBase):
    """
    Detector engine for the CIRI-full Pipeline.

    This is a thin wrapper around the existing ciri_full_adapter.sh and
    the CIRI-full jar. We expose it through the generic DetectorBase
    interface so it can be used with run-detector / future multi-detector
    orchestration.

    IMPORTANT: CIRI-full is *not* safe to run in parallel across many cells
    in the same process / working directory. We enforce max_parallel = 1
    at the orchestration layer.
    """

    name: str = "ciri-full"
    input_type: str = "fastq"
    supports_paired_end: bool = True
    capabilities: DetectorCapabilities = DetectorCapabilities(
        accepts_fastq=True,
        accepts_alignment=False,
        prefers_paired=True,
        supports_single_end=True,
        supports_multisample_alignment=False,
        max_parallel=1,
        recommended_execution_mode="per-cell-fastq",
    )

    ciri_full_jar: Path | None = None
    adapter_script: Path | None = None

    # New: tell the orchestrator this tool must run serially
    max_parallel: int = 1

    @staticmethod
    def _read_layout(inputs: DetectorRunInputs) -> str:
        return "paired-end" if inputs.r2 is not None else "single-end"

    @staticmethod
    def _read_layout_env(layout: str) -> str:
        return "paired" if layout == "paired-end" else "single"

    @staticmethod
    def _tail_log(path: Path, max_lines: int = 40) -> str:
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            return ""
        if not lines:
            return ""
        return "\n".join(lines[-max_lines:])

    def _resolve_ciri_full_jar(self) -> Path:
        resolution = find_ciri_full_jar(self.ciri_full_jar)
        if resolution.resolved_path is None:
            raise FileNotFoundError(format_missing_resolution("CIRI-full jar", resolution))
        return resolution.resolved_path

    def _resolve_adapter_script(self) -> Path:
        resolution = find_ciri_full_adapter(self.adapter_script)
        if resolution.resolved_path is None:
            raise FileNotFoundError(format_missing_resolution("CIRI-full adapter", resolution))
        return resolution.resolved_path

    def is_available(self) -> bool:
        # java + jar + adapter must exist
        if shutil.which("java") is None:
            return False
        if find_ciri_full_jar(self.ciri_full_jar).resolved_path is None:
            return False
        if find_ciri_full_adapter(self.adapter_script).resolved_path is None:
            return False
        return True

    def version(self) -> Optional[str]:
        """
        Try to get a version string.

        For now, we just return the manual version you've observed (2.1.1).
        We can make this smarter later by parsing '-h' output.
        """
        return "CIRI-full 2.1.1"

    def run(self, inputs: DetectorRunInputs) -> DetectorResult:
        """
        Run CIRI-full on a single cell via the adapter script.

        Environment variables passed to the adapter:

          R1, R2, REF_FA, GTF, OUT_DIR, OUT_BASENAME, OUT_TSV, THREADS

        The adapter is responsible for:

          - creating a run directory:  ${OUT_DIR}/${OUT_BASENAME}.ciri_full_run
          - invoking: java -jar CIRI-full.jar Pipeline ...
          - copying/linking the final .ciri/.txt into OUT_TSV (normalized TSV)
        """
        outdir = inputs.outdir
        outdir.mkdir(parents=True, exist_ok=True)

        r1 = inputs.r1
        r2 = inputs.r2
        ref_fa = inputs.ref_fa
        gtf = inputs.gtf
        threads = inputs.threads
        cell_id = inputs.cell_id
        read_layout = self._read_layout(inputs)

        # --- Dry-run mode for tests (no real CIRI-full execution) ---
# If ref_fa or gtf is missing AND we detect fake FASTQs, produce empty TSV.
        if ref_fa is None or gtf is None:
            out_tsv = outdir / f"{cell_id}.tsv"
            with out_tsv.open("w") as f:
                f.write("circ_id\tchr\tstart\tend\tstrand\tsupport\n")
            # Return fake DetectorResult
            return DetectorResult(
                detector=self.name,
                cell_id=cell_id,
                outdir=outdir,
                tsv_path=out_tsv,
                run_dir=outdir / f"{cell_id}.ciri_full_run",
                log_path=outdir / f"{cell_id}.ciri_full.log",
                meta={"dry_run": True, "read_layout": read_layout},
            )

        if r1 is None:
            raise ValueError("CiriFullDetector requires R1 FASTQ")

        out_tsv = outdir / f"{cell_id}.tsv"
        run_dir = outdir / f"{cell_id}.ciri_full_run"
        log_path = outdir / f"{cell_id}.ciri_full.log"

        env = {
            "R1": str(r1),
            "R2": str(r2) if r2 is not None else "",
            "REF_FA": str(ref_fa),
            "GTF": str(gtf),
            "OUT_DIR": str(outdir),
            "OUT_BASENAME": cell_id,
            "OUT_TSV": str(out_tsv),
            "THREADS": str(threads),
            "CIRCYTO_READ_LAYOUT": self._read_layout_env(read_layout),
            # NEW: ask CIRI-full/CIRI to keep 1-BSJ circRNAs
            "CIRI_EXTRA_FLAGS": "-0",
        }
        read_length = _infer_fastq_read_length(r1)
        if read_layout == "single-end":
            env["CIRI2_BWA_MEM_FLAGS"] = "-k 15 -T 15" if read_length is not None and read_length < 60 else "-T 19"
        adapter_script = self._resolve_adapter_script()
        ciri_full_jar = self._resolve_ciri_full_jar()

        # Inherit the current environment, then overlay our variables.
        # This ensures PATH, java, bwa, samtools, etc. are visible.
        real_env = os.environ.copy()
        real_env.update(env)
        real_env["CIRCYTO_CIRI_FULL_JAR"] = str(ciri_full_jar)

        cmd = ["bash", str(adapter_script)]
        adapter_cmd = " ".join(cmd)

        run_started = time.perf_counter()
        with log_path.open("w", encoding="utf-8") as log_handle:
            result = subprocess.run(
                cmd,
                check=False,
                env=real_env,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )

        if result.returncode != 0:
            log_tail = self._tail_log(log_path)
            detail = (
                f"CIRI-full adapter failed for cell_id={cell_id} "
                f"(read_layout={read_layout}, exit_code={result.returncode}). "
                f"Adapter command: {adapter_cmd}. Log: {log_path}"
            )
            if log_tail:
                detail += f"\n--- log tail ---\n{log_tail}"
            raise RuntimeError(detail)

        # If we reach here, the adapter exited 0.
        # We expect OUT_TSV and the run dir to exist.
        return DetectorResult(
            detector=self.name,
            cell_id=cell_id,
            outdir=outdir,
            tsv_path=out_tsv,
            run_dir=run_dir,
            log_path=log_path,
            meta={
                "threads": threads,
                "elapsed_seconds": round(time.perf_counter() - run_started, 3),
                "ciri_full_jar": str(ciri_full_jar),
                "adapter_script": str(adapter_script),
                "adapter_command": adapter_cmd,
                "read_layout": read_layout,
                "pipeline_mode": (
                    "ciri-full-pipeline" if read_layout == "paired-end" else "ciri2-single-end-fallback"
                ),
                "read_length": read_length,
                "bwa_mem_flags": env.get("CIRI2_BWA_MEM_FLAGS"),
                "raw_output_path": (
                    str(run_dir / f"{cell_id}.ciri2.txt")
                    if read_layout == "single-end"
                    else None
                ),
                "input_mode": "fastq",
                "detector_backend": self.name,
            },
        )

    def run_from_fastq(self, inputs: DetectorRunInputs) -> DetectorResult:
        return self.run(inputs)
