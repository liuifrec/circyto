from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from circyto.cli.circyto import app
from circyto.pipeline.scomatic_full_length_adapter import HAS_ANNDATA, HAS_MUDATA


runner = CliRunner()


def _write_fake_samtools(tmp_path: Path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir()
    samtools = bindir / "samtools"
    samtools.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            from pathlib import Path
            import sys

            args = sys.argv[1:]
            if not args:
                raise SystemExit("missing command")
            cmd = args[0]
            rest = args[1:]
            if cmd == "view":
                if "-b" in rest:
                    sys.stdout.write(sys.stdin.read())
                else:
                    Path(rest[-1]).read_text(encoding="utf-8")
                    sys.stdout.write(Path(rest[-1]).read_text(encoding="utf-8"))
            elif cmd == "sort":
                out = rest[rest.index("-o") + 1]
                Path(out).write_text(sys.stdin.read(), encoding="utf-8")
            elif cmd == "merge":
                out = rest[rest.index("-f") + 1]
                start = rest.index("-f") + 2
                inputs = [item for item in rest[start:] if not item.startswith("-")]
                Path(out).write_text(
                    "".join(Path(item).read_text(encoding="utf-8") for item in inputs),
                    encoding="utf-8",
                )
            elif cmd == "index":
                Path(rest[0] + ".bai").write_text("index\\n", encoding="utf-8")
            else:
                raise SystemExit(f"unsupported fake samtools command: {cmd}")
            """
        ),
        encoding="utf-8",
    )
    samtools.chmod(0o755)
    return bindir


def _write_sam(path: Path, *, read_name: str, cb_tag: str = "") -> None:
    tag = f"\t{cb_tag}" if cb_tag else ""
    path.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        f"{read_name}\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tnM:i:0\tNH:i:1{tag}\n",
        encoding="utf-8",
    )


def test_prepare_scomatic_input_one_cell_per_sam(tmp_path: Path) -> None:
    cell_a = tmp_path / "cellA.sam"
    cell_b = tmp_path / "cellB.sam"
    _write_sam(cell_a, read_name="readA")
    _write_sam(cell_b, read_name="readB", cb_tag="CB:Z:old")
    manifest = tmp_path / "alignment_manifest.tsv"
    manifest.write_text(
        "cell_id\tsam\tcell_type\n"
        f"cellA\t{cell_a}\tTypeA\n"
        f"cellB\t{cell_b}\tTypeB\n",
        encoding="utf-8",
    )
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "prepared"
    env = os.environ.copy()
    env["PATH"] = f"{_write_fake_samtools(tmp_path)}{os.pathsep}{env['PATH']}"

    result = runner.invoke(
        app,
        [
            "prepare-scomatic-input",
            "--alignment-manifest",
            str(manifest),
            "--outdir",
            str(outdir),
            "--reference-fasta",
            str(ref),
            "--protocol",
            "Smart-seq3",
            "--sample-id",
            "sample1",
            "--cell-type-column",
            "cell_type",
        ],
        env=env,
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["protocol"] == "smartseq3"
    assert payload["mode"] == "one-cell-per-bam"
    assert payload["runs_scomatic"] is False
    assert payload["scomatic_result_terminology"] == "RNA-derived candidate variant signals"
    assert payload["tag_summary"]["missing_cb_before_injection"] == 1
    assert payload["tag_summary"]["replaced_cb_not_matching_cell_id"] == 1

    merged = outdir / "merged" / "sample1.scomatic.bam"
    assert merged.exists()
    assert Path(str(merged) + ".bai").exists()
    assert "CB:Z:cellA" in (outdir / "per_cell" / "cellA.scomatic.sorted.bam").read_text(encoding="utf-8")
    assert "CB:Z:cellB" in (outdir / "per_cell" / "cellB.scomatic.sorted.bam").read_text(encoding="utf-8")
    annotation = (outdir / "cell_annotations.tsv").read_text(encoding="utf-8").splitlines()
    assert annotation == ["Index\tCell_type", "cellA\tTypeA", "cellB\tTypeB"]


def test_prepare_scomatic_input_merged_bam_mode(tmp_path: Path) -> None:
    merged = tmp_path / "merged.bam"
    merged.write_text("already merged\n", encoding="utf-8")
    metadata = tmp_path / "cells.tsv"
    metadata.write_text("cell_id\tcell_type\ncellA\tTypeA\ncellB\tTypeB\n", encoding="utf-8")
    outdir = tmp_path / "prepared"

    result = runner.invoke(
        app,
        [
            "prepare-scomatic-input",
            "--merged-bam",
            str(merged),
            "--cell-metadata",
            str(metadata),
            "--outdir",
            str(outdir),
            "--protocol",
            "scrr-rna",
            "--sample-id",
            "merged_sample",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["mode"] == "merged-bam"
    assert payload["protocol"] == "scrr-rna"
    assert payload["merged_bam"] == str(merged.resolve())
    prepared_bams = pd.read_csv(outdir / "scomatic_prepared_bams.tsv", sep="\t")
    assert prepared_bams.loc[0, "bam"] == str(merged.resolve())


def test_run_scomatic_writes_dry_run_plan(tmp_path: Path) -> None:
    merged = tmp_path / "merged.bam"
    merged.write_text("already merged\n", encoding="utf-8")
    metadata = tmp_path / "cells.tsv"
    metadata.write_text("cell_id\tcell_type\ncellA\tTypeA\n", encoding="utf-8")
    prepared = tmp_path / "prepared"
    prepare_result = runner.invoke(
        app,
        [
            "prepare-scomatic-input",
            "--merged-bam",
            str(merged),
            "--cell-metadata",
            str(metadata),
            "--outdir",
            str(prepared),
            "--protocol",
            "ramda",
        ],
    )
    assert prepare_result.exit_code == 0, prepare_result.output

    scomatic_dir = tmp_path / "SComatic"
    (scomatic_dir / "scripts" / "BaseCellCounter").mkdir(parents=True)
    (scomatic_dir / "scripts" / "BaseCellCalling").mkdir(parents=True)
    (scomatic_dir / "scripts" / "BaseCellCounter" / "BaseCellCounter.py").write_text("# base\n", encoding="utf-8")
    (scomatic_dir / "scripts" / "BaseCellCalling" / "BaseCellCalling.step1.py").write_text("# step1\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    outdir = tmp_path / "scomatic_run"

    result = runner.invoke(
        app,
        [
            "run-scomatic",
            "--prepared-dir",
            str(prepared),
            "--outdir",
            str(outdir),
            "--scomatic-dir",
            str(scomatic_dir),
            "--reference-fasta",
            str(ref),
            "--run-step1",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["execute"] is False
    assert [command["stage"] for command in payload["commands"]] == ["BaseCellCounter", "Step1"]
    assert Path(payload["command_plan"]).exists()
    plan = Path(payload["command_plan"]).read_text(encoding="utf-8")
    assert "BaseCellCounter.py" in plan
    assert "RNA-derived candidate variant signals" in plan


def test_run_scomatic_without_execute_does_not_shell_out(tmp_path: Path) -> None:
    merged = tmp_path / "merged.bam"
    merged.write_text("already merged\n", encoding="utf-8")
    metadata = tmp_path / "cells.tsv"
    metadata.write_text("cell_id\tcell_type\ncellA\tTypeA\n", encoding="utf-8")
    prepared = tmp_path / "prepared"
    prepare_result = runner.invoke(
        app,
        [
            "prepare-scomatic-input",
            "--merged-bam",
            str(merged),
            "--cell-metadata",
            str(metadata),
            "--outdir",
            str(prepared),
            "--protocol",
            "ramda",
        ],
    )
    assert prepare_result.exit_code == 0, prepare_result.output

    marker = tmp_path / "executed.marker"
    fake_python = tmp_path / "fake_python.sh"
    fake_python.write_text(
        "#!/usr/bin/env bash\n"
        f"touch {marker}\n"
        "exit 0\n",
        encoding="utf-8",
    )
    fake_python.chmod(0o755)
    scomatic_dir = tmp_path / "SComatic"
    script_dir = scomatic_dir / "scripts" / "BaseCellCounter"
    script_dir.mkdir(parents=True)
    (script_dir / "BaseCellCounter.py").write_text("raise SystemExit('should not run')\n", encoding="utf-8")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")

    result = runner.invoke(
        app,
        [
            "run-scomatic",
            "--prepared-dir",
            str(prepared),
            "--outdir",
            str(tmp_path / "scomatic_run"),
            "--scomatic-dir",
            str(scomatic_dir),
            "--reference-fasta",
            str(ref),
            "--python-executable",
            str(fake_python),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["execute"] is False
    assert payload["execution_results"] == []
    assert not marker.exists()


def test_import_scomatic_alias_writes_candidate_summary(tmp_path: Path) -> None:
    raw = tmp_path / "scomatic.tsv"
    raw.write_text(
        "cell_id\tchrom\tpos\tref\talt\tgene\tfilter_status\tread_support\tvaf\n"
        "cellA\tchr1\t10\tA\tG\tGENE1\tPASS\t4\t0.2\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "imported"

    result = runner.invoke(
        app,
        [
            "import-scomatic",
            "--scomatic-output",
            str(raw),
            "--outdir",
            str(outdir),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["command_name"] == "circyto import-scomatic"
    assert payload["n_candidates"] == 1
    assert (outdir / "scomatic_candidate_summary.tsv").exists()
    assert (outdir / "import_scomatic_summary.json").exists()


@pytest.mark.skipif(not (HAS_ANNDATA and HAS_MUDATA), reason="anndata or mudata not installed")
def test_merge_scomatic_adds_candidate_snv_modality(tmp_path: Path) -> None:
    import anndata as ad
    import mudata as mu

    obs = pd.DataFrame(index=["cellA", "cellB"])
    rna = ad.AnnData(X=np.array([[1], [2]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["G1"]))
    circ = ad.AnnData(X=np.array([[0], [3]], dtype=np.int32), obs=obs.copy(), var=pd.DataFrame(index=["C1"]))
    input_h5mu = tmp_path / "rna_circ.h5mu"
    mu.MuData({"rna": rna, "circ": circ}).write_h5mu(input_h5mu)
    candidates = tmp_path / "scomatic_candidate_summary.tsv"
    candidates.write_text(
        "\t".join(
            [
                "variant_id",
                "cell_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "gene",
                "filter_status",
                "candidate_variant_class",
                "read_support",
                "vaf",
                "caller",
            ]
        )
        + "\n"
        + "var1\tcellA\tchr1\t10\tA\tG\tGENE1\tPASS\tRNA-derived candidate variant signal\t4\t0.2\tSComatic\n",
        encoding="utf-8",
    )
    output = tmp_path / "with_scomatic.h5mu"

    result = runner.invoke(
        app,
        [
            "merge-scomatic",
            "--input",
            str(input_h5mu),
            "--scomatic-candidates",
            str(candidates),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert "candidate_snv" in payload["modalities"]
    merged = mu.read_h5mu(output)
    assert "candidate_snv" in merged.mod
    candidate = merged.mod["candidate_snv"]
    assert list(candidate.obs_names) == ["cellA", "cellB"]
    assert list(candidate.var_names) == ["chr1:10:A>G"]
    assert candidate.X.tolist() == [[4.0], [0.0]]
    assert candidate.layers["vaf"].tolist() == [[0.2], [0.0]]
