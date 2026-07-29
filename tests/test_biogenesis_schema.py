from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from circyto.biogenesis import (
    BIOGENESIS_SCHEMA_VERSION,
    BiogenesisSchemaError,
    export_biogenesis_bundle,
    validate_biogenesis_bundle,
    validate_candidate_records,
    validate_cell_context_records,
    validate_observation_records,
)
from circyto.cli.circyto import app


FIXTURE_DIR = Path("testdata/biogenesis")
CANDIDATES = FIXTURE_DIR / "circRNA_candidates.v1.tsv"
CELL_CONTEXTS = FIXTURE_DIR / "cell_contexts.v1.tsv"
OBSERVATIONS = FIXTURE_DIR / "cell_circ_observations.v1.tsv"
runner = CliRunner()


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _fixture_frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return _read(CANDIDATES), _read(CELL_CONTEXTS), _read(OBSERVATIONS)


def test_synthetic_biogenesis_bundle_validates_positive_unlabelled_semantics() -> None:
    bundle = validate_biogenesis_bundle(*_fixture_frames())

    assert set(bundle.candidates["schema_version"]) == {BIOGENESIS_SCHEMA_VERSION}
    assert set(bundle.candidates["label_status"]) == {"positive", "unlabelled"}
    assert set(bundle.candidates["known_status"]) == {"known", "unknown"}
    undetected_unlabelled = bundle.observations[
        (bundle.observations["cell_id"] == "cellA")
        & (bundle.observations["circ_id"] == "circB")
    ].iloc[0]
    assert not bool(undetected_unlabelled["detected"])
    assert undetected_unlabelled["candidate_label_status"] == "unlabelled"


def test_candidate_validation_reports_missing_fields_and_invalid_coordinates() -> None:
    candidates, _, _ = _fixture_frames()
    with pytest.raises(
        BiogenesisSchemaError,
        match="missing required columns: reference_genome",
    ):
        validate_candidate_records(candidates.drop(columns=["reference_genome"]))

    invalid = candidates.copy()
    invalid.loc[0, "end"] = invalid.loc[0, "start"]
    with pytest.raises(BiogenesisSchemaError, match="end must be greater than start"):
        validate_candidate_records(invalid)

    invalid_one_based = candidates.copy()
    invalid_one_based.loc[0, ["coordinate_system", "start"]] = [
        "1-based-closed",
        0,
    ]
    with pytest.raises(BiogenesisSchemaError, match="require start >= 1"):
        validate_candidate_records(invalid_one_based)


def test_candidate_validation_rejects_negative_label_and_unsupported_version() -> None:
    candidates, _, _ = _fixture_frames()
    invalid_label = candidates.copy()
    invalid_label.loc[0, "label_status"] = "negative"
    with pytest.raises(BiogenesisSchemaError, match="label_status.*positive.*unlabelled"):
        validate_candidate_records(invalid_label)

    invalid_version = candidates.copy()
    invalid_version["schema_version"] = "2.0"
    with pytest.raises(BiogenesisSchemaError, match="expected '1.0'"):
        validate_candidate_records(invalid_version)


def test_context_validation_checks_qc_and_extensible_numeric_features() -> None:
    _, contexts, _ = _fixture_frames()
    invalid_qc = contexts.copy()
    invalid_qc.loc[0, "mitochondrial_fraction"] = 1.5
    with pytest.raises(BiogenesisSchemaError, match="mitochondrial_fraction"):
        validate_cell_context_records(invalid_qc)

    invalid_program = contexts.astype({"rbp_program__QKI": object}).copy()
    invalid_program.loc[0, "rbp_program__QKI"] = "not-a-number"
    with pytest.raises(BiogenesisSchemaError, match="rbp_program__QKI"):
        validate_cell_context_records(invalid_program)


def test_observation_validation_checks_detection_and_pair_uniqueness() -> None:
    _, _, observations = _fixture_frames()
    inconsistent = observations.copy()
    inconsistent.loc[0, "detected"] = False
    with pytest.raises(BiogenesisSchemaError, match="must equal whether count or bsj_support"):
        validate_observation_records(inconsistent)

    duplicated = pd.concat([observations, observations.iloc[[0]]], ignore_index=True)
    with pytest.raises(BiogenesisSchemaError, match="duplicate key"):
        validate_observation_records(duplicated)


def test_bundle_validation_checks_foreign_keys_and_provenance() -> None:
    candidates, contexts, observations = _fixture_frames()
    unknown = observations.copy()
    unknown.loc[0, "circ_id"] = "missing_circ"
    with pytest.raises(BiogenesisSchemaError, match="unknown circ_id"):
        validate_biogenesis_bundle(candidates, contexts, unknown)

    mismatched = observations.copy()
    mismatched.loc[0, "donor_id"] = "wrong_donor"
    with pytest.raises(
        BiogenesisSchemaError,
        match="observations.donor_id does not match cell_contexts.donor_id",
    ):
        validate_biogenesis_bundle(candidates, contexts, mismatched)


def test_export_biogenesis_bundle_round_trip_and_provenance(tmp_path: Path) -> None:
    outdir = tmp_path / "biogenesis"
    summary = export_biogenesis_bundle(
        candidates_path=CANDIDATES,
        cell_contexts_path=CELL_CONTEXTS,
        observations_path=OBSERVATIONS,
        outdir=outdir,
    )

    assert summary["schema_version"] == BIOGENESIS_SCHEMA_VERSION
    assert summary["record_counts"] == {
        "candidates": 2,
        "cell_contexts": 2,
        "observations": 4,
    }
    exported = validate_biogenesis_bundle(
        _read(outdir / "circRNA_candidates.v1.tsv"),
        _read(outdir / "cell_contexts.v1.tsv"),
        _read(outdir / "cell_circ_observations.v1.tsv"),
    )
    assert exported.observations.shape[0] == 4

    provenance = json.loads(
        (outdir / "biogenesis_provenance.v1.json").read_text(encoding="utf-8")
    )
    assert provenance["schema_versions"]["candidates"] == "1.0"
    assert provenance["provenance_values"] == {
        "datasets": ["dataset_synthetic"],
        "detectors": ["ciri3"],
        "donors": ["donor1", "donor2"],
        "protocols": ["ramda"],
        "reference_annotations": ["gencode.v47"],
        "reference_genomes": ["GRCh38"],
        "workflow_uuids": ["workflow-synthetic-001"],
    }
    assert provenance["label_semantics"]["non_detection_is_negative"] is False


def test_export_biogenesis_cli_and_overwrite_guard(tmp_path: Path) -> None:
    outdir = tmp_path / "cli_bundle"
    args = [
        "export-biogenesis",
        "--candidates",
        str(CANDIDATES),
        "--cell-contexts",
        str(CELL_CONTEXTS),
        "--observations",
        str(OBSERVATIONS),
        "--outdir",
        str(outdir),
    ]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["record_counts"]["observations"] == 4

    second = runner.invoke(app, args)
    assert second.exit_code != 0
    assert "Use --overwrite" in second.output

    overwritten = runner.invoke(app, [*args, "--overwrite"])
    assert overwritten.exit_code == 0, overwritten.output
