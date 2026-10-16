"""Retention/completion bookkeeping, not extra scientific configurations."""

import gzip
import json

from handoff import combined_prevalence, grid_complete, grid_rows, pending_outcome, prevalence_rows, sanitize


def test_complete_requires_record_and_successful_stage(tmp_path):
    records = tmp_path / "records"
    records.mkdir()
    plan = {"cases": [{"case": "control"}], "arms": ["circyto", "reference", "unsafe"]}
    observed = [{"case": "control", "arm": arm} for arm in plan["arms"]]
    assert not grid_complete(tmp_path, observed, plan)
    (records / "a04_complete.json").write_text(json.dumps({"passed": True}))
    (records / "a04-grid.json").write_text(json.dumps({"launcher_exit_code": 0}))
    assert grid_complete(tmp_path, observed, plan)
    (records / "a04_stop.json").write_text("{}")
    assert not grid_complete(tmp_path, observed, plan)


def test_prevalence_retains_and_flags_truncated_stream(tmp_path):
    p = tmp_path / "table.tsv.gz"
    data = gzip.compress(b"candidate\tnumerator\nA\t1\n")
    p.write_bytes(data)
    state = {}
    assert list(prevalence_rows(p, state)) == [{"candidate": "A", "numerator": "1"}]
    assert state["gzip_integrity"]
    truncated = data[:-5]
    p.write_bytes(truncated)
    list(prevalence_rows(p, state))
    assert not state["gzip_integrity"]
    assert p.read_bytes() == truncated


def test_path_aliases_are_applied_to_keys_and_values(tmp_path):
    root = tmp_path / "A02" / "evaluation"
    repo = tmp_path / "worktree"
    x = {str(root / "a"): [str(root.parent / "original"), str(repo / "script.py")]}
    assert sanitize(x, root, repo) == {
        "${EVAL_ROOT}/a": ["${A02_ROOT}/original", "${ANALYSIS_WORKTREE}/script.py"]
    }


def test_interrupted_attempt_is_not_called_not_started(tmp_path):
    assert pending_outcome(tmp_path, "case", "circyto")[0] == "NOT_STARTED"
    directory = tmp_path / "a04/cases/case/circyto"
    directory.mkdir(parents=True)
    assert pending_outcome(tmp_path, "case", "circyto")[0] == "STARTED_NO_COMPLETION_RECORD"
    (directory / "execution.json").write_text(json.dumps({"exit_code": 0}))
    assert pending_outcome(tmp_path, "case", "circyto")[0] == "EXECUTED_CHECK_INCOMPLETE"


def test_continued_case_retains_prior_attempt_and_merges_by_identity(tmp_path):
    for segment in ["a04", "a04_continued"]:
        folder = tmp_path / segment
        (folder / "cases/retried").mkdir(parents=True)
        (folder / "robustness_outcomes.tsv").write_text("case\tarm\noutside\treference\n" if segment == "a04" else "case\tarm\nretried\tcircyto\n")
    old = tmp_path / "a04/robustness_outcomes.tsv"
    old.write_text(old.read_text()+"retried\tcircyto\n")
    selected, prior, preferred = grid_rows(tmp_path)
    assert {(r["case"], r["segment"]) for r in selected} == {("outside", "a04"), ("retried", "a04_continued")}
    assert prior[0]["segment"] == "a04"
    for segment, rows in [("a04", "outside\tA\t2\nretried\tA\t99\n"), ("a04_continued", "retried\tA\t3\n")]:
        with gzip.open(tmp_path / segment / "completed_cohort_prevalence.tsv.gz", "wt") as f:
            f.write("case\tcandidate\tnumerator\n"+rows)
    integrity = {"error": "file absent"}
    merged = list(combined_prevalence(tmp_path, preferred, integrity))
    assert [(r["case"],r["numerator"]) for r in merged] == [("outside","2"),("retried","3")]
    assert integrity["gzip_integrity"]
    assert integrity["error"] is None
    assert integrity["sources"][0]["prior_attempt_rows"] == 1


def test_final_segment_retains_both_earlier_attempts_and_uses_final_completion(tmp_path):
    records = tmp_path / "records"
    records.mkdir()
    for segment in ["a04", "a04_continued", "a04_final"]:
        folder = tmp_path / segment
        (folder / "cases/retried").mkdir(parents=True)
        (folder / "robustness_outcomes.tsv").write_text("case\tarm\nretried\tcircyto\nretried\treference\nretried\tunsafe\n")
    selected, prior, preferred = grid_rows(tmp_path)
    assert len(selected) == 3 and len(prior) == 6
    assert preferred == {"retried": "a04_final"}
    assert {r["segment"] for r in selected} == {"a04_final"}
    planned = {"cases": [{"case": "retried"}], "arms": ["circyto", "reference", "unsafe"]}
    assert not grid_complete(tmp_path, selected, planned)
    (records / "a04_final_complete.json").write_text(json.dumps({"passed": True}))
    (records / "a04-final-grid.json").write_text(json.dumps({"launcher_exit_code": 0}))
    assert grid_complete(tmp_path, selected, planned)
    assert not grid_complete(tmp_path, [selected[0]]*3, planned)
    foreign = [dict(r, case="foreign") for r in selected]
    assert not grid_complete(tmp_path, foreign, planned)
    (records / "a04_continued_stop.json").write_text("{}")
    assert not grid_complete(tmp_path, selected, planned)
