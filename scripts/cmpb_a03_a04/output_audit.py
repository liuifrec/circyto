"""Additional output-only audit; independent CSV/Matrix Market/HTML reader.

This does not change either arm, the perturbations or frozen grid outcomes.
It checks summary paragraphs, declared counts and the separate exported matrix.
"""

import argparse
from collections import Counter
import csv
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import tarfile

import h5py

from run_stage import digest, dump, stamp


def rows(data):
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))


def names(data):
    return data.decode().splitlines()


def exact_matrix(data, cells, features):
    """Read coordinate integer Matrix Market without scientific-library coercion."""
    lines = data.decode().splitlines()
    if lines[0].lower() != "%%matrixmarket matrix coordinate integer general":
        raise ValueError("Expected integer coordinate matrix")
    content = [s for s in lines[1:] if s.strip() and not s.startswith("%")]
    nr, nc, nnz = map(int, content[0].split())
    if (nr, nc, nnz) != (len(features), len(cells), len(content) - 1):
        raise ValueError("Matrix shape or record count mismatch")
    values = {}
    seen = set()
    for line in content[1:]:
        r, c, value = map(int, line.split())
        if not 1 <= r <= nr or not 1 <= c <= nc or value < 0:
            raise ValueError("Invalid coordinate/support")
        key = (cells[c - 1], features[r - 1])
        if key in seen:
            raise ValueError("Duplicate matrix record")
        seen.add(key)
        if value:
            values[key] = value
    return values


def truth(config, source):
    metadata = rows(source(config["cells"]))
    cells = [r["cell_id"] for r in metadata]
    states = {r["cell_id"]: r.get("detector_status", "unknown") for r in metadata}
    values, features = {}, []
    if config["mode"] == "matrix":
        folder = Path(config["input"])
        source_cells = names(source(folder / "cell_index.txt"))
        features = names(source(folder / "circ_index.txt"))
        values = exact_matrix(source(folder / "circ_counts.mtx"), source_cells, features)
    elif config["mode"] == "native":
        manifest = Path(config["input"])
        for mapping in rows(source(manifest)):
            path = manifest.parent / mapping["path"]
            for row in rows(source(path)):
                feature = row["circRNA_ID"]
                if feature not in features:
                    features.append(feature)
                value = int(row["#junction_reads"])
                if value:
                    values[mapping["cell_id"], feature] = value
    else:
        for row in rows(source(config["input"])):
            feature = row["circ_id"]
            if feature not in features:
                features.append(feature)
            if int(row["support"]):
                values[row["cell_id"], feature] = int(row["support"])
    rna_ids = []
    if config.get("rna"):
        with h5py.File(io.BytesIO(source(config["rna"])), "r") as h:
            key = h["obs"].attrs["_index"]
            rna_ids = [v.decode() if isinstance(v, bytes) else str(v) for v in h["obs"][key][:]]
    return cells, states, features, values, rna_ids


class Paragraphs(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.paragraphs = []
        self.current = None
        self.remote_resources = []

    def handle_starttag(self, tag, attrs):
        if tag == "p":
            self.current = []
        for key, value in attrs:
            if key in {"src", "href"} and value and value.startswith(("http:", "https:", "//")):
                self.remote_resources.append(value)

    def handle_data(self, data):
        if self.current is not None:
            self.current.append(data)

    def handle_endtag(self, tag):
        if tag == "p" and self.current is not None:
            self.paragraphs.append(" ".join("".join(self.current).split()))
            self.current = None


def check(config, source, output, arm):
    cells, states, features, expected_values, rna = truth(config, source)
    checks = {}
    exported_cells = names(output("matrix/cell_index.txt"))
    exported_features = names(output("matrix/circ_index.txt"))
    got = exact_matrix(output("matrix/circ_counts.mtx"), exported_cells, exported_features)
    checks["exported_matrix_axes"] = exported_cells == cells and exported_features == features
    checks["exported_matrix_all_integer_values"] = got == expected_values
    counts = Counter(c for c, f in expected_values)
    completed = [c for c in cells if states[c] == "completed"]
    den = len(completed)
    positive = sum(counts[c] > 0 for c in completed)
    metrics = json.loads(output("qc/metrics.json"))
    wanted = dict(n_cells=len(cells), n_candidates=len(features), nnz=len(expected_values), total_support=sum(expected_values.values()))
    checks["json_summary_counts"] = all(metrics.get(k) == v for k, v in wanted.items())
    observations = Counter(
        ("completed_detected" if counts[c] else "completed_zero")
        if states[c] == "completed" else states[c] for c in cells
    )
    checks["json_observation_states"] = metrics.get("observation_status_counts") == dict(observations)
    for name, numerator in [("detection_fraction_completed", positive), ("zero_fraction_completed", den-positive)]:
        checks["json_"+name] = metrics.get(name) == dict(numerator=numerator, denominator=den, value=numerator/den if den else None)
    declarations = {k: config[k] for k in ["genome_build", "coordinate_convention", "support_definition", "detector"]}
    checks["json_declarations"] = metrics.get("declarations") == declarations
    overlap = dict(n_shared=len(set(cells)&set(rna)), n_circ_only=len(set(cells)-set(rna)), n_rna_only=len(set(rna)-set(cells)))
    if config.get("rna"):
        checks["json_RNA_overlap"] = all(metrics.get("rna_overlap", {}).get(k) == v for k, v in overlap.items())
    document = output("qc/report.html").decode()
    html = Paragraphs()
    html.feed(document)
    checks["offline_no_remote_resource"] = not html.remote_resources and not re.search(r"(?:url\(|@import\s+)[\s\"']*(?:https?:)?//", document)
    if arm == "circyto":
        summary = f"{len(cells)} metadata cells; {len(features)} candidates; {sum(expected_values.values())} total support."
        rate = f"{positive} / {den}" if den else "unavailable"
        checks["html_summary_counts"] = summary in html.paragraphs
        checks["html_completed_fraction"] = f"Cells with candidates / explicitly completed cells: {rate}." in html.paragraphs
        checks["html_declarations"] = "Genome and support declarations (user supplied): "+json.dumps(declarations,sort_keys=True) in html.paragraphs
        overlap_paragraphs = [p[len("RNA overlap: "):] for p in html.paragraphs if p.startswith("RNA overlap: ")]
        checks["html_RNA_overlap"] = len(overlap_paragraphs)==1 and all(json.loads(overlap_paragraphs[0]).get(k)==v for k,v in overlap.items())
    else:
        summary = f"{len(cells)} cells; {len(features)} candidates; {sum(expected_values.values())} support. Completed denominator {den}."
        checks["html_summary_counts_and_denominator"] = summary in html.paragraphs
    return checks


def archive_payload(path):
    """One streaming pass; skip large output HDF5 without copying it."""
    payload = {}
    with tarfile.open(path, "r|gz") as t:
        for member in t:
            name = member.name
            if member.isfile() and (Path(name).suffix in {".json", ".tsv", ".txt", ".mtx", ".html"} or name == "inputs/rna.h5ad"):
                payload[name] = t.extractfile(member).read()
    return payload


def audit(root):
    record = root / "records/additional_output_audit.json"
    if record.exists():
        raise FileExistsError(record)
    dump(root / "records/additional_output_audit_source.json", dict(recorded_utc=stamp(), source_sha256=digest(__file__)))
    result = []
    a03_config = json.loads((root / "a03/config.json").read_text())
    for directory in sorted((root / "a03").iterdir()):
        if not (directory / "execution.json").exists():
            continue
        arm = directory.name.rsplit("-",1)[1]
        checks = check(a03_config, lambda p: Path(p).read_bytes(), lambda p: (directory / "output" / p).read_bytes(), arm)
        result.append(dict(case="A03-"+directory.name, arm=arm, checks=checks, passed=all(checks.values())))
        if not all(checks.values()):
            dump(record, dict(passed=False, stopped=True, outcomes=result, reason="A03 retained-output inconsistency; review required"))
            raise RuntimeError("STOP: additional A03 output audit failure")
    audited_cases = 0
    for segment in ["a04", "a04_continued", "a04_final"]:
        if not (root / segment / "robustness_outcomes.tsv").exists():
            continue
        outcomes = rows((root / segment / "robustness_outcomes.tsv").read_bytes())
        by_case = {}
        for row in outcomes:
            if row["outcome"] == "accepted_valid":
                by_case.setdefault(row["case"], []).append(row["arm"])
        for case, arms in by_case.items():
            folder = root / segment / "cases" / case
            archive = folder.with_suffix(".tar.gz")
            if archive.exists() and folder.with_suffix(".tar.manifest.json").exists():
                payload = archive_payload(archive)
                config = json.loads(payload["inputs/config.json"])
                source_root = str(Path(config["cells"]).parent)
                source = lambda p: payload["inputs/" + str(Path(p).relative_to(source_root))]
                output = lambda arm, p: payload[arm+"/output/"+p]
            elif folder.is_dir():
                config = json.loads((folder / "inputs/config.json").read_text())
                source = lambda p: Path(p).read_bytes()
                output = lambda arm, p: (folder / arm / "output" / p).read_bytes()
            else:
                raise RuntimeError("Missing retained case: "+case)
            for arm in arms:
                checks = check(config, source, lambda p: output(arm, p), arm)
                result.append(dict(case=segment+":"+case, arm=arm, checks=checks, passed=all(checks.values())))
                if arm != "unsafe" and not all(checks.values()):
                    dump(record, dict(passed=False, stopped=True, outcomes=result, reason="Unexpected evaluated-arm output inconsistency; review required"))
                    raise RuntimeError("STOP: additional output audit failure")
            audited_cases += 1
            if audited_cases % 25 == 0:
                progress = dict(recorded_utc=stamp(), archive_cases_checked=audited_cases,
                                accepted_outputs_checked=len(result), current_case=segment+":"+case)
                dump(root / "records/additional_output_progress.json", progress)
                dump(root / "records/additional_output_audit_partial.json", dict(
                    complete=False, progress=progress, outcomes=result))
                print(json.dumps(progress), flush=True)
    dump(record, dict(passed=True, stopped=False, outcomes=result, finished_utc=stamp(), scope="Additional check of retained outputs; frozen main outcomes remain unchanged"))
    print(json.dumps(dict(audited=len(result), evaluated_arm_failures=sum(not r["passed"] for r in result if r["arm"]!="unsafe"))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve())
