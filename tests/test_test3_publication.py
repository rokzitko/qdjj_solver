"""Offline, stdlib-only public dataset checks; no raw campaign or backend required."""

from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

from NRG_comparisons.test3 import generate


CASE = Path(__file__).resolve().parents[1] / "NRG_comparisons/test3"


@pytest.fixture(scope="module")
def report():
    return json.loads((CASE / "output/results.json").read_text())


def test_public_projection_inventory_science_and_hashes(report):
    assert generate.public_projection(report) == report
    assert len(report["tasks"]) == 27
    assert sum(t["stage"] == "controls" for t in report["tasks"]) == 17
    assert sum(len(t["result"]["states"]) for t in report["tasks"]) == 330
    assert len(report["comparisons"]) == 858
    assert len(report["controls"]["checks"]) == 42
    assert all(check["passed"] for check in report["controls"]["checks"])
    assert report["inputs"] == json.loads((CASE / "input/physical.json").read_text())
    assert report["gates"]["production"] == "deferred"
    assert not report["gates"]["continuum_target_established"]
    assert not report["gates"]["spectral_completeness_established"]
    assert report["gates"]["uncertainty_kind"] == "not_established"
    # Pinned against the verified original archive, excluding only provenance
    # and private execution fields. This includes every nested solver diagnostic.
    science = {key: report[key] for key in ("inputs", "controls", "gates", "costs", "limitations", "comparisons")}
    science["tasks"] = [{key: value for key, value in task.items() if key != "provenance"}
                        for task in report["tasks"]]
    assert generate.digest(science) == "29830e98aa3ee0c1b3f02750c621649c21cc3f9d30228e1b5c95ac2425671c13"
    for task in report["tasks"]:
        assert task["status"] == "completed"
        for key in ("request", "result"):
            assert task["provenance"][f"public_{key}_sha256"] == generate.digest(task[key])
        assert task["provenance"]["request_sha256"] != task["provenance"]["public_request_sha256"]
        for row in task["result"]["states"]:
            v = row["values"]
            assert set(v) == set(generate.OBSERVABLES)
            assert v["P0"] + v["P1"] + v["P2"] == pytest.approx(1, abs=1e-8)
            assert v["m_dot"] + v["m_bath"] == pytest.approx(row["twice_sz"] / 2, abs=1e-8)
            assert v["P_up"] == pytest.approx(v["P1"] / 2 + v["m_dot"], abs=1e-8)
            assert v["Czz_connected"] == pytest.approx(v["m_dot"] ** 2 - v["P1"] / 4, abs=1e-8)
            if task["backend"] == "nrg":
                assert row["residual"] is None
                assert row["analysis"]["finite_problem_qualified"] is None
        if task["backend"] == "nrg":
            chain = task["result"]["diagnostics"]["wilson_chain"]
            assert len(chain["xi_over_gap"]) == len(chain["delta_over_gap"])
    provenance = json.loads((CASE / "reference/provenance.json").read_text())
    for name, digest in provenance["published_file_sha256"].items():
        assert hashlib.sha256((CASE / name).read_bytes()).hexdigest() == digest
    assert provenance["recorded_original_file_sha256"]["output/results.json"] != provenance[
        "published_file_sha256"]["output/results.json"]
    assert provenance["recorded_source_hashes"]["NRG_comparisons/test3/analyze.py"] != provenance[
        "published_file_sha256"]["analyze.py"]


def test_standalone_stdlib_regeneration_and_read_only_outputs(tmp_path):
    # Ship only the generator and data into an otherwise empty tree. Even the
    # physics helper and analyzer are absent, as are ignored runs and backends.
    isolated = tmp_path / "standalone"
    isolated.mkdir()
    (isolated / "output").mkdir()
    shutil.copy2(CASE / "generate.py", isolated / "generate.py")
    shutil.copy2(CASE / "output/results.json", isolated / "output/results.json")
    command = [sys.executable, "-B", "-S", str(isolated / "generate.py"), "--output-dir"]
    destination = tmp_path / "regenerated"
    completed = subprocess.run([*command, str(destination)], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert completed.returncode == 0, completed.stderr
    assert {p.name for p in destination.iterdir()} == set(generate.OUTPUT_NAMES)
    for name in generate.OUTPUT_NAMES:
        assert (destination / name).read_bytes() == (CASE / "output" / name).read_bytes()
    alias = tmp_path / "alias"
    alias.symlink_to(isolated / "output", target_is_directory=True)
    file_alias = tmp_path / "file-alias"
    file_alias.mkdir()
    (file_alias / "report.md").symlink_to(isolated / "output/report.md")
    for path in (isolated / "output", alias, file_alias):
        completed = subprocess.run([*command, str(path)], cwd=tmp_path, capture_output=True, text=True, timeout=30)
        assert completed.returncode == 2 and "published outputs are read-only" in completed.stderr


def test_csvs_retain_exact_operands_all_states_and_four_state_report(report):
    with (CASE / "output/observables.csv").open(newline="") as stream:
        observations = list(csv.DictReader(stream))
    assert len(observations) == 330
    by_id = {task["id"]: task for task in report["tasks"]}
    for observation in observations:
        state = next(row for row in by_id[observation["task"]]["result"]["states"] if row["id"] == observation["state"])
        assert {key: float(observation[key]) for key in generate.OBSERVABLES} == state["values"]
    with (CASE / "output/comparisons.csv").open(newline="") as stream:
        comparisons = list(csv.DictReader(stream))
    for actual, expected in zip(comparisons, report["comparisons"], strict=True):
        assert actual == {key: "" if value is None else str(value) for key, value in expected.items()}
    text = (CASE / "output/report.md").read_text()
    assert text == generate.report_markdown(report)
    upfront = text.split("## Supplemental material")[0]
    assert "PILOT ONLY" in upfront and "not a reference or z average" in upfront
    for taskid in ("cosh8-qp-qfull", "cosh8-dmrg-chi128", "nrg-k2000"):
        for state in ("odd_down", "odd_up", "even", "sz+0-r1"):
            assert f"| {taskid} | {state} |" in upfront
    for row in by_id["cosh8-qp-qfull"]["result"]["states"]:
        if row["id"] in ("odd_down", "odd_up", "even", "sz+0-r1"):
            assert row["analysis"]["finite_problem_qualified"] and not row["degenerate"]
    for row in by_id["cosh8-dmrg-chi128"]["result"]["states"]:
        if row["id"] in ("odd_down", "odd_up", "even", "sz+0-r1"):
            assert not row["analysis"]["finite_problem_qualified"] and row["residual"] > 1e-7


def test_allowlist_drops_recursive_execution_records_without_changing_science(report):
    dirty = deepcopy(report)
    secret = dict(username="private-operator", hostname="private-node", uid=99999,
                  session="2030-01-02", execution={"path": "/home/private-operator/work"})
    dirty["execution"] = secret
    dirty["provenance"]["identity"] = secret
    dirty["tasks"][0]["request"]["tools"] = secret
    dirty["tasks"][0]["result"]["diagnostics"]["session"] = secret
    dirty["tasks"][0]["result"]["states"][0]["execution"] = secret
    dirty["controls"]["checks"][0]["identity"] = secret
    assert generate.public_projection(dirty) == report


@pytest.mark.parametrize("payload", ["/home/private-operator/work/result.json", "/project/private-work/results",
                                     "/srv/private-operator/receipts", "C:\\private-operator\\receipts",
                                     "private-operator@private-node", "started 2030-01-02"])
def test_allowed_text_cannot_smuggle_private_metadata(report, tmp_path, payload):
    dirty = deepcopy(report)
    dirty["tasks"][0]["provenance"]["issues"].append(payload)
    destination = tmp_path / "rejected"
    with pytest.raises(ValueError, match="private execution metadata"):
        generate.write_report(dirty, destination)
    assert not destination.exists()


def test_public_assets_privacy_headings_and_local_links(report):
    generate.privacy_guard(report)
    for path in [*(CASE / "output").iterdir(), CASE / "README.md", CASE / "reference/provenance.json"]:
        text = path.read_text()
        assert not re.search(r"/(?:project|home|Users|opt)/|qdjj_solver_old|\b\d{4}-\d{2}-\d{2}\b|"
                             r'"(?:hostname|uid|username|session|execution|identity|started_unix)"\s*:', text), path
    for path in (CASE / "README.md", CASE / "output/report.md"):
        text = path.read_text()
        headings = ["## TL;DR", "## Model", "## Energies and Gap", "## Probabilities and Moment",
                    "## Achieved Accuracy", "## Data", "## Supplemental material"]
        assert [text.index(h) for h in headings] == sorted(text.index(h) for h in headings)
        for link in re.findall(r"\[[^\]]*\]\(([^\s)]+)\)", text):
            if not re.match(r"[a-zA-Z][\w+.-]*:", link) and not link.startswith("#"):
                assert (path.parent / link.split("#", 1)[0]).exists(), link
