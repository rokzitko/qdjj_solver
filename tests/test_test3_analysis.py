"""Offline receipt/report regressions; no NRG, Mathematica, or pilot execution."""

from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

pytest.importorskip("resource", reason="synthetic task inventory imports the POSIX-only runner")

from NRG_comparisons.test3 import analyze as analysis
from NRG_comparisons.test3 import generate
from NRG_comparisons.test3.manual_run.physics import OBSERVABLES, PHYSICAL, finish_states, state_values
from NRG_comparisons.test3.run import digest, tasks


BATH = dict(xi=[.23], weights=[.019], delta=1., bandwidth=100.)
CHAIN = dict(xi_over_gap=[.2, .3, 999.], zeta_over_gap=[0., 0., 0.],
             delta_over_gap=[1., 1., 1.], kappa_over_gap=[0., 0., 0.],
             unused_last_bond_coefficient=True)
TARGETS = dict(qp=dict(absolute=1e-9, relative=1e-6), dmrg=dict(absolute=1e-9, relative=1e-6),
                nrg=dict(absolute=1e-6, relative=.001))


def test_reviewed_pilot_snapshot_records_unexpected_state_and_failed_dmrg_residuals():
    output = Path(__file__).resolve().parents[1] / "NRG_comparisons/test3/output"
    report = json.loads((output / "results.json").read_text())
    assert len(report["tasks"]) == 27
    assert all(task["status"] == "completed" for task in report["tasks"])
    assert report["gates"]["controls"] == "passed"
    assert report["gates"]["production"] == "deferred"
    assert not report["gates"]["continuum_target_established"]
    by_id = {task["id"]: task for task in report["tasks"]}
    for taskid in ("cosh8-qp-qfull", "nrg-k2000"):
        result = by_id[taskid]["result"]
        states = {row["id"]: row for row in result["states"]}
        assert set(result["coverage"]["subgap_ids"]) == {"even", "odd_down", "odd_up", "sz+0-r1"}
        assert not result["coverage"]["expected_three_subgap_states"]
        assert .93 < states["sz+0-r1"]["excitation"] < .94
        assert result["coverage"]["ground_ids"] == ["odd_down"]
        assert states["even"]["values"]["m_dot"] == pytest.approx(-states["even"]["values"]["m_bath"], abs=1e-10)
    assert by_id["cosh8-qp-qfull"]["qualification"]["finite_problem_qualified"]
    dmrg = by_id["cosh8-dmrg-chi128"]
    assert not dmrg["qualification"]["finite_problem_qualified"]
    assert all(row["residual"] > 1e-7 for row in dmrg["result"]["states"] if row["analysis"]["primary"])
    assert (output / "report.md").read_text() == analysis.report_markdown(report)
    assert "Additional Subgap States" in (output / "report.md").read_text()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def synthetic_result(backend, *, shift=0., moment=.1):
    rows = []
    for spin, root, energy in ((0, 0, .2), (-1, 0, 0.), (1, 0, .03), (0, 1, 1.2)):
        values = state_values(.8, .1, moment, .15, spin/2)
        rows.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=dict(values, energy=energy+shift),
                         degenerate=False, residual=None if backend == "nrg" else 1e-12,
                         **({} if backend == "nrg" else dict(finite_converged=True))))
    result = dict(states=rows, coverage=finish_states(rows), status="completed", monitored_seconds=2., peak_rss_gib=.25,
                  diagnostics=dict(residual_scope="projected_cutoff_hamiltonian"))
    if backend == "nrg":
        result["diagnostics"].update(wilson_chain=CHAIN, all_sector_spectra={
            "1": dict(twice_sz=0, indices_in_energy_order=list(range(3000)),
                      ground_relative_energy_over_gap=[.2, 1.2, 1.4] + [2.+i for i in range(2997)]),
            "0": dict(twice_sz=-1, indices_in_energy_order=[0, 1], ground_relative_energy_over_gap=[0., 1.5]),
        })
        result["numerical"] = {"lambda": 2., "z": 1., "nmax": 42}
    else:
        result["qualified"] = True
    return result


@pytest.fixture
def campaign(tmp_path):
    directory = tmp_path / "campaign"
    directory.mkdir()
    inventory = {stage: tasks(stage) for stage in ("controls", "pilot")}
    manifest = dict(identity=dict(
        inputs=dict(physical=PHYSICAL, targets=TARGETS,
                    limits=dict(cpus=4, memory_gib=16, task_seconds=1200, campaign_seconds=7200,
                                mathematica_instances=1)),
        tasks=inventory, source_hashes={"source.py": "frozen-source-hash"},
        tools={"nrg": dict(path="/not/executed", sha256="tool-hash")}, baths_sha256="bath-hash"),
        environment=dict(python="fixture"), solver=dict(version="fixture"))
    save_json(directory / "campaign.json", manifest)

    def emit(taskid, result=None):
        task = next(task for stage in inventory.values() for task in stage if task["id"] == taskid)
        work = directory / taskid
        work.mkdir(exist_ok=True)
        request = dict(task=task, physical=PHYSICAL | task.get("physical_changes", {}), tools={"nrg": "/not/executed"})
        if task["backend"] != "nrg":
            request["bath"] = dict(BATH, metadata=dict(chain=CHAIN)) if "bath_source" in task else BATH
        save_json(work / "request.json", request)
        result = deepcopy(synthetic_result(task["backend"]) if result is None else result)
        result["request_sha256"] = digest(request)
        save_json(work / "result.json", result)

    return directory, inventory, emit


def complete(campaign):
    directory, inventory, emit = campaign
    for stage in inventory.values():
        for task in stage:
            emit(task["id"])
    save_json(directory / "controls.json", dict(passed=True, checks=[
        dict(name="synthetic control", maximum_error=1e-12, tolerance=1e-8, passed=True)]))


def test_complete_snapshot_compact_provenance_and_all_saved_states(campaign, tmp_path):
    complete(campaign)
    directory, inventory, _ = campaign
    original = (directory / "nrg-k2000/result.json").read_bytes()
    report = analysis.analyze(directory)
    assert report["gates"]["controls"] == "passed"
    assert report["gates"]["pilot_escalation_allowed"]
    assert report["gates"]["production"] == "deferred"
    assert not report["gates"]["continuum_target_established"]
    assert not report["gates"]["spectral_completeness_established"]
    assert report["gates"]["uncertainty_kind"] == "not_established"
    assert len(report["tasks"]) == sum(map(len, inventory.values()))
    task = next(task for task in report["tasks"] if task["id"] == "nrg-k2000")
    assert task["provenance"]["result_sha256"] == hashlib.sha256(original).hexdigest()
    request = json.loads((directory / "nrg-k2000/request.json").read_text())
    assert task["provenance"]["request_sha256"] == digest(request)
    assert task["provenance"]["request_file_sha256"] != digest(request)
    assert task["provenance"]["source_hashes"] == {"source.py": "frozen-source-hash"}
    diagnostics = task["result"]["diagnostics"]
    assert "all_sector_spectra" not in diagnostics
    assert diagnostics["wilson_chain"] == CHAIN
    summary = diagnostics["sector_spectrum_summary"]["1"]
    assert summary["root_count"] == 3000
    assert summary["first_two_excitations"] == [.2, 1.2]
    assert summary["first_omitted_root"] == dict(root=2, excitation=1.4)
    assert summary["nearest_excluded_root"] == dict(root=1, excitation=1.2)
    assert summary["nearest_gap_edge_root"]["distance"] == pytest.approx(.2)
    assert task["qualification"]["finite_problem_qualified"] is None
    assert all(row["residual"] is None for row in task["result"]["states"])
    assert task["result"]["coverage"]["expected_three_subgap_states"]
    assert report["costs"]["total_monitored_seconds"] == 2*len(report["tasks"])
    output = tmp_path / "presentation"
    analysis.write_report(report, output)
    assert {path.name for path in output.iterdir()} == set(analysis.OUTPUT_NAMES)
    published = json.loads((output / "results.json").read_text())
    assert published == generate.public_projection(report)
    assert "campaign" not in published
    assert published["inputs"]["targets"] == TARGETS
    assert "tools" not in published["tasks"][0]["request"]
    assert published["tasks"][0]["provenance"]["public_request_sha256"] == generate.digest(
        published["tasks"][0]["request"])
    with (output / "observables.csv").open() as stream:
        observations = list(csv.DictReader(stream))
    assert len(observations) == 4*len(report["tasks"])
    assert len([row for row in observations if row["primary"] == "False"]) == len(report["tasks"])
    assert set(OBSERVABLES) <= observations[0].keys()
    markdown = (output / "report.md").read_text()
    assert "Czz_connected" in markdown and "P_up" in markdown
    assert "unmatched-bath diagnostic" in markdown and "not a reference" in markdown
    assert (directory / "nrg-k2000/result.json").read_bytes() == original


def test_comparisons_all_axes_new_backend_normalization_and_controls(campaign):
    complete(campaign)
    directory, _, emit = campaign
    emit("cosh8-dmrg-chi128", synthetic_result("dmrg", shift=.001, moment=0.))
    emit("nrg-k2000", synthetic_result("nrg", shift=.002, moment=.2))
    emit("atomic-qp", synthetic_result("qp", shift=5e-9))
    report = analysis.analyze(directory)
    rows = report["comparisons"]
    assert {row["kind"] for row in rows} == {
        "finite_control", "qp_cutoff", "qp_full_finite_bath", "dmrg_bond", "dmrg_vs_full_qp",
        "nrg_retention", "unmatched_bath_single_z"}
    for kind in {row["kind"] for row in rows}:
        assert {row["observable"] for row in rows if row["kind"] == kind} == set(analysis.COMPARISON_OBSERVABLES)
    row = next(row for row in rows if row["kind"] == "dmrg_vs_full_qp" and row["new_task"] == "cosh8-dmrg-chi128"
               and row["observable"] == "m_dot" and row["state"] == "even")
    assert row["new"] == 0 and row["old"] == .1
    assert row["target_backend"] == "dmrg" and row["target_tolerance"] == 1e-9
    assert row["normalized_change"] == pytest.approx(1e8)
    assert not row["within_target"] and row["finite_control_tolerance"] is None
    cross = next(row for row in rows if row["kind"] == "unmatched_bath_single_z" and row["observable"] == "m_dot")
    assert cross["new_task"] == "nrg-k2000" and cross["target_backend"] == "nrg"
    assert cross["target_tolerance"] == pytest.approx(1e-6+.001*.2)
    assert not cross["matched_bath"] and "not a reference or uncertainty" in cross["reason"]
    assert all(row["finite_control_tolerance"] == 1e-8 for row in rows if row["kind"] == "finite_control")
    control = next(row for row in rows if row["new_task"] == "atomic-qp" and row["observable"] == "energy"
                   and row["state"] == "odd_down")
    assert control["finite_control_passed"] and not control["within_target"]
    # ED is allowed as the new comparison endpoint, with explicitly QP diagnostic targets.
    by_id = {task["id"]: task for task in report["tasks"]}
    ed_rows = analysis.comparison_rows("finite_control", by_id["chain-qp"], by_id["chain-ed"], TARGETS, True)
    assert all(row["target_backend"] == "qp" and row["finite_control_passed"] for row in ed_rows)


def test_failure_incomplete_unstarted_and_missing_controls(campaign, tmp_path):
    directory, inventory, emit = campaign
    emit("atomic-nrg", dict(status="failed", error="resource ceiling", monitored_seconds=1200., peak_rss_gib=15.))
    (directory / "atomic-ed").mkdir()
    report = analysis.analyze(directory)
    by_id = {task["id"]: task for task in report["tasks"]}
    assert by_id["atomic-nrg"]["status"] == "failed"
    assert by_id["atomic-ed"]["status"] == "incomplete"
    assert by_id["nrg-k2000"]["status"] == "unstarted"
    assert len(by_id) == sum(map(len, inventory.values()))
    assert report["gates"]["controls"] == "blocked_controls"
    assert not report["gates"]["pilot_escalation_allowed"]
    assert all(row["status"] == "unavailable" and row["within_target"] is None for row in report["comparisons"])
    assert report["costs"]["total_monitored_seconds"] == 1200.
    analysis.write_report(report, tmp_path / "partial")
    markdown = (tmp_path / "partial/report.md").read_text()
    assert "blocked_controls" in markdown and "resource ceiling" in markdown and "unstarted" in markdown


def test_empty_campaign_and_bad_json_are_reportable(campaign):
    directory, _, _ = campaign
    empty = analysis.analyze(directory)
    assert empty["costs"]["total_monitored_seconds"] is None
    assert empty["costs"]["maximum_peak_rss_gib"] is None
    work = directory / "atomic-nrg"
    work.mkdir()
    (work / "request.json").write_text("{")
    (work / "result.json").write_text("{")
    (directory / "controls.json").write_text("{")
    report = analysis.analyze(directory)
    assert len(report["tasks"][0]["provenance"]["issues"]) == 2
    assert "unreadable" in report["gates"]["controls_issues"][0]


@pytest.mark.parametrize("corruption", ["hash", "request", "missing", "control"])
def test_integrity_and_control_receipt_block_escalation(campaign, corruption):
    complete(campaign)
    directory, _, _ = campaign
    work = directory / "atomic-nrg"
    if corruption == "hash":
        result = json.loads((work / "result.json").read_text())
        result["request_sha256"] = "wrong"
        save_json(work / "result.json", result)
    elif corruption == "request":
        request = json.loads((work / "request.json").read_text())
        request["physical"]["field"] = .04
        save_json(work / "request.json", request)
        result = json.loads((work / "result.json").read_text())
        result["request_sha256"] = digest(request)
        save_json(work / "result.json", result)
    elif corruption == "missing":
        (work / "request.json").unlink()
    else:
        save_json(directory / "controls.json", dict(passed=True, checks=[dict(name="failed", passed=False)]))
    report = analysis.analyze(directory)
    assert report["gates"]["controls"] == "blocked_controls"
    if corruption != "control":
        rows = [row for row in report["comparisons"] if row["new_task"] == "atomic-nrg"]
        assert all(row["status"] == "unavailable" and row["finite_control_passed"] is None for row in rows)


def test_primary_and_exclusion_qualification_and_symmetry_ties(campaign):
    directory, _, emit = campaign
    result = synthetic_result("dmrg")
    result["states"][-1]["finite_converged"] = False
    result["states"][-1]["residual"] = .2
    # A symmetry-resolved cross-sector tie does not make either vector ambiguous.
    result["states"][2]["values"]["energy"] = 0.
    result["coverage"] = finish_states(result["states"])
    emit("cosh8-dmrg-chi32", result)
    emit("cosh8-qp-qfull")
    result["states"][0]["degenerate"] = True
    result["qualified"] = False
    emit("cosh8-dmrg-chi64", result)
    result["states"][0]["degenerate"] = False
    result["states"][0]["finite_converged"] = False
    result["states"][0]["residual"] = .3
    emit("cosh8-dmrg-chi128", result)
    report = analysis.analyze(directory)
    by_id = {task["id"]: task for task in report["tasks"]}
    first = by_id["cosh8-dmrg-chi32"]
    assert first["qualification"]["finite_problem_qualified"]
    assert first["result"]["states"][-1]["analysis"]["qualification"] == "not_qualified"
    assert not first["result"]["states"][2]["degenerate"]
    assert not by_id["cosh8-dmrg-chi64"]["qualification"]["finite_problem_qualified"]
    rows = [row for row in report["comparisons"] if row["kind"] == "dmrg_vs_full_qp" and row["state"] == "even"]
    ambiguous = next(row for row in rows if row["new_task"] == "cosh8-dmrg-chi64" and row["observable"] == "m_dot")
    assert ambiguous["status"] == "basis_dependent" and ambiguous["within_target"] is None
    energy = next(row for row in rows if row["new_task"] == "cosh8-dmrg-chi64" and row["observable"] == "energy")
    assert energy["status"] == "compared"  # Degeneracy does not make energies basis dependent.
    unqualified = next(row for row in rows if row["new_task"] == "cosh8-dmrg-chi128" and row["observable"] == "energy")
    assert unqualified["status"] == "unqualified" and not unqualified["qualified_pair"]


@pytest.mark.parametrize("taskid,coefficient", [("chain-qp", "weights"), ("chain-nrg", "xi_over_gap")])
def test_claimed_matched_controls_require_actual_bath_evidence(campaign, taskid, coefficient):
    complete(campaign)
    directory, _, _ = campaign
    work = directory / taskid
    request = json.loads((work / "request.json").read_text())
    result = json.loads((work / "result.json").read_text())
    if coefficient == "weights":
        request["bath"][coefficient][0] *= 2
        save_json(work / "request.json", request)
        result["request_sha256"] = digest(request)
    else:
        result["diagnostics"]["wilson_chain"][coefficient][0] *= 2
    save_json(work / "result.json", result)
    report = analysis.analyze(directory)
    rows = [row for row in report["comparisons"] if row["new_task"] == taskid]
    assert all(row["status"] == "unavailable" and not row["matched_bath"] for row in rows)
    assert report["gates"]["controls"] == "blocked_controls"


def test_highest_available_nrg_is_not_assumed_highest_planned(campaign):
    directory, _, emit = campaign
    emit("cosh8-qp-qfull")
    emit("nrg-k500")
    emit("nrg-k1000")
    emit("nrg-k2000", dict(status="failed", error="timeout"))
    report = analysis.analyze(directory)
    rows = [row for row in report["comparisons"] if row["kind"] == "unmatched_bath_single_z"]
    assert rows and {row["new_task"] for row in rows} == {"nrg-k1000"}


def test_output_cli_safety_and_stdlib_execution(campaign, tmp_path, monkeypatch):
    directory, _, emit = campaign
    emit("atomic-ed")
    fake_here = tmp_path / "test3"
    monkeypatch.setattr(analysis, "HERE", fake_here)
    analysis.main([str(directory)])
    output = fake_here / "runs/presentation"
    before = (output / "results.json").read_bytes()
    with pytest.raises(SystemExit, match="2"):
        analysis.main([str(directory)])
    assert (output / "results.json").read_bytes() == before
    analysis.main([str(directory), "--output-dir", str(output)])
    with pytest.raises(ValueError, match="outside"):
        analysis.main([str(directory), "--output-dir", str(directory / "nested")])
    assert not (directory / "nested").exists()
    # -S removes installed dependencies: the analysis path must still work.
    completed = subprocess.run([sys.executable, "-B", "-S", str(Path(analysis.__file__)), str(directory),
                                "--output-dir", str(tmp_path / "stdlib")], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert "production deferred" in completed.stdout


@pytest.mark.parametrize("change", ["duplicate", "path"])
def test_invalid_inventory_is_rejected(campaign, change):
    directory, _, _ = campaign
    manifest = json.loads((directory / "campaign.json").read_text())
    task = manifest["identity"]["tasks"]["controls"][0]
    if change == "duplicate":
        manifest["identity"]["tasks"]["controls"].append(task)
    else:
        task["id"] = "../escape"
    save_json(directory / "campaign.json", manifest)
    with pytest.raises(ValueError, match="duplicate|basename"):
        analysis.analyze(directory)


def test_compiled_finite_receipt_when_available(campaign, tmp_path):
    pytest.importorskip("qdjj_solver")
    from NRG_comparisons.test3.finite import run_finite

    directory, inventory, emit = campaign
    task = next(task for task in inventory["controls"] if task["id"] == "atomic-qp")
    # One explicit bath level, tiny dense sectors, and no external process.
    result = run_finite(task, dict(PHYSICAL, gamma=0.), BATH)
    emit(task["id"], dict(result, status="completed"))
    report = analysis.analyze(directory)
    actual = next(task for task in report["tasks"] if task["id"] == "atomic-qp")
    assert actual["qualification"]["finite_problem_qualified"]
    even = next(row for row in actual["result"]["states"] if row["id"] == "even")
    assert even["values"]["energy"] == pytest.approx(.5, abs=1e-12)
    assert even["analysis"]["finite_problem_qualified"]
    analysis.write_report(report, tmp_path / "compiled-receipt")
