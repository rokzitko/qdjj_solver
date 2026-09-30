"""Bounded checks of retained inputs and fresh-directory reproduction commands."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

from NRG_comparisons.test1 import control, focused_internal, nrg, run, targeted_nrg


HERE = Path(__file__).resolve().parents[1] / "NRG_comparisons/test1"


def test_final_inputs_cover_every_published_endpoint():
    manifest, _ = run.load_inputs()
    published = json.loads((HERE / "output/results.json").read_text())
    spec = json.loads((HERE / "input/qp-dmrg-focused.json").read_text())
    baths = json.loads((HERE / "input/final-baths.json").read_text())
    assert manifest["physical"] == published["model"]
    assert manifest["targets"] == published["targets"]
    assert [e["id"] for e in spec["tasks"]] == list(published["finite_results"])
    assert len(baths) == 7
    assert "runs/" not in json.dumps(spec)
    for entry in spec["tasks"]:
        task = entry["calculation"]
        bath = baths[entry["bath_input"]]
        focused_internal.validate_bath(manifest["physical"], task["bath"], bath)
        focused_internal.preflight(task)
        assert bath["xi"] == [-x for x in reversed(bath["xi"])]
        assert bath["weights"] == list(reversed(bath["weights"]))
        assert min(bath["weights"]) > 0
        expected = published["finite_results"][entry["id"]]["settings"]
        assert expected == dict(backend=task["backend"], bath_family=task["bath"]["kind"],
                               bath_levels=task["bath"]["levels"], qp_cutoff=task.get("cutoff"),
                               bond_dimension=task.get("chi"), frequency_window=task["bath"].get("frequency_cutoff"))
    # The window-100 surrogate is deliberately not normalized to unity.
    assert sum(baths["surrogate24-window100"]["weights"]) == pytest.approx(1 - 5.293473890199962e-5, abs=1e-15)


def test_final_nrg_plan_contains_all_checks_without_baseline():
    spec = json.loads((HERE / "input/nrg-targeted.json").read_text())
    plan = targeted_nrg.plan_for("retention2", spec)
    tasks, anchors = targeted_nrg.tasks_for(plan)
    assert len(anchors) == 160
    assert len({run.digest(t) for t in tasks}) == len(tasks)
    clean = [t for t in tasks if t["clean"]]
    assert {(t["numerical"]["nmax"], t["units"]) for t in clean} == {
        (n, u) for n in (2, 4) for u in ("gap", "bandwidth")}
    physical = run.load_inputs()[0]["physical"]
    for task in tasks:
        assert "[param]" in nrg.render_param(physical, task["numerical"], task["units"], task["clean"])
    assert {a["lambda"]: a["nmax"] for a in plan["anchors"]} == {4.: 26, 3.: 30, 2.5: 34, 2.: 42, 1.8: 48}
    assert all(a["nz"] == 32 and a["keep"] == 4000 for a in plan["anchors"])
    assert all(a["checks"]["keep"] == [2000, 4000, 8000] for a in plan["anchors"] if a["lambda"] <= 2)


@pytest.fixture
def fresh_case(tmp_path, monkeypatch):
    shutil.copytree(HERE / "input", tmp_path / "input")
    shutil.copyfile(HERE / "twist_analysis.py", tmp_path / "twist_analysis.py")
    for module in (run, focused_internal, targeted_nrg):
        monkeypatch.setattr(module, "HERE", tmp_path)
    monkeypatch.setattr(run, "provenance", lambda *a, **kw: {})
    monkeypatch.setattr(targeted_nrg, "implementation", lambda: {})
    return tmp_path


@pytest.mark.skipif(sys.platform != "linux", reason="original research supervisor uses Linux locks")
def test_finite_prepare_and_queue_need_no_historical_states(fresh_case, monkeypatch):
    focused_internal.main(["prepare"])
    spec = json.loads((fresh_case / "input/qp-dmrg-focused.json").read_text())
    calls = []

    def execute(directory, limits, remaining):
        request = json.loads((directory / "request.json").read_text())
        assert "initial" not in request
        assert (directory / "bath.json").is_file()
        assert remaining > 0
        calls.append(request["task"])
        result = dict(status="completed", task=request["task"], physical=request["physical"],
                      monitored_seconds=0.01, finite_converged=False, values={"signed_gap": -0.17})
        run.write_json(directory / "result.json", result)
        return result

    monkeypatch.setattr(run, "execute", execute)
    focused_internal.main(["run"])
    assert calls == [e["calculation"] for e in spec["tasks"]]
    focused_internal.main(["run"])
    assert len(calls) == 12  # Completed tasks are not rerun.
    report = json.loads((fresh_case / "runs/finite-reproduction/summary.json").read_text())
    assert report["complete"] and not report["continuum_converged"]
    assert any(c["records"][0] != c["records"][1] for c in report["comparisons"])


@pytest.mark.skipif(sys.platform != "linux", reason="original research supervisor uses Linux locks")
def test_finite_interrupted_attempt_can_retry_the_same_retained_bath(fresh_case, monkeypatch):
    attempts = []

    def execute(directory, limits, remaining):
        request = json.loads((directory / "request.json").read_text())
        attempts.append(directory)
        result = dict(task=request["task"], physical=request["physical"], monitored_seconds=.01,
                      status="failed" if len(attempts) == 1 else "completed",
                      finite_converged=False, values={"signed_gap": -.17})
        run.write_json(directory / "result.json", result)
        if len(attempts) == 1:
            raise KeyboardInterrupt
        return result

    monkeypatch.setattr(run, "execute", execute)
    with pytest.raises(KeyboardInterrupt):
        focused_internal.main(["run", "--max-tasks", "1"])
    with pytest.raises(RuntimeError, match="explicit review/retry"):
        focused_internal.main(["run", "--max-tasks", "1"])
    focused_internal.main(["run", "--max-tasks", "1", "--retry-failed"])
    assert [p.name for p in attempts] == ["attempt-001", "attempt-002"]
    assert attempts[0].parent == attempts[1].parent
    assert (attempts[0] / "bath.json").read_bytes() == (attempts[1] / "bath.json").read_bytes()
    for directory in attempts:
        request = json.loads((directory / "request.json").read_text())
        assert request["bath_sha256"] == hashlib.sha256((directory / "bath.json").read_bytes()).hexdigest()
    campaign = fresh_case / "runs/finite-reproduction"
    assert not list(campaign.glob("bath-*.json"))
    assert not (campaign / "bath-inputs").exists()
    summary = json.loads((campaign / "summary.json").read_text())
    assert len(summary["records"]) == len(summary["failures"]) == 1


@pytest.mark.skipif(sys.platform != "linux", reason="original research supervisor uses Linux locks")
def test_controls_nrg_campaign_can_supply_a_completed_reuse_record(fresh_case, monkeypatch):
    run.main(["prepare", "--profile", "controls", "--backend", "nrg", "--campaign", "controls",
              "--max-tasks", "1"])
    campaign = fresh_case / "runs/controls"
    identity = json.loads((campaign / "manifest.json").read_text())
    directory = next((campaign / "tasks").glob("nrg-*/attempt-001"))
    request = json.loads((directory / "request.json").read_text())
    assert request["task"] in identity["tasks"]
    # Synthetic external artifacts; test campaign/task validation and replay wiring,
    # not an installed NRG executable or its numerical accuracy.
    for name in ("data", "nrginit.log", "nrg.log", "raw.h5"):
        (directory / name).write_bytes(b"synthetic NRG artifact")
    dataset = targeted_nrg.nrg_dataset
    payload = dict(values={"signed_gap": -.17}, finite_converged=True, diagnostics={},
                   finite_convergence_meaning="synthetic finite flow")
    artifacts = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
                 for name in dataset.REQUIRED_ARTIFACTS}
    run.write_json(directory / "result.json", dict(request, **payload, status="completed",
                                                   monitored_seconds=.01, artifacts=artifacts))
    monkeypatch.setattr(run, "postprocess_nrg", lambda *args: payload)
    monkeypatch.setattr(run, "nrg_runtime", lambda _: dict(environment={}, nrg_sha256="test"))
    monkeypatch.setattr(dataset, "fingerprint", lambda: dict(environment={}, nrg_sha256="test"))
    targeted_nrg.main(["prepare", "--campaign", "reuse", "--source", "controls"])
    baseline = json.loads((fresh_case / "runs/reuse/baseline.json").read_text())
    assert len(baseline["records"]) == 1
    assert baseline["records"][0]["task"] == request["task"]
    assert baseline["records"][0]["values"] == payload["values"]


@pytest.mark.skipif(sys.platform != "linux", reason="original research supervisor uses Linux locks")
def test_targeted_prepare_and_single_control_need_no_source(fresh_case, monkeypatch):
    targeted_nrg.main(["prepare"])
    campaign = fresh_case / "runs/nrg-reproduction"
    identity = json.loads((campaign / "manifest.json").read_text())
    assert identity["source_campaign"] is None
    assert "baseline" not in identity
    assert not (campaign / "baseline.json").exists()

    # A newly prepared campaign is also a valid optional reuse source.
    monkeypatch.setattr(run, "nrg_runtime", lambda _: dict(environment={}, nrg_sha256="test"))
    monkeypatch.setattr(targeted_nrg.nrg_dataset, "fingerprint", lambda: dict(environment={}, nrg_sha256="test"))
    targeted_nrg.main(["prepare", "--campaign", "reuse", "--source", "nrg-reproduction"])
    reused = json.loads((fresh_case / "runs/reuse/baseline.json").read_text())
    assert reused["records"] == []

    def execute(directory, limits, remaining):
        request = json.loads((directory / "request.json").read_text())
        assert request["task"]["backend"] == "nrg"
        assert (directory / "param").read_text().startswith("[param]")
        result = dict(status="completed", task=request["task"], physical=request["physical"],
                      monitored_seconds=0.01, finite_converged=True, values={})
        run.write_json(directory / "result.json", result)
        return result

    monkeypatch.setattr(run, "execute", execute)
    targeted_nrg.main(["run", "--max-tasks", "1"])
    targeted_nrg.main(["report"])
    report = json.loads((campaign / "summary.json").read_text())
    assert report["new_completed"] == 1
    assert not report["prerequisites_passed"]
    assert all(not r["accepted"] for r in report["final_reference"]["nrg"].values())


def test_small_runner_qp_agrees_with_independent_electron_ed(tmp_path):
    physical = run.load_inputs()[0]["physical"]
    task = dict(backend="qp", cutoff=None, bath=dict(kind="cosh", levels=2))
    record = run.construct_bath(physical, task["bath"])
    run.write_json(tmp_path / "bath.json", record)
    result = run.run_internal(task, physical, tmp_path)
    expected = control.electron_ed(physical, record)
    for key, value in run.flatten(expected["branches"]).items():
        assert result["values"][key] == pytest.approx(value, abs=2e-11, rel=0)
    assert result["finite_converged"]


def test_bath_validation_rejects_changed_physics():
    physical = run.load_inputs()[0]["physical"]
    bath = json.loads((HERE / "input/final-baths.json").read_text())["cosh12"]
    changed = deepcopy(bath)
    changed["weights"][0] *= 1.01
    with pytest.raises(ValueError, match="coefficients"):
        focused_internal.validate_bath(physical, dict(kind="cosh", levels=12), changed)


def test_extrapolation_of_constant_synthetic_data_and_missing_checks():
    from NRG_comparisons.test1.analysis import OBSERVABLES
    from NRG_comparisons.test1.twist_analysis import assess

    manifest, _ = run.load_inputs()
    spec = json.loads((HERE / "input/nrg-targeted.json").read_text())
    plan = targeted_nrg.plan_for("retention2", spec)
    tasks, _ = targeted_nrg.tasks_for(plan)
    values = dict.fromkeys(OBSERVABLES, .25)
    values.update({"singlet.energy": -1., "doublet.energy": -1.25, "signed_gap": -.25,
                   "singlet.moment": 0., "singlet.P1": .5, "doublet.P1": .5})
    records = [dict(task=t, physical=manifest["physical"], status="completed", finite_converged=True,
                    values={"clean.vacuum_error": 0.} if t["clean"] else values) for t in tasks]
    report = assess(records, plan, manifest["targets"])
    assert report["prerequisites_passed"]
    for name in set(OBSERVABLES) - {"singlet.moment"}:
        fitted = report["final_reference"]["nrg"][name]
        assert fitted["accepted"]
        assert fitted["value"] == pytest.approx(values[name], abs=1e-13, rel=0)
        assert fitted["empirical_uncertainty"] < 1e-12
    report = assess([r for r in records if not r["task"]["clean"]], plan, manifest["targets"])
    assert not report["prerequisites_passed"]
    assert all(not r["accepted"] for r in report["final_reference"]["nrg"].values())


def test_presentation_reproduces_tables_without_touching_publication(tmp_path):
    pytest.importorskip("matplotlib")
    from NRG_comparisons.test1 import generate

    before = {p.name: p.read_bytes() for p in (HERE / "output").iterdir() if p.is_file()}
    generate.main(["--output-dir", str(tmp_path)])
    for name in ("results.json", "observables.csv", "comparisons.csv"):
        assert (tmp_path / name).read_bytes() == before[name]
    assert all(p.read_bytes() == before[p.name] for p in (HERE / "output").iterdir() if p.is_file())
    for name in ("comparison", "convergence"):
        assert (tmp_path / f"{name}.svg").is_file()
        assert np.std(generate.plt.imread(tmp_path / f"{name}.png")) > .05
    with pytest.raises(SystemExit):
        generate.main(["--output-dir", str(HERE / "output")])
