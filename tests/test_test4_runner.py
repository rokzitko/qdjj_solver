"""Offline runner contracts: synthetic spectra, mocked processes and signals only."""

from contextlib import nullcontext
from copy import deepcopy
import json
import math
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import Mock, call

import pytest

resource = pytest.importorskip("resource", reason="pilot supervisor is POSIX-only")

from NRG_comparisons.test4 import run as runner
from NRG_comparisons.test4.physics import PHYSICAL, paper_offset, state_values


LIMITS = dict(cpus=4, memory_gib=16, task_seconds=1200, campaign_seconds=7200,
              mathematica_instances=1)


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("unexpected process, signal, or numerical calculation in offline test")

    # Replace runner references, not stdlib objects used by pytest or other tests.
    monkeypatch.setattr(runner, "os", SimpleNamespace(
        environ={}, getpid=Mock(return_value=4100), getppid=Mock(return_value=4000),
        getpgrp=Mock(return_value=4100), getuid=Mock(return_value=1234),
        killpg=Mock(side_effect=forbidden), sched_setaffinity=Mock(),
        sched_getaffinity=Mock(return_value={2, 4, 7, 9, 11})))
    monkeypatch.setattr(runner, "subprocess", SimpleNamespace(
        Popen=Mock(side_effect=forbidden), run=Mock(side_effect=forbidden),
        STDOUT=subprocess.STDOUT, TimeoutExpired=subprocess.TimeoutExpired))
    monkeypatch.setattr(runner, "signal", SimpleNamespace(
        SIGTERM=signal.SIGTERM, SIGKILL=signal.SIGKILL, SIG_IGN=signal.SIG_IGN, signal=Mock()))
    monkeypatch.setattr(runner, "ctypes", SimpleNamespace(CDLL=Mock(side_effect=forbidden)))
    monkeypatch.setattr(runner, "time", SimpleNamespace(
        time=Mock(return_value=1000.), monotonic=Mock(return_value=100.), sleep=Mock()))
    monkeypatch.setattr(runner, "shutil", SimpleNamespace(
        which=Mock(side_effect=forbidden), copyfile=Mock(wraps=shutil.copyfile)))
    monkeypatch.setattr(runner, "lock", Mock(side_effect=lambda path: nullcontext()))
    for name in ("process_memory", "chain_bath", "extract_result"):
        monkeypatch.setattr(runner, name, Mock(side_effect=forbidden))
    limits = SimpleNamespace(RLIMIT_AS=resource.RLIMIT_AS, setrlimit=Mock())
    finite = SimpleNamespace(run_finite=Mock(side_effect=forbidden))
    monkeypatch.setitem(sys.modules, "resource", limits)
    monkeypatch.setitem(sys.modules, "NRG_comparisons.test4.finite", finite)
    monkeypatch.setitem(sys.modules, "qdjj_solver.common.io", SimpleNamespace(
        numerical_environment=Mock(return_value={"source": "offline"}),
        source_provenance=Mock(return_value={"version": "offline"})))
    return SimpleNamespace(finite=finite, resource=limits)


@pytest.fixture
def controls():
    """Use actual task declarations, not spectra or solver calls, to build receipts."""
    inventory, results = runner.tasks("controls"), {}
    for task in inventory:
        if task["backend"] == "nrg" or task["id"].startswith("hf-"):
            continue
        label = task["id"].rsplit("-", 1)[0]
        physical = PHYSICAL | task.get("physical_changes", {})
        states = []
        for spin in task["spins"]:
            original = -spin if label == "reverse" else spin
            count = min(task["roots_by_spin"].get(str(spin), task["roots"]), math.comb(10, 5+spin))
            for root in range(count):
                md = original/20 + .02
                charge, double, z2 = .8, .1, .15
                if abs(original) == 5:
                    md, charge, double, z2 = original/10, 1., 0., .25
                values = state_values(charge, double, md, z2, original/2)
                energy = .2 + .3*abs(original) + root + .03*md
                if label == "atomic":
                    energy = {(0, 0): .5, (0, 1): 1.5, (-1, 0): -.015, (1, 0): .015}.get((spin, root), energy)
                current = 0. if label == "atomic" else .2 + .01*original + .02*root
                if label == "reverse":
                    for name in ("m_dot", "m_bath"):
                        values[name] *= -1
                    values["P_up"], values["P_down"] = values["P_down"], values["P_up"]
                if label == "phase-reverse":
                    current *= -1
                values.update(energy=energy, energy_paper=energy+paper_offset(physical), current=current)
                states.append(dict(twice_sz=spin, root=root, parity=spin % 2, values=values,
                                   degenerate=False, residual=0., finite_converged=True))
        results[task["id"]] = dict(states=states, qualified=True,
                                   diagnostics=dict(all_returned_states_converged=True))
    for task in inventory:
        if task["backend"] == "nrg":
            label = "junction" if task["id"] in ("junction-bandwidth", "junction-sds") else task["id"].rsplit("-", 1)[0]
            results[task["id"]] = deepcopy(results[label+"-ed"])
        elif task["id"].startswith("hf-"):
            result = deepcopy(results[task["bath_source"].replace("-nrg", "-ed")])
            result["states"] = [r for r in result["states"] if r["twice_sz"] in task["spins"] and r["root"] < task["roots"]]
            # Artificially linear E(phi) gives an exact, nonzero branch current.
            for row in result["states"]:
                for name in ("energy", "energy_paper"):
                    row["values"][name] += (task["physical_changes"]["phi"]-PHYSICAL["phi"])*row["values"]["current"]
            results[task["id"]] = result
    return results


def state(results, task, spin, root=0):
    return next(r for r in results[task]["states"] if (r["twice_sz"], r["root"]) == (spin, root))


def checks(results):
    gate = runner.check_controls(results)
    return gate["passed"], {row["name"]: row for row in gate["checks"]}


def test_frozen_inventory_limits_and_pilot_zero_bridges():
    expected = []
    for label, changes, backends in (
        ("atomic", {"gamma": 0.}, ("ed", "qp")), ("junction", {}, ("ed", "qp", "dmrg")),
        ("normal", {"wn": 20.}, ("ed",)), ("spin", {"ws": 10.}, ("ed",)),
        ("bridge", {"wn": 20., "ws": 10.}, ("ed", "qp")),
        ("direct", {"gamma": 0., "wn": 20., "ws": 10.}, ("ed",)),
        ("free", {"u": 0.}, ("ed", "qp")), ("reverse", {"wn": 20., "ws": -10., "field": -.03}, ("ed",)),
        ("phase-reverse", {"wn": 20., "ws": -10., "phi": -3*math.pi/5}, ("ed",)),
    ):
        numerical = dict(nmax=1, keep=1024, keepmin=0, keepenergy=-1., untruncated=True, z=1., **{"lambda": 4.})
        expected.append(dict(id=label+"-nrg", backend="nrg", physical_changes=changes, numerical=numerical))
        expected.extend(dict(id=label+"-"+backend, backend=backend, physical_changes=changes,
                             bath_source=label+"-nrg", chi=32, cutoff=None, roots=2,
                             roots_by_spin={"0": 3}, spins=list(range(-5, 6))) for backend in backends)
    expected += [dict(id="junction-bandwidth", backend="nrg", units="bandwidth", numerical=numerical),
                 dict(id="junction-sds", backend="nrg", representation="sds", numerical=numerical)]
    for label, changes in (("junction", {}), ("bridge", {"wn": 20., "ws": 10.}),
                           ("direct", {"gamma": 0., "wn": 20., "ws": 10.})):
        expected.extend(dict(id=f"hf-{label}-{step:g}-{sign:+d}", backend="ed", bath_source=label+"-nrg",
                             roots=2, spins=[0, -1, 1], physical_changes=changes | {"phi": 3*math.pi/5+sign*step})
                        for step in (2e-4, 1e-4) for sign in (-1, 1))
    assert runner.tasks("controls") == expected
    pilot = [dict(id=f"nrg-k{k}", backend="nrg", numerical=dict(
        nmax=20, keep=k, keepmin=k//5, keepenergy=14., z=1., **{"lambda": 4.})) for k in (200, 500, 1000)]
    pilot += [dict(id="nrg-length22", backend="nrg", numerical=dict(
        nmax=22, keep=500, keepmin=100, keepenergy=14., z=1., **{"lambda": 4.})),
        dict(id="cosh2-qp-full", backend="qp", bath="cosh2", cutoff=None, roots=3)]
    pilot += [dict(id=f"cosh4-qp-q{q if q is not None else 'full'}", backend="qp", bath="cosh4", cutoff=q,
                   roots=2, roots_by_spin={"0": 3, "-2": 1, "2": 1}) for q in (2, 4, 6, None)]
    pilot += [dict(id=f"cosh4-dmrg-chi{chi}", backend="dmrg", bath="cosh4", chi=chi,
                   roots=2, roots_by_spin={"0": 3, "-2": 1, "2": 1}) for chi in (32, 64, 128)]
    assert runner.tasks("pilot") == pilot
    inputs = json.loads((runner.HERE / "input/physical.json").read_text())
    assert inputs["limits"] == runner.LIMITS == LIMITS and inputs["minimum_lambda"] == 4
    assert inputs["physical"] == PHYSICAL
    for task in pilot:
        physical = inputs["physical"] | task.get("physical_changes", {})
        assert physical["wn"] == physical["ws"] == 0
    expected[0]["numerical"]["lambda"] = 2
    assert runner.tasks("controls")[0]["numerical"]["lambda"] == 4


@pytest.mark.parametrize("stage", ["controls", "pilot"])
def test_lambda_floor_is_hard_even_for_untruncated_controls(stage):
    for task in runner.tasks(stage):
        if task["backend"] != "nrg":
            continue
        physical = PHYSICAL | task.get("physical_changes", {})
        args = task.get("units", "gap"), task.get("representation", "paper")
        assert "Lambda=4\n" in runner.render_param(physical, task["numerical"], *args)
        with pytest.raises(ValueError, match="Lambda must be >= 4"):
            runner.render_param(physical, task["numerical"] | {"lambda": 3.999}, *args)


def test_complete_controls_and_all_comparisons(controls):
    passed, rows = checks(controls)
    assert passed
    assert set(controls) == {task["id"] for task in runner.tasks("controls")}
    for task in runner.tasks("controls"):
        if task["backend"] == "nrg" or task["id"].startswith("hf-"):
            continue
        name = task["id"]
        new, old = (task["bath_source"], name) if task["backend"] == "ed" else (name, name.rsplit("-", 1)[0]+"-ed")
        assert rows[f"{name} declared inventory"]["missing"] == []
        assert rows[f"{name} all-root residual/orthogonality"]["passed"]
        assert rows[f"{new} covers {old}"]["missing"] == []
        assert rows[f"{new} vs {old}"]["maximum_error"] == 0
    for name in ("junction-bandwidth", "junction-sds"):
        assert rows[f"{name} vs junction-nrg"]["maximum_error"] == 0
    assert len([r for name, r in rows.items() if "Hellmann-Feynman" in name and r["maximum_error"] < 1e-10]) == 6
    assert rows["bridge carries current with Gamma=0"]["measured"] > 0


@pytest.mark.parametrize("task,spin,root", [
    ("junction-ed", -5, 0), ("junction-ed", 5, 0), ("junction-ed", 0, 2),
    ("normal-ed", 2, 1), ("junction-qp", 0, 2), ("junction-dmrg", -5, 0),
    ("junction-nrg", 5, 0), ("junction-sds", 0, 2), ("junction-bandwidth", 0, 2),
    ("hf-direct-0.0001-+1", 0, 1),
])
def test_missing_declared_roots_must_not_shrink_gate(controls, task, spin, root):
    controls[task]["states"].remove(state(controls, task, spin, root))
    assert not runner.check_controls(controls)["passed"], f"accepted missing {task} root {(spin, root)}"


@pytest.mark.parametrize("task", [t["id"] for t in runner.tasks("controls") if t["backend"] != "nrg"])
def test_every_finite_controls_all_root_residual_flag_is_required(controls, task):
    state(controls, task, 0, 1).update(residual=.1, finite_converged=False)
    controls[task]["diagnostics"]["all_returned_states_converged"] = False
    assert controls[task]["qualified"]  # The primary-state flag alone is insufficient.
    assert not runner.check_controls(controls)["passed"], f"accepted unconverged {task} excitation"


@pytest.mark.parametrize("endpoint", ["junction-ed", "junction-qp"])
@pytest.mark.parametrize("energy", ["energy", "energy_paper"])
def test_same_sector_degeneracy_compares_both_energies_only(controls, endpoint, energy):
    for name in ("junction-ed", "junction-qp", "junction-dmrg", "junction-nrg", "junction-sds", "junction-bandwidth"):
        for key in ("energy", "energy_paper"):
            state(controls, name, 2, 1)["values"][key] = state(controls, name, 2)["values"][key]
    for root in (0, 1):
        state(controls, endpoint, 2, root)["degenerate"] = True
    values = state(controls, "junction-qp", 2, 1)["values"]
    for name in runner.OBSERVABLES:
        if name not in ("energy", "energy_paper"):
            values[name] += .1
    assert runner.check_controls(controls)["passed"]
    values[energy] += 1e-5
    passed, rows = checks(controls)
    assert not passed and not rows["junction-qp vs junction-ed"]["passed"]


@pytest.mark.parametrize("observable", runner.OBSERVABLES)
def test_nondegenerate_comparisons_do_not_ignore_any_observable(controls, observable):
    state(controls, "junction-qp", 0, 2)["values"][observable] += 1e-5
    assert not runner.check_controls(controls)["passed"]


@pytest.mark.parametrize("label", ["junction", "bridge", "direct"])
@pytest.mark.parametrize("step", [2e-4, 1e-4])
def test_hf_checks_each_phase_step_and_even_branch_current(controls, label, step):
    current = state(controls, label+"-nrg", 0)["values"]["current"]
    assert current != 0
    state(controls, f"hf-{label}-{step:g}-+1", 0)["values"]["energy"] -= 2*step*current
    passed, rows = checks(controls)
    assert not passed and not rows[f"{label} current Hellmann-Feynman h={step:g}"]["passed"]
    other = 1e-4 if step == 2e-4 else 2e-4
    assert rows[f"{label} current Hellmann-Feynman h={other:g}"]["passed"]


def test_gamma_zero_direct_current_gate_is_independent_of_backend_agreement(controls):
    for name, result in controls.items():
        if name.startswith(("direct-", "hf-direct-")):
            for row in result["states"]:
                if name.startswith("hf-"):
                    reference = state(controls, "direct-ed", row["twice_sz"], row["root"])
                    for key in ("energy", "energy_paper"):
                        row["values"][key] = reference["values"][key]
                row["values"]["current"] = 0.
    passed, rows = checks(controls)
    assert not passed
    assert [name for name, row in rows.items() if not row["passed"]] == ["bridge carries current with Gamma=0"]


@pytest.mark.parametrize("label,observable", [("atomic", "energy_paper"), ("atomic", "current"),
                                               ("reverse", "m_bath"), ("phase-reverse", "current")])
def test_absolute_and_reversal_controls_survive_common_backend_errors(controls, label, observable):
    for name in controls:
        if name.startswith(label+"-"):
            state(controls, name, 1)["values"][observable] += .01
    passed, rows = checks(controls)
    assert not passed and rows[f"{label}-nrg vs {label}-ed"]["passed"]
    name = "atomic energies, paper constant, Zeeman sign and zero current" if label == "atomic" else label+" symmetry"
    assert not rows[name]["passed"]


@pytest.fixture
def campaign(tmp_path, monkeypatch, controls):
    root = tmp_path / "repo"
    here = root / "NRG_comparisons/test4"
    (here / "input").mkdir(parents=True)
    for relative in [*(f"NRG_comparisons/test4/{name}" for name in
                       ("run.py", "finite.py", "nrg.py", "physics.py", *runner.MODELS)),
                     "NRG_comparisons/test1/run.py", "tools/nrginit", "tools/nrg", "tools/math", "tools/model.m"]:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"offline identity fixture: {relative}\n")
    monkeypatch.setattr(runner, "ROOT", root)
    monkeypatch.setattr(runner, "HERE", here)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "home"))
    inputs = dict(format_version=1, physical=deepcopy(PHYSICAL), limits=deepcopy(LIMITS), minimum_lambda=4)
    runner.write_json(here / "input/physical.json", inputs)
    bath = dict(xi=[.25], weights=[.03], delta=1., bandwidth=100.)
    runner.write_json(here / "input/baths.json", {"cosh2": bath, "cosh4": bath})
    runner.shutil.which.side_effect = lambda name: str(root / "tools" / name)
    runner.subprocess.run.side_effect = None
    runner.subprocess.run.return_value = SimpleNamespace(returncode=1)
    runner.chain_bath.side_effect = None
    runner.chain_bath.return_value = bath

    def execute(work, remaining):
        request = json.loads((work / "request.json").read_text())
        result = deepcopy(controls.get(request["task"]["id"], {}))
        result.update(status="completed", request_sha256=runner.digest(request), monitored_seconds=.25, peak_rss_gib=.01)
        if request["task"]["backend"] == "nrg":
            result["artifacts"] = {name: runner.sha256(work / name) for name in (*runner.MODELS, "param")}
        runner.write_json(work / "result.json", result)
        return result

    execute_mock = Mock(side_effect=execute)
    monkeypatch.setattr(runner, "execute", execute_mock)
    return SimpleNamespace(directory=here / "runs/offline", root=root, inputs=inputs, execute=execute_mock)


def test_campaign_snapshot_replay_and_sequential_pilot_gate(campaign, capsys):
    directory = campaign.directory
    shared = campaign.root / "tools/model.m"
    original = shared.read_bytes()
    results = runner.campaign("controls", directory)
    assert list(results) == [t["id"] for t in runner.tasks("controls")]
    assert campaign.execute.call_count == 37
    manifest = json.loads((directory / "campaign.json").read_text())
    assert manifest["identity"]["tasks"] == {s: runner.tasks(s) for s in ("controls", "pilot")}
    assert len(manifest["identity"]["source_hashes"]) == 7
    assert manifest["identity"]["models"] == {name: runner.sha256(runner.HERE / name) for name in runner.MODELS}
    assert manifest["identity"]["external_sources"] == {str(shared): runner.sha256(shared)}
    for task in runner.tasks("controls"):
        if task["backend"] == "nrg":
            request = json.loads((directory / task["id"] / "request.json").read_text())
            assert request["models"] == manifest["identity"]["models"]
            for name in runner.MODELS:
                snapshot = directory / task["id"] / name
                assert snapshot.read_bytes() == (runner.HERE / name).read_bytes()
                assert not snapshot.samefile(runner.HERE / name)
    before = {p: p.read_bytes() for p in directory.rglob("*") if p.is_file()}
    campaign.execute.reset_mock()
    assert runner.campaign("controls", directory) == results
    campaign.execute.assert_not_called()
    assert all(p.read_bytes() == content for p, content in before.items())
    results = runner.campaign("pilot", directory)
    assert list(results) == [t["id"] for s in ("controls", "pilot") for t in runner.tasks(s)]
    assert [c.args[0].name for c in campaign.execute.call_args_list] == [t["id"] for t in runner.tasks("pilot")]
    assert all(c.args[1] == 7200 for c in campaign.execute.call_args_list)
    assert json.loads((directory / "cosh4-qp-q2/request.json").read_text())["bath"]["weights"] == [.03]
    assert shared.read_bytes() == original
    assert all(Path(c.args[1]).is_relative_to(directory) for c in runner.shutil.copyfile.call_args_list)
    assert runner.lock.call_args_list == [call(directory.parent / ".session.lock")]*3
    assert "Production not launched." in capsys.readouterr().out
    runner.subprocess.Popen.assert_not_called()


@pytest.mark.parametrize("setting", [*LIMITS, "minimum_lambda", "outside"])
def test_unreviewed_limits_and_outside_paths_rejected_before_tools(campaign, setting):
    directory = campaign.directory
    if setting == "outside":
        directory = runner.HERE / "runs/../outside"
    else:
        target = campaign.inputs if setting == "minimum_lambda" else campaign.inputs["limits"]
        target[setting] += 1
        runner.write_json(runner.HERE / "input/physical.json", campaign.inputs)
    with pytest.raises(ValueError):
        runner.campaign("pilot", directory)
    runner.shutil.which.assert_not_called()
    campaign.execute.assert_not_called()


@pytest.mark.parametrize("problem", ["missing-tool", "active-process", "pgrep-error"])
def test_environment_failure_does_not_stop_external_work(campaign, problem):
    if problem == "missing-tool":
        runner.shutil.which.side_effect = None
        runner.shutil.which.return_value = None
    else:
        runner.subprocess.run.return_value.returncode = 0 if problem == "active-process" else 2
    with pytest.raises(RuntimeError, match="load NRG/Mathematica|nothing stopped"):
        runner.campaign("controls", campaign.directory)
    campaign.execute.assert_not_called()
    runner.os.killpg.assert_not_called()


@pytest.mark.parametrize("relative", ["NRG_comparisons/test4/run.py", "NRG_comparisons/test4/model.m",
                                       "NRG_comparisons/test4/modeloperators.m", "NRG_comparisons/test4/input/baths.json",
                                       "NRG_comparisons/test1/run.py", "tools/nrg", "tools/model.m"])
def test_identity_changes_invalidate_replay_and_preserve_old_files(campaign, relative):
    runner.campaign("controls", campaign.directory)
    before = {p: p.read_bytes() for p in campaign.directory.rglob("*") if p.is_file()}
    path = campaign.root / relative
    path.write_text(path.read_text()+"\n")
    campaign.execute.reset_mock()
    with pytest.raises(ValueError, match="code/inputs/tools changed"):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert all(p.read_bytes() == content for p, content in before.items())


@pytest.mark.parametrize("corruption", ["request", "digest", "failed", "missing", "param", "model.m", "modeloperators.m"])
def test_replay_never_retries_changed_or_incomplete_work(campaign, corruption):
    runner.campaign("controls", campaign.directory)
    work = campaign.directory / "atomic-nrg"
    if corruption in ("request", "digest", "failed"):
        path = work / ("request.json" if corruption == "request" else "result.json")
        value = json.loads(path.read_text())
        value[{"request": "physical", "digest": "request_sha256", "failed": "status"}[corruption]] = "changed"
        runner.write_json(path, value)
    elif corruption == "missing":
        (work / "result.json").unlink()
    else:
        (work / corruption).write_text("changed snapshot\n")
    before = {p: p.read_bytes() for p in work.iterdir()}
    campaign.execute.reset_mock()
    with pytest.raises((ValueError, FileNotFoundError)):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert {p: p.read_bytes() for p in work.iterdir()} == before
    assert not (campaign.directory / "nrg-k200").exists()


def test_failed_controls_block_every_pilot_and_worker_failure_stops_immediately(campaign):
    runner.campaign("controls", campaign.directory)
    path = campaign.directory / "junction-dmrg/result.json"
    result = json.loads(path.read_text())
    result["diagnostics"]["all_returned_states_converged"] = False
    runner.write_json(path, result)
    campaign.execute.reset_mock()
    with pytest.raises(RuntimeError, match="model controls failed; cost pilots blocked"):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert not json.loads((campaign.directory / "controls.json").read_text())["passed"]
    assert not (campaign.directory / "nrg-k200").exists()
    campaign.execute.side_effect = None
    campaign.execute.return_value = dict(status="failed", error="synthetic", monitored_seconds=.25, peak_rss_gib=.01)
    fresh = campaign.directory.parent / "fresh"
    with pytest.raises(RuntimeError, match="atomic-nrg failed: synthetic"):
        runner.campaign("pilot", fresh)
    campaign.execute.assert_called_once()
    assert not (fresh / "atomic-ed").exists()


def test_session_budget_is_shared_by_test4_siblings_but_independent_of_test3(campaign):
    old = campaign.root / "NRG_comparisons/test3/runs/old/campaign.json"
    old.parent.mkdir(parents=True)
    runner.write_json(old, {"started_unix": -100000})
    runner.campaign("controls", campaign.directory)
    assert json.loads((campaign.directory / "campaign.json").read_text())["started_unix"] == 1000
    runner.time.time.return_value = 1010
    campaign.execute.reset_mock()
    fresh = campaign.directory.parent / "sibling"
    runner.campaign("controls", fresh)
    assert json.loads((fresh / "campaign.json").read_text())["started_unix"] == 1000
    assert all(c.args[1] == 7190 for c in campaign.execute.call_args_list)
    assert runner.lock.call_args_list == [call(campaign.directory.parent / ".session.lock")]*2
    assert json.loads(old.read_text()) == {"started_unix": -100000}


@pytest.mark.parametrize("cause", ["timeout", "memory"])
def test_execute_bounds_and_owned_group_cleanup(tmp_path, offline, cause):
    child = Mock(pid=9001, returncode=-signal.SIGKILL)
    child.poll.return_value = None
    child.wait.side_effect = [subprocess.TimeoutExpired("mock worker", 3), -signal.SIGKILL, -signal.SIGKILL]
    runner.subprocess.Popen.side_effect = lambda *a, **kw: (kw["preexec_fn"](), child)[1]
    runner.process_memory.side_effect = None
    runner.process_memory.return_value = 17*1024**3 if cause == "memory" else 1024**3
    runner.time.monotonic.side_effect = [100., 105., 107.] if cause == "timeout" else [100., 101.]
    runner.os.killpg.side_effect = None
    events = Mock()
    events.attach_mock(runner.os.killpg, "killpg")
    events.attach_mock(child.wait, "wait")
    result = runner.execute(tmp_path, remaining=5)
    assert result["status"] == "failed" and result["error"] == "memory/time resource limit reached"
    assert result["peak_rss_gib"] == (17 if cause == "memory" else 1)
    assert result["monitored_seconds"] == (7 if cause == "timeout" else 1)
    assert json.loads((tmp_path / "result.json").read_text()) == result
    assert json.loads((tmp_path / "launch.json").read_text()) == dict(supervisor_pid=4100, timeout_seconds=5)
    args, kwargs = runner.subprocess.Popen.call_args
    assert args[0] == [sys.executable, "-B", str(runner.HERE / "run.py"), "worker", str(tmp_path)]
    assert kwargs["start_new_session"] is True and kwargs["cwd"] == runner.ROOT
    assert all(kwargs["env"][k] == "1" for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE"))
    offline.resource.setrlimit.assert_called_once_with(resource.RLIMIT_AS, (16*1024**3, 16*1024**3))
    runner.os.sched_setaffinity.assert_called_once_with(0, [2, 4, 7, 9])
    runner.process_memory.assert_called_once_with(9001)
    assert events.mock_calls == [call.killpg(9001, signal.SIGTERM), call.wait(timeout=3),
                                 call.killpg(9001, signal.SIGKILL), call.wait(timeout=3), call.wait()]


@pytest.mark.parametrize("exitcode,receipt", [(0, {"status": "completed"}), (7, {"status": "failed", "error": "specific"}), (0, None)])
def test_execute_finished_worker_and_absent_group(tmp_path, exitcode, receipt):
    if receipt is not None:
        runner.write_json(tmp_path / "result.json", receipt)
    child = Mock(pid=9001, returncode=exitcode)
    child.poll.return_value = exitcode
    runner.subprocess.Popen.side_effect = None
    runner.subprocess.Popen.return_value = child
    runner.os.killpg.side_effect = ProcessLookupError
    result = runner.execute(tmp_path, remaining=7200)
    assert result["status"] == ("completed" if receipt and exitcode == 0 else "failed")
    if exitcode:
        assert result["error"] == "specific"
    assert result["peak_rss_gib"] == 0
    assert json.loads((tmp_path / "launch.json").read_text())["timeout_seconds"] == 1200
    assert runner.os.killpg.call_args_list == [call(9001, signal.SIGTERM), call(9001, signal.SIGKILL)]
    runner.process_memory.assert_not_called()


@pytest.mark.parametrize("remaining", [0, -1])
def test_execute_exhausted_budget_never_launches_or_signals(tmp_path, remaining):
    with pytest.raises(RuntimeError, match="budget exhausted"):
        runner.execute(tmp_path, remaining)
    runner.subprocess.Popen.assert_not_called()
    runner.os.killpg.assert_not_called()
    assert list(tmp_path.iterdir()) == []


@pytest.fixture
def worker(tmp_path):
    task = runner.tasks("controls")[0]
    request = dict(task=task, physical=PHYSICAL | task["physical_changes"],
                   tools={name: f"/mock/{name}" for name in ("nrginit", "math", "nrg")})
    for name in runner.MODELS:
        (tmp_path / name).write_text(f"offline {name}\n")
    request["models"] = {name: runner.sha256(tmp_path / name) for name in runner.MODELS}
    (tmp_path / "param").write_text(runner.render_param(request["physical"], task["numerical"]))
    runner.write_json(tmp_path / "request.json", request)
    runner.write_json(tmp_path / "launch.json", {"supervisor_pid": 4000})
    libc = Mock()
    libc.prctl.return_value = 0
    runner.ctypes.CDLL.side_effect = None
    runner.ctypes.CDLL.return_value = libc
    return tmp_path, request, libc


@pytest.mark.parametrize("changed", ["param", "model.m", "modeloperators.m"])
def test_worker_rejects_changed_prepared_inputs_before_external_work(worker, changed):
    directory, _, _ = worker
    (directory / changed).write_text("changed\n")
    with pytest.raises(ValueError, match="prepared model or NRG input changed"):
        runner.worker(directory)
    assert json.loads((directory / "result.json").read_text())["status"] == "failed"
    runner.subprocess.run.assert_not_called()
    assert not (directory / "scratch").exists()


@pytest.mark.parametrize("returncode", [0, 2])
def test_worker_refuses_active_mathematica_or_failed_probe(worker, returncode):
    directory, _, _ = worker
    runner.subprocess.run.side_effect = None
    runner.subprocess.run.return_value = SimpleNamespace(returncode=returncode)
    with pytest.raises(RuntimeError, match="Mathematica active or process check failed"):
        runner.worker(directory)
    assert runner.subprocess.run.call_count == 1
    assert runner.subprocess.run.call_args.args[0][0] == "pgrep"
    runner.os.killpg.assert_not_called()


@pytest.mark.parametrize("cause", ["wrong-group", "prctl-failed", "parent-changed"])
def test_worker_owns_group_and_parentdeath_cleanup_is_local(worker, offline, cause):
    directory, _, libc = worker
    if cause == "wrong-group":
        runner.os.getpgrp.return_value = 4000
        with pytest.raises(RuntimeError, match="must own its process group"):
            runner.worker(directory)
        libc.prctl.assert_not_called()
        runner.os.killpg.assert_not_called()
    else:
        if cause == "prctl-failed":
            libc.prctl.return_value = -1
        else:
            runner.os.getppid.return_value = 1
        runner.os.killpg.side_effect = [None, SystemExit("mock SIGKILL")]
        with pytest.raises(SystemExit, match="mock SIGKILL"):
            runner.worker(directory)
        libc.prctl.assert_called_once_with(1, signal.SIGTERM, 0, 0, 0)
        assert runner.os.killpg.call_args_list == [call(4100, signal.SIGTERM), call(4100, signal.SIGKILL)]
        runner.time.sleep.assert_called_once_with(2)
        runner.signal.signal.assert_any_call(signal.SIGTERM, signal.SIG_IGN)
    offline.finite.run_finite.assert_not_called()
    runner.subprocess.run.assert_not_called()


def test_production_cli_rejected_and_manual_export_is_not_execution(tmp_path, monkeypatch):
    from NRG_comparisons.test4 import manual

    campaign = Mock()
    monkeypatch.setattr(runner, "campaign", campaign)
    monkeypatch.setattr(sys, "argv", ["run.py", "production", str(tmp_path)])
    with pytest.raises(ValueError, match="production is not authorized"):
        runner.tasks("production")
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
    campaign.assert_not_called()
    monkeypatch.setattr(manual, "subprocess", runner.subprocess)
    prepare, execute = Mock(), Mock()
    monkeypatch.setattr(manual, "prepare", prepare)
    monkeypatch.setattr(manual, "execution_tools", execute)
    destination = tmp_path / "manual"
    manual.main(["export", str(destination)])
    assert (destination / "bundle.json").is_file()
    assert not (destination / "sweep.json").exists()
    assert "#AUTOLOOP: true" in (destination / "param.loop").read_text()
    prepare.assert_not_called()
    manual.main(["run"])
    prepare.assert_called_once_with(manual.HERE, execute=False)
    execute.assert_not_called()
    runner.subprocess.run.assert_not_called()
    runner.subprocess.Popen.assert_not_called()
