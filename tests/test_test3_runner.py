"""Offline runner contracts: synthetic receipts, no solvers, processes, or signals."""

from contextlib import nullcontext
from copy import deepcopy
import json
from pathlib import Path
import signal
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import Mock, call

import pytest

pytest.importorskip("resource", reason="pilot supervisor is POSIX-only")

from NRG_comparisons.test3 import run as runner
from NRG_comparisons.test3.manual_run.physics import PHYSICAL, state_values


LIMITS = dict(cpus=4, memory_gib=16, task_seconds=1200, campaign_seconds=7200,
              mathematica_instances=1)


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    """Replace module references, not the stdlib objects used by pytest itself."""
    def forbidden(*args, **kwargs):
        pytest.fail("unexpected external or numerical operation in offline runner test")

    monkeypatch.setattr(runner, "os", SimpleNamespace(
        environ={}, getpid=Mock(return_value=4100), getppid=Mock(return_value=4000),
        getpgrp=Mock(return_value=4100), getuid=Mock(return_value=1234),
        killpg=Mock(side_effect=forbidden), sched_setaffinity=Mock(),
        sched_getaffinity=Mock(return_value={2, 4, 7, 9, 11, 12})))
    monkeypatch.setattr(runner, "subprocess", SimpleNamespace(
        Popen=Mock(side_effect=forbidden), run=Mock(side_effect=forbidden),
        STDOUT=subprocess.STDOUT, TimeoutExpired=subprocess.TimeoutExpired))
    monkeypatch.setattr(runner, "signal", SimpleNamespace(
        SIGTERM=signal.SIGTERM, SIGKILL=signal.SIGKILL, SIG_IGN=signal.SIG_IGN,
        signal=Mock()))
    monkeypatch.setattr(runner, "resource", SimpleNamespace(
        RLIMIT_AS=runner.resource.RLIMIT_AS, setrlimit=Mock()))
    monkeypatch.setattr(runner, "ctypes", SimpleNamespace(CDLL=Mock(side_effect=forbidden)))
    monkeypatch.setattr(runner, "time", SimpleNamespace(
        time=Mock(return_value=1000.), monotonic=Mock(return_value=100.), sleep=Mock()))
    monkeypatch.setattr(runner, "shutil", SimpleNamespace(which=Mock(side_effect=forbidden)))
    monkeypatch.setattr(runner, "lock", Mock(side_effect=lambda path: nullcontext()))
    for name in ("process_memory", "chain_bath", "extract_result"):
        monkeypatch.setattr(runner, name, Mock(side_effect=forbidden))
    finite = SimpleNamespace(run_finite=Mock(side_effect=forbidden))
    monkeypatch.setitem(sys.modules, "NRG_comparisons.test3.finite", finite)
    monkeypatch.setitem(sys.modules, "qdjj_solver.common.io", SimpleNamespace(
        numerical_environment=Mock(return_value={"source": "offline fixture"}),
        source_provenance=Mock(return_value={"version": "offline fixture"})))
    return finite


@pytest.fixture
def controls():
    results = {}
    for label, maximum in (("atomic", 3), ("chain", 4), ("zero", 4)):
        states = []
        for spin in range(-maximum, maximum+1):
            for root in range(1 if abs(spin) == maximum else 2):
                md = spin/20 + (0 if label == "zero" else .02)
                charge, double, z2 = .8, .1, .15
                if abs(spin) == maximum:
                    md, charge, double, z2 = (.5 if spin > 0 else -.5), 1., 0., .25
                energy = .2 + .3*abs(spin) + root + .03*md
                if label == "atomic" and (spin, root) in {(0, 0), (0, 1), (-1, 0), (1, 0)}:
                    energy = {(0, 0): .5, (0, 1): 1.5, (-1, 0): -.015, (1, 0): .015}[spin, root]
                    md, charge, double, z2 = spin/2, float(abs(spin)), 0., abs(spin)/4
                    if (spin, root) == (0, 1):
                        charge, double = 2., 1.
                values = dict(state_values(charge, double, md, z2, spin/2), energy=energy)
                states.append(dict(twice_sz=spin, root=root, parity=spin % 2, values=values,
                                   degenerate=False, residual=0., finite_converged=True))
        for backend in ("nrg", "ed", "qp", "dmrg") if label != "zero" else ("nrg", "ed"):
            results[f"{label}-{backend}"] = dict(
                states=deepcopy(states), qualified=True,
                diagnostics=dict(all_returned_states_converged=True))
    reverse = deepcopy(results["chain-nrg"])
    for row in reverse["states"]:
        row["twice_sz"] *= -1
        values = row["values"]
        for name in ("m_dot", "m_bath"):
            values[name] *= -1
        values["P_up"], values["P_down"] = values["P_down"], values["P_up"]
    for backend in ("nrg", "ed"):
        results[f"reverse-{backend}"] = deepcopy(reverse)
    results["chain-bandwidth"] = deepcopy(results["chain-nrg"])
    for step in (2e-4, 1e-4):
        for sign in (-1, 1):
            result = deepcopy(results["chain-ed"])
            result["states"] = [row for row in result["states"] if row["twice_sz"] in runner.PRIMARY]
            for row in result["states"]:
                row["values"]["energy"] += sign*step*row["values"]["m_dot"]
            results[f"hf-{step:g}-{sign:+d}"] = result
    return results


def state(results, task, spin, root=0):
    return next(row for row in results[task]["states"] if (row["twice_sz"], row["root"]) == (spin, root))


def checks(results):
    gate = runner.check_controls(results)
    return gate, {row["name"]: row for row in gate["checks"]}


def test_frozen_inventory_and_no_production(monkeypatch):
    controls = runner.tasks("controls")
    assert [task["id"] for task in controls] == [
        "atomic-nrg", "atomic-ed", "atomic-qp", "atomic-dmrg",
        "chain-nrg", "chain-ed", "chain-qp", "chain-dmrg", "zero-nrg", "zero-ed",
        "reverse-nrg", "reverse-ed", "chain-bandwidth",
        "hf-0.0002--1", "hf-0.0002-+1", "hf-0.0001--1", "hf-0.0001-+1"]
    for task in controls:
        if task["backend"] == "nrg":
            assert task["numerical"] == dict(
                nmax=1 if task["id"] == "atomic-nrg" else 2, keep=512, keepmin=0,
                keepenergy=-1., untruncated=True, z=1., **{"lambda": 2.})
        elif not task["id"].startswith("hf-"):
            maximum = 3 if task["id"].startswith("atomic-") else 4
            assert task["spins"] == list(range(-maximum, maximum+1))
            assert task["roots"] == 2 and task["cutoff"] is None and task["chi"] == 32
    by_id = {task["id"]: task for task in controls}
    assert by_id["atomic-nrg"]["physical_changes"] == {"gamma": 0.}
    assert by_id["zero-nrg"]["physical_changes"] == {"field": 0.}
    assert by_id["reverse-nrg"]["physical_changes"] == {"field": -.03}
    assert by_id["chain-bandwidth"]["units"] == "bandwidth"
    for step in (2e-4, 1e-4):
        for sign in (-1, 1):
            assert by_id[f"hf-{step:g}-{sign:+d}"] == dict(
                id=f"hf-{step:g}-{sign:+d}", backend="ed", bath_source="chain-nrg",
                physical_changes={"field": .03+sign*step}, roots=2, spins=[0, -1, 1])
    expected = [dict(id=f"cosh8-qp-q{q if q is not None else 'full'}", backend="qp",
                     bath="cosh8", cutoff=q, roots=2) for q in (2, 4, 6, None)]
    expected += [dict(id=f"cosh8-dmrg-chi{chi}", backend="dmrg", bath="cosh8", chi=chi, roots=2)
                 for chi in (32, 64, 128)]
    expected += [dict(id=f"nrg-k{k}", backend="nrg", numerical=dict(
        nmax=42, keep=k, keepmin=k//5, keepenergy=14., z=1., **{"lambda": 2.}))
        for k in (500, 1000, 2000)]
    assert runner.tasks("pilot") == expected
    expected[0]["cutoff"] = 100
    assert runner.tasks("pilot")[0]["cutoff"] == 2
    with pytest.raises(ValueError, match="controls or pilot"):
        runner.tasks("production")
    campaign = Mock()
    monkeypatch.setattr(runner, "campaign", campaign)
    monkeypatch.setattr(sys, "argv", ["run.py", "production", "unused"])
    with pytest.raises(SystemExit) as error:
        runner.main()
    assert error.value.code == 2
    campaign.assert_not_called()


def test_complete_synthetic_controls_cover_every_comparison(controls):
    gate, rows = checks(controls)
    assert gate["passed"]
    pairs = [(f"{label}-{backend}", f"{label}-ed") for label in ("atomic", "chain")
             for backend in ("nrg", "qp", "dmrg")]
    pairs += [("zero-nrg", "zero-ed"), ("reverse-nrg", "reverse-ed"), ("chain-bandwidth", "chain-nrg")]
    for new, old in pairs:
        assert rows[f"{new} covers all {old} control roots"]["missing"] == []
        assert rows[f"{new} vs {old}"]["maximum_error"] == 0
    for label in ("atomic", "chain"):
        for backend in ("ed", "qp", "dmrg"):
            assert rows[f"{label}-{backend} primary residual/uniqueness qualification"]["passed"]
            assert rows[f"{label}-{backend} all control-root residuals and orthogonality"]["passed"]
    assert rows["atomic energies and Zeeman sign"]["passed"] and rows["NRG spin reversal"]["passed"]
    for step in (2e-4, 1e-4):
        assert rows[f"Hellmann-Feynman h={step:g}"]["maximum_error"] < 1e-12


@pytest.mark.parametrize("task,spin,root,reference", [
    ("atomic-nrg", -3, 0, "atomic-ed"), ("chain-qp", 4, 0, "chain-ed"),
    ("chain-dmrg", -2, 0, "chain-ed"), ("chain-dmrg", 1, 1, "chain-ed"),
    ("chain-bandwidth", 0, 1, "chain-nrg"),
])
def test_missing_polarized_exclusion_and_excited_roots_fail(controls, task, spin, root, reference):
    controls[task]["states"].remove(state(controls, task, spin, root))
    gate, rows = checks(controls)
    coverage = rows[f"{task} covers all {reference} control roots"]
    assert not gate["passed"] and not coverage["passed"]
    assert coverage["missing"] == [[spin, root]]


def test_missing_primary_is_a_hard_failure(controls):
    controls["chain-qp"]["states"].remove(state(controls, "chain-qp", -1))
    with pytest.raises(ValueError, match="lacks a primary state"):
        runner.check_controls(controls)


@pytest.mark.parametrize("spin,root", [(-4, 0), (2, 1)])
def test_missing_ed_reference_root_must_not_shrink_required_coverage(controls, spin, root):
    # ED's declared spins/roots remain required even if its receipt drops one.
    controls["chain-ed"]["states"].remove(state(controls, "chain-ed", spin, root))
    assert not runner.check_controls(controls)["passed"]


@pytest.mark.parametrize("task", ["zero-ed", "reverse-ed"])
def test_unqualified_field_control_excitations_must_block_gate(controls, task):
    state(controls, task, 0, 1).update(residual=.1, finite_converged=False)
    controls[task]["diagnostics"]["all_returned_states_converged"] = False
    assert controls[task]["qualified"]
    assert not runner.check_controls(controls)["passed"]


@pytest.mark.parametrize("task", ["atomic-ed", "chain-qp", "chain-dmrg"])
def test_unqualified_primary_and_excitation_fail_independently(controls, task):
    controls[task]["qualified"] = False
    gate, rows = checks(controls)
    assert not gate["passed"]
    assert not rows[f"{task} primary residual/uniqueness qualification"]["passed"]
    controls[task]["qualified"] = True
    state(controls, task, 1, 1).update(residual=.1, finite_converged=False)
    controls[task]["diagnostics"]["all_returned_states_converged"] = False
    gate, rows = checks(controls)
    assert not gate["passed"]
    assert rows[f"{task} primary residual/uniqueness qualification"]["passed"]
    assert not rows[f"{task} all control-root residuals and orthogonality"]["passed"]


@pytest.mark.parametrize("degenerate_endpoint", ["chain-qp", "chain-ed"])
def test_same_sector_degeneracy_compares_energies_only(controls, degenerate_endpoint):
    for task in ("chain-nrg", "chain-ed", "chain-qp", "chain-dmrg", "chain-bandwidth"):
        state(controls, task, 2, 1)["values"]["energy"] = state(controls, task, 2)["values"]["energy"]
    for root in (0, 1):
        state(controls, degenerate_endpoint, 2, root)["degenerate"] = True
    values = state(controls, "chain-qp", 2, 1)["values"]
    for name in runner.OBSERVABLES:
        if name != "energy":
            values[name] += .1
    assert runner.check_controls(controls)["passed"]
    values["energy"] += 1e-5
    gate, rows = checks(controls)
    assert not gate["passed"] and not rows["chain-qp vs chain-ed"]["passed"]


def test_cross_sector_tie_does_not_hide_observable_errors(controls):
    for task in ("chain-nrg", "chain-ed", "chain-qp", "chain-dmrg", "chain-bandwidth"):
        state(controls, task, 2)["values"]["energy"] = state(controls, task, -2)["values"]["energy"]
    assert runner.check_controls(controls)["passed"]
    state(controls, "chain-qp", 2)["values"]["P_up"] += 1e-5
    gate, rows = checks(controls)
    assert not gate["passed"] and not rows["chain-qp vs chain-ed"]["passed"]


def test_absolute_atomic_sign_and_reversal_controls_are_independent(controls):
    for task in ("atomic-nrg", "atomic-ed", "atomic-qp", "atomic-dmrg"):
        for spin in (-1, 1):
            state(controls, task, spin)["values"]["energy"] *= -1
    for task in ("reverse-nrg", "reverse-ed"):
        state(controls, task, -1)["values"]["m_bath"] *= -1
    gate, rows = checks(controls)
    assert not gate["passed"]
    assert rows["atomic-nrg vs atomic-ed"]["passed"] and rows["reverse-nrg vs reverse-ed"]["passed"]
    assert not rows["atomic energies and Zeeman sign"]["passed"]
    assert not rows["NRG spin reversal"]["passed"]


@pytest.mark.parametrize("step", [2e-4, 1e-4])
def test_hf_checks_each_step_and_nonzero_even_magnetization(controls, step):
    md = state(controls, "chain-qp", 0)["values"]["m_dot"]
    assert md != 0
    state(controls, f"hf-{step:g}-+1", 0)["values"]["energy"] -= 2*step*md
    gate, rows = checks(controls)
    assert not gate["passed"]
    assert not rows[f"Hellmann-Feynman h={step:g}"]["passed"]
    other = 1e-4 if step == 2e-4 else 2e-4
    assert rows[f"Hellmann-Feynman h={other:g}"]["passed"]


@pytest.fixture
def campaign(tmp_path, monkeypatch, controls):
    root = tmp_path / "repo"
    here = root / "NRG_comparisons/test3"
    (here / "input").mkdir(parents=True)
    for relative in ("NRG_comparisons/test3/run.py", "NRG_comparisons/test3/finite.py",
                     "NRG_comparisons/test3/manual_run/nrg.py", "NRG_comparisons/test3/manual_run/physics.py",
                     "NRG_comparisons/test1/run.py", "tools/nrginit", "tools/nrg", "tools/math", "tools/model.m"):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"offline identity fixture: {relative}\n")
    monkeypatch.setattr(runner, "ROOT", root)
    monkeypatch.setattr(runner, "HERE", here)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "home"))
    inputs = dict(format_version=1, physical=deepcopy(PHYSICAL), limits=deepcopy(LIMITS))
    runner.write_json(here / "input/physical.json", inputs)
    bath = dict(xi=[.25], weights=[.03], delta=1., bandwidth=100.)
    runner.write_json(here / "input/baths.json", {"cosh8": bath})
    runner.shutil.which.side_effect = lambda name: str(root / "tools" / name)
    runner.subprocess.run.side_effect = None
    runner.subprocess.run.return_value = SimpleNamespace(returncode=1)
    runner.chain_bath.side_effect = None
    runner.chain_bath.return_value = bath

    def execute(work, limits, remaining):
        assert limits == LIMITS and remaining == 7200
        request = json.loads((work / "request.json").read_text())
        result = deepcopy(controls.get(request["task"]["id"], {}))
        result.update(status="completed", request_sha256=runner.digest(request),
                      monitored_seconds=.25, peak_rss_gib=.01)
        if request["task"]["backend"] == "nrg":
            result["artifacts"] = {"param": runner.sha256(work / "param")}
        runner.write_json(work / "result.json", result)
        return result

    execute_mock = Mock(side_effect=execute)
    monkeypatch.setattr(runner, "execute", execute_mock)
    return SimpleNamespace(directory=here / "runs/offline", inputs=inputs, execute=execute_mock, root=root)


def test_campaign_manifest_replay_and_sequential_pilot_gate(campaign, capsys):
    directory = campaign.directory
    controls = runner.campaign("controls", directory)
    assert list(controls) == [task["id"] for task in runner.tasks("controls")]
    assert campaign.execute.call_count == 17
    manifest_bytes = (directory / "campaign.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    assert manifest["identity"]["tasks"] == {stage: runner.tasks(stage) for stage in ("controls", "pilot")}
    assert manifest["identity"]["inputs"] == campaign.inputs
    assert len(manifest["identity"]["source_hashes"]) == 5
    tool = campaign.root / "tools/nrg"
    assert manifest["identity"]["tools"]["nrg"] == dict(path=str(tool), sha256=runner.sha256(tool))
    assert manifest["identity"]["external_sources"] == {
        str(campaign.root / "tools/model.m"): runner.sha256(campaign.root / "tools/model.m")}
    receipts = {path: path.read_bytes() for path in directory.glob("*/result.json")}
    campaign.execute.reset_mock()
    assert runner.campaign("controls", directory) == controls
    campaign.execute.assert_not_called()
    assert (directory / "campaign.json").read_bytes() == manifest_bytes
    assert all(path.read_bytes() == content for path, content in receipts.items())
    results = runner.campaign("pilot", directory)
    assert list(results) == list(controls) + [task["id"] for task in runner.tasks("pilot")]
    assert [entry.args[0].name for entry in campaign.execute.call_args_list] == [
        task["id"] for task in runner.tasks("pilot")]
    request = json.loads((directory / "cosh8-qp-q2/request.json").read_text())
    assert request["bath"]["weights"] == [.03]
    assert json.loads((directory / "controls.json").read_text())["passed"]
    assert "Production has not been launched." in capsys.readouterr().out
    runner.subprocess.Popen.assert_not_called()


@pytest.mark.parametrize("key,value", [("cpus", 5), ("memory_gib", 17), ("task_seconds", 1201),
                                       ("campaign_seconds", 7201), ("mathematica_instances", 2)])
def test_campaign_rejects_unreviewed_resource_caps(campaign, key, value):
    campaign.inputs["limits"][key] = value
    runner.write_json(runner.HERE / "input/physical.json", campaign.inputs)
    with pytest.raises(ValueError, match="resource ceilings"):
        runner.campaign("pilot", campaign.directory)
    assert not campaign.directory.exists()
    runner.shutil.which.assert_not_called()
    campaign.execute.assert_not_called()


def test_campaign_rejects_outside_directory_before_tool_discovery(campaign):
    with pytest.raises(ValueError, match="inside test3/runs"):
        runner.campaign("controls", runner.HERE / "runs/../outside")
    runner.shutil.which.assert_not_called()
    campaign.execute.assert_not_called()


@pytest.mark.parametrize("problem", ["missing-tool", "active-process", "pgrep-error"])
def test_campaign_environment_failure_never_stops_external_work(campaign, problem):
    if problem == "missing-tool":
        runner.shutil.which.side_effect = None
        runner.shutil.which.return_value = None
        message = "load NRG/Mathematica"
    else:
        runner.subprocess.run.return_value.returncode = 0 if problem == "active-process" else 2
        message = "nothing stopped"
    with pytest.raises(RuntimeError, match=message):
        runner.campaign("controls", campaign.directory)
    campaign.execute.assert_not_called()
    runner.os.killpg.assert_not_called()
    assert not (campaign.directory / "campaign.json").exists()


@pytest.mark.parametrize("relative", ["NRG_comparisons/test3/run.py", "NRG_comparisons/test3/input/baths.json",
                                      "tools/nrg", "tools/model.m"])
def test_campaign_identity_digest_changes_reject_replay(campaign, relative):
    runner.campaign("controls", campaign.directory)
    manifest = (campaign.directory / "campaign.json").read_bytes()
    path = campaign.root / relative
    path.write_text(path.read_text()+"\n")
    campaign.execute.reset_mock()
    with pytest.raises(ValueError, match="code/inputs/tools changed"):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert (campaign.directory / "campaign.json").read_bytes() == manifest


@pytest.mark.parametrize("corruption,exception,match", [
    ("request", ValueError, "prepared request changed"),
    ("digest", ValueError, "no automatic retry"),
    ("failed", ValueError, "no automatic retry"),
    ("artifact", ValueError, "raw NRG artifact changed"),
    ("missing", FileNotFoundError, "result.json"),
])
def test_replay_rejects_changed_or_incomplete_receipts_without_retry(campaign, corruption, exception, match):
    runner.campaign("controls", campaign.directory)
    work = campaign.directory / "atomic-nrg"
    if corruption == "request":
        request = json.loads((work / "request.json").read_text())
        request["physical"]["field"] = .04
        runner.write_json(work / "request.json", request)
    elif corruption in {"digest", "failed"}:
        result = json.loads((work / "result.json").read_text())
        result["request_sha256" if corruption == "digest" else "status"] = "failed"
        runner.write_json(work / "result.json", result)
    elif corruption == "artifact":
        (work / "param").write_text("changed deck\n")
    else:
        (work / "result.json").unlink()
    original = {path: path.read_bytes() for path in work.iterdir()}
    campaign.execute.reset_mock()
    with pytest.raises(exception, match=match):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert {path: path.read_bytes() for path in work.iterdir()} == original
    assert not (campaign.directory / "cosh8-qp-q2").exists()


def test_failed_control_gate_blocks_all_pilot_work(campaign):
    runner.campaign("controls", campaign.directory)
    path = campaign.directory / "chain-dmrg/result.json"
    result = json.loads(path.read_text())
    result["diagnostics"]["all_returned_states_converged"] = False
    runner.write_json(path, result)
    campaign.execute.reset_mock()
    with pytest.raises(RuntimeError, match="controls failed; pilot escalation blocked"):
        runner.campaign("pilot", campaign.directory)
    campaign.execute.assert_not_called()
    assert not json.loads((campaign.directory / "controls.json").read_text())["passed"]
    assert not (campaign.directory / "cosh8-qp-q2").exists()


def test_campaign_stops_at_first_worker_failure_without_retry(campaign):
    campaign.execute.side_effect = None
    campaign.execute.return_value = dict(status="failed", error="synthetic failure",
                                         monitored_seconds=.25, peak_rss_gib=.01)
    with pytest.raises(RuntimeError, match="atomic-nrg failed: synthetic failure"):
        runner.campaign("pilot", campaign.directory)
    assert campaign.execute.call_count == 1
    assert not (campaign.directory / "atomic-ed").exists()


@pytest.mark.parametrize("cause", ["timeout", "memory"])
def test_execute_enforces_limits_and_cleans_only_owned_group(tmp_path, cause):
    child = Mock(pid=9001, returncode=-signal.SIGKILL)
    child.poll.return_value = None
    child.wait.side_effect = [subprocess.TimeoutExpired("mock worker", 3), -signal.SIGKILL, -signal.SIGKILL]
    runner.subprocess.Popen.side_effect = lambda *args, **kwargs: (kwargs["preexec_fn"](), child)[1]
    runner.process_memory.side_effect = None
    runner.process_memory.return_value = 17*1024**3 if cause == "memory" else 1024**3
    runner.time.monotonic.side_effect = [100., 106., 107.] if cause == "timeout" else [100., 101.]
    runner.os.killpg.side_effect = None
    events = Mock()
    events.attach_mock(runner.os.killpg, "killpg")
    events.attach_mock(child.wait, "wait")
    result = runner.execute(tmp_path, LIMITS, remaining=5)
    assert result["status"] == "failed" and result["error"] == "memory/time resource limit reached"
    assert result["peak_rss_gib"] == (17 if cause == "memory" else 1)
    assert result["monitored_seconds"] == (7 if cause == "timeout" else 1)
    assert json.loads((tmp_path / "result.json").read_text()) == result
    assert json.loads((tmp_path / "launch.json").read_text()) == dict(supervisor_pid=4100, timeout_seconds=5)
    args, kwargs = runner.subprocess.Popen.call_args
    assert args[0] == [sys.executable, "-B", str(runner.HERE / "run.py"), "worker", str(tmp_path)]
    assert kwargs["start_new_session"] is True and kwargs["cwd"] == runner.ROOT
    assert all(kwargs["env"][name] == "1" for name in
               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONDONTWRITEBYTECODE"))
    runner.resource.setrlimit.assert_called_once_with(runner.resource.RLIMIT_AS, (16*1024**3, 16*1024**3))
    runner.os.sched_setaffinity.assert_called_once_with(0, [2, 4, 7, 9])
    assert events.mock_calls == [call.killpg(9001, signal.SIGTERM), call.wait(timeout=3),
                                 call.killpg(9001, signal.SIGKILL), call.wait(timeout=3), call.wait()]


@pytest.mark.parametrize("exitcode,receipt,error", [
    (0, {"status": "completed", "answer": "synthetic"}, None),
    (7, {"status": "failed", "error": "specific worker failure"}, "specific worker failure"),
    (0, None, "worker exit 0"),
])
def test_execute_finished_worker_receipt_and_absent_group(tmp_path, exitcode, receipt, error):
    if receipt is not None:
        runner.write_json(tmp_path / "result.json", receipt)
    child = Mock(pid=9001, returncode=exitcode)
    child.poll.return_value = exitcode
    runner.subprocess.Popen.side_effect = None
    runner.subprocess.Popen.return_value = child
    runner.os.killpg.side_effect = ProcessLookupError
    result = runner.execute(tmp_path, LIMITS, remaining=7200)
    assert result["status"] == ("completed" if error is None else "failed")
    if error is None:
        assert result["answer"] == "synthetic"
    else:
        assert result["error"] == error
    assert result["peak_rss_gib"] == 0
    assert json.loads((tmp_path / "launch.json").read_text())["timeout_seconds"] == 1200
    assert runner.os.killpg.call_args_list == [call(9001, signal.SIGTERM), call(9001, signal.SIGKILL)]
    assert child.wait.call_count == 3
    runner.process_memory.assert_not_called()


def test_execute_exhausted_budget_does_not_launch_or_signal(tmp_path):
    with pytest.raises(RuntimeError, match="campaign budget exhausted"):
        runner.execute(tmp_path, LIMITS, remaining=0)
    runner.subprocess.Popen.assert_not_called()
    runner.os.killpg.assert_not_called()
    assert not (tmp_path / "launch.json").exists()


@pytest.fixture
def worker(tmp_path):
    request = dict(task=dict(backend="ed"), physical=PHYSICAL, bath={"synthetic": True})
    runner.write_json(tmp_path / "request.json", request)
    runner.write_json(tmp_path / "launch.json", dict(supervisor_pid=4000))
    libc = Mock()
    libc.prctl.return_value = 0
    runner.ctypes.CDLL.side_effect = None
    runner.ctypes.CDLL.return_value = libc
    return tmp_path, request, libc


def test_worker_requires_own_group_before_installing_parentdeath(worker):
    directory, _, libc = worker
    runner.os.getpgrp.return_value = 4000
    with pytest.raises(RuntimeError, match="must own its process group"):
        runner.worker(directory)
    runner.signal.signal.assert_not_called()
    libc.prctl.assert_not_called()
    runner.os.killpg.assert_not_called()


@pytest.mark.parametrize("cause", ["prctl-failed", "parent-changed", "handler"])
def test_worker_parentdeath_term_grace_kill_without_real_signals(worker, offline, cause):
    directory, request, libc = worker
    if cause == "prctl-failed":
        libc.prctl.return_value = -1
    elif cause == "parent-changed":
        runner.os.getppid.return_value = 1
    else:
        offline.run_finite.side_effect = None
        offline.run_finite.return_value = {"synthetic": True}
        runner.worker(directory)
        result = json.loads((directory / "result.json").read_text())
        assert result == dict(synthetic=True, status="completed", request_sha256=runner.digest(request))
        offline.run_finite.assert_called_once_with(request["task"], request["physical"], request["bath"])

    def killpg(group, signum):
        if signum == signal.SIGKILL:
            raise SystemExit("mock SIGKILL")

    runner.os.killpg.side_effect = killpg
    events = Mock()
    events.attach_mock(runner.signal.signal, "signal")
    events.attach_mock(runner.os.killpg, "killpg")
    events.attach_mock(runner.time.sleep, "sleep")
    with pytest.raises(SystemExit, match="mock SIGKILL"):
        if cause == "handler":
            handler = runner.signal.signal.call_args_list[0].args[1]
            handler(signal.SIGTERM, None)
        else:
            runner.worker(directory)
    libc.prctl.assert_called_once_with(1, signal.SIGTERM, 0, 0, 0)
    assert events.mock_calls[-4:] == [call.signal(signal.SIGTERM, signal.SIG_IGN),
                                     call.killpg(4100, signal.SIGTERM), call.sleep(2),
                                     call.killpg(4100, signal.SIGKILL)]
    if cause != "handler":
        offline.run_finite.assert_not_called()
        assert not (directory / "result.json").exists()


def test_worker_failure_receipt_preserves_exception(worker, offline):
    directory, _, _ = worker
    offline.run_finite.side_effect = ValueError("synthetic solver rejection")
    with pytest.raises(ValueError, match="synthetic solver rejection"):
        runner.worker(directory)
    assert json.loads((directory / "result.json").read_text()) == dict(
        status="failed", error="ValueError: synthetic solver rejection")
    runner.os.killpg.assert_not_called()
