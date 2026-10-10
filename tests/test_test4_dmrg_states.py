"""Tiny three-site checkpoints only: no NRG, campaign, bath fit, or pilot."""

from copy import deepcopy
from dataclasses import asdict
import hashlib
import inspect
import json
import logging
from pathlib import Path
import shutil

import numpy as np
from numpy.testing import assert_allclose
import pytest

from NRG_comparisons.test4 import dmrg_states
from NRG_comparisons.test4.finite import run_finite
from NRG_comparisons.test4.physics import OBSERVABLES, PHYSICAL, PRIMARY, finish_states
from qdjj_solver import load_result
from qdjj_solver.common.io import complex_array, json_value


pytestmark = pytest.mark.dmrg
BATH = dict(xi=[.23], weights=[.019], delta=1., bandwidth=100.)
COMPLEX = dict(PHYSICAL, phi=.83, wn=7., ws=5.)
TASK = dict(backend="dmrg", twice_sz=0, roots=1, chi=32, max_sweeps=10)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def initial_for(directory, report):
    return dict(path=str(directory/report["checkpoint"]["path"]), sha256=report["checkpoint"]["sha256"])


@pytest.fixture(scope="module")
def small_runs(tmp_path_factory):
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    directory = tmp_path_factory.mktemp("test4-three-sites")
    runs = {}
    original_solve = dmrg_states.solve
    live = []

    def capture(*args, **kwargs):
        result = original_solve(*args, **kwargs)
        live.append(result)
        return result

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(dmrg_states, "solve", capture)
        for label, physical, spin, roots in (
            ("even", PHYSICAL, 0, 1), ("odd_down", PHYSICAL, -1, 1),
            ("odd_up", PHYSICAL, 1, 1), ("complex", COMPLEX, 0, 2),
        ):
            target = directory/label
            target.mkdir()
            task = dict(TASK, twice_sz=spin, roots=roots)
            before = deepcopy((task, physical, BATH))
            report = dmrg_states.run_state(task, physical, BATH, target)
            assert (task, physical, BATH) == before
            runs[label] = (target, report, live[-1], physical)
    return runs


def test_exact_public_apis_and_bath_builder():
    assert str(inspect.signature(dmrg_states.build_model)) == "(physical, bath)"
    assert str(inspect.signature(dmrg_states.rows_from_result)) == "(result, physical)"
    assert str(inspect.signature(dmrg_states.run_state)) == "(task, physical, bath_record, directory, initial=None) -> dict"
    from qdjj_solver import DiscreteBath

    h = dmrg_states.build_model(PHYSICAL, BATH)
    assert h.fingerprint() == dmrg_states.build_model(PHYSICAL, DiscreteBath.from_record(BATH)).fingerprint()
    assert h.metadata["gamma"] == .4
    assert h.metadata["detuning"] == .5
    assert h.metadata["phi"] == 3*np.pi/5
    assert h.metadata["field"] == .03
    assert h.metadata["rho_wn"] == h.metadata["rho_ws"] == 0
    assert not h.metadata["paired_mode_compression"]
    assert h.metadata["coordinates"]["restriction"] == "none"
    assert [bath["weights"] for bath in h.metadata["bath"]] == [[.019], [.019]]
    assert "rho_ws_derivative" not in h.observables
    for bath in (dict(BATH, delta=2.), dict(BATH, bandwidth=50.), dict(BATH, weights=[0.]),
                 dict(BATH, xi=[]), dict(BATH, xi=[float("nan")])):
        with pytest.raises(ValueError):
            dmrg_states.build_model(PHYSICAL, bath)


def test_primary_complex_ed_and_frozen_full_qp(small_runs):
    for label, (_, report, live, physical) in small_runs.items():
        spin, roots = live.prepared.sector.twice_sz, len(live.energies)
        for backend in ("ed", "qp"):
            expected = run_finite(dict(backend=backend, spins=[spin], roots=roots, cutoff=None), physical, BATH)
            for actual, reference in zip(report["states"], expected["states"], strict=True):
                assert set(actual["values"]) == set(OBSERVABLES)
                # Full tiny Hilbert space; tolerance allows Lanczos roundoff, not bath error.
                assert_allclose([actual["values"][key] for key in OBSERVABLES],
                                [reference["values"][key] for key in OBSERVABLES], rtol=0, atol=2e-8)
                assert actual["residual"] < 1e-8
                assert actual["finite_converged"] and not actual["degenerate"]
        assert report["qualified"]
        assert any(abs(c.imag) > 1e-5 for c in live.hamiltonian.operator.terms.values())
        assert abs(report["states"][0]["values"]["current"]) > 1e-5
        if label != "complex":
            assert report["states"][0]["id"] == label


def test_live_and_loaded_rows_checkpoint_receipt_and_scope(small_runs):
    for directory, report, live, physical in small_runs.values():
        receipt = report["checkpoint"]
        path = directory/receipt["path"]
        assert receipt["path"] == "eigenstates.json"
        assert receipt["sha256"] == digest(path)
        artifact = receipt["state_artifact"]
        assert Path(artifact["file"]).name == artifact["file"]
        assert artifact["sha256"] == digest(directory/artifact["file"])
        loaded = load_result(path)
        assert isinstance(live.metadata["gram"], np.ndarray)
        assert set(loaded.metadata["gram"]) == {"real", "imag"}
        assert dmrg_states.rows_from_result(live, physical) == report["states"]
        assert dmrg_states.rows_from_result(loaded, physical) == report["states"]
        assert receipt["hamiltonian_sha256"] == live.hamiltonian.fingerprint() == loaded.hamiltonian.fingerprint()
        assert receipt["coordinate_id"] == loaded.prepared.coordinate_id
        assert receipt["twice_sz"] == loaded.prepared.sector.twice_sz
        assert receipt["roots"] == len(loaded.states) and receipt["chi"] == 32
        assert loaded.metadata["test4_dmrg_states"]["root_prefix"] == list(range(receipt["roots"]))
        diagnostics = report["diagnostics"]
        assert diagnostics["physical"] == physical and diagnostics["bath"] == BATH
        assert diagnostics["mode_order"] == [2, 3, 0, 1, 4, 5]
        assert loaded.prepared.groups == ((2, 3), (0, 1), (4, 5))
        assert not diagnostics["weights_renormalized"]
        assert diagnostics["bath_weight_sum_per_lead"] == .019
        assert diagnostics["coverage_scope"].startswith("single-sector")
        for key in ("spectral_completeness_established", "lowest_root_completeness_established",
                    "bath_convergence_checked", "chi_convergence_checked", "full_qp_comparison_checked"):
            assert diagnostics[key] is False
        assert all("excitation" not in row and "subgap" not in row for row in report["states"])
        assert not {"coverage", "status", "request_sha256", "monitored_seconds"} & report.keys()
        assert diagnostics["norm_error_tolerance"] == 1e-8
        assert max(diagnostics["norm_error"]) < 1e-8
        assert diagnostics["timings"]["total"] > 0
        assert all(info["trial_count"] == 2 for info in diagnostics["root_selection_seed_independence"])
        for info in diagnostics["solver_metadata"]["roots"]:
            assert len(info["seed_trials"]) == 2 and info["sweep_statistics"]
        json.dumps(report, allow_nan=False)
    combined = [deepcopy(small_runs[name][1]["states"][0]) for name in PRIMARY.values()]
    coverage = finish_states(combined)
    assert len(coverage["ground_ids"]) == 1
    assert sum(row["excitation"] > 0 for row in combined) == 2


@pytest.mark.parametrize("roots", [1, 2])
def test_warm_chi64_continues_intact_prefix_and_reoptimizes_root_zero(small_runs, tmp_path, monkeypatch, caplog, roots):
    directory, old, _, physical = small_runs["even"]
    before = {path.name: path.read_bytes() for path in directory.iterdir()}
    original_solve = dmrg_states.solve
    seen = []

    def checked_solve(h, **kwargs):
        assert kwargs["options"]["eigenpairs"] == roots
        assert kwargs["options"]["chi_max"] == 64
        assert len(kwargs["initial"].states) == 1
        assert kwargs["initial"].metadata["test4_dmrg_states"]["root_prefix"] == [0]
        result = original_solve(h, **kwargs)
        seen.append(result)
        return result

    monkeypatch.setattr(dmrg_states, "solve", checked_solve)
    logger = logging.getLogger("qdjj_solver.dmrg_solver.solver")
    old_level = logger.level
    report = dmrg_states.run_state(dict(TASK, chi=64, roots=roots), physical, BATH, tmp_path,
                                   initial_for(directory, old))
    assert logger.level == old_level
    assert any("root=0 trial=0" in record.message for record in caplog.records)
    assert report["qualified"] and len(report["states"]) == roots
    assert seen[0].metadata["roots"][0]["seed_trials"][0]["seed"] == "continuation"
    assert all(path.read_bytes() == before[path.name] for path in directory.iterdir())
    assert {path.name for path in directory.iterdir()} == set(before)
    loaded = load_result(tmp_path/"eigenstates.json")
    assert loaded.metadata["roots"][0]["seed_trials"][0]["seed"] == "continuation"
    assert_allclose(complex_array(loaded.metadata["gram"]), np.eye(roots), atol=1e-8, rtol=0)
    expected = run_finite(dict(backend="ed", spins=[0], roots=roots), physical, BATH)
    for row, reference in zip(report["states"], expected["states"], strict=True):
        assert_allclose([row["values"][name] for name in OBSERVABLES],
                        [reference["values"][name] for name in OBSERVABLES], atol=2e-8, rtol=0)
    if roots == 2:
        assert [row["root"] for row in report["states"]] == [0, 1]
        assert report["states"][1]["id"] == "sz+0-r1"
        assert loaded.metadata["roots"][1]["seed_trials"][0]["seed"] != "continuation"
        assert any("root=1 trial=0" in record.message for record in caplog.records)


@pytest.mark.parametrize("changes", [dict(backend="qp"), dict(twice_sz=True), dict(twice_sz=.0),
    dict(roots=True), dict(roots=0), dict(roots=3), dict(chi=16), dict(chi=64.),
    dict(seed_trials=True), dict(seed_trials=0), dict(seed_trials=2.), dict(seed=-1), dict(seed=1.),
    dict(min_sweeps=7), dict(min_sweeps=True), dict(max_sweeps=41), dict(max_sweeps=7),
    dict(max_sweeps=10.), dict(min_sweeps=11), dict(layout="chain"), dict(cutoff=2),
    dict(threads=4), dict(group_size=1), dict(mode_order=[0, 1]), dict(orthogonal_to=[]),
    dict(eigenpairs=1), dict(chi_max=64), dict(twice_sz=4), dict(twice_sz=3, roots=2)])
def test_task_preflight_before_solve(tmp_path, monkeypatch, changes):
    monkeypatch.setattr(dmrg_states, "solve", lambda *a, **kw: pytest.fail("solve before validation"))
    with pytest.raises(ValueError):
        dmrg_states.run_state(dict(TASK, **changes), PHYSICAL, BATH, tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name", ["eigenstates.json", "eigenstates.old.mps.h5", "result.json", "other.json", "other"])
def test_refuse_existing_or_unrelated_outputs(tmp_path, monkeypatch, name):
    monkeypatch.setattr(dmrg_states, "solve", lambda *a, **kw: pytest.fail("solve in occupied directory"))
    path = tmp_path/name
    path.write_text("preserve")
    with pytest.raises(ValueError, match="existing checkpoint or unrelated"):
        dmrg_states.run_state(TASK, PHYSICAL, BATH, tmp_path)
    assert path.read_text() == "preserve"


def test_directory_must_exist_and_no_symlinks(tmp_path):
    with pytest.raises(ValueError, match="already exist"):
        dmrg_states.run_state(TASK, PHYSICAL, BATH, tmp_path/"missing")
    link = tmp_path/"link"
    link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="ordinary directory"):
        dmrg_states.run_state(TASK, PHYSICAL, BATH, link)
    link.unlink()
    (tmp_path/"request.json").symlink_to(tmp_path/"missing")
    with pytest.raises(ValueError, match="existing checkpoint or unrelated"):
        dmrg_states.run_state(TASK, PHYSICAL, BATH, tmp_path)


@pytest.mark.parametrize("kind", ["hash", "relative", "label", "extra", "physical", "bath", "sector"])
def test_initial_rejections_before_load(small_runs, tmp_path, monkeypatch, kind):
    directory, old, _, physical = small_runs["even"]
    initial = initial_for(directory, old)
    task, bath = dict(TASK, chi=64), deepcopy(BATH)
    if kind == "hash":
        initial["sha256"] = "0"*64
    elif kind == "relative":
        initial["path"] = "eigenstates.json"
    elif kind == "label":
        initial = "even-chi32"
    elif kind == "extra":
        initial["root"] = 1
    elif kind == "physical":
        physical = dict(physical, gamma=.41)
    elif kind == "bath":
        bath["weights"] = [.02]
    else:
        task["twice_sz"] = -1
    monkeypatch.setattr(dmrg_states, "load_result", lambda *a, **kw: pytest.fail("loaded rejected initial"))
    monkeypatch.setattr(dmrg_states, "solve", lambda *a, **kw: pytest.fail("solved rejected initial"))
    with pytest.raises(ValueError):
        dmrg_states.run_state(task, physical, bath, tmp_path, initial)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("kind", ["layout", "order", "option_order", "group", "groups", "coordinate",
                                      "prefix", "missing_prefix", "count", "physical", "bath", "hamiltonian"])
def test_initial_metadata_checked_before_loading_hdf5(small_runs, tmp_path, monkeypatch, kind):
    directory, report, _, physical = small_runs["even"]
    record = json.loads((directory/"eigenstates.json").read_text())
    metadata = record["calculation"]
    workflow = metadata["test4_dmrg_states"]
    if kind == "layout":
        workflow["layout"] = "chain"
    elif kind == "order":
        metadata["mode_order"] = list(range(6))
    elif kind == "option_order":
        metadata["options"]["mode_order"] = list(range(6))
    elif kind == "group":
        metadata["options"]["group_size"] = 1
    elif kind == "groups":
        metadata["groups"] = [[0, 1], [2, 3], [4, 5]]
    elif kind == "coordinate":
        metadata["coordinate_id"] = "0"*64
    elif kind == "prefix":
        workflow["root_prefix"] = [1]
    elif kind == "missing_prefix":
        del metadata["test4_dmrg_states"]
    elif kind == "count":
        metadata["options"]["eigenpairs"] = 2
    elif kind == "physical":
        workflow["physical"]["field"] = -.03
    elif kind == "bath":
        workflow["bath"]["xi"] = [.24]
    else:
        record["hamiltonian"]["operator"][0]["real"] += 1.
    path = tmp_path/"wrong.json"
    path.write_text(json.dumps(record))
    target = tmp_path/"target"
    target.mkdir()
    monkeypatch.setattr(dmrg_states, "load_result", lambda *a, **kw: pytest.fail("loaded bad metadata"))
    with pytest.raises(ValueError):
        dmrg_states.run_state(dict(TASK, chi=64), physical, BATH, target,
                             dict(path=str(path), sha256=digest(path)))
    assert list(target.iterdir()) == []
    assert digest(directory/"eigenstates.json") == report["checkpoint"]["sha256"]


def test_companion_hash_checked_and_no_root_one_relabel(small_runs, tmp_path):
    directory, report, _, physical = small_runs["complex"]
    with pytest.raises(ValueError, match="complete root-0 prefix"):
        dmrg_states.run_state(dict(TASK, chi=64), physical, BATH, tmp_path, initial_for(directory, report))
    source = tmp_path/"source"
    shutil.copytree(directory, source)
    artifact = source/report["checkpoint"]["state_artifact"]["file"]
    artifact.write_bytes(artifact.read_bytes()+b"broken")
    target = tmp_path/"target"
    target.mkdir()
    with pytest.raises(ValueError, match="MPS checkpoint checksum"):
        dmrg_states.run_state(dict(TASK, chi=64, roots=2), physical, BATH, target, initial_for(source, report))
    assert not (target/"eigenstates.json").exists()


def test_rows_reject_different_physical_and_hamiltonian(small_runs):
    _, _, live, physical = small_runs["even"]
    with pytest.raises(ValueError, match="Hamiltonian differs"):
        dmrg_states.rows_from_result(live, dict(physical, nu=.74))
    changed = deepcopy(live)
    changed.hamiltonian.observables["phase_derivative"] *= 2
    with pytest.raises(ValueError, match="Hamiltonian differs"):
        dmrg_states.rows_from_result(changed, physical)


@pytest.mark.parametrize("failure", ["residual", "sweeps", "norm", "overlap", "gram", "order", "degenerate", "no_residual"])
def test_unqualified_return_still_saved(small_runs, tmp_path, monkeypatch, failure):
    _, _, live, physical = small_runs["complex"]
    result = deepcopy(live)
    if failure == "residual":
        result.residuals[0] = 2e-7
    elif failure == "sweeps":
        result.metadata["roots"][0]["sweep_energy_converged"] = False
    elif failure == "norm":
        result.metadata["roots"][0]["norm_error"] = 2e-8
    elif failure == "overlap":
        result.metadata["roots"][1]["lower_root_overlap"] = 2e-8
    elif failure == "gram":
        result.metadata["gram"][0, 1] = 2e-8j
    elif failure == "order":
        result.energies = result.energies[::-1]
    elif failure == "degenerate":
        result.energies[:] = result.energies[0]
    else:
        result.residuals = None
    monkeypatch.setattr(dmrg_states, "solve", lambda *a, **kw: result)
    report = dmrg_states.run_state(dict(TASK, roots=2), physical, BATH, tmp_path)
    assert not report["qualified"]
    assert len(load_result(tmp_path/"eigenstates.json").states) == 2
    if failure == "degenerate":
        assert all(row["degenerate"] for row in report["states"])
    else:
        assert any(not row["finite_converged"] for row in report["states"])
    assert all("excitation" not in row and "subgap" not in row for row in report["states"])
    assert report["checkpoint"]["sha256"] == digest(tmp_path/"eigenstates.json")


def test_checkpoint_precedes_conversion_and_is_not_rewritten(small_runs, tmp_path, monkeypatch):
    _, _, live, physical = small_runs["even"]
    monkeypatch.setattr(dmrg_states, "solve", lambda *a, **kw: deepcopy(live))
    saved = {}

    def fail_conversion(result, physical):
        path = tmp_path/"eigenstates.json"
        saved["bytes"] = path.read_bytes()
        assert len(load_result(path).states) == 1
        raise ValueError("injected conversion failure")

    monkeypatch.setattr(dmrg_states, "rows_from_result", fail_conversion)
    with pytest.raises(ValueError, match="injected conversion"):
        dmrg_states.run_state(TASK, physical, BATH, tmp_path)
    assert (tmp_path/"eigenstates.json").read_bytes() == saved["bytes"]
    assert not (tmp_path/"result.json").exists()


@pytest.mark.parametrize("chi", [32, 64, 128, 256])
def test_allowed_options_and_parent_files_without_new_numerics(small_runs, tmp_path, monkeypatch, chi):
    from qdjj_solver.dmrg_solver.solver import SolverOptions

    _, _, live, physical = small_runs["even"]
    for name in ("request.json", "launch.json", "worker.log"):
        (tmp_path/name).write_text("parent-owned")

    def fake_solve(h, **kwargs):
        options = kwargs["options"]
        assert options["chi_max"] == chi and options["eigenpairs"] == 1
        assert options["max_sweeps"] == 40 and options["min_sweeps"] == 8
        assert options["seed_trials"] == 2 and options["seed"] == 1729
        assert options["threads"] == 1 and options["group_size"] == 2
        assert options["mode_order"] == [2, 3, 0, 1, 4, 5]
        assert options["energy_tolerance"] == 1e-12 and options["residual_tolerance"] == 1e-7
        assert options["lanczos_probability_tolerance"] == 1e-22
        assert options["calculate_residuals"] and not options["require_convergence"]
        assert kwargs["initial"] is None and kwargs["cutoff"] is None
        result = deepcopy(live)
        result.metadata["options"] = asdict(SolverOptions(**options))
        return result

    monkeypatch.setattr(dmrg_states, "solve", fake_solve)
    task = {key: value for key, value in dict(TASK, chi=chi, id="parent-label").items() if key != "max_sweeps"}
    report = dmrg_states.run_state(task, physical, BATH, tmp_path)
    assert report["checkpoint"]["chi"] == report["diagnostics"]["chi"] == chi
    for name in ("request.json", "launch.json", "worker.log"):
        assert (tmp_path/name).read_text() == "parent-owned"
    assert load_result(tmp_path/"eigenstates.json").metadata["options"]["chi_max"] == chi


def test_interrupted_solve_has_no_checkpoint(tmp_path, monkeypatch):
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    logger = logging.getLogger("qdjj_solver.dmrg_solver.solver")
    before = logger.level, list(logger.handlers)
    monkeypatch.setattr(dmrg_states, "solve", interrupted)
    with pytest.raises(KeyboardInterrupt):
        dmrg_states.run_state(TASK, PHYSICAL, BATH, tmp_path)
    assert list(tmp_path.iterdir()) == []
    assert (logger.level, logger.handlers) == before


@pytest.mark.parametrize("trials", [1, 3])
def test_positive_seed_trials_forwarded_without_extra_dmrg(tmp_path, monkeypatch, trials):
    class StopBeforeNumerics(Exception):
        pass

    def stop(h, **kwargs):
        assert kwargs["options"]["seed_trials"] == trials
        raise StopBeforeNumerics

    monkeypatch.setattr(dmrg_states, "solve", stop)
    with pytest.raises(StopBeforeNumerics):
        dmrg_states.run_state(dict(TASK, seed_trials=trials), PHYSICAL, BATH, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_gap_units_and_loaded_gram_decoder_without_extra_dmrg(small_runs):
    _, report, live, physical = small_runs["complex"]
    dimensional = {"gap", "bandwidth", "u", "gamma", "field", "wn", "ws"}
    scaled = {key: value*2 if key in dimensional else value for key, value in physical.items()}
    bath = dict(BATH, xi=[2*BATH["xi"][0]], delta=2., bandwidth=200.)
    result = deepcopy(live)
    result.hamiltonian = dmrg_states.build_model(scaled, bath)
    result.energies *= 2
    result.residuals *= 2
    result.observables["phase_derivative"] *= 2
    result.metadata["gram"] = json_value(result.metadata["gram"])
    rows = dmrg_states.rows_from_result(result, scaled)
    for row, original in zip(rows, report["states"], strict=True):
        assert_allclose([row["values"][key] for key in OBSERVABLES],
                        [original["values"][key] for key in OBSERVABLES], atol=1e-14, rtol=0)
