"""Bounded three-site numerics and checkpoint guards; never launch a campaign."""

from copy import deepcopy
from dataclasses import asdict
import hashlib
import inspect
import json
import logging
import shutil
import subprocess
import sys
import weakref

import numpy as np
from numpy.testing import assert_allclose
import pytest

from NRG_comparisons.test4 import continuum_finite as worker, dmrg_states
from NRG_comparisons.test4.finite import run_finite
from NRG_comparisons.test4.physics import OBSERVABLES, PHYSICAL
from qdjj_solver import load_result
from qdjj_solver.common.io import complex_array, write_json


BATH = dict(xi=[.23], weights=[.019], delta=1., bandwidth=100.)
COMPLEX = dict(PHYSICAL, phi=.83, wn=7., ws=5.)
QP = dict(backend="qp", cutoff=4, method="sparse", spins=[0, -1, 1, -2, 2],
          roots=2, roots_by_spin={"0": 3, "-2": 1, "2": 1})
DMRG = dict(backend="dmrg", twice_sz=0, roots=1, chi=128, max_sweeps=10)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def initial_for(directory, report):
    return dict(path=str(directory/report["checkpoint"]["path"]), sha256=report["checkpoint"]["sha256"])


def compare_rows(actual, expected):
    for row, reference in zip(actual, expected, strict=True):
        assert set(row["values"]) == set(OBSERVABLES)
        # Exact tiny finite Hilbert space; this allows solver roundoff, not bath error.
        assert_allclose([row["values"][key] for key in OBSERVABLES],
                        [reference["values"][key] for key in OBSERVABLES], rtol=0, atol=2e-8)
        for key in ("twice_sz", "parity", "root", "id", "finite_converged", "degenerate"):
            assert row[key] == reference[key]
        assert row["residual"] < 1e-8


def verify_output_artifacts(directory, task, physical, bath, report):
    assert not {"status", "request_sha256", "monitored_seconds", "peak_rss_gib", "artifacts"} & report.keys()
    records = ([report["checkpoint"]] if task["backend"] == "dmrg" else
               [sector["scalar_receipt"] for sector in report["diagnostics"]["sectors"]])
    for receipt in records:
        path = directory / receipt["path"]
        assert digest(path) == receipt["sha256"]
        record = json.loads(path.read_text())
        workflow = record["calculation"]["test4_dmrg_states" if task["backend"] == "dmrg" else "test4_continuum_finite"]
        assert workflow["task"] == task and workflow["physical"] == physical and workflow["bath"] == bath
        if task["backend"] == "dmrg":
            companion = record["state_artifact"]
            assert digest(directory / companion["file"]) == companion["sha256"]


def test_api_and_lazy_import():
    assert str(inspect.signature(worker.run_task)) == "(task, physical, bath_record, directory, initial=None) -> dict"
    assert worker.FORMAT_VERSION == 1
    assert worker.build_model is dmrg_states.build_model
    assert worker.rows_from_result is dmrg_states.rows_from_result
    subprocess.run([sys.executable, "-B", "-c",
                    "import sys; from NRG_comparisons.test4 import continuum_finite; "
                    "assert not {'tenpy', 'h5py'} & sys.modules.keys()"], check=True)


@pytest.mark.parametrize("scale", [1., 2.])
def test_qp_matches_frozen_sparse_and_independent_ed(tmp_path, monkeypatch, scale):
    dimensional = {"gap", "bandwidth", "u", "gamma", "field", "wn", "ws"}
    physical = {key: value*scale if key in dimensional else value for key, value in COMPLEX.items()}
    bath = dict(BATH, xi=[.23*scale], delta=scale, bandwidth=100*scale)
    before = deepcopy((QP, physical, bath))
    original = worker.solve
    live = []

    def checked_solve(h, **kwargs):
        assert all(ref() is None for ref in live), "previous sector vectors retained across solves"
        options = kwargs["options"]
        assert options == dict(method="sparse", eigenpairs={0: 4, -1: 3, 1: 3, -2: 2, 2: 2}[
            kwargs["sector"].twice_sz], residual_tolerance=1e-8*scale, tolerance=1e-12,
            ncv=40, seed=1729, threads=1, max_dimension=10_000_000, max_memory_gib=48., max_nnz=500_000_000)
        result = original(h, **kwargs)
        live.append(weakref.ref(result))
        return result

    monkeypatch.setattr(worker, "solve", checked_solve)
    report = worker.run_task(QP, physical, bath, tmp_path)
    assert (QP, physical, bath) == before
    assert len(live) == 5 and all(ref() is None for ref in live)
    expected = run_finite(QP, physical, bath)
    compare_rows(report["states"], expected["states"])
    compare_rows(report["states"], run_finite(dict(QP, backend="ed", cutoff=None), physical, bath)["states"])
    assert len(report["states"]) == 9 and report["qualified"] == expected["qualified"]
    for row, reference in zip(report["states"], expected["states"], strict=True):
        assert row["subgap"] == reference["subgap"]
        assert row["excitation"] == pytest.approx(reference["excitation"], abs=1e-12)
    for key, value in expected["coverage"].items():
        assert report["coverage"][key] == (pytest.approx(value) if isinstance(value, float) else value)
    for key, value in expected["diagnostics"].items():
        if key != "sectors":
            assert report["diagnostics"][key] == value
    assert not report["diagnostics"]["weights_renormalized"]
    assert report["diagnostics"]["bath_weight_sum_per_lead"] == .019
    assert len(list(tmp_path.iterdir())) == 5
    for info, reference in zip(report["diagnostics"]["sectors"], expected["diagnostics"]["sectors"], strict=True):
        for key in ("dimension", "requested_roots", "returned_roots", "solved_roots", "complete_sector_diagonalized",
                    "degeneracy_scope"):
            assert info[key] == reference[key]
        path = tmp_path/info["scalar_receipt"]["path"]
        assert digest(path) == info["scalar_receipt"]["sha256"]
        scalar = json.loads(path.read_text())
        assert scalar["format"] == "qdjj-eigenstates" and scalar["backend"] == "qp"
        assert not {"state_artifact", "state_amplitudes", "coverage", "qualified", "states"} & scalar.keys()
        assert scalar["qp_weights"] == info["qp_weights"]
        assert_allclose(np.sum(info["qp_weights"], axis=1), 1., rtol=0, atol=1e-12)
        assert scalar["calculation"] == info["solver_metadata"]
        assert scalar["timings_seconds"] == info["timings"]
        assert scalar["hamiltonian_sha256"] == worker.build_model(physical, bath).fingerprint()
        assert info["timings"]["total"] > 0 and info["solver_metadata"]["nnz"] > 0
        assert info["solver_metadata"]["method"] == "sparse"
        assert scalar["calculation"]["test4_continuum_finite"]["root_prefix"] == list(range(info["solved_roots"]))
    json.dumps(report, allow_nan=False)
    verify_output_artifacts(tmp_path, QP, physical, bath, report)


@pytest.mark.parametrize("changes", [
    dict(cutoff=None), dict(cutoff=2), dict(cutoff=5), dict(cutoff=12), dict(cutoff=True), dict(cutoff=4.),
    dict(method="dense"), dict(method="auto"), dict(method="matrix-free"), dict(method=None),
    dict(spins=[0, -1, 1]), dict(spins=[False, -1, 1, -2, 2]), dict(spins=[0, 1, -1, -2, 2]),
    dict(roots=0), dict(roots=True), dict(roots=2.), dict(roots_by_spin={}), dict(roots_by_spin=None),
    dict(roots_by_spin={0: 3}), dict(roots_by_spin={"3": 1}), dict(roots_by_spin={"0": True}),
    dict(max_dimension=10_000_001), dict(max_dimension=0), dict(max_dimension=True), dict(max_dimension=10.),
    dict(max_nnz=500_000_001), dict(max_nnz=0), dict(max_nnz=True), dict(max_nnz=100.),
    dict(max_memory_gib=48.1), dict(max_memory_gib=0), dict(max_memory_gib=True),
    dict(max_memory_gib=float("nan")), dict(max_memory_gib=float("inf")),
    dict(seed=17), dict(seed=1729.), dict(ncv=80), dict(threads=4), dict(residual_tolerance=1e-6),
    dict(options={}), dict(eigenpairs=1),
])
def test_qp_guards_before_solve(tmp_path, monkeypatch, changes):
    monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solve before validation"))
    with pytest.raises(ValueError):
        worker.run_task(dict(QP, **changes), PHYSICAL, BATH, tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("changes", [
    dict(twice_sz=True), dict(twice_sz=0.), dict(twice_sz=-2), dict(twice_sz=2),
    dict(roots=0), dict(roots=3), dict(roots=True), dict(roots=1.),
    dict(chi=32), dict(chi=64), dict(chi=2048), dict(chi=128.), dict(chi=True),
    dict(seed_trials=1), dict(seed_trials=3), dict(seed_trials=2.), dict(seed=-1), dict(seed=True),
    dict(min_sweeps=7), dict(min_sweeps=True), dict(max_sweeps=41), dict(max_sweeps=7),
    dict(max_sweeps=10.), dict(min_sweeps=11), dict(layout="chain"), dict(cutoff=4),
    dict(mode_order=[0, 1]), dict(group_size=1), dict(threads=4), dict(orthogonal_to=[]),
    dict(eigenpairs=1), dict(chi_max=512), dict(solver={}), dict(sector={}), dict(initial={}),
])
def test_dmrg_guards_before_solve(tmp_path, monkeypatch, changes):
    monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solve before validation"))
    with pytest.raises(ValueError):
        worker.run_task(dict(DMRG, **changes), PHYSICAL, BATH, tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("task", [None, {}, {"backend": "ed"}, {"backend": "nrg"}])
def test_backend_guard(tmp_path, task):
    with pytest.raises(ValueError, match="backend"):
        worker.run_task(task, PHYSICAL, BATH, tmp_path)


@pytest.mark.parametrize("task", [QP, DMRG])
@pytest.mark.parametrize("name", ["eigenstates.json", "old.mps.h5", "qp-sz+0.json", "result.json", "other"])
def test_no_overwrite_or_unrelated_files(tmp_path, monkeypatch, task, name):
    path = tmp_path/name
    path.write_text("preserve")
    monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solve in occupied directory"))
    with pytest.raises(ValueError, match="fresh directory"):
        worker.run_task(task, PHYSICAL, BATH, tmp_path)
    assert path.read_text() == "preserve"


@pytest.mark.parametrize("task", [QP, DMRG])
def test_ordinary_directory_and_parent_files_only(tmp_path, monkeypatch, task):
    with pytest.raises(ValueError, match="already exist"):
        worker.run_task(task, PHYSICAL, BATH, tmp_path/"missing")
    link = tmp_path/"link"
    link.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="ordinary directory"):
        worker.run_task(task, PHYSICAL, BATH, link)
    link.unlink()
    (tmp_path/"request.json").symlink_to(tmp_path/"missing")
    with pytest.raises(ValueError, match="fresh directory"):
        worker.run_task(task, PHYSICAL, BATH, tmp_path)
    (tmp_path/"request.json").unlink()
    (tmp_path/"launch.json").mkdir()
    with pytest.raises(ValueError, match="fresh directory"):
        worker.run_task(task, PHYSICAL, BATH, tmp_path)
    (tmp_path/"launch.json").rmdir()
    for name in ("request.json", "launch.json", "worker.log"):
        (tmp_path/name).write_text("parent-owned")
    logger = logging.getLogger("qdjj_solver.dmrg_solver.solver")
    previous = logger.level, list(logger.handlers)

    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(worker, "solve", interrupted)
    with pytest.raises(KeyboardInterrupt):
        worker.run_task(task, PHYSICAL, BATH, tmp_path)
    assert (logger.level, logger.handlers) == previous
    assert {path.name: path.read_text() for path in tmp_path.iterdir()} == {
        name: "parent-owned" for name in ("request.json", "launch.json", "worker.log")}


@pytest.mark.parametrize("task", [QP, DMRG])
@pytest.mark.parametrize("change", [dict(delta=2.), dict(bandwidth=50.), dict(xi=[]),
                                   dict(xi=[float("nan")]), dict(weights=[0.]), dict(weights=[.1, .2])])
def test_bath_guards(tmp_path, monkeypatch, task, change):
    monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solve with invalid bath"))
    with pytest.raises(ValueError):
        worker.run_task(task, PHYSICAL, dict(BATH, **change), tmp_path)


@pytest.mark.parametrize("cutoff", [4, 6, 8, 10])
def test_qp_cutoffs_and_lower_limits_forwarded_without_solving(tmp_path, monkeypatch, cutoff):
    bath = dict(BATH, xi=[-.4, .23, .5], weights=[.01, .019, .02])

    def stop(h, **kwargs):
        assert kwargs["cutoff"] == cutoff and kwargs["backend"] == "qp"
        assert kwargs["options"]["max_dimension"] == 1_000_000
        assert kwargs["options"]["max_memory_gib"] == 8.
        assert kwargs["options"]["max_nnz"] == 50_000_000
        raise InterruptedError("before numerics")

    monkeypatch.setattr(worker, "solve", stop)
    with pytest.raises(InterruptedError, match="before numerics"):
        worker.run_task(dict(QP, cutoff=cutoff, max_dimension=1_000_000, max_memory_gib=8., max_nnz=50_000_000),
                        PHYSICAL, bath, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_qp_preflight_dimension_initial_and_oversized_cutoff(tmp_path, monkeypatch):
    monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solve before preflight"))
    with pytest.raises(MemoryError, match="dimension"):
        worker.run_task(dict(QP, max_dimension=1), PHYSICAL, BATH, tmp_path)
    with pytest.raises(ValueError, match="initial"):
        worker.run_task(QP, PHYSICAL, BATH, tmp_path, initial={})
    with pytest.raises(ValueError, match="bath mode count"):
        worker.run_task(dict(QP, cutoff=6), PHYSICAL, BATH, tmp_path)


@pytest.mark.parametrize("failure", ["next_sector", "conversion", "incomplete_prefix"])
def test_qp_failure_keeps_only_scalar_sectors_and_never_globalizes(tmp_path, monkeypatch, failure):
    original = worker.solve
    calls = []

    def solve(h, **kwargs):
        calls.append(kwargs["sector"].twice_sz)
        if len(calls) == 2:
            raise RuntimeError("injected next sector failure")
        result = original(h, **kwargs)
        if failure == "incomplete_prefix":
            result.energies = result.energies[:-1]
        return result

    def fail_conversion(*args, **kwargs):
        assert (tmp_path/"qp-sz+0.json").is_file()
        raise ValueError("injected conversion failure")

    monkeypatch.setattr(worker, "solve", solve)
    monkeypatch.setattr(worker, "finish_states", lambda *a: pytest.fail("global reference on incomplete task"))
    if failure == "conversion":
        monkeypatch.setattr(worker, "state_values", fail_conversion)
    with pytest.raises((ValueError, RuntimeError), match="injected|complete requested prefix"):
        worker.run_task(QP, COMPLEX, BATH, tmp_path)
    paths = list(tmp_path.iterdir())
    if failure == "incomplete_prefix":
        assert not paths
    else:
        assert paths == [tmp_path/"qp-sz+0.json"]
        scalar = json.loads(paths[0].read_text())
        assert len(scalar["energies"]) == 4
        assert not {"coverage", "qualified", "status", "state_artifact"} & scalar.keys()
        with pytest.raises(ValueError, match="fresh directory"):
            worker.run_task(QP, COMPLEX, BATH, tmp_path)


@pytest.mark.parametrize("failure", ["primary_residual", "extra_residual", "lookahead_tie"])
def test_qp_qualification_and_lookahead_match_frozen_scope(tmp_path, monkeypatch, failure):
    original = worker.solve

    def altered(h, **kwargs):
        result = original(h, **kwargs)
        if kwargs["sector"].twice_sz == 0:
            if failure == "lookahead_tie":
                result.energies[3] = result.energies[2]
            else:
                result.residuals[0 if failure == "primary_residual" else 2] = 2e-8
        return result

    monkeypatch.setattr(worker, "solve", altered)
    report = worker.run_task(QP, COMPLEX, BATH, tmp_path)
    assert report["qualified"] == (failure != "primary_residual")
    assert report["diagnostics"]["all_returned_states_converged"] == (failure == "lookahead_tie")
    assert report["states"][2]["degenerate"] == (failure == "lookahead_tie")
    assert len(report["states"]) == 9


@pytest.fixture(scope="module")
def small_dmrg(tmp_path_factory):
    pytest.importorskip("tenpy")
    pytest.importorskip("h5py")
    base = tmp_path_factory.mktemp("continuum-three-sites")
    runs = {}
    original = worker.solve
    live = []

    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        live.append(result)
        return result

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(worker, "solve", capture)
        for spin, roots in ((0, 2), (-1, 1), (1, 1)):
            task = dict(DMRG, twice_sz=spin, roots=roots)
            target, old_target = base/f"new{spin}", base/f"old{spin}"
            target.mkdir()
            old_target.mkdir()
            before = deepcopy((task, COMPLEX, BATH))
            report = worker.run_task(task, COMPLEX, BATH, target)
            old = dmrg_states.run_state(task, COMPLEX, BATH, old_target)
            assert (task, COMPLEX, BATH) == before
            runs[spin] = (target, report, live[-1], old_target, old, task)
    return runs


@pytest.mark.dmrg
class TestDMRG:
    def test_frozen_equivalence_checkpoint_fields_and_parent_validation(self, small_dmrg, tmp_path):
        for spin, (directory, report, live, _, old, task) in small_dmrg.items():
            compare_rows(report["states"], old["states"])
            expected = run_finite(dict(backend="ed", spins=[spin], roots=task["roots"]), COMPLEX, BATH)
            compare_rows(report["states"], expected["states"])
            assert report["qualified"] == old["qualified"] is True
            assert set(report) == set(old)
            assert set(report["checkpoint"]) == set(old["checkpoint"])
            assert set(report["diagnostics"]) == set(old["diagnostics"]) | {"continuum_finite_version"}
            checkpoint, diagnostics = report["checkpoint"], report["diagnostics"]
            path = directory/checkpoint["path"]
            assert digest(path) == checkpoint["sha256"]
            assert digest(directory/checkpoint["state_artifact"]["file"]) == checkpoint["state_artifact"]["sha256"]
            loaded = load_result(path)
            assert worker.rows_from_result(loaded, COMPLEX) == report["states"]
            assert checkpoint["hamiltonian_sha256"] == live.hamiltonian.fingerprint()
            assert checkpoint["coordinate_id"] == loaded.prepared.coordinate_id
            assert loaded.metadata["test4_dmrg_states"]["continuum_finite_version"] == 1
            assert diagnostics["mode_order"] == [2, 3, 0, 1, 4, 5]
            assert diagnostics["bath_weight_sum_per_lead"] == .019
            for key in ("weights_renormalized", "bath_convergence_checked", "chi_convergence_checked",
                        "spectral_completeness_established", "lowest_root_completeness_established",
                        "full_qp_comparison_checked"):
                assert diagnostics[key] is False
            assert diagnostics["timings"]["total"] > 0 and diagnostics["checkpoint_seconds"] > 0
            assert all(info["trial_count"] == 2 for info in diagnostics["root_selection_seed_independence"])
            assert all(info["sweep_statistics"] for info in diagnostics["solver_metadata"]["roots"])
            assert all(not {"excitation", "subgap"} & row.keys() for row in report["states"])
            assert "coverage" not in report
            target = tmp_path/f"parent{spin}"
            shutil.copytree(directory, target)
            verify_output_artifacts(target, task, COMPLEX, BATH, report)
            json.dumps(report, allow_nan=False)

    @pytest.mark.parametrize("spin,legacy,roots,chi", [(0, True, 2, 256), (-1, False, 2, 512), (1, False, 1, 1024)])
    def test_same_bath_continuation_reoptimizes_intact_prefix(self, small_dmrg, tmp_path, monkeypatch,
                                                            spin, legacy, roots, chi):
        directory, report, _, old_directory, old, _ = small_dmrg[spin]
        if legacy:
            directory, report = old_directory, old
        before = {path.name: path.read_bytes() for path in directory.iterdir()}
        original = worker.solve

        def checked(h, **kwargs):
            assert len(kwargs["initial"].states) == len(report["states"])
            assert kwargs["options"]["chi_max"] == chi and kwargs["options"]["eigenpairs"] == roots
            return original(h, **kwargs)

        monkeypatch.setattr(worker, "solve", checked)
        new = worker.run_task(dict(DMRG, twice_sz=spin, roots=roots, chi=chi), COMPLEX, BATH, tmp_path,
                              initial_for(directory, report))
        assert {path.name: path.read_bytes() for path in directory.iterdir()} == before
        loaded = load_result(tmp_path/"eigenstates.json")
        assert len(loaded.states) == roots
        assert loaded.metadata["test4_dmrg_states"]["root_prefix"] == list(range(roots))
        for root, info in enumerate(loaded.metadata["roots"]):
            assert (info["seed_trials"][0]["seed"] == "continuation") == (root < len(report["states"]))
        assert_allclose(complex_array(loaded.metadata["gram"]), np.eye(roots), rtol=0, atol=1e-8)
        expected = run_finite(dict(backend="ed", spins=[spin], roots=roots), COMPLEX, BATH)
        compare_rows(new["states"], expected["states"])
        assert new["qualified"]

    @pytest.mark.parametrize("chi", [128, 256, 512, 1024])
    def test_approved_options_without_extra_numerics(self, small_dmrg, tmp_path, monkeypatch, chi):
        from qdjj_solver.dmrg_solver.solver import SolverOptions

        live = small_dmrg[-1][2]

        def fake(h, **kwargs):
            options = kwargs["options"]
            assert options == dict(eigenpairs=1, chi_max=chi, min_sweeps=8, max_sweeps=40,
                energy_tolerance=1e-12, residual_tolerance=1e-7, lanczos_probability_tolerance=1e-22,
                seed_trials=2, seed=1729, calculate_residuals=True, require_convergence=False, threads=1,
                group_size=2, mode_order=[2, 3, 0, 1, 4, 5])
            assert kwargs["initial"] is None and kwargs["cutoff"] is None
            result = deepcopy(live)
            result.metadata["options"] = asdict(SolverOptions(**options))
            return result

        monkeypatch.setattr(worker, "solve", fake)
        report = worker.run_task(dict(backend="dmrg", twice_sz=-1, roots=1, chi=chi), COMPLEX, BATH, tmp_path)
        assert report["checkpoint"]["chi"] == chi
        assert load_result(tmp_path/"eigenstates.json").metadata["options"]["chi_max"] == chi

    @pytest.mark.parametrize("kind", ["hash", "relative", "label", "extra", "physical", "bath", "sector", "slice"])
    def test_initial_rejections_before_load(self, small_dmrg, tmp_path, monkeypatch, kind):
        directory, report, _, _, _, _ = small_dmrg[0]
        initial = initial_for(directory, report)
        task, physical, bath = dict(DMRG, roots=2), COMPLEX, BATH
        if kind == "hash":
            initial["sha256"] = "0"*64
        elif kind == "relative":
            initial["path"] = "eigenstates.json"
        elif kind == "label":
            initial = "old-even"
        elif kind == "extra":
            initial["root"] = 1
        elif kind == "physical":
            physical = dict(COMPLEX, gamma=.41)
        elif kind == "bath":
            bath = dict(BATH, weights=[.02])
        elif kind == "sector":
            task["twice_sz"] = -1
        else:
            task["roots"] = 1
        monkeypatch.setattr(worker, "load_result", lambda *a: pytest.fail("loaded rejected initial"))
        monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solved rejected initial"))
        with pytest.raises(ValueError):
            worker.run_task(task, physical, bath, tmp_path, initial)
        assert list(tmp_path.iterdir()) == []

    @pytest.mark.parametrize("kind", ["format", "hamiltonian_hash", "hamiltonian", "observable", "prefix",
        "missing_prefix", "count", "root_count", "physical", "bath", "layout", "sector", "backend",
        "order", "option_order", "group", "groups", "coordinate", "cutoff"])
    def test_metadata_rejection_before_loading_mps(self, small_dmrg, tmp_path, monkeypatch, kind):
        directory, report, _, _, _, _ = small_dmrg[0]
        record = json.loads((directory/"eigenstates.json").read_text())
        metadata = record["calculation"]
        workflow = metadata["test4_dmrg_states"]
        if kind == "format":
            record["format_version"] = 2
        elif kind == "hamiltonian_hash":
            record["hamiltonian_sha256"] = "0"*64
        elif kind == "hamiltonian":
            record["hamiltonian"]["operator"][0]["real"] += 1.
        elif kind == "observable":
            record["hamiltonian"]["observables"]["phase_derivative"][0]["real"] += 1.
        elif kind == "prefix":
            workflow["root_prefix"] = [1, 2]
        elif kind == "missing_prefix":
            del metadata["test4_dmrg_states"]
        elif kind == "count":
            metadata["options"]["eigenpairs"] = 1
        elif kind == "root_count":
            metadata["roots"].pop()
        elif kind == "physical":
            workflow["physical"]["field"] *= -1
        elif kind == "bath":
            workflow["bath"]["weights"] = [.02]
        elif kind == "layout":
            workflow["layout"] = "chain"
        elif kind == "sector":
            metadata["sector"]["eta"] = 1
        elif kind == "backend":
            metadata["backend"] = "qp"
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
        else:
            metadata["qp_cutoff"] = 4
        path, target = tmp_path/"wrong.json", tmp_path/"target"
        write_json(path, record)
        target.mkdir()
        monkeypatch.setattr(worker, "load_result", lambda *a: pytest.fail("loaded inconsistent metadata"))
        monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solved inconsistent metadata"))
        with pytest.raises(ValueError):
            worker.run_task(dict(DMRG, roots=2), COMPLEX, BATH, target, dict(path=str(path), sha256=digest(path)))
        assert list(target.iterdir()) == []
        assert digest(directory/"eigenstates.json") == report["checkpoint"]["sha256"]

    @pytest.mark.parametrize("kind", ["mps_bytes", "mps_symlink", "json_symlink", "path", "binding", "load_race"])
    def test_json_mps_hash_binding_and_ordinary_files(self, small_dmrg, tmp_path, monkeypatch, kind):
        directory, report, _, _, _, _ = small_dmrg[0]
        source, target = tmp_path/"source", tmp_path/"target"
        shutil.copytree(directory, source)
        target.mkdir()
        path = source/"eigenstates.json"
        artifact = source/report["checkpoint"]["state_artifact"]["file"]
        if kind == "mps_bytes":
            artifact.write_bytes(artifact.read_bytes()+b"broken")
        elif kind == "mps_symlink":
            artifact.unlink()
            artifact.symlink_to(directory/artifact.name)
        elif kind == "json_symlink":
            path.unlink()
            path.symlink_to(directory/path.name)
        elif kind in ("path", "binding"):
            record = json.loads(path.read_text())
            if kind == "path":
                record["state_artifact"]["file"] = "../"+artifact.name
            else:
                record["observables"]["impurity_charge"][0] += .01
            write_json(path, record)
        else:
            original = worker.load_result

            def changed(path):
                result = original(path)
                path.write_text(path.read_text()+"\n")
                return result

            monkeypatch.setattr(worker, "load_result", changed)
        initial = dict(path=str(path), sha256=digest(path))
        monkeypatch.setattr(worker, "solve", lambda *a, **kw: pytest.fail("solved bad checkpoint"))
        with pytest.raises(ValueError, match="checksum|ordinary file|local to|do not match|changed"):
            worker.run_task(dict(DMRG, roots=2), COMPLEX, BATH, target, initial)
        assert list(target.iterdir()) == []

    @pytest.mark.parametrize("failure", ["residual", "sweeps", "norm", "overlap", "gram", "order", "degenerate", "no_residual"])
    def test_unqualified_prefix_saved_and_never_globalized(self, small_dmrg, tmp_path, monkeypatch, failure):
        result = deepcopy(small_dmrg[0][2])
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
        monkeypatch.setattr(worker, "solve", lambda *a, **kw: result)
        report = worker.run_task(dict(DMRG, roots=2), COMPLEX, BATH, tmp_path)
        assert not report["qualified"]
        assert all(not {"excitation", "subgap"} & row.keys() for row in report["states"])
        assert len(load_result(tmp_path/"eigenstates.json").states) == 2
        assert report["checkpoint"]["sha256"] == digest(tmp_path/"eigenstates.json")

    def test_checkpoint_precedes_conversion_and_parent_owns_failure(self, small_dmrg, tmp_path, monkeypatch):
        monkeypatch.setattr(worker, "solve", lambda *a, **kw: deepcopy(small_dmrg[0][2]))
        saved = {}

        def failed_conversion(*args):
            path = tmp_path/"eigenstates.json"
            saved["bytes"] = path.read_bytes()
            assert len(load_result(path).states) == 2
            raise ValueError("injected conversion failure")

        monkeypatch.setattr(worker, "rows_from_result", failed_conversion)
        with pytest.raises(ValueError, match="injected conversion"):
            worker.run_task(dict(DMRG, roots=2), COMPLEX, BATH, tmp_path)
        assert (tmp_path/"eigenstates.json").read_bytes() == saved["bytes"]
        assert not (tmp_path/"result.json").exists()
        companion = json.loads((tmp_path/"eigenstates.json").read_text())["state_artifact"]
        assert digest(tmp_path / companion["file"]) == companion["sha256"]
