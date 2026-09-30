"""Independent physics checks for the documentation benchmark, not timing gates."""

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_limits

from qdjj_solver import DiscreteBath, Sector, cosh_grid, reference_model, solve
from tools import benchmark_baths as bench
from tools.plot_bath_benchmarks import compare_interacting
from tools.refine_bath_chain import chain_hamiltonian


@pytest.mark.parametrize("bath", [
    cosh_grid(1, delta=.7, bandwidth=5.),
    DiscreteBath([-2., 0., 2.], [.2, .1, .2], delta=.7, bandwidth=5.),
])
def test_quadratic_integrals_and_poles_against_physical_bdg_and_fock_ed(bath):
    """Nonunit contact measures, detuning, phase and finite bandwidth all matter."""
    parameters = dict(gamma=.63, phi=1.1, detuning=.37)
    integral = bench.quadratic_observables(bandwidth=bath.bandwidth, delta=bath.delta,
                                          bath=bath, **parameters)
    with threadpool_limits(limits=1):
        bdg = bench.quadratic_bdg(bath, **parameters)
    for name in ("excitation", "current", "charge"):
        assert_allclose(integral[name], bdg[name], rtol=0, atol=2e-11)
    # Keep Fock ED bounded: the two-level junction has 10 electron modes.
    if bath.levels == 2:
        h = reference_model(bath, u=0., compress=False, **parameters)
        states = [solve(h, cutoff=None, sector=Sector(p, p)) for p in (0, 1)]
        assert_allclose(states[0].energies[0], bdg["energy"], atol=3e-12, rtol=0)
        assert_allclose(states[1].energies[0]-states[0].energies[0], bdg["excitation"], atol=3e-12, rtol=0)
        assert_allclose(states[0].observables["phase_derivative"][0], bdg["current"], atol=3e-12, rtol=0)
        assert_allclose(states[0].observables["impurity_charge"][0], bdg["charge"], atol=3e-12, rtol=0)
    step = 1e-5
    with threadpool_limits(limits=1):
        plus = bench.quadratic_bdg(bath, **dict(parameters, phi=parameters["phi"]+step))
        minus = bench.quadratic_bdg(bath, **dict(parameters, phi=parameters["phi"]-step))
    assert_allclose((plus["energy"]-minus["energy"])/(2*step), bdg["current"], atol=2e-9, rtol=0)


@pytest.mark.parametrize("bandwidth", [5., 100., 2000.])
def test_finite_band_continuum_oracle_against_refined_independent_bdg(bandwidth):
    parameters = dict(gamma=2., phi=1.7, detuning=.6)
    reference = bench.quadratic_observables(bandwidth=bandwidth, **parameters)
    # No Fock/MPS or scalar discrete-g evaluation is used by the BdG reference.
    with threadpool_limits(limits=1):
        finite = bench.quadratic_bdg(cosh_grid(48, bandwidth=bandwidth), **parameters)
    for name in ("excitation", "current", "charge"):
        assert_allclose(finite[name], reference[name], atol=2e-10, rtol=0)
    assert reference["current_quadrature_error"] < 2e-11
    assert reference["charge_quadrature_error"] < 2e-11


def test_direct_gauss_measure_replay_and_weak_coupling_reference():
    bath = bench.construct_bath(bench.bath_spec("linear-gl", 4, 10.))
    for power in (0, 2, 4, 6):
        assert_allclose(np.dot(bath.weights, (bath.xi/10)**power), 1/(power+1), atol=2e-15, rtol=0)
    replay = DiscreteBath.from_record(bench.bath_record(bath))
    assert_allclose(replay.weights, bath.weights, atol=0, rtol=0)
    assert replay.paired
    refined = bench.kernel_metrics(cosh_grid(32, bandwidth=100.))
    assert refined["weak_spin_slope_error"] < 2e-13
    # Absolute objective fitting does not imply a normalized bath measure.
    coarse = bench.kernel_metrics(cosh_grid(1, bandwidth=100.))
    assert abs(coarse["moment_0_error"]) > .1


@pytest.mark.parametrize("bath,u,gamma", [
    (cosh_grid(1, bandwidth=10.), 15., 20.),
    (DiscreteBath([-2., 0., 2.], [.3, .05, .3], delta=.7, bandwidth=5.), 3., .6),
])
def test_independent_chain_reference_matches_unrestricted_star(bath, u, gamma):
    parameters = dict(u=u, gamma=gamma, geometry="single", compress=False)
    chain = chain_hamiltonian(bath, **parameters)
    star = reference_model(bath, **parameters)
    for parity in (0, 1):
        left = solve(chain, cutoff=None, sector=Sector(parity, parity))
        right = solve(star, cutoff=None, sector=Sector(parity, parity))
        assert_allclose(left.energies, right.energies, rtol=0, atol=2e-11)
        for name in ("impurity_spin_z", "impurity_charge", "double_occupancy_0"):
            assert_allclose(left.observables[name], right.observables[name], rtol=0, atol=2e-11)


def test_invalid_quadratic_reference_and_timing_inputs():
    with pytest.raises(ValueError, match="crossing"):
        bench.quadratic_observables(bandwidth=10., gamma=.3, phi=np.pi)
    with pytest.raises(ValueError, match="symmetric"):
        bench.quadratic_observables(bandwidth=10., gamma=.3, phi=.8,
                                   bath=DiscreteBath([1.], [1.], bandwidth=10.))
    with pytest.raises(ValueError, match="even"):
        bench.construct_bath(bench.bath_spec("cosh-grid", 3, 10.))
    with pytest.raises(ValueError, match="timing"):
        bench.timing_summary([float("nan")])


def test_report_excludes_unconverged_mps_and_truncated_reference_candidates():
    """A larger declared bath must not masquerade as a better finite reference."""
    specs = {str(levels): bench.bath_spec("cosh-grid", levels, 100.) for levels in (4, 8, 16)}

    def record(levels, value, solver, converged=True):
        return dict(status="ok", converged=converged,
                    task=dict(type="solve", study="knight", case="point", bath_key=str(levels), solver=solver),
                    values={"kappa": value}, states=[dict(residual=1e-12)], peak_process_rss_bytes=100,
                    timings_seconds={"reused_bath_total": bench.timing_summary([1., 2., 3.])})

    records = [record(4, .13, {"backend": "qp", "cutoff": None}),
               record(8, .14, {"backend": "qp", "cutoff": None}),
               record(16, .99, {"backend": "qp", "cutoff": 2}),
               record(16, .91, {"backend": "dmrg", "options": {"chi_max": 128}}, converged=False)]
    rows, references = compare_interacting(records, specs, {"fitting": {}})
    assert len(references) == 1 and len(rows) == 2
    assert references[0]["value"] == .14
    assert references[0]["label"] == "largest finite ED bath"
    assert references[0]["stability_indicators"]["bond_step"] is None
    assert_allclose(rows[0]["difference"], .01, atol=1e-16, rtol=0)


def test_benchmark_smoke_resume_and_portable_export(tmp_path):
    output, archive = tmp_path/"run", tmp_path/"archive"
    env = dict(os.environ, **{k: "1" for k in bench.THREAD_VARIABLES})
    command = [sys.executable, "-B", str(Path(bench.__file__)), "--profile", "smoke", "--output", str(output)]
    subprocess.run(command, env=env, check=True, capture_output=True, text=True)
    manifest = bench.read_json(output/"manifest.json")
    before = {p.name: p.read_bytes() for p in (output/"jobs").glob("*.json")}
    subprocess.run(command, env=env, check=True, capture_output=True, text=True)
    assert bench.read_json(output/"manifest.json")["completed"] == manifest["completed"]
    assert before == {p.name: p.read_bytes() for p in (output/"jobs").glob("*.json")}
    bench.export_results(output, archive)
    assert bench.read_json(archive/"failures.json") == []
    records = bench.read_json(archive/"measurements.json")
    baths = bench.read_json(archive/"baths.json")
    assert all(r["task"]["bath_key"] in baths for r in records)
    assert len([r for r in records if r["task"]["type"] == "solve"]) == 2
    for record in records:
        assert record["status"] == "ok"
        if record["task"]["type"] == "solve":
            assert record["converged"]
            assert record["states"][0]["residual"] < 2e-8
    # A changed input must never silently reuse timings/observables from another run.
    changed = dict(manifest["input"], solve_repetitions=99)
    config = tmp_path/"changed.json"
    config.write_text(json.dumps(changed), encoding="utf-8")
    failed = subprocess.run(command+["--input", str(config)], env=env, capture_output=True, text=True)
    assert failed.returncode != 0 and "resume requires identical" in failed.stderr


def test_reviewed_evidence_reconstructs_recorded_task_checksums():
    """Detect mismatched bath coefficients, dropped rows, or broken archive joins."""
    root = bench.ROOT/"docs/benchmarks/baths"
    for directory in (root, root/"refinement", root/"chain", root/"controls", root/"edge"):
        manifest = bench.read_json(directory/"manifest.json")
        records = bench.read_json(directory/"measurements.json")
        baths = bench.read_json(directory/"baths.json")
        assert len(records) == len(manifest["completed"])
        for key, record in zip(manifest["completed"], records, strict=True):
            task = dict(record["task"])
            if task["type"] == "bath":
                assert bench.task_id(task["spec"]) == task["bath_key"]
            else:
                task["bath"] = baths[task["bath_key"]]
            assert bench.task_id(task) == key
