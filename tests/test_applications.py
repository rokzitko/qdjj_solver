"""Offline scientific regressions with archived baths and independent identities."""

import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import DiscreteBath, Sector
from qdjj_solver.qp_solver import ProjectedOperator
from applications._common import eigenstate, parity_states, read_csv, read_json, scalar
from applications.zitko_2023_knight_shift import run as knight
from applications.zonda_2023_double_dot import run as double
from applications.zalom_2024_multiterminal import run as multi

ROOT = Path(__file__).resolve().parents[1]
CASES = (knight, double, multi)


def archived(module, *, direct=False, closure=False):
    config = read_json(module.CASE/"input"/"parameters.json")
    settings = config["direct_settings"] if direct else config["profiles"]["quick"]["settings"]
    if closure:
        settings = config["closure_settings"]
    manifest = read_json(module.CASE/"output"/"quick"/"manifest.json")
    records = [b for b in manifest["baths"].values() if len(b["xi"]) == settings["levels"]]
    assert len(records) == 1
    return config, settings, DiscreteBath.from_record(records[0])


def test_knight_archived_curve_and_zeeman_identity():
    config, settings, bath = archived(knight)
    u = config["model"]["u"]
    for row in read_csv(knight.CASE/"output"/"quick"/"phase.csv"):
        p = float(row["phi_over_pi"])*np.pi
        value, state = knight.knight(bath, u, .1*u, p, settings)
        assert_allclose(value, float(row["kappa"]), atol=2e-9, rtol=0)
        assert np.max(state.residuals) < 1e-8
    field = min(config["finite_fields"])
    phi = np.pi*config["finite_field_phi_over_pi"]
    zero, _ = knight.knight(bath, u, .1*u, phi, settings)
    up = knight.knight(bath, u, .1*u, phi, settings, field, 1)[1]
    down = knight.knight(bath, u, .1*u, phi, settings, field, -1)[1]
    assert_allclose(1-(up.energies[0]-down.energies[0])/field, zero, atol=2e-8, rtol=0)


def test_knight_perturbative_orders():
    config, settings, bath = archived(knight)
    u = config["model"]["u"]
    # Independent second-order perturbation expression for this exact finite bath.
    slope = np.sum(bath.weights/(np.pi*bath.rho)/(u/2+bath.energies)**2)
    gamma = 1e-5
    value, _ = knight.knight(bath, u, gamma, .7, settings)
    assert_allclose(value/gamma, slope, rtol=2e-5)
    # Fourth order in hopping means second order in Gamma for the phase difference.
    differences = []
    for gamma in (.01, .02):
        values = [knight.knight(bath, u, gamma, p, settings)[0] for p in (0., np.pi)]
        differences.append(values[1]-values[0])
    assert differences[0] > 0
    assert_allclose(differences[1]/differences[0], 4., rtol=.02)
    assert_allclose(knight.leading_coefficient(10., 100.), .11415375251228316, atol=1e-12)


@pytest.mark.parametrize("phi", [0., .8, np.pi])
def test_single_dot_compression_preserves_bound_doublet(phi):
    config, settings, bath = archived(multi, direct=True)
    settings = settings | {"solver": {"method": "dense"}}
    models = [multi.effective_model(bath, 3., 1., np.cos(phi/2), compress) for compress in (False, True)]
    results = [[eigenstate(h, Sector(p, p), settings, roots=2 if phi == np.pi and p == 0 else 1)
                for p in (0, 1)] for h in models]
    for p in (0, 1):
        assert_allclose(results[0][p].energies, results[1][p].energies, atol=1e-10, rtol=0)
        if phi == np.pi and p == 0:
            # At the singlet crossing compare the operator on the whole
            # degenerate space; individual eigensolver eigenvectors are arbitrary.
            eigenvalues = []
            for h, states in zip(models, results, strict=True):
                state = states[p]
                operator = ProjectedOperator(h.observables["phase_derivative"], state.basis)
                matrix = np.column_stack([state.vectors.conj().T @ operator.action(v) for v in state.vectors.T])
                eigenvalues.append(np.linalg.eigvalsh(matrix))
            assert_allclose(*eigenvalues, atol=1e-10, rtol=0)
        else:
            assert_allclose(results[0][p].observables["phase_derivative"],
                            results[1][p].observables["phase_derivative"], atol=1e-10, rtol=0)


@pytest.mark.parametrize("label", ["weak", "strong"])
def test_double_dot_current_baseline_and_derivative(label):
    config, settings, bath = archived(double)
    curve = next(c for c in config["curves"] if c["label"] == label)
    rows = [r for r in read_csv(double.CASE/"output"/"quick"/"current.csv") if r["curve"] == label]
    for row in rows:
        result, _ = double.point(double.model(bath, 4., curve, float(row["phi_over_pi"])*np.pi), settings)
        for key in ("signed_gap", "even_current", "odd_current"):
            assert_allclose(result[key], float(row[key]), atol=2e-9, rtol=0)
        assert result["ground_parity"] == int(row["ground_parity"])
    phi, step = np.pi/2, config["finite_difference_step"]
    center = double.point(double.model(bath, 4., curve, phi), settings)[0]
    minus = double.point(double.model(bath, 4., curve, phi-step), settings)[0]
    plus = double.point(double.model(bath, 4., curve, phi+step), settings)[0]
    for parity in ("even", "odd"):
        current = (plus[f"{parity}_energy"]-minus[f"{parity}_energy"])/step
        assert_allclose(current, center[f"{parity}_current"], atol=2e-8, rtol=0)


def test_double_dot_gal_analytical_phase_boundaries():
    # Independent Eq. (24), not a diagonalization-generated reference.
    u, gamma = 4., 1.4
    boundaries = (gamma-u/(2*(1+gamma)),
                  gamma/3+np.sqrt((2*gamma/3)**2-u**2/(12*(1+gamma)**2)))
    for hopping in boundaries:
        curve = dict(gamma_left=gamma, gamma_right=gamma, hopping=hopping)
        result, _ = double.point(double.gal_model(u, curve, np.pi), {})
        assert abs(result["signed_gap"]) < 1e-10


def test_multiterminal_mapping_currents_and_gauge():
    config, settings, bath = archived(multi, direct=True)
    for pair in config["direct_phase_pairs_over_pi"]:
        phases = np.r_[0., np.pi*np.array(pair)]
        rows, states, _ = multi.compare_mapping(bath, config["model"], config["relative_couplings"], phases, settings)
        for row in rows:
            assert abs(row["energy_error"]) < 1e-9
            assert row["current_mapping_error"] < 1e-9
            assert abs(row["current_sum"]) < 1e-9
        shifted = parity_states(multi.direct_model(bath, 3., 1., config["relative_couplings"], phases+.37), settings)
        for state, other in zip(states, shifted, strict=True):
            assert_allclose(state.energies, other.energies, atol=1e-9, rtol=0)
            for j in range(3):
                assert_allclose(scalar(state, f"lead_current_{j}"),
                                scalar(other, f"lead_current_{j}"), atol=1e-9, rtol=0)


def test_multiterminal_phase_crossing_and_excited_degeneracy():
    config, settings, bath = archived(multi)
    reference = read_json(multi.CASE/"output"/"quick"/"summary.json")
    critical = reference["critical"]["chi_critical"]
    signs = []
    for chi in (critical-.01, critical, critical+.01):
        states = parity_states(multi.effective_model(bath, 3., 1., chi), settings)
        signs.append(states[1].energies[0]-states[0].energies[0])
    assert signs[0] < 0 < signs[2]
    assert abs(signs[1]) < 2e-7
    _, direct_settings, direct_bath = archived(multi, closure=True)
    theta = np.arccos(-3/7)
    h = multi.direct_model(direct_bath, 3., 1., config["relative_couplings"], [0., theta, -theta])
    even = eigenstate(h, Sector(0, 0), direct_settings, roots=2)
    odd = eigenstate(h, Sector(1, 1), direct_settings)
    assert_allclose(even.energies[0], even.energies[1], atol=1e-10, rtol=0)
    assert 0 < even.energies[0]-odd.energies[0] < 1


@pytest.mark.parametrize("module", CASES, ids=[m.CASE.name for m in CASES])
def test_application_runner_offline_outputs(module, tmp_path):
    output = tmp_path/"calculation with spaces ž"
    env = os.environ | {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                       "VECLIB_MAXIMUM_THREADS": "1"}
    command = [sys.executable, "-B", "-m", module.__name__, "--profile", "quick", "--output", str(output)]
    if module is knight:
        command.append("--raw")
    completed = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True, timeout=90)
    assert completed.returncode == 0, completed.stdout+completed.stderr
    manifest = read_json(output/"manifest.json")
    assert manifest["case"] == module.CASE.name
    assert manifest["state_records"] > 0
    states = read_json(output/"eigenstates.json")
    assert len(states) == manifest["state_records"]
    assert max(max(s["residuals"]) for s in states) < 1e-8
    assert manifest["baths"] and manifest["implementation"]["sha256"]
    if module is knight:
        record = read_json(output/"raw"/"00001.json")
        assert record["format"] == "qdjj-eigenstates" and record["backend"] == "qp"
        assert "hamiltonian" in record


@pytest.mark.slow
@pytest.mark.convergence
@pytest.mark.dmrg
def test_double_dot_independent_dmrg_checkpoint():
    pytest.importorskip("tenpy")
    config, settings, bath = archived(double)
    h = double.model(bath, 4., config["curves"][1], np.pi/2)
    qp = parity_states(h, settings)
    dmrg = parity_states(h, dict(backend="dmrg", solver=dict(chi_max=64, seed_trials=1,
                                                            residual_tolerance=1e-7)))
    for exact, mps in zip(qp, dmrg, strict=True):
        assert_allclose(exact.energies, mps.energies, atol=1e-8, rtol=0)
        assert_allclose(exact.observables["phase_derivative"],
                        mps.observables["phase_derivative"], atol=1e-7, rtol=0)


def test_reference_extractors_parse_independent_formats():
    extractor = importlib.import_module("applications.zitko_2023_knight_shift.extract_reference")
    sets = extractor.grace_sets("# comment\n@target G0.S0\n@type xy\n0 1\n1 2\n&\n@type xy\n0 3\n&\n")
    assert_allclose(sets[0], [[0, 1], [1, 2]])
    assert_allclose(sets[1], [[0, 3]])
    # The committed NRG table retains original SVG coordinates and explicit uncertainty.
    rows = read_csv(double.CASE/"reference"/"figure9a.csv")
    for row in rows:
        assert_allclose(float(row["current"]), -.1+(float(row["svg_y"])-5188)/1666.5*.5, atol=1e-12)
    provenance = json.loads((double.CASE/"reference"/"provenance.json").read_text())
    assert provenance["current_extraction_uncertainty"] > 0
