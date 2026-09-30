"""Unequal reservoirs, exact coupling asymmetry and dynamical spectral oracles."""

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

from qdjj_solver import DiscreteBath, Impurity
from applications._common import bath_for, read_csv, read_json
from applications.zonda_2016_unequal_gaps import run as unequal
from applications.kadlecova_2017_asymmetry import run as asymmetry
from applications.hecht_2008_gap_edge import run as edge
from applications.hecht_2008_gap_edge import spectral
from electron_oracle import electron_hamiltonian, electron_operators, sector_indices

ROOT = Path(__file__).resolve().parents[1]


def archived_bath(module, delta=1., profile="quick", levels=2):
    manifest = read_json(module.CASE/"output"/profile/"manifest.json")
    return DiscreteBath.from_record(next(b for b in manifest["baths"].values()
                                        if b["delta"] == delta and len(b["xi"]) == levels))


def test_unequal_gaps_against_physical_electron_ed_and_current():
    baths = [DiscreteBath([-.4, .4], [.2, .2], d, 4.) for d in (1., .25)]
    parameters = {"u": 1., "gammas": [.5, .25]}
    detuning, phi, step = .3, 1.2, 1e-4
    row, states = unequal.point(baths, parameters, detuning, phi, {})
    contacts = [np.sqrt(g/(np.pi*b.rho))*np.eye(2)
                for g, b in zip(parameters["gammas"], baths, strict=True)]
    matrix = electron_hamiltonian(Impurity.anderson(1., detuning=detuning), baths,
                                   contacts, [-phi/2, phi/2])
    for p, state in enumerate(states):
        indices = sector_indices(10, p, p)
        expected = eigh(matrix[indices][:, indices].toarray(), eigvals_only=True)[0]
        assert_allclose(state.energies[0], expected, atol=2e-11)
    minus = unequal.point(baths, parameters, detuning, phi-step, {})[0]
    plus = unequal.point(baths, parameters, detuning, phi+step, {})[0]
    key = ("even_energy", "odd_energy")[row["ground_parity"]]
    assert_allclose((plus[key]-minus[key])/step, row["current"], atol=2e-8)
    assert abs(row["current_sum"]) < 2e-10
    reflected = unequal.point(baths, parameters, -detuning, phi, {})[0]
    assert_allclose(row["charge"]+reflected["charge"], 2., atol=2e-10)


def test_variable_gap_bath_replay_and_normal_measure():
    bath = archived_bath(unequal, .25)
    record = {"bath_record": bath.record()}
    assert bath_for(record, 100., .25).delta == .25
    with pytest.raises(ValueError, match="physical energy scales"):
        bath_for(record, 100.)
    grid = bath_for({"bath_kind": "cosh", "pairs": 3}, 7., .25)
    assert grid.delta == .25
    assert abs(sum(grid.weights)-1) > 1e-6  # No silent coarse-grid normalization.


def test_unequal_gap_published_current_checkpoint():
    parameters = read_json(unequal.CASE/"input"/"parameters.json")["model"]
    baths = [archived_bath(unequal, d, "paper", 4) for d in (1., .5)]
    row, _ = unequal.point(baths, parameters, 0., np.pi/2, {})
    ref = min((r for r in read_csv(unequal.CASE/"reference"/"figure8.csv") if r["method"] == "NRG"),
              key=lambda r: abs(float(r["detuning"])))
    # Separate ~4e-5 vector precision from ~1e-3 bath/continuum-method accuracy.
    assert_allclose(row["current"], float(ref["current"]), atol=.0015, rtol=0)


@pytest.mark.parametrize("a", [1., 4., 11., 1/11])
def test_asymmetry_energy_charge_and_current_identity(a):
    bath = archived_bath(asymmetry)
    rows, _, _ = asymmetry.compare_mapping(bath, 4., .8, .35, .73*np.pi, a, {})
    for row in rows:
        assert_allclose([row[k] for k in ("energy_error", "charge_error", "current_error")], 0., atol=2e-10)
    chi, phase, jacobian = asymmetry.mapping(.73*np.pi, a)
    assert_allclose(chi, abs((a*np.exp(-.365j*np.pi)+np.exp(.365j*np.pi))/(a+1))**2, atol=2e-15)
    assert_allclose(asymmetry.inverse_phase(phase, a), .73*np.pi, atol=2e-14)
    step = 1e-5
    derivative = (asymmetry.mapping(.73*np.pi+step, a)[1]-asymmetry.mapping(.73*np.pi-step, a)[1])/(2*step)
    assert_allclose(jacobian, derivative, atol=2e-10)


def test_asymmetry_endpoints_and_inaccessible_transition():
    assert asymmetry.inverse_phase(.8*np.pi, 11.) is None
    assert asymmetry.mapping(np.pi, 1.) == (0., np.pi, 1.)
    assert abs(asymmetry.mapping(np.pi, 11.)[2]) < 1e-15
    assert asymmetry.boundary(lambda x: {"signed_gap": 1+x}, 1e-6) is None
    with pytest.raises(ValueError, match="not bracketed"):
        asymmetry.boundary(lambda x: {"signed_gap": -1-x}, 1e-6)
    with pytest.raises(ValueError, match="asymmetry"):
        asymmetry.mapping(.5, 0.)


def test_asymmetry_published_gate_boundary_bracket():
    bath = archived_bath(asymmetry, profile="paper", levels=4)
    ref = min((r for r in read_csv(asymmetry.CASE/"reference"/"figure1.csv") if float(r["u_meV"]) == 3.2),
              key=lambda r: abs(float(r["phi_over_pi"])-.6))
    x, phi = float(ref["tilde_epsilon"]), float(ref["phi_over_pi"])*np.pi
    # Published gate-coordinate precision is .0018; the finite-bath/wide-band
    # discrepancy is larger (~.02), so test a declared physical bracket.
    a = asymmetry.point(bath, 3.2/.17, .44/.17, x-.03, phi, {})[0]["signed_gap"]
    b = asymmetry.point(bath, 3.2/.17, .44/.17, x+.03, phi, {})[0]["signed_gap"]
    assert a < 0 < b


@pytest.mark.parametrize("u,gamma", [(0., .7), (.8, .7), (3., .12)])
def test_spectral_lehmann_resolvent_against_independent_electron_ed(u, gamma):
    bath = DiscreteBath([-.4, .4], [.17, .17], 1., 4.)
    poles, states, _ = spectral.fock_measure(bath, u, gamma, .2, {}, 128)
    z = np.array([-.8, .2, 1.4, 3.1])+.12j
    computed = sum(p["weight"]/(z-p["omega"]) for p in poles)
    contacts = [np.sqrt(gamma/(np.pi*bath.rho))*np.eye(2)]
    matrix = electron_hamiltonian(Impurity.anderson(u, detuning=.2), [bath], contacts, [0.])
    energies, vectors = eigh(matrix.toarray())
    ground = np.flatnonzero(energies-energies[0] < 1e-10)
    expected = np.zeros_like(z)
    for j in ground:
        for op in electron_operators(6)[:2]:
            for sign, operator in ((1, op.getH()), (-1, op)):
                weights = abs(vectors.conj().T@(operator@vectors[:, j]))**2
                expected += np.sum(weights[:, None]/(z-sign*(energies-energies[0])[:, None]), axis=0)/(2*len(ground))
    assert_allclose(computed, expected, atol=5e-12)
    assert_allclose(sum(p["weight"] for p in poles), 1., atol=2e-14)
    assert min(p["weight"] for p in poles) >= 0
    assert_allclose(min(s.energies[0] for s in states), energies[0], atol=1e-12)
    if u == 0:
        assert_allclose(computed, spectral.quadratic_green(z, gamma, .2, bath=bath), atol=5e-12)


def test_gap_edge_quadrature_broadening_and_published_formula():
    bath = spectral.logarithmic_bath(1024, 1e4)
    assert_allclose(sum(bath.weights), 1., atol=5e-15)
    gamma = 80.
    offsets = np.geomspace(1e-6, .1, 25)
    z = 1+offsets+.04j*offsets
    discrete = spectral.quadratic_green(z, gamma, bath=bath)
    exact = spectral.quadratic_green(z, gamma, bandwidth=1e4)
    assert_allclose(discrete, exact, atol=.0001, rtol=.0002)
    ref = [r for r in read_csv(edge.CASE/"reference"/"spectra.csv") if r["method"] == "analytic"
           and float(r["omega_prime_over_D"]) < 1e-5]
    for row in ref[::20]:
        offset = float(row["omega_prime_over_D"])/1e-4
        epsilon = float(row["epsilon_over_D"])/1e-4
        value = spectral.analytic_continuum(offset, gamma, epsilon)*np.pi*gamma
        assert_allclose(value, float(row["pi_gamma_A"]), rtol=.003)
    peak = minimize_scalar(lambda x: -spectral.analytic_continuum(np.exp(x), gamma),
                           bounds=(-12, -3), method="bounded")
    assert_allclose(-peak.fun*np.pi*gamma, 20., rtol=.001)
    poles = [{"omega": 1.4, "weight": .3}, {"omega": .5, "weight": .7}]
    integral = quad(lambda x: float(spectral.broaden_continuum(np.exp(x), poles, .3))*np.exp(x), -20, 20)[0]
    assert_allclose(integral, .3, atol=1e-10)  # ABS weight is not smeared into the continuum.


@pytest.mark.parametrize("module", (unequal, asymmetry, edge), ids=lambda m: m.CASE.name)
def test_third_set_offline_runner(module, tmp_path):
    output = tmp_path/"résultats with spaces"
    env = os.environ | {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    result = subprocess.run([sys.executable, "-B", "-m", module.__name__, "--output", str(output)],
                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout+result.stderr
    manifest = read_json(output/"manifest.json")
    states = read_json(output/"eigenstates.json")
    assert manifest["state_records"] == len(states) > 0
    assert max(max(r["residuals"]) for r in states) < 2e-8
    assert manifest["baths"] and manifest["scripts_sha256"]


@pytest.mark.slow
@pytest.mark.convergence
def test_gap_edge_refined_quadrature_and_zero_broadening_limit():
    offsets = np.geomspace(1e-8, .1, 61)
    gamma, epsilon = 80., 320.
    bath = spectral.logarithmic_bath(8192, 1e4)
    z = 1+offsets+.01j*offsets
    computed = -spectral.quadratic_green(z, gamma, epsilon, bath=bath).imag
    finite = -spectral.quadratic_green(z, gamma, epsilon, bandwidth=1e4).imag
    limit = spectral.analytic_continuum(offsets, gamma, epsilon)*np.pi
    assert max(abs(computed-finite))/max(limit) < 5e-5
    assert max(abs(finite-limit))/max(limit) < .012
