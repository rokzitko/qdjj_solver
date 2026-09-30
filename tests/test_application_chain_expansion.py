"""ChE/GAL literature regressions and independent multiorbital spin checks."""

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import eigh

from qdjj_solver import Impurity, Sector, chain_expansion, number
from qdjj_solver.common.baths import discrete_g, hybridization_g
from applications._common import eigenstate, parity_states, read_csv, read_json, scalar
from applications.bobok_2025_chain_expansion import models, run
from applications.bobok_2025_chain_expansion.chain_check import chain_model
from electron_oracle import electron_hamiltonian, electron_operators, evaluate

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("u", [1., 3., 10.])
def test_shared_lead_spin_and_observables_against_electron_ed(u):
    bath = chain_expansion(2, wide_band=True)
    row, _, _ = run.shared_point(bath, u, {"roots": 6})
    charges = [number(2*j)+number(2*j+1) for j in (0, 1)]
    impurity = Impurity(2, sum(u/2*(n-1)*(n-1) for n in charges))
    contacts = [np.sqrt(.5/(np.pi*bath.rho))*np.tile(np.eye(2), (1, 2))]
    physical = electron_hamiltonian(impurity, [bath], contacts, [0.])
    energies, vectors = eigh(physical.toarray())
    cs = electron_operators(8)
    sz = [(cs[2*j].getH()@cs[2*j]-cs[2*j+1].getH()@cs[2*j+1])/2 for j in range(4)]
    raising = [cs[2*j].getH()@cs[2*j+1] for j in range(4)]
    total_z, total_plus = sum(sz), sum(raising)
    spin_square = total_z@total_z+(total_plus@total_plus.getH()+total_plus.getH()@total_plus)/2
    spin = np.real(np.sum(vectors.conj()*(spin_square@vectors), axis=0))
    for s, name in ((0., "singlet"), (.5, "doublet"), (1., "triplet")):
        expected = min(energies[abs(spin-s*(s+1)) < 1e-8])
        assert_allclose(row[f"{name}_energy"], expected, atol=3e-12)
    ground = vectors[:, energies-energies[0] < 1e-9]
    observables = dict(spin_correlation=sz[0]@sz[1]+(raising[0]@raising[1].getH()+raising[0].getH()@raising[1])/2,
                       pairing=sum(cs[2*j].getH()@cs[2*j+1].getH()+cs[2*j+1]@cs[2*j] for j in (0, 1))/4)
    for key, op in observables.items():
        expected = np.trace(ground.conj().T@(op@ground)).real/ground.shape[1]
        assert_allclose(row[key], expected, atol=4e-12)
    # Full finite-space chain/star spectra, including high-energy states, match.
    chain = evaluate(chain_model(bath, u, .5).operator, cs).toarray()
    assert_allclose(eigh(chain, eigvals_only=True), energies, atol=5e-12)


def test_triplet_is_a_same_parity_transition_and_zbw_misses_it():
    rows = [run.shared_point(chain_expansion(l, wide_band=True), 10., {"roots": 4})[0] for l in (1, 2)]
    assert rows[0]["ground_spin"] == 0. and rows[1]["ground_spin"] == 1.
    assert rows[0]["triplet_singlet_gap"] > 0 > rows[1]["triplet_singlet_gap"]
    assert all(r["doublet_singlet_gap"] > 0 for r in rows)
    assert max(r["triplet_projection_error"] for r in rows) < 1e-10


def test_spin_resolution_at_exact_degeneracy():
    # U=0 has a free dark impurity orbital and an exactly degenerate S/T
    # excited subspace. Dense ED retains the entire sector before rotation.
    h = models.double_dot(chain_expansion(1, wide_band=True), 0., .5)
    state = eigenstate(h, Sector(0, 0), {"solver": {"method": "dense"}}, roots=20)
    spin, error = run.resolve_spin(state)
    assert error < 2e-12
    assert set(spin) >= {0., 1.}
    assert max(state.residuals) < 1e-11
    assert_allclose(state.vectors.T.conj()@state.vectors, np.eye(20), atol=5e-14)


@pytest.mark.parametrize("cutoff", [2, None])
def test_exchange_symmetry_keeps_both_interacting_impurity_orbitals(cutoff):
    bath = chain_expansion(2, wide_band=True)
    direct, _, _ = run.shared_point(bath, 5., {"cutoff": cutoff})
    reduced, _, states = run.shared_point(bath, 5., {"cutoff": cutoff, "exchange_basis": True})
    for key in ("singlet_energy", "doublet_energy", "triplet_energy", "pairing", "spin_correlation"):
        assert_allclose(direct[key], reduced[key], atol=2e-11)
    assert {s.basis.sector.eta for s in states} == {-1, 1}
    h = states[0].hamiltonian
    original_double = h.physical_operator(number(0)*number(1))
    assert max((abs(c) for c in (original_double-h.observables["double_occupancy_0"]).terms.values()), default=0.) < 1e-13


def test_gal_and_che_currents_against_published_vectors():
    references = read_csv(run.CASE/"reference"/"figure8.csv")
    for u in (2., 4., 8.):
        row, states = run.current_point(chain_expansion(4), u, .5*np.pi, {})
        ref = next(r for r in references if r["method"] == "NRG" and float(r["u"]) == u
                   and float(r["phi_over_pi"]) == .5)
        assert_allclose(row["current"], float(ref["current"]), atol=.003, rtol=0)
        assert max(s.residuals[0] for s in states) < 2e-8
        gal, _ = run.current_point(None, u, .9*np.pi, {})
        assert gal["current"] == 0.  # Dot-only GAL has no doublet current.
    for u in (2., 4., 8.):
        gal, _ = run.current_point(None, u, .3*np.pi, {})
        assert_allclose(gal["even_current"], 2/3*np.sin(.15*np.pi), atol=1e-14)
    phi, step = .37*np.pi, 1e-4
    bath = chain_expansion(2)
    center, _ = run.current_point(bath, 4., phi, {})
    minus, _ = run.current_point(bath, 4., phi-step, {})
    plus, _ = run.current_point(bath, 4., phi+step, {})
    for branch in ("even", "odd"):
        fd = (plus[f"{branch}_energy"]-minus[f"{branch}_energy"])/step
        assert_allclose(fd, center[f"{branch}_current"], atol=3e-8)


def test_finite_band_tail_and_pade_order_are_distinct():
    bath = chain_expansion(4, bandwidth=10.)
    errors = [abs(discrete_g(bath, x)-hybridization_g(x, bandwidth=10.)) for x in (.2, .1)]
    assert 200 < errors[0]/errors[1] < 300  # First omitted term is omega^8.
    assert abs(bath.weights.sum()-1) > .1  # Low-frequency matching does not fix the tail.


def test_chain_expansion_application_offline(tmp_path):
    output = tmp_path/"résultats with spaces"
    env = os.environ | {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    result = subprocess.run([sys.executable, "-B", "-m", run.__name__, "--output", str(output)],
                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout+result.stderr
    summary = read_json(output/"summary.json")
    assert summary["current_points"] == 18 and summary["shared_points"] == 6
    assert summary["maximum_current_derivative_error"] < 2e-8
    assert summary["maximum_spin_quantization_error"] < 1e-8
    manifest = read_json(output/"manifest.json")
    assert manifest["baths"] and manifest["state_records"] > 0


@pytest.mark.slow
@pytest.mark.dmrg
def test_che_chain_dmrg_and_star_qp_agree():
    pytest.importorskip("tenpy")
    bath = chain_expansion(4, wide_band=True)
    star, chain = models.double_dot(bath, 7., .5), chain_model(bath, 7., .5)
    settings = read_json(run.CASE/"input"/"dmrg_check.json")["dmrg"]
    for sector in (Sector(0, 0), Sector(0, 2), Sector(1, 1)):
        exact, mps = eigenstate(star, sector, {}), eigenstate(chain, sector, settings)
        assert_allclose(exact.energies, mps.energies, atol=1e-9, rtol=0)
        for key in ("pairing", "spin_correlation", "total_spin_squared"):
            assert_allclose(scalar(exact, key), scalar(mps, key), atol=2e-7)


def test_uncompressed_junction_spectrum_agrees_for_bound_ground_states():
    bath = chain_expansion(2)
    full = parity_states(models.junction(bath, 4., 1., .6*np.pi, compress=False), {})
    reduced = parity_states(models.junction(bath, 4., 1., .6*np.pi), {})
    for a, b in zip(full, reduced, strict=True):
        assert_allclose(a.energies, b.energies, atol=2e-12)
        assert_allclose(scalar(a, "phase_derivative"), scalar(b, "phase_derivative"), atol=2e-12)
