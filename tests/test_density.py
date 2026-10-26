"""Independent partial traces, analytic coherences, and backend/basis checks."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from qdjj_solver import (DiscreteBath, Hamiltonian, Impurity, Sector, annihilate, create,
                         hermitian_pair, make_model, number, orbital_density_matrix,
                         orbital_density_operators, orbital_density_properties, solve)


def local_pair(phase):
    return Hamiltonian(-hermitian_pair(np.exp(1j*phase)*create(0)*create(1)), 2, (1, -1))


@pytest.mark.parametrize("phase", [0., np.pi, np.pi/2, .37])
def test_analytic_complex_pair_and_incoherent_mixture(phase):
    result = solve(local_pair(phase), sector=Sector(0, None), options={"method": "dense", "eigenpairs": 2})
    rho = orbital_density_matrix(result, physical=False)
    ket = np.array([1, 0, 0, np.exp(1j*phase)])/np.sqrt(2)
    assert_allclose(rho, np.outer(ket, ket.conj()), atol=1e-14)
    props = orbital_density_properties(rho)
    assert_allclose([props["occupation"], props["empty"], props["single"], props["double"]], [1, .5, 0, .5])
    assert_allclose(props["pair_amplitude"], .5*np.exp(1j*phase), atol=1e-14)
    assert_allclose([props["plus"], props["minus"]], [(1+np.cos(phase))/2, (1-np.cos(phase))/2], atol=1e-14)
    assert_allclose([props["phase_optimized_plus"], props["phase_optimized_minus"], props["purity"]], [1, 0, 1], atol=1e-14)
    all_rho = orbital_density_matrix(result, state=None, physical=False)
    assert all_rho.shape == (2, 4, 4)
    mixed = orbital_density_matrix(result, state=None, weights=[2, 2], physical=False)
    assert_allclose(mixed, np.diag([.5, 0, 0, .5]), atol=1e-14)
    assert_allclose(mixed, all_rho.mean(axis=0), atol=1e-14)
    assert orbital_density_properties(mixed)["phase_optimized_plus"] == pytest.approx(.5)


def test_spin_coherence_and_singly_occupied_state():
    h = Hamiltonian(hermitian_pair(.5j*create(0)*annihilate(1)), 2, (1, -1))
    result = solve(h, sector=Sector(1, None))
    rho = orbital_density_matrix(result, physical=False)
    ket = np.array([0, 1, 1j, 0])/np.sqrt(2)
    assert_allclose(rho, np.outer(ket, ket.conj()), atol=1e-14)
    props = orbital_density_properties(rho)
    assert_allclose([props["occupation"], props["single"], props["even_weight"], props["plus"], props["minus"]],
                    [1, 1, 0, 0, 0], atol=1e-14)


def test_two_orbitals_against_independent_bit_basis_partial_trace():
    rng = np.random.default_rng(740)
    one_body = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    one_body += one_body.conj().T
    impurity = Impurity.from_integrals(one_body)
    op = impurity.operator + .7*number(0)*number(1)+.4*number(2)*number(3)
    op += hermitian_pair((.3+.6j)*create(0)*create(1)+(.6-.2j)*create(2)*create(3))
    h = Hamiltonian(op, 4, (1, -1, 1, -1))
    for parity in (0, 1):
        result = solve(h, sector=Sector(parity, None), options={"method": "dense", "eigenpairs": 8})
        occupations = result.basis.occupations()[:, 0]
        for orbital in (0, 1):
            actual = orbital_density_matrix(result, orbital, state=None, physical=False)
            for index, vector in enumerate(result.vectors.T):
                amplitudes = np.zeros((4, 4), complex)
                for state, coefficient in zip(occupations, vector, strict=True):
                    state = int(state)
                    local = (state >> (2*orbital)) & 3
                    rest = (state >> (2*(1-orbital))) & 3
                    amplitudes[local, rest] = coefficient
                expected = amplitudes @ amplitudes.conj().T
                assert_allclose(actual[index], expected, atol=3e-14, rtol=0)
                assert np.linalg.eigvalsh(actual[index]).min() >= -1e-14


@pytest.mark.parametrize("cutoff", [0, 2, None])
def test_electron_eta_and_compressed_coordinates(cutoff):
    bath = DiscreteBath([0.], [1.], 1., 2.)
    impurity = Impurity.anderson(u=.5)
    params = dict(impurity=impurity, baths=[bath, bath], tunneling=[.3*np.eye(2)]*2, phases=[.2, .8])
    options = {"eigenpairs": 2, "method": "dense"}
    physical = solve(make_model(**params), cutoff=cutoff, sector=Sector(0, None), options=options)
    # At cutoff zero the even impurity pair is degenerate. Compare its complete
    # equally weighted eigenspace rather than coordinate-dependent eigenvectors.
    weights = [1, 1] if cutoff == 0 else [1, 0]
    expected = orbital_density_matrix(physical, state=None, weights=weights)
    eta = solve(make_model(**params, eta_basis=True), cutoff=cutoff, sector=Sector(0, None), options=options)
    assert_allclose(orbital_density_matrix(eta, state=None, weights=weights), expected, atol=2e-13, rtol=0)
    compressed = make_model(**params, eta_basis=True, compress_pairs=True)
    compressed_cutoff = None if cutoff is None else min(cutoff, compressed.bath_modes)
    compressed_result = solve(compressed, cutoff=compressed_cutoff, sector=Sector(0, None), options=options)
    assert_allclose(orbital_density_matrix(compressed_result, state=None, weights=weights), expected, atol=2e-13, rtol=0)


@pytest.mark.dmrg
def test_dmrg_complex_pair_density_matches_qp():
    pytest.importorskip("tenpy")
    h = local_pair(.71)
    sector = Sector(0, None)
    qp = solve(h, sector=sector, options={"eigenpairs": 2})
    dmrg = solve(h, sector=sector, backend="dmrg", options={"eigenpairs": 2, "require_convergence": True})
    assert_allclose(orbital_density_matrix(dmrg, state=None, physical=False),
                    orbital_density_matrix(qp, state=None, physical=False), atol=1e-12)
    assert_allclose(orbital_density_matrix(dmrg, state=None, weights=[1, 1], physical=False),
                    np.diag([.5, 0, 0, .5]), atol=1e-12)


@pytest.mark.parametrize("orbital", [-1, 1.2, True, np.bool_(False)])
def test_invalid_orbitals(orbital):
    with pytest.raises(ValueError, match="nonnegative integer"):
        orbital_density_operators(orbital)


def test_missing_or_compressed_invalid_coordinate_maps():
    h = local_pair(.2)
    result = solve(h, sector=Sector(0, None))
    with pytest.raises(ValueError, match="uncompressed declared electron map"):
        orbital_density_matrix(result)
    with pytest.raises(ValueError, match="outside the declared impurity"):
        orbital_density_matrix(result, 1, physical=False)
    malformed = replace(h, metadata={"coordinates": {"restriction": "dark-vacuum"}})
    with pytest.raises(ValueError, match="verified impurity-only eta map"):
        orbital_density_operators(0, hamiltonian=malformed)
    with pytest.raises(ValueError, match="QP or DMRG"):
        orbital_density_matrix(object())
    with pytest.raises(ValueError, match="result Hamiltonian"):
        orbital_density_matrix(SimpleNamespace(backend="qp", hamiltonian=None))


@pytest.mark.parametrize("kwargs", [
    {"state": -1}, {"state": 2}, {"state": True}, {"state": 1.5},
    {"weights": [1, 1]}, {"state": None, "weights": [1j, 1]},
    {"state": None, "weights": [1]}, {"state": None, "weights": [0, 0]},
    {"state": None, "weights": [1, -1]}, {"state": None, "weights": [np.inf, 1]},
])
def test_invalid_state_or_mixture(kwargs):
    result = solve(local_pair(.2), sector=Sector(0, None), options={"eigenpairs": 2})
    with pytest.raises(ValueError):
        orbital_density_matrix(result, physical=False, **kwargs)


@pytest.mark.parametrize("density", [np.eye(2), np.full((4, 4), np.nan), np.eye(4),
                                    np.diag([-.1, .5, .5, .1]),
                                    np.array([[1, 1j, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])])
def test_invalid_density_properties(density):
    with pytest.raises(ValueError):
        orbital_density_properties(density)
