"""Physical-coordinate equivalence, bath measures, and native-core regressions."""

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import (DiscreteBath, FermionOperator, Impurity, Sector, annihilate, cosh_grid,
                         create, fit_surrogate, hermitian_pair, make_model, number, reference_model, solve)
from qdjj_solver.common.baths import discrete_g, hybridization_g
from qdjj_solver.qp_solver import FockBasis, ProjectedOperator, SolverOptions
from electron_oracle import electron_hamiltonian, electron_operators, evaluate, sector_indices


@pytest.mark.parametrize("backend", ["qp", "dmrg"])
@pytest.mark.parametrize("spin_mixing", [False, True])
def test_physical_electron_oracle(backend, spin_mixing):
    if backend == "dmrg":
        pytest.importorskip("tenpy")
    impurity = Impurity.anderson(u=1.73, detuning=.19, field=.08)
    bath = DiscreteBath([.37], [.81], delta=.93, bandwidth=3.7)
    hopping = [np.array([[.31, .09j if spin_mixing else 0.], [0., .27]]), np.eye(2)*.23]
    phases = [-.41, .68]
    direct = {(0, 1): np.array([[.08+.11j, .04j if spin_mixing else 0.], [0., .08-.11j]])}
    h = make_model(impurity, [bath, bath], hopping, phases=phases, direct=direct,
                   phase_velocities=[-.5, .5])
    sector = Sector(1, None if spin_mixing else 1)
    result = solve(h, sector=sector, backend=backend, options={"eigenpairs": 3})
    matrix = electron_hamiltonian(impurity, [bath, bath], hopping, phases, direct)
    indices = sector_indices(6, sector.parity, sector.twice_sz)
    energies = np.linalg.eigvalsh(matrix[indices][:, indices].toarray())[:3]
    assert_allclose(result.energies, energies, atol=2e-10)
    assert np.max(result.residuals) < 2e-8


def test_multi_orbital_exchange_and_pair_hopping():
    op = sum((-.25+.13*i)*number(i) for i in range(4))
    op += 1.1*number(0)*number(1)+.7*number(2)*number(3)
    op += hermitian_pair(.13j*create(0)*annihilate(3))
    op += hermitian_pair(.08*create(0)*annihilate(1)*create(3)*annihilate(2))
    op += hermitian_pair(-.07*create(0)*create(1)*annihilate(3)*annihilate(2))
    impurity = Impurity(2, op)
    bath = DiscreteBath([.3], [.9], .7, 2.4)
    hopping = [np.array([[.2, .03j, .13, 0.], [0., .22, .06j, .11]])]
    h = make_model(impurity, [bath], hopping, [.32])
    matrix = electron_hamiltonian(impurity, [bath], hopping, [.32])
    selected = sector_indices(6, 1)
    expected = np.linalg.eigvalsh(matrix[selected][:, selected].toarray())[:3]
    qp = solve(h, sector=Sector(1, None), options={"eigenpairs": 3})
    assert_allclose(qp.energies, expected, atol=2e-12)
    pytest.importorskip("tenpy")
    mps = solve(h, sector=Sector(1, None), backend="dmrg", options={"eigenpairs": 3})
    assert_allclose(mps.energies, expected, atol=3e-10)


def test_normal_ordering_and_projection_boundary():
    cs = electron_operators(4)
    rng = np.random.default_rng(123)
    for length in range(1, 9):
        for _ in range(10):
            key = tuple(int(x) for x in rng.choice([-4, -3, -2, -1, 1, 2, 3, 4], size=length))
            expected = np.eye(16)
            for op in key:
                expected = expected @ (cs[op-1].getH() if op > 0 else cs[-op-1])
            assert_allclose(evaluate(FermionOperator({key: 1}), cs).toarray(), expected, atol=1e-14)
    basis = FockBasis(2, (1, -1)*2, 0, Sector(0, 0))
    correct = ProjectedOperator(annihilate(2)*create(2), basis).sparse_matrix()
    separately = ProjectedOperator(annihilate(2), basis).sparse_matrix() @ ProjectedOperator(create(2), basis).sparse_matrix()
    assert_allclose(correct.toarray(), np.eye(basis.dimension))
    assert_allclose(separately.toarray(), 0.)


def test_word_boundary_exchange_and_number_sector():
    basis = FockBasis(0, (1, -1)*40, 2, Sector(0, 2, particle_number=2))
    occupations = basis.occupations()
    source = np.flatnonzero((occupations[:, 0] == 1) & (occupations[:, 1] == 1))[0]
    target = np.flatnonzero((occupations[:, 0] == 0) & (occupations[:, 1] == 1+(1 << 14)))[0]
    op = ProjectedOperator(create(78)*annihilate(0), basis)
    vector = np.zeros(basis.dimension)
    vector[source] = 1.
    assert_allclose(op.action(vector)[target], -1.)
    assert basis.dimension == 40*39//2


def test_bath_measure_and_surrogate():
    bath = cosh_grid(24)
    assert bath.paired and not bath.weights.flags.writeable
    assert_allclose(bath.weights.sum(), 1., atol=2e-14)
    omega = np.geomspace(1e-4, 1e4, 100)
    assert_allclose(discrete_g(bath, omega), hybridization_g(omega), rtol=3e-10)
    assert abs(cosh_grid(1).metadata["measure_error"]) > .01
    fitted = fit_surrogate(2, bandwidth=10., frequency_cutoff=10., starts=2)
    # Published rounded coefficients, Baran et al., PRB 108, L220506 (2023).
    assert_allclose(fitted.xi, [-1.3099, 1.3099], atol=5e-5, rtol=0.)
    assert_allclose(fitted.weights/(np.pi*fitted.rho), [1.2642, 1.2642], atol=5e-5, rtol=0.)
    assert not DiscreteBath([-1., 1.+1e-12], [.5, .5]).paired


def test_eigensolvers_cutoffs_and_coordinate_derivative():
    bath = cosh_grid(1, bandwidth=4.)
    h = reference_model(bath, phi=.81, rho_ws=.04)
    sector = Sector(1, 1, 1)
    energies = [solve(h, cutoff=q, sector=sector).energies[0] for q in range(5)]
    assert np.all(np.diff(energies) <= 1e-12)
    results = [solve(h, cutoff=3, sector=sector, options=SolverOptions(method=method, eigenpairs=2))
               for method in ("dense", "sparse", "matrix-free")]
    for result in results[1:]:
        assert_allclose(result.energies, results[0].energies, atol=2e-11)
    step = 1e-5
    energies = [solve(reference_model(bath, phi=.81+s*step, rho_ws=.04), cutoff=3, sector=sector).energies[0]
                for s in (-1, 1)]
    assert_allclose((energies[1]-energies[0])/(2*step), results[0].observables["phase_derivative"][0], atol=3e-8)
    for vacuum in ("isolated", "coupled"):
        model = reference_model(bath, phi=.81, bath_reference=vacuum, compress=False)
        reduced = solve(model, cutoff=2, options={"method": "schur"})
        full = solve(model, cutoff=2, options={"method": "sparse"})
        assert_allclose(reduced.energies, full.energies, atol=2e-11)
