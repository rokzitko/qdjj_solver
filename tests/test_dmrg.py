"""Independent tensor-product ED checks, including genuinely positive target roots."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_limits

pytest.importorskip("tenpy")

from qdjj_solver.common import Hamiltonian, Sector, annihilate, create, hermitian_pair, number
from qdjj_solver.dmrg_solver import SolverOptions, solve


def dense_operator(operator, modes):
    # Mode zero is the leftmost tensor factor; no QP basis/native-core code is used.
    identity = np.eye(2)
    z = np.diag([1., -1.])
    a = np.array([[0., 1.], [0., 0.]])
    annihilators = []
    for i in range(modes):
        matrix = np.ones((1, 1))
        for j in range(modes):
            matrix = np.kron(matrix, z if j < i else (a if j == i else identity))
        annihilators.append(matrix)
    result = np.zeros((2**modes, 2**modes), complex)
    for key, value in operator.terms.items():
        matrix = value*np.eye(2**modes, dtype=complex)
        for op in key:
            a = annihilators[abs(op)-1]
            matrix = matrix @ (a.T if op > 0 else a)
        result += matrix
    return result


def exact(hamiltonian, sector, roots):
    n = len(hamiltonian.spins)
    occ = np.array([[(mask >> (n-1-i)) & 1 for i in range(n)] for mask in range(2**n)])
    selected = occ.sum(axis=1) % 2 == sector.parity
    if sector.twice_sz is not None:
        selected &= occ @ hamiltonian.spins == sector.twice_sz
    if sector.eta is not None:
        selected &= np.prod(np.where(occ, hamiltonian.eta_labels, 1), axis=1) == sector.eta
    if sector.particle_number is not None:
        selected &= occ.sum(axis=1) == sector.particle_number
    indices = np.flatnonzero(selected)
    with threadpool_limits(limits=1):
        full = dense_operator(hamiltonian.operator, n)
        energies, vectors = np.linalg.eigh(full[np.ix_(indices, indices)])
    embedded = np.zeros((2**n, roots), complex)
    embedded[indices] = vectors[:, :roots]
    return energies[:roots], embedded


def general_problem(spin_mixing=False):
    op = sum((.13+.09*i)*number(i) for i in range(6)) + 5.
    op += .8*number(0)*number(1) + .23*number(0)*number(4)
    op += hermitian_pair((.21+.13j)*create(0)*annihilate(4))
    op += hermitian_pair((-.17+.09j)*create(1)*annihilate(5))
    op += hermitian_pair((.14-.05j)*create(2)*create(3))
    op += hermitian_pair(.11j*create(0)*create(3)*annihilate(2)*annihilate(1))
    if spin_mixing:
        op += hermitian_pair(.19j*create(0)*annihilate(3))
    return Hamiltonian(op, 2, (1, -1)*3, (1, 1, -1, -1, 1, 1),
                       {"density": number(0), "pair": hermitian_pair(create(2)*create(3))},
                       {"coordinates": {"kind": "declared-mode-order", "labels": list(range(6))}})


@pytest.mark.parametrize("group_size,order,spin_mixing,eta", [
    (2, None, True, None),
    (1, (4, 2, 0, 5, 3, 1), False, None),
    (2, (5, 0, 2, 4, 1, 3), False, -1),
])
def test_excited_roots_against_independent_ed(group_size, order, spin_mixing, eta):
    h = general_problem(spin_mixing)
    sector = Sector(1, None if spin_mixing else 1, eta)
    energies, vectors = exact(h, sector, 3)
    result = solve(h, sector=sector, options=SolverOptions(
        eigenpairs=3, chi_max=32, group_size=group_size, mode_order=order,
        residual_tolerance=2e-8, require_convergence=True))
    assert np.min(energies) > 0
    assert_allclose(result.energies, energies, atol=2e-10)
    assert np.max(result.residuals) < 2e-8
    assert_allclose(result.overlaps(), np.eye(3), atol=2e-9)
    for name, operator in h.observables.items():
        with threadpool_limits(limits=1):
            dense = vectors.conj().T @ dense_operator(operator, 6) @ vectors
        assert_allclose(result.observables[name], np.diag(dense), atol=2e-8)
        assert_allclose(abs(result.matrix_elements(operator)), abs(dense), atol=2e-8)


def test_parity_and_spin_changing_transitions():
    h = general_problem()
    options = SolverOptions(eigenpairs=2, chi_max=32, require_convergence=True)
    odd = solve(h, sector=Sector(1, 1), options=options)
    even = solve(h, sector=Sector(0, 0), options=options)
    _, v_odd = exact(h, Sector(1, 1), 2)
    _, v_even = exact(h, Sector(0, 0), 2)
    operator = create(4)+.3j*create(0)
    with threadpool_limits(limits=1):
        expected = v_odd.conj().T @ dense_operator(operator, 6) @ v_even
    assert_allclose(abs(odd.matrix_elements(operator, even)), abs(expected), atol=3e-8)
    assert_allclose(even.matrix_elements(operator.dagger(), odd),
                    odd.matrix_elements(operator, even).conj().T, atol=2e-12)


def test_single_site_and_zero_operator():
    op = 2.+.3*number(0)+.7*number(1)+hermitian_pair(.2j*create(0)*create(1))
    h = Hamiltonian(op, 2, (1, -1))
    result = solve(h, sector=Sector(0, 0), options=SolverOptions(eigenpairs=2, require_convergence=True))
    energies, _ = exact(h, Sector(0, 0), 2)
    assert_allclose(result.energies, energies, atol=2e-14)
    assert np.max(result.residuals) < 2e-14
    zero = Hamiltonian(0.*op, 2, (1, -1)*2)
    result = solve(zero, options=SolverOptions(eigenpairs=2, require_convergence=True))
    assert_allclose(result.energies, 0.)
    assert_allclose(result.overlaps(), np.eye(2), atol=1e-15)


def test_cutoff_and_incompatible_continuation_are_rejected():
    h = general_problem()
    with pytest.raises(ValueError, match="no QP cutoff"):
        solve(h, cutoff=2)
    with pytest.raises(ValueError, match="eta reduction is not defined"):
        solve(Hamiltonian(0.*h.operator, h.nimp, h.spins), sector=Sector(1, 1, 1))
    result = solve(h)
    with pytest.raises(ValueError, match="incompatible coordinates"):
        solve(h, initial=result, options=SolverOptions(group_size=1))


def test_particle_number_sector_shared_with_native_backend():
    from qdjj_solver import solve as dispatch
    from qdjj_solver.common import FermionOperator
    h = general_problem(spin_mixing=True)
    h.operator = FermionOperator({key: value for key, value in h.operator.terms.items()
                                  if sum(1 if op > 0 else -1 for op in key) == 0})
    sector = Sector(1, None, particle_number=3)
    energies, _ = exact(h, sector, 3)
    qp = dispatch(h, backend="qp", sector=sector, options={"eigenpairs": 3})
    mps = solve(h, sector=sector, options=SolverOptions(eigenpairs=3, require_convergence=True,
                                                      excitation_operator="density"))
    assert qp.basis.dimension == 20
    assert_allclose(qp.energies, energies, atol=2e-12)
    assert_allclose(mps.energies, energies, atol=2e-10)
    with pytest.raises(ValueError, match="does not conserve particle number"):
        solve(general_problem(), sector=Sector(1, 1, particle_number=3))


def test_long_heterogeneous_layout_is_a_finite_chain():
    modes = 40
    energies = 1.+.03*np.arange(modes)
    op = sum(energy*number(i) for i, energy in enumerate(energies))
    op += sum(hermitian_pair(.07j*create(i)*annihilate(i+1)) for i in range(modes-1))
    h = Hamiltonian(op, 0, (1, -1)*(modes//2))
    result = solve(h, sector=Sector(1, None, particle_number=1),
                   options=SolverOptions(eigenpairs=2, chi_max=8, seed_trials=1, require_convergence=True))
    exact_matrix = np.diag(energies).astype(complex)
    exact_matrix += np.diag(np.full(modes-1, .07j), 1)+np.diag(np.full(modes-1, -.07j), -1)
    assert_allclose(result.energies, np.linalg.eigvalsh(exact_matrix)[:2], atol=2e-10)


def test_uncompressed_residual_of_an_intentionally_truncated_mps():
    from tenpy.tools.misc import TenpyInconsistencyWarning
    h = general_problem(spin_mixing=True)
    with pytest.warns(TenpyInconsistencyWarning, match="Maximum truncation error"):
        result = solve(h, sector=Sector(1, None), options=SolverOptions(
            chi_max=1, group_size=1, seed_trials=1, max_sweeps=12))
    psi = result.states[0]
    tensor = psi.get_theta(0, n=psi.L).to_ndarray().reshape([2]*psi.L)
    for i, site in enumerate(psi.sites):
        tensor = np.take(tensor, [site.state_index("empty"), site.state_index("full")], axis=i)
    vector = tensor.reshape(-1)
    with threadpool_limits(limits=1):
        matrix = dense_operator(h.operator, 6)
        expected = np.linalg.norm(matrix @ vector-result.energies[0]*vector)
    assert expected > .01
    assert_allclose(result.residuals[0], expected, atol=2e-13)
    assert not result.metadata["finite_problem_converged"]


def test_bond_ladder_against_a_small_exact_problem():
    from qdjj_solver import cosh_grid, reference_model
    from qdjj_solver.dmrg_solver.convergence import converge_bond_dimension
    h = reference_model(cosh_grid(1, bandwidth=3.), geometry="single", detuning=.13)
    archived = []
    result, assessment = converge_bond_dimension(
        h, [4, 8, 16], options=SolverOptions(seed_trials=1),
        callback=lambda state, entry: archived.append(entry["chi_max"]))
    expected, _ = exact(h, Sector(), 1)
    assert archived == [4, 8, 16]
    assert assessment["empirical_convergence"]
    assert not assessment["bath_convergence_checked"]
    assert_allclose(result.energies, expected, atol=2e-11)


def test_exact_mpo_channel_merging_for_a_centered_star():
    from qdjj_solver.dmrg_solver import prepare
    from tenpy.networks.mps import MPS
    modes = 26
    op = sum((i+1)*number(i) for i in range(modes))
    op += sum(hermitian_pair((.1+.01j*i)*create(0)*annihilate(i)) for i in range(2, modes, 2))
    h = Hamiltonian(op, 2, (1, -1)*(modes//2))
    order = tuple(range(2, 14))+(0, 1)+tuple(range(14, modes))
    prepared = prepare(h, Sector(1, 1, particle_number=1), mode_order=order)
    unreduced = prepared.compile(op, merge_parallel=False)
    reduced = prepared.compile(op)
    assert max(reduced.chi) < max(unreduced.chi)
    # A delocalized charged MPS probes off-diagonal as well as onsite entries.
    psi = MPS.from_product_state(prepared.sites, ["up"]+["empty"]*(len(prepared.sites)-1),
                                bc="finite", dtype=complex, unit_cell_width=1)
    from qdjj_solver.dmrg_solver.solver import _randomize_seed
    _randomize_seed(psi, np.random.default_rng(77), 12, 8)
    a, b = psi.copy(), psi.copy()
    unreduced.apply_naively(a)
    reduced.apply_naively(b)
    a.canonical_form(renormalize=False)
    b.canonical_form(renormalize=False)
    assert_allclose(a.overlap(a), b.overlap(b), atol=3e-12)
    assert_allclose(a.overlap(b), a.overlap(a), atol=3e-12)
