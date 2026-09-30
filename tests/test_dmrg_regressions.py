"""Scientific and API regressions supplementing the baseline DMRG tests."""

from dataclasses import replace
from itertools import product

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from threadpoolctl import threadpool_limits

pytest.importorskip("tenpy")

from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS

from qdjj_solver.common import (
    FermionOperator, Hamiltonian, Sector, annihilate, create, hermitian_pair, number,
)
from qdjj_solver.dmrg_solver import SolverOptions, solve, solver
from qdjj_solver.dmrg_solver.convergence import compare_states, converge_bond_dimension
from qdjj_solver.dmrg_solver.mpo import prepare, residual_norm
from qdjj_solver.dmrg_solver.solver import _randomize_seed


pytestmark = pytest.mark.dmrg


@pytest.fixture(autouse=True)
def single_blas_thread():
    with threadpool_limits(limits=1):
        yield


def jw_matrix(operator, nmodes):
    """Independent full Jordan--Wigner matrix, with mode zero leftmost."""
    annihilators = []
    for mode in range(nmodes):
        matrix = np.ones((1, 1))
        for index in range(nmodes):
            local = (np.diag([1., -1.]) if index < mode else
                     np.array([[0., 1.], [0., 0.]]) if index == mode else np.eye(2))
            matrix = np.kron(matrix, local)
        annihilators.append(matrix)
    result = np.zeros((2**nmodes, 2**nmodes), complex)
    for key, coefficient in operator.terms.items():
        term = coefficient*np.eye(2**nmodes, dtype=complex)
        for op in key:
            local = annihilators[abs(op)-1]
            term = term @ (local.conj().T if op > 0 else local)
        result += term
    return result


def sector_indices(hamiltonian, sector):
    occ = np.array(list(product((0, 1), repeat=len(hamiltonian.spins))))
    selected = occ.sum(axis=1) % 2 == sector.parity
    if sector.twice_sz is not None:
        selected &= occ @ hamiltonian.spins == sector.twice_sz
    if sector.eta is not None:
        selected &= np.prod(np.where(occ, hamiltonian.eta_labels, 1), axis=1) == sector.eta
    if sector.particle_number is not None:
        selected &= occ.sum(axis=1) == sector.particle_number
    return np.flatnonzero(selected)


def fock_coordinates(prepared):
    """Local tensor-basis kets in physical Fock coordinates, including signs."""
    nmodes = len(prepared.order)
    transform = np.zeros((2**nmodes, 2**nmodes))
    labels = {(0,): "empty", (1,): "full", (0, 0): "empty",
              (1, 0): "up", (0, 1): "down", (1, 1): "full"}
    for row, occ in enumerate(product((0, 1), repeat=nmodes)):
        occupied_order = [mode for mode in prepared.order if occ[mode]]
        inversions = sum(a > b for i, a in enumerate(occupied_order)
                         for b in occupied_order[i+1:])
        indices = tuple(site.state_index(labels[tuple(occ[m] for m in group)])
                        for site, group in zip(prepared.sites, prepared.groups, strict=True))
        column = np.ravel_multi_index(indices, tuple(site.dim for site in prepared.sites))
        transform[row, column] = (-1)**inversions
    return transform


def mpo_matrix(mpo):
    """Contract virtual indices directly, independently of MPS expectations."""
    labels = ["wL", "wR", "p", "p*"]
    full = mpo.get_W(0).transpose(labels).to_ndarray()[mpo.get_IdL(0)]
    for index in range(1, mpo.L):
        tensor = mpo.get_W(index).transpose(labels).to_ndarray()
        full = np.einsum("aij,abkl->bikjl", full, tensor)
        full = full.reshape(tensor.shape[1], full.shape[1]*full.shape[2], -1)
    return full[mpo.get_IdR(mpo.L-1)]


def state_from_vector(prepared, vector):
    local = fock_coordinates(prepared).T @ vector
    if len(prepared.sites) == 1:
        return MPS.from_product_state(prepared.sites, [local], bc="finite", dtype=complex,
                                      permute=False, unit_cell_width=1)
    tensor = npc.Array.from_ndarray(
        local.reshape(tuple(site.dim for site in prepared.sites)),
        [site.leg for site in prepared.sites],
        labels=[f"p{i}" for i in range(len(prepared.sites))],
        qtotal=prepared.target_charge, cutoff=0.)
    return MPS.from_full(prepared.sites, tensor, cutoff=0., unit_cell_width=1)


def state_vector(prepared, state):
    tensor = state.get_theta(0, n=state.L).to_ndarray().reshape(-1)
    return fock_coordinates(prepared) @ tensor * state.norm


def random_states(prepared, count=2, seed=1909):
    indices = sector_indices(prepared.hamiltonian, prepared.sector)
    count = min(count, len(indices))
    rng = np.random.default_rng(seed)
    vectors, _ = np.linalg.qr(rng.normal(size=(len(indices), count))
                             + 1j*rng.normal(size=(len(indices), count)))
    full = np.zeros((2**len(prepared.order), count), complex)
    full[indices] = vectors
    return full, [state_from_vector(prepared, vector) for vector in full.T]


def diagonal_problem(nmodes=4):
    op = sum((.2+.31*i)*number(i) for i in range(nmodes)) + 2.
    return Hamiltonian(op, min(2, nmodes), tuple(1 if i % 2 == 0 else -1 for i in range(nmodes)),
                       tuple(1 if i % 3 == 0 else -1 for i in range(nmodes)),
                       {"density": number(0)}, {"coordinates": {"labels": list(range(nmodes))}})


@pytest.mark.parametrize("group_size,order", [
    (1, None), (1, (4, 2, 0, 3, 1)), (2, None), (2, (3, 0, 4, 2, 1)),
])
def test_full_mpo_in_signed_reordered_and_ragged_coordinates(group_size, order):
    op = .43 + sum((-.21+.11*i)*number(i) for i in range(5))
    op += .37*number(0)*number(3) + .06*number(0)*number(2)*number(4)
    op += hermitian_pair((.27+.19j)*create(0)*annihilate(4))
    op += hermitian_pair((-.17+.07j)*create(1)*create(3))
    op += hermitian_pair((.13-.09j)*create(0)*create(2)*annihilate(3)*annihilate(1))
    h = Hamiltonian(op, 2, (1, -1, 1, -1, 1))
    prepared = prepare(h, Sector(0, None), group_size=group_size, mode_order=order)
    transform = fock_coordinates(prepared)
    assert_array_equal(transform.T @ transform, np.eye(32))
    if order is not None:
        assert np.any(transform == -1)  # A bosonic permutation is not sufficient.
    for merge in (False, True):
        actual = transform @ mpo_matrix(prepared.compile(op, merge_parallel=merge)) @ transform.T
        assert_allclose(actual, jw_matrix(op, 5), atol=1e-12, rtol=0.)
    assert prepared.groups[-1] == (prepared.order[-1],)


@pytest.mark.parametrize("coefficient", [1e-100, (1.+2j)*1e-35])
def test_tiny_mpo_coefficients_are_retained(coefficient):
    prepared = prepare(diagonal_problem(3), Sector(1, None), mode_order=(2, 0, 1))
    operator = hermitian_pair(coefficient*create(0)*annihilate(2)) + abs(coefficient)*number(1)
    transform = fock_coordinates(prepared)
    for merge in (False, True):
        actual = transform @ mpo_matrix(prepared.compile(operator, merge_parallel=merge)) @ transform.T
        # No absolute floor: silently discarding these coefficients must fail.
        assert_allclose(actual, jw_matrix(operator, 3), atol=0., rtol=1e-12)


@pytest.mark.parametrize("bra_sector,ket_sector,operator", [
    (Sector(1, None), Sector(0, None), create(4)+.31j*create(0)*number(2)),
    (Sector(0, None, particle_number=2), Sector(0, None, particle_number=0),
     (.7+.2j)*create(1)*create(4)),
    (Sector(1, 1), Sector(1, -1), (.2-.8j)*create(0)*annihilate(1)),
    (Sector(1, None, -1), Sector(1, None, 1), (.3+.7j)*create(2)*annihilate(0)),
    (Sector(1, None), Sector(1, None), .4j*number(1)+(.2-.3j)*create(4)*annihilate(0)+.17),
])
def test_complex_transition_amplitudes(bra_sector, ket_sector, operator):
    h = diagonal_problem(5)
    bras = prepare(h, bra_sector, mode_order=(4, 1, 0, 3, 2))
    kets = prepare(h, ket_sector, mode_order=(4, 1, 0, 3, 2))
    left, bra_states = random_states(bras, seed=811)
    right, ket_states = random_states(kets, seed=97)
    expected = left.conj().T @ jw_matrix(operator, 5) @ right
    assert np.max(abs(expected)) > .01
    actual = bras.matrix_elements(operator, bra_states, ket_states)
    assert_allclose(actual, expected, atol=1e-12, rtol=0.)
    assert_allclose(kets.matrix_elements(operator.dagger(), ket_states, bra_states),
                    actual.conj().T, atol=1e-12, rtol=0.)


def test_single_site_randomization_regression():
    prepared = prepare(diagonal_problem(2), Sector(1, None))
    vector = np.array([0., 1., 1j, 0.])/np.sqrt(2.)
    state = state_from_vector(prepared, vector)
    # TeNPy's finite canonicalization sweep asserts L > 1. This must be a no-op.
    _randomize_seed(state, np.random.default_rng(83), 3, 4)
    assert_allclose(state_vector(prepared, state), vector, atol=1e-12, rtol=0.)


@pytest.mark.parametrize("group_size,nmodes", [(2, 2), (1, 3), (2, 3)])
def test_residual_of_inexact_states_and_wrong_energies(group_size, nmodes):
    h = diagonal_problem(nmodes)
    h.operator += hermitian_pair((.17+.23j)*create(0)*annihilate(nmodes-1))
    prepared = prepare(h, Sector(1, None), group_size=group_size,
                       mode_order=tuple(reversed(range(nmodes))))
    vectors, states = random_states(prepared, count=1)
    matrix = jw_matrix(h.operator, nmodes)
    energy = -1.7
    expected = np.linalg.norm((matrix-energy*np.eye(2**nmodes)) @ vectors[:, 0])
    assert expected > 1.
    assert_allclose(residual_norm(prepared, states[0], energy), expected, atol=1e-12, rtol=0.)
    indices = sector_indices(h, prepared.sector)
    energies, eigenvectors = np.linalg.eigh(matrix[np.ix_(indices, indices)])
    vector = np.zeros(2**nmodes, complex)
    vector[indices] = eigenvectors[:, 0]
    eigenstate = state_from_vector(prepared, vector)
    assert_allclose(residual_norm(prepared, eigenstate, energies[0]+.125),
                    .125, atol=1e-12, rtol=0.)


def test_convergence_refuses_a_valid_mps_with_nonzero_physical_residual(monkeypatch):
    h = diagonal_problem(3)
    sector = Sector(1, None)
    prepared = prepare(h, sector, group_size=1)
    vector = np.zeros(8, complex)
    vector[[1, 2]] = np.array([1., 1j])/np.sqrt(2.)
    state = state_from_vector(prepared, vector)

    class InexactEngine:
        sweep_stats = {"E": [-1., -1., -1.]}
        sweeps = 4

        def __init__(self, *args, **kwargs):
            pass

        def run(self):
            return -1., state.copy()

        def is_converged(self):
            return True

    monkeypatch.setattr(solver.dmrg, "TwoSiteDMRGEngine", InexactEngine)
    options = SolverOptions(group_size=1, seed_trials=1, mixer=False, min_sweeps=2, max_sweeps=4)
    result = solve(h, sector=sector, options=options)
    expected = np.linalg.norm(jw_matrix(h.operator, 3) @ vector-result.energies[0]*vector)
    assert expected > .1
    assert_allclose(result.residuals, [expected], atol=1e-12, rtol=0.)
    assert not result.metadata["finite_problem_converged"]
    with pytest.raises(RuntimeError, match="finite-problem convergence failed"):
        solve(h, sector=sector, options=replace(options, require_convergence=True))


def test_constant_shift_covariance_of_multiple_roots():
    base = .2*number(0)+.6*number(1)+.9*number(2)
    base += hermitian_pair((.17+.11j)*create(0)*annihilate(2))
    base += hermitian_pair((-.08+.13j)*create(1)*annihilate(2))
    h = Hamiltonian(base, 1, (1, -1, 1), observables={"density": number(0)},
                    metadata={"coordinates": {"labels": ["a", "b", "c"]}})
    options = SolverOptions(eigenpairs=2, group_size=1, chi_max=4, seed_trials=1,
                            mixer=False, min_sweeps=4, max_sweeps=12, require_convergence=True,
                            mode_order=(2, 0, 1), residual_tolerance=2e-10)
    sector = Sector(1, None, particle_number=1)
    initial = solve(h, sector=sector, options=options)
    shifted = solve(replace(h, operator=base+250.), sector=sector, options=options, initial=initial)
    indices = sector_indices(h, sector)
    energies = np.linalg.eigvalsh(jw_matrix(base, 3)[np.ix_(indices, indices)])[:2]
    assert_allclose(initial.energies, energies, atol=1e-12, rtol=0.)
    assert_allclose(shifted.energies, energies+250., atol=1e-11, rtol=0.)
    assert_allclose(abs(initial.overlaps(shifted)), np.eye(2), atol=2e-11, rtol=0.)
    assert_allclose(initial.observables["density"], shifted.observables["density"],
                    atol=2e-11, rtol=0.)


def test_lanczos_probability_tolerance_default_and_positional_options():
    assert SolverOptions().lanczos_probability_tolerance == 1e-14
    # The new field must not consume the existing positional threads argument.
    options = SolverOptions(2, 16, 16, 10, 1e-11, 1e-8, 1e-7, 1e-8, 1e-14, 100,
                            1, 1729, 4, None, True, 1e-5, 6, 1, None, None, True, False, 3)
    assert options.threads == 3
    assert options.lanczos_probability_tolerance == 1e-14


@pytest.mark.parametrize("value", [0., -1e-14, np.nan, np.inf, -np.inf])
def test_invalid_lanczos_probability_tolerance(value):
    with pytest.raises(ValueError, match="lanczos_probability_tolerance must be positive and finite"):
        SolverOptions(lanczos_probability_tolerance=value)


def test_lanczos_probability_tolerance_improves_independent_state_residuals(monkeypatch):
    # Physical electrons: an interacting spin-mixed impurity and two paired bath
    # levels. Six modes allow full ED; the high level exposes local stopping error.
    op = sum(e*number(i) for i, e in enumerate((-.63, -.77, .23, .23, 1000., 1000.)))
    op += 1.9*number(0)*number(1) + hermitian_pair(.13j*create(0)*annihilate(1))
    for mode, hopping, phase in ((2, .37, .4), (4, .29, -.7)):
        op += hermitian_pair(.8*np.exp(1j*phase)*create(mode)*create(mode+1))
        for spin in (0, 1):
            op += hermitian_pair(hopping*np.exp(.17j*(mode+spin))*create(spin)*annihilate(mode+spin))
    h = Hamiltonian(op, 2, (1, -1)*3)
    sector = Sector(0, None)
    matrix = jw_matrix(op, 6)
    indices = sector_indices(h, sector)
    energies, eigenvectors = np.linalg.eigh(matrix[np.ix_(indices, indices)])

    probabilities = []
    engine_init = solver.dmrg.TwoSiteDMRGEngine.__init__

    def track_engine(self, *args, **kwargs):
        engine_init(self, *args, **kwargs)
        probabilities.append(self.lanczos_params["P_tol"])
        assert self.options["P_tol_to_trunc"] is None

    monkeypatch.setattr(solver.dmrg.TwoSiteDMRGEngine, "__init__", track_engine)
    options = SolverOptions(eigenpairs=2, group_size=1, chi_max=16, seed_trials=1,
                            max_sweeps=16, seed=1729)
    residuals = []
    # The deliberately inexact default run need only meet its configured overlap
    # tolerance; retain the stricter orthogonality check for the tightened solve.
    for current, overlap_tolerance in (
        (options, options.orthogonality_tolerance),
        (replace(options, lanczos_probability_tolerance=1e-22), 2e-10),
    ):
        result = solve(h, sector=sector, options=current)
        vectors = np.column_stack([state_vector(result.prepared, state) for state in result.states])
        independent = np.linalg.norm(matrix @ vectors-vectors*result.energies, axis=0)
        residuals.append(independent)
        assert result.metadata["options"]["lanczos_probability_tolerance"] == current.lanczos_probability_tolerance
        assert_allclose(result.energies, energies[:2], atol=2e-10, rtol=0.)
        assert_allclose(vectors.conj().T @ vectors, np.eye(2), atol=overlap_tolerance, rtol=0.)
        assert_allclose(abs(eigenvectors[:, :2].conj().T @ vectors[indices]), np.eye(2),
                        atol=1e-8, rtol=0.)
        assert_allclose(result.residuals, independent, atol=2e-12, rtol=0.)
    assert probabilities == [1e-14, 1e-14, 1e-22, 1e-22]
    assert max(residuals[1]) < 1e-8
    # Compare the worst root, not roundoff-limited individual residuals. Tenfold
    # improvement leaves ample margin without fixing platform-specific values.
    assert max(residuals[1]) < .1*max(residuals[0])


@pytest.mark.parametrize("kwargs", [
    {"eigenpairs": 0}, {"chi_max": True}, {"max_sweeps": 1.5},
    {"seed": -1}, {"seed": True}, {"seed_randomization_steps": True},
    {"residual_tolerance": 0.}, {"energy_tolerance": np.nan}, {"orthogonality_tolerance": np.inf},
    {"svd_min": -1.}, {"svd_min": np.nan}, {"excitation_operator": 7},
    {"min_sweeps": 61}, {"min_sweeps": 6},
    {"group_size": 3}, {"group_size": True}, {"group_size": 1.0},
    {"require_convergence": True, "calculate_residuals": False},
])
def test_invalid_solver_options(kwargs):
    with pytest.raises(ValueError):
        SolverOptions(**kwargs)


@pytest.mark.parametrize("kwargs", [
    {"group_size": 3}, {"group_size": True}, {"group_size": 1.0},
    {"mode_order": (0, 1, 2)}, {"mode_order": (0, 1, 2, 2)}, {"mode_order": (0, 1, 2, 4)},
])
def test_invalid_prepared_layouts(kwargs):
    with pytest.raises(ValueError):
        prepare(diagonal_problem(), **kwargs)


def test_invalid_mode_and_sector_dimensions():
    with pytest.raises(ValueError, match="at least one fermionic mode"):
        prepare(Hamiltonian(FermionOperator(), 0, ()))
    with pytest.raises(ValueError, match="more eigenpairs"):
        solve(diagonal_problem(2), options=SolverOptions(eigenpairs=2))
    with pytest.raises(ValueError, match="more eigenpairs"):
        solve(diagonal_problem(), sector=Sector(1, 7))
    with pytest.raises(ValueError, match="undeclared mode"):
        prepare(diagonal_problem()).matrix_elements(number(7), [], [])


def test_incompatible_coordinates_layouts_and_charge_schemas():
    h = diagonal_problem(2)
    sector = Sector(1, None)
    reference = solve(h, sector=sector)
    preparations = [
        prepare(replace(h, metadata={"coordinates": {"labels": ["different"]}}), sector),
        prepare(h, sector, group_size=1),
        prepare(h, sector, mode_order=(1, 0)),
        prepare(h, Sector(1, 1)),
    ]
    for prepared in preparations:
        # The same physical one-particle eigenstate, in each declared layout.
        state = state_from_vector(prepared, np.array([0., 0., 1., 0.]))
        other = replace(reference, prepared=prepared, hamiltonian=prepared.hamiltonian, states=[state])
        with pytest.raises(ValueError, match="different coordinates, layouts, or conserved charges"):
            reference.overlaps(other)
        with pytest.raises(ValueError, match="different coordinates, layouts, or conserved charges"):
            reference.matrix_elements(number(0), other)


def test_invalid_continuation_and_observable_requests():
    h = diagonal_problem(2)
    initial = solve(h)
    with pytest.raises(ValueError, match="incompatible coordinates"):
        solve(h, initial=object())
    with pytest.raises(ValueError, match="different sector"):
        solve(h, sector=Sector(1, -1), initial=initial)
    with pytest.raises(ValueError, match="unknown requested observable"):
        solve(h, options=SolverOptions(observables=("missing",)))
    with pytest.raises(ValueError, match="name a model observable"):
        solve(h, options=SolverOptions(excitation_operator="missing"))
    h.observables["charge_changing"] = create(0)
    with pytest.raises(ValueError, match="preserve the target sector charges"):
        solve(h, options=SolverOptions(excitation_operator="charge_changing"))


def test_state_comparison_matches_permuted_and_phased_roots():
    previous = solve(diagonal_problem(2), sector=Sector(1, None), options=SolverOptions(eigenpairs=2))
    permutation = [1, 0]
    phases = np.exp(1j*np.array([.3, -1.4]))
    states = [state_from_vector(previous.prepared, phase*state_vector(previous.prepared, previous.states[i]))
              for i, phase in zip(permutation, phases, strict=True)]
    current = replace(previous, states=states, energies=previous.energies[permutation],
                      observables={name: values[permutation] for name, values in previous.observables.items()})
    comparison = compare_states(previous, current)
    assert comparison["assignment"] == [(0, 1), (1, 0)]
    assert_allclose(comparison["energy_changes"], 0., atol=1e-12, rtol=0.)
    assert_allclose(comparison["observable_changes"]["density"], 0., atol=1e-12, rtol=0.)
    assert_allclose(comparison["subspace_singular_values"], 1., atol=1e-12, rtol=0.)


def test_state_comparison_recognizes_a_rotated_degenerate_subspace():
    h = replace(diagonal_problem(2), operator=number(0)+number(1))
    previous = solve(h, sector=Sector(1, None), options=SolverOptions(eigenpairs=2))
    vectors = np.column_stack([state_vector(previous.prepared, state) for state in previous.states])
    unitary = np.array([[1., 1j], [1j, 1.]])/np.sqrt(2.)
    states = [state_from_vector(previous.prepared, vector) for vector in (vectors @ unitary).T]
    current = replace(previous, states=states, observables={})
    comparison = compare_states(previous, current)
    assert_allclose(comparison["overlap_squared"], .5, atol=1e-12, rtol=0.)
    assert_allclose(comparison["subspace_singular_values"], 1., atol=1e-12, rtol=0.)
    assert_allclose(comparison["energy_changes"], 0., atol=1e-12, rtol=0.)


@pytest.mark.parametrize("kwargs", [
    {"chi_values": []}, {"chi_values": [0]}, {"chi_values": [True]},
    {"chi_values": [2, 2]}, {"chi_values": [4, 2]},
    {"stable_steps": 0}, {"stable_steps": 1.5}, {"stable_steps": True}, {"stable_steps": np.nan},
    {"energy_tolerance": np.nan}, {"energy_tolerance": 0.}, {"observable_tolerance": 0.},
])
def test_invalid_bond_ladder_targets(kwargs):
    with pytest.raises(ValueError):
        converge_bond_dimension(diagonal_problem(2), **({"chi_values": [2, 4]} | kwargs))
