"""QP solver contracts, independent occupation oracles, and numerical failures."""

import itertools
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.sparse import csc_matrix, issparse
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator

from qdjj_solver.common.algebra import FermionOperator, annihilate, create, number
from qdjj_solver.common.problem import Hamiltonian, Sector
from qdjj_solver.qp_solver import basis as basis_module
from qdjj_solver.qp_solver import solver as solver_module
from qdjj_solver.qp_solver.basis import FockBasis, dimension
from qdjj_solver.qp_solver.solver import ProjectedOperator, SolverOptions, _two_qp_schur, solve


def enumerate_occupations(nimp, spins, cutoff, sector, eta=None):
    result = {}
    for occupied in itertools.product((0, 1), repeat=len(spins)):
        particles = sum(occupied)
        q = sum(occupied[nimp:])
        if q > cutoff or particles % 2 != sector.parity:
            continue
        if sector.twice_sz is not None and sum(s * n for s, n in zip(spins, occupied, strict=True)) != sector.twice_sz:
            continue
        if sector.eta is not None and np.prod([e for e, n in zip(eta, occupied, strict=True) if n]) != sector.eta:
            continue
        if sector.particle_number is not None and particles != sector.particle_number:
            continue
        result[sum(n << i for i, n in enumerate(occupied))] = q
    return result


def states(basis):
    return [sum(int(word) << (64 * i) for i, word in enumerate(row)) for row in basis.occupations()]


def direct_matrix(operator, source, destination=None):
    """Bit-string CAR oracle: no normal-ordering, native action, or sparse assembly."""
    destination = source if destination is None else destination
    target_rows = {state: row for row, state in enumerate(states(destination))}
    matrix = np.zeros((destination.dimension, source.dimension), dtype=complex)
    for column, occupation in enumerate(states(source)):
        for string, coefficient in operator.terms.items():
            ket, sign = occupation, 1
            for op in reversed(string):
                mask = 1 << (abs(op) - 1)
                if bool(ket & mask) == (op > 0):
                    break
                sign *= (-1) ** (ket & (mask - 1)).bit_count()
                ket ^= mask
            else:
                if ket in target_rows:
                    matrix[target_rows[ket], column] += coefficient * sign
    return matrix


@pytest.mark.native
@pytest.mark.parametrize("spins", [(1, -1) * 3, (1,) * 6, (-1, 1, 1, -1, -1, 1)])
@pytest.mark.parametrize("nimp", [0, 2, 6])
def test_dimensions_and_basis_against_independent_enumeration(spins, nimp):
    labels = (-1, 1, -1, 1, 1, -1)
    sectors = [Sector(0, None), Sector(1, None), Sector(0, 0), Sector(1, -1),
               Sector(0, 2, 1), Sector(0, None, -1), Sector(1, 1, -1),
               Sector(0, None, particle_number=0), Sector(0, 0, particle_number=2),
               Sector(1, None, 1, particle_number=3), Sector(0, None, particle_number=8)]
    for cutoff, sector in itertools.product(range(len(spins) - nimp + 1), sectors):
        expected = enumerate_occupations(nimp, spins, cutoff, sector, labels)
        assert dimension(nimp, spins, np.int64(cutoff), sector, labels) == len(expected)
        if expected:
            basis = FockBasis(nimp, spins, cutoff, sector, labels)
            assert dict(zip(states(basis), basis.qp_counts(), strict=True)) == expected
            assert basis.dimension == len(expected)
            assert basis.bytes >= basis.occupations().nbytes + basis.qp_counts().nbytes
        else:
            with pytest.raises(ValueError, match="empty symmetry sector"):
                FockBasis(nimp, spins, cutoff, sector, labels)


@pytest.mark.parametrize("nimp, spins, cutoff, sector, labels, message", [
    (-1, (1, -1), 0, Sector(), None, "impurity"),
    (3, (1, -1), 0, Sector(), None, "impurity"),
    (0.5, (1, -1), 0, Sector(), None, "impurity"),
    (0, (1, 0), 0, Sector(), None, "spin labels"),
    (0, (1, -1), -1, Sector(), None, "cutoff"),
    (1, (1, -1), 2, Sector(), None, "cutoff"),
    (0, (1, -1), 0.5, Sector(), None, "cutoff"),
    (0, (1, -1), "1", Sector(), None, "cutoff"),
    (0, (1, -1), 0, Sector(0, None, 1), None, "no diagonal eta"),
    (0, (1, -1), 0, Sector(0, None, 1), (1,), "eta labels"),
    (0, (1, -1), 0, Sector(0, None), (1, 0), "eta labels"),
])
def test_dimension_input_validation(nimp, spins, cutoff, sector, labels, message):
    with pytest.raises(ValueError, match=message):
        dimension(nimp, spins, cutoff, sector, labels)


@pytest.mark.native
def test_basis_iterable_labels_are_materialized_once():
    basis = FockBasis(np.int64(2), iter((1, -1, 1, -1)), 2,
                      Sector(0, None, -1), iter((1, -1, 1, -1)))
    expected = enumerate_occupations(2, (1, -1, 1, -1), 2, Sector(0, None, -1), (1, -1, 1, -1))
    assert set(states(basis)) == set(expected)


def test_arbitrary_precision_dimension_and_preallocation_guard(monkeypatch):
    assert dimension(100, (1,) * 100, 0, Sector(0, None)) == 2**99
    def forbidden(*args, **kwargs):
        pytest.fail("native enumeration must not run beyond the dimension limit")
    monkeypatch.setattr(basis_module._core, "Basis", forbidden)
    with pytest.raises(MemoryError, match="max_dimension=8"):
        FockBasis(100, (1,) * 100, 0, Sector(0, None), max_dimension=8)


def test_basis_counting_disagreement_is_controlled(monkeypatch):
    monkeypatch.setattr(basis_module._core, "Basis", lambda *args: SimpleNamespace(dimension=999))
    with pytest.raises(RuntimeError, match="counting and enumeration disagree"):
        FockBasis(1, (1, -1), 1, Sector(0, None))


@pytest.mark.native
@pytest.mark.parametrize("complex_operator", [False, True])
def test_projected_action_sparse_linear_diagonal_and_expectations(complex_operator):
    basis = FockBasis(2, (1, -1, 1, -1), 1, Sector(0, None))
    operator = FermionOperator({(): .3, (1, -1): -.7, (1, 3): .2,
                               (1, -4): .4j if complex_operator else .4,
                               (1, 2, -2, -1): .9})
    projected = ProjectedOperator(operator, basis)
    matrix = direct_matrix(operator, basis)
    vector = np.arange(basis.dimension, dtype=float) / 3
    ket = vector + 1j * vector[::-1]
    assert_allclose(projected.action(vector), matrix @ vector, rtol=0, atol=1e-14)
    assert_allclose(projected.action(ket), matrix @ ket, rtol=0, atol=1e-14)
    assert projected.dtype == np.dtype(complex if complex_operator else float)
    assert projected.sparse_matrix().dtype == projected.dtype
    assert_allclose(projected.sparse_matrix().toarray(), matrix, rtol=0, atol=1e-14)
    assert_allclose(projected.diagonal(), np.diag(matrix), rtol=0, atol=1e-14)
    assert_allclose(projected.expectation(ket), np.vdot(ket, matrix @ ket), rtol=0, atol=1e-13)
    assert_allclose(projected.expectation(vector, ket), np.vdot(vector, matrix @ ket), rtol=0, atol=1e-13)
    linear = projected.linear_operator()
    assert linear.dtype == projected.dtype
    assert_allclose(linear @ ket, matrix @ ket, rtol=0, atol=1e-14)
    assert_allclose(linear @ vector[:, None], (matrix @ vector)[:, None], rtol=0, atol=1e-14)
    assert_allclose(projected.linear_operator(complex) @ np.eye(basis.dimension), matrix, rtol=0, atol=1e-14)
    # Pure imaginary diagonals must retain their complex dtype.
    imag = ProjectedOperator(1j * number(0), basis).diagonal()
    assert np.iscomplexobj(imag)
    assert_allclose(imag, np.diag(direct_matrix(1j * number(0), basis)), rtol=0, atol=0)


@pytest.mark.native
@pytest.mark.parametrize("coefficient", [1., .3 + .4j])
def test_cross_sector_and_cutoff_transitions(coefficient):
    source = FockBasis(2, (1, -1, 1, -1), 1, Sector(0, 0, 1), (1, -1, 1, -1))
    destination = FockBasis(2, source.spins, 2, Sector(1, 1, -1), (-1, 1, -1, 1))
    operator = coefficient * (create(0) + create(2))
    projected = ProjectedOperator(operator, source)
    expected = direct_matrix(operator, source, destination)
    assert np.count_nonzero(expected) > 0
    ket = np.arange(1, source.dimension + 1) * (1 - .2j)
    bra = np.arange(1, destination.dimension + 1) * (.1 + .4j)
    assert_allclose(projected.action_between(ket, destination), expected @ ket, rtol=0, atol=1e-14)
    assert_allclose(projected.transition(bra, destination, ket), np.vdot(bra, expected @ ket), rtol=0, atol=1e-13)
    # Projection into the original parity sector is zero, unlike the transition.
    assert_allclose(projected.action(ket), 0, rtol=0, atol=0)


@pytest.mark.native
def test_projected_operator_input_errors():
    basis = FockBasis(2, (1, -1, 1, -1), 1, Sector(0, None))
    with pytest.raises(ValueError, match="mode outside"):
        ProjectedOperator(create(4), basis)
    projected = ProjectedOperator(number(0), basis)
    vector = np.ones(basis.dimension)
    for invalid in (vector[:-1], vector[:, None]):
        with pytest.raises(ValueError, match="state-vector dimension"):
            projected.action(invalid)
        with pytest.raises(ValueError, match="transition vectors"):
            projected.action_between(invalid, basis)
        with pytest.raises(ValueError, match="transition vectors"):
            projected.transition(invalid, basis, vector)
    for other in (FockBasis(1, basis.spins, 1, Sector(0, None)),
                  FockBasis(2, (-1, 1, 1, -1), 1, Sector(0, None))):
        with pytest.raises(ValueError, match="transition vectors"):
            projected.action_between(vector, other)
        with pytest.raises(ValueError, match="transition vectors"):
            projected.transition(np.ones(other.dimension), other, vector)


@pytest.mark.parametrize("method", ["auto", "dense", "sparse", "matrix-free", "schur"])
def test_solver_options_valid_methods(method):
    assert SolverOptions(method=method).method == method


def test_solver_options_invalid_method():
    with pytest.raises(ValueError, match="method must"):
        SolverOptions(method="lanczos")


@pytest.mark.parametrize("name, value", [
    ("eigenpairs", 0), ("eigenpairs", 1.5), ("maxiter", -1), ("ncv", 0),
    ("max_dimension", 0), ("max_nnz", -1), ("threads", 0), ("threads", "2"),
])
def test_solver_options_positive_integers(name, value):
    with pytest.raises(ValueError, match=name):
        SolverOptions(**{name: value})


@pytest.mark.parametrize("name, value", [
    ("tolerance", 0.), ("tolerance", np.nan),
    ("residual_tolerance", -1.), ("residual_tolerance", np.inf),
    ("max_memory_gib", np.nan), ("max_memory_gib", 0.),
])
def test_solver_options_positive_finite_numbers(name, value):
    with pytest.raises(ValueError, match=name):
        SolverOptions(**{name: value})


def one_particle_problem(complex_operator=False, nmodes=6):
    diagonal = np.linspace(-.8, 1.4, nmodes)
    matrix = np.diag(diagonal).astype(complex if complex_operator else float)
    for i in range(nmodes - 1):
        matrix[i, i + 1] = .23 + (.17j if complex_operator else 0)
        matrix[i + 1, i] = matrix[i, i + 1].conjugate()
    terms = {(i + 1, -(j + 1)): matrix[i, j] for i in range(nmodes) for j in range(nmodes) if matrix[i, j]}
    h = Hamiltonian(FermionOperator(terms), 2, (1,) * nmodes,
                    observables={"first": number(0), "imaginary": 1j * number(1), "zero": FermionOperator()},
                    metadata={"test_label": "analytic one-particle chain"})
    return h, matrix, Sector(1, 1, particle_number=1)


@pytest.mark.native
@pytest.mark.parametrize("method, complex_operator, initial_kind", [
    ("dense", False, None), ("dense", True, None),
    ("sparse", False, None), ("sparse", True, None),
    ("matrix-free", False, None), ("matrix-free", True, None),
    ("sparse", False, "real"), ("sparse", False, "complex"),
    ("matrix-free", True, "complex"),
])
def test_true_solver_dispatch_and_analytic_one_particle_spectrum(monkeypatch, method, complex_operator, initial_kind):
    h, analytic, sector = one_particle_problem(complex_operator)
    expected_energies = np.linalg.eigvalsh(analytic)[:2]
    calls = []
    actual_eigsh = solver_module.eigsh
    def recording_eigsh(matrix, **kwargs):
        calls.append(matrix)
        assert (isinstance(matrix, LinearOperator) if method == "matrix-free" else issparse(matrix))
        return actual_eigsh(matrix, **kwargs)
    monkeypatch.setattr(solver_module, "eigsh", recording_eigsh)
    if method == "matrix-free":
        def forbidden_sparse(*args, **kwargs):
            pytest.fail("matrix-free dispatch must not materialize a sparse matrix")
        monkeypatch.setattr(ProjectedOperator, "sparse_matrix", forbidden_sparse)
    initial = None if initial_kind is None else np.linspace(.4, 1.2, 6)
    if initial_kind == "complex":
        initial = initial + 1j * initial[::-1]
    result = solve(h, cutoff=1, sector=sector,
                   options=SolverOptions(method=method, eigenpairs=2, ncv=3), initial=initial)
    assert len(calls) == (method != "dense")
    assert result.metadata["method"] == method
    assert result.metadata["real_arithmetic"] is (not complex_operator and initial_kind != "complex")
    assert result.metadata["test_label"] == h.metadata["test_label"]
    assert result.hamiltonian is h
    assert_allclose(result.energies, expected_energies, rtol=0, atol=2e-13)
    mode_order = [state.bit_length() - 1 for state in states(result.basis)]
    matrix = analytic[np.ix_(mode_order, mode_order)]
    residual = matrix @ result.vectors - result.vectors * result.energies
    assert_allclose(np.linalg.norm(residual, axis=0), result.residuals, rtol=0, atol=1e-15)
    assert np.max(result.residuals) < 1e-12
    assert_allclose(result.vectors.conj().T @ result.vectors, np.eye(2), rtol=0, atol=1e-13)
    counts = np.array([int(mode >= 2) for mode in mode_order])
    weights = np.array([[sum(abs(v[counts == q])**2) for q in (0, 1)] for v in result.vectors.T])
    assert_allclose(result.qp_weights, weights, rtol=0, atol=1e-14)
    assert_allclose(result.qp_weights.sum(axis=1), 1., rtol=0, atol=1e-14)
    assert np.all(result.qp_weights >= 0)
    assert_allclose(result.observables["first"], abs(result.vectors[mode_order.index(0)])**2, rtol=0, atol=1e-14)
    assert_allclose(result.observables["imaginary"], 1j * abs(result.vectors[mode_order.index(1)])**2, rtol=0, atol=1e-14)
    assert_allclose(result.observables["zero"], 0., rtol=0, atol=0)
    assert result.metadata["nnz"] == (np.count_nonzero(matrix) if method == "sparse" else None)


@pytest.mark.native
@pytest.mark.parametrize("nmodes, expected_method", [(8, "dense"), (9, "sparse")])
def test_auto_dispatch_with_actual_eigensolutions(nmodes, expected_method):
    operator = .125 + sum((i + 1) * number(i) for i in range(nmodes))
    h = Hamiltonian(operator, 0, (1,) * nmodes)
    result = solve(h, sector=Sector(0, None), options=SolverOptions(eigenpairs=2))
    assert result.metadata["method"] == expected_method
    assert_allclose(result.energies, [.125, 3.125], rtol=0, atol=2e-13)
    assert result.qp_weights.shape == (2, nmodes + 1)
    assert_allclose(result.qp_weights.sum(axis=1), 1., rtol=0, atol=2e-14)
    assert_allclose(result.qp_weights[0, 0], 1., rtol=0, atol=1e-13)
    assert_allclose(result.qp_weights[1, 2], 1., rtol=0, atol=1e-13)


@pytest.mark.native
@pytest.mark.parametrize("requested", ["sparse", "matrix-free"])
@pytest.mark.parametrize("k", [1, 2])
def test_small_problem_dense_fallback_is_explicit(requested, k):
    h = Hamiltonian(number(0) + 2 * number(1), 0, (1, -1))
    result = solve(h, sector=Sector(1, None), options=SolverOptions(method=requested, eigenpairs=k))
    assert result.metadata["method"] == "dense"
    assert_allclose(result.energies, [1., 2.][:k], rtol=0, atol=1e-14)


def test_solve_preallocation_errors(monkeypatch):
    h, _, sector = one_particle_problem()
    def forbidden(*args, **kwargs):
        pytest.fail("basis allocation must follow input/resource validation")
    monkeypatch.setattr(solver_module, "FockBasis", forbidden)
    for options, selected, error, message in [
        (SolverOptions(), Sector(0, -2), ValueError, "empty symmetry"),
        (SolverOptions(max_dimension=5), sector, MemoryError, "permitted maximum"),
        (SolverOptions(eigenpairs=7), sector, ValueError, "more eigenpairs"),
        (SolverOptions(max_memory_gib=1e-12), sector, MemoryError, "estimated basis/eigensolver"),
        (SolverOptions(method="schur"), sector, ValueError, "QP cutoff 2"),
    ]:
        with pytest.raises(error, match=message):
            solve(h, cutoff=1, sector=selected, options=options)
    with pytest.raises(ValueError, match="lowest state"):
        solve(h, cutoff=2, sector=sector, options=SolverOptions(method="schur", eigenpairs=2))


@pytest.mark.native
def test_dense_memory_guard_precedes_matrix_materialization(monkeypatch):
    h = Hamiltonian(FermionOperator(), 0, (1,) * 8)
    def forbidden(*args, **kwargs):
        pytest.fail("dense guard must run before matrix allocation")
    monkeypatch.setattr(ProjectedOperator, "sparse_matrix", forbidden)
    with pytest.raises(MemoryError, match="dense eigensolution"):
        solve(h, sector=Sector(0, None), options=SolverOptions(method="dense", max_memory_gib=200_000 / 1024**3))


@pytest.mark.native
@pytest.mark.parametrize("nimp, limit", [(8, 0.0005), (16, 0.125)])
@pytest.mark.parametrize("complex_operator", [False, True])
def test_schur_memory_guard_precedes_basis_allocation(monkeypatch, nimp, limit, complex_operator):
    operator = sum(number(i) for i in range(nimp + 2))
    if complex_operator:
        hopping = 1j * create(0) * annihilate(nimp)
        operator += hopping + hopping.dagger()
    h = Hamiltonian(operator, nimp, (1, -1) * (nimp // 2 + 1))
    def forbidden(*args, **kwargs):
        pytest.fail("Schur storage must be checked before allocating the basis or matrix")
    monkeypatch.setattr(solver_module, "FockBasis", forbidden)
    with pytest.raises(MemoryError, match="estimated basis/eigensolver"):
        solve(h, cutoff=2, sector=Sector(0, None),
              options=SolverOptions(method="schur", max_memory_gib=limit))


@pytest.mark.native
def test_sparse_storage_limits_without_large_allocations():
    h, _, sector = one_particle_problem()
    with pytest.raises(ValueError, match="max_nnz"):
        solve(h, sector=sector, options=SolverOptions(method="sparse", max_nnz=1))
    # The six-dimensional basis/Krylov estimate is 1296 bytes for k=1, ncv=6.
    estimate = 6 * (8 + 40 + 8 * (2 * 6 + 1 + 8))
    with pytest.raises(ValueError, match="max_nnz"):
        solve(h, sector=sector, options=SolverOptions(method="sparse", max_memory_gib=(estimate + 32) / 1024**3))


@pytest.mark.native
@pytest.mark.parametrize("method, initial", [
    ("sparse", np.zeros(6)), ("matrix-free", np.ones(5)), ("sparse", np.ones((6, 1))),
    ("matrix-free", np.array([1., 1., np.nan, 1., 1., 1.])),
    ("sparse", np.array([1., 1., np.inf, 1., 1., 1.])),
    ("matrix-free", np.full(6, complex(0, np.inf))),
])
def test_invalid_iterative_initial_vectors(method, initial):
    h, _, sector = one_particle_problem()
    with pytest.raises(ValueError, match="invalid starting vector"):
        solve(h, sector=sector, options=SolverOptions(method=method), initial=initial)


@pytest.mark.native
def test_arpack_nonconvergence_preserves_cause_and_partial_count(monkeypatch):
    h, _, sector = one_particle_problem()
    error = ArpackNoConvergence("injected iteration exhaustion", np.array([0.]), np.ones((6, 1)))
    def fail(*args, **kwargs):
        raise error
    monkeypatch.setattr(solver_module, "eigsh", fail)
    with pytest.raises(RuntimeError, match=r"did not converge \(1/2 eigenpairs\)") as captured:
        solve(h, sector=sector, options=SolverOptions(method="sparse", eigenpairs=2))
    assert captured.value.__cause__ is error


@pytest.mark.native
def test_iterative_eigenpairs_are_sorted_together(monkeypatch):
    h, _, sector = one_particle_problem()
    basis = FockBasis(h.nimp, h.spins, h.bath_modes, sector)
    energies, vectors = np.linalg.eigh(direct_matrix(h.operator, basis))
    def unordered(*args, **kwargs):
        return energies[[1, 0]], vectors[:, [1, 0]]
    monkeypatch.setattr(solver_module, "eigsh", unordered)
    result = solve(h, sector=sector, options=SolverOptions(method="sparse", eigenpairs=2))
    assert_allclose(result.energies, energies[:2], rtol=0, atol=2e-14)
    assert_allclose(abs(vectors[:, :2].conj().T @ result.vectors), np.eye(2), rtol=0, atol=2e-14)
    assert np.max(result.residuals) < 1e-13


@pytest.mark.native
@pytest.mark.parametrize("method,nmodes", [("dense", 6), ("sparse", 6), ("matrix-free", 6), ("auto", 130)])
@pytest.mark.parametrize("splitting", [0., 1e-10, 1e-6])
def test_complex_degenerate_and_near_degenerate_ritz_subspaces(method, nmodes, splitting):
    analytic = np.diag(2. + np.arange(nmodes)/nmodes).astype(complex)
    analytic[:4, :4] = 0.
    analytic[0, 1] = analytic[2, 3] = 1j
    analytic[1, 0] = analytic[3, 2] = -1j
    analytic[2, 2] = analytic[3, 3] = splitting
    terms = {(i+1, -(j+1)): analytic[i, j] for i, j in zip(*np.nonzero(analytic), strict=True)}
    h = Hamiltonian(FermionOperator(terms), 0, (1, -1)*(nmodes//2))
    result = solve(h, cutoff=1, sector=Sector(1, None, particle_number=1),
                   options=SolverOptions(method=method, eigenpairs=2))
    expected, expected_vectors = np.linalg.eigh(analytic)
    mode_order = [state.bit_length()-1 for state in states(result.basis)]
    vectors = result.vectors[np.argsort(mode_order)]
    assert result.metadata["method"] == ("sparse" if method == "auto" else method)
    assert_allclose(result.energies, expected[:2], rtol=0, atol=2e-13)
    assert_allclose(vectors.conj().T @ vectors, np.eye(2), rtol=0, atol=2e-14)
    # Compare the entire retained subspace, independently of rotations/phases.
    singular_values = np.linalg.svd(expected_vectors[:, :2].conj().T @ vectors, compute_uv=False)
    assert_allclose(singular_values, 1., rtol=0, atol=2e-14)
    assert_allclose(analytic @ vectors, vectors*result.energies, rtol=0, atol=2e-13)
    assert np.max(result.residuals) < 2e-13


@pytest.mark.native
def test_iterative_rank_deficiency_cannot_be_repaired_into_extra_roots(monkeypatch):
    h, _, sector = one_particle_problem(complex_operator=True)
    basis = FockBasis(h.nimp, h.spins, h.bath_modes, sector)
    energies, vectors = np.linalg.eigh(direct_matrix(h.operator, basis))
    monkeypatch.setattr(solver_module, "eigsh", lambda *a, **kw: (energies[[0, 0]], vectors[:, [0, 0]]))
    with pytest.raises(RuntimeError, match="linearly dependent"):
        solve(h, sector=sector, options=SolverOptions(method="sparse", eigenpairs=2))


@pytest.mark.native
@pytest.mark.parametrize("fault, message", [
    ("energy_nan", "non-finite"), ("energy_inf", "non-finite"),
    ("vector_nan", "non-finite"), ("vector_inf", "non-finite"),
    ("residual", "Ritz residual"), ("residual_nonfinite", "Ritz residual"),
    ("nonorthogonal", "not orthonormal"), ("unnormalized", "not orthonormal"),
    ("wrong_rows", "invalid dimensions"), ("missing_pairs", "invalid dimensions"),
])
def test_corrupt_eigensystems_are_rejected(monkeypatch, fault, message):
    h = Hamiltonian(FermionOperator(), 0, (1,) * 4)
    energies, vectors = np.zeros(4), np.eye(4)
    if fault.startswith("energy_"):
        energies[0] = np.nan if fault.endswith("nan") else np.inf
    elif fault.startswith("vector_"):
        vectors[0, 0] = np.nan if fault.endswith("nan") else np.inf
    elif fault == "residual":
        energies[0] = 1.
    elif fault == "residual_nonfinite":
        monkeypatch.setattr(ProjectedOperator, "action", lambda self, v: np.full_like(v, np.inf))
    elif fault == "nonorthogonal":
        vectors[:, 1] = vectors[:, 0]
    elif fault == "unnormalized":
        vectors[:, 0] *= 2
    elif fault == "wrong_rows":
        vectors = vectors[:-1]
    elif fault == "missing_pairs":
        energies, vectors = energies[:1], vectors[:, :1]
    monkeypatch.setattr(np.linalg, "eigh", lambda matrix: (energies, vectors))
    with pytest.raises(RuntimeError, match=message):
        solve(h, sector=Sector(1, None, particle_number=1),
              options=SolverOptions(method="dense", eigenpairs=2))


def schur_basis(counts, cutoff=2):
    return SimpleNamespace(cutoff=cutoff, dimension=len(counts), qp_counts=lambda: np.array(counts))


@pytest.mark.parametrize("complex_blocks", [False, True])
def test_schur_analytic_blocks_and_full_residual(complex_blocks):
    ll = np.array([[-.4, .13], [.13, .7]], dtype=complex if complex_blocks else float)
    coupling = np.array([[.3, .1], [-.2, .4]], dtype=ll.dtype)
    if complex_blocks:
        ll[0, 1] += .2j
        ll[1, 0] -= .2j
        coupling[0, 1] += .17j
        coupling[1, 0] -= .21j
    diagonal = np.array([1.7, 2.4])
    matrix = np.block([[ll, coupling.conj().T], [coupling, np.diag(diagonal)]])
    energy, vectors, metadata = _two_qp_schur(csc_matrix(matrix), schur_basis([0, 1, 2, 2]), 1e-13)
    expected, expected_vectors = np.linalg.eigh(matrix)
    assert_allclose(energy, expected[:1], rtol=0, atol=1e-14)
    assert_allclose(abs(np.vdot(vectors[:, 0], expected_vectors[:, 0])), 1., rtol=0, atol=1e-14)
    assert_allclose(matrix @ vectors, vectors * energy, rtol=0, atol=1e-13)
    assert_allclose(vectors.conj().T @ vectors, [[1]], rtol=0, atol=1e-14)
    low, high = vectors[:2, 0], vectors[2:, 0]
    assert_allclose(high, coupling @ low / (energy[0] - diagonal), rtol=0, atol=1e-14)
    effective = ll + coupling.conj().T @ (coupling / (energy[0] - diagonal)[:, None])
    assert_allclose(effective @ low, energy[0] * low, rtol=0, atol=1e-13)
    assert metadata["schur_dimension"] == 2
    assert metadata["schur_evaluations"] >= 3


@pytest.mark.parametrize("pole, phase, strength, tolerance", [
    (0., 1., 1e-6, 1e-11), (0., 1j, 1e-9, 1e-20),
    (1., 1., 1e-6, 1e-11), (-1., 1j, 1e-6, 1e-20),
    (1., 1j, 1e-9, 1e-11), (-1., 1., 1e-9, 1e-11),
])
def test_schur_valid_root_closer_than_old_absolute_pole_margin(pole, phase, strength, tolerance):
    coupling = strength * phase
    matrix = np.array([[pole + 1., np.conj(coupling)], [coupling, pole]])
    energy, vectors, _ = _two_qp_schur(csc_matrix(matrix), schur_basis([0, 2]), tolerance)
    # Stable analytic two-level root; both shifts were formerly excluded.
    shift = 2 * abs(coupling)**2 / (1 + np.hypot(1., 2 * abs(coupling)))
    assert 0 <= pole - energy[0] < 1e-10
    if pole == 0:
        assert_allclose(energy, [-shift], rtol=2e-14, atol=0)
    assert_allclose(energy, [pole - shift], atol=5e-16, rtol=0)
    assert_allclose(vectors[0, 0] / vectors[1, 0], -np.conj(coupling) / (1 + shift), rtol=2e-14, atol=0)
    assert_allclose(matrix @ vectors, vectors * energy, atol=5e-16, rtol=0)
    assert_allclose(np.linalg.norm(vectors), 1., rtol=0, atol=1e-14)


@pytest.mark.parametrize("counts, cutoff, matrix, message", [
    ([0, 2], 1, np.eye(2), "cutoff 2"),
    ([0, 1], 2, np.eye(2), "no two-QP"),
    ([2, 2], 2, np.eye(2), "below two QPs"),
    ([0, 2, 2], 2, np.array([[0., .1, .2], [.1, 1., .3], [.2, .3, 2.]]), "diagonal two-QP"),
    ([0, 2], 2, np.diag([2., 0.]), "regular Schur branch"),
    ([0, 2], 2, np.diag([0., 1. + 1j]), "diagonal two-QP"),
])
def test_schur_invalid_blocks_fail_cleanly(counts, cutoff, matrix, message):
    with pytest.raises(ValueError, match=message):
        _two_qp_schur(csc_matrix(matrix), schur_basis(counts, cutoff), 1e-12)


def test_schur_tolerates_only_roundoff_in_high_block():
    matrix = np.array([[-1., .1, .2], [.1, 1., 1e-15], [.2, 1e-15, 2.]])
    energy, vector, _ = _two_qp_schur(csc_matrix(matrix), schur_basis([0, 2, 2]), 1e-12)
    assert_allclose(matrix @ vector, vector * energy, rtol=0, atol=1e-13)


@pytest.mark.native
def test_solve_schur_empty_low_sector_is_a_controlled_failure():
    h = Hamiltonian(number(0) + 2 * number(1) + 3 * number(2), 0, (1,) * 3)
    with pytest.raises(ValueError, match="below two QPs"):
        solve(h, cutoff=2, sector=Sector(0, None, particle_number=2), options=SolverOptions(method="schur"))


@pytest.mark.native
@pytest.mark.parametrize("complex_operator", [False, True])
def test_solve_schur_true_dispatch_and_qp_probabilities(complex_operator):
    pairing = FermionOperator({(1, 2): .2, (1, 3): .17j if complex_operator else .17, (2, 3): -.13})
    operator = .1 + number(0) + 1.3 * number(1) + 1.7 * number(2) + pairing + pairing.dagger()
    h = Hamiltonian(operator, 0, (1,) * 3, observables={"number": sum(number(i) for i in range(3))})
    result = solve(h, cutoff=2, sector=Sector(0, None), options=SolverOptions(method="schur"))
    matrix = direct_matrix(operator, result.basis)
    assert result.metadata["method"] == "schur"
    assert result.metadata["schur_dimension"] == 1
    assert_allclose(result.energies, np.linalg.eigvalsh(matrix)[:1], rtol=0, atol=1e-13)
    assert_allclose(matrix @ result.vectors, result.vectors * result.energies, rtol=0, atol=1e-12)
    assert_allclose(result.qp_weights.sum(axis=1), 1., rtol=0, atol=1e-14)
    assert_allclose(result.qp_weights[:, 1], 0., rtol=0, atol=0)
    assert_allclose(result.observables["number"], 2 * result.qp_weights[:, 2], rtol=0, atol=1e-14)
