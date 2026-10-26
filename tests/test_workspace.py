"""Prepared execution checked against fresh solves and independent matrix spectra."""

from dataclasses import replace
import gc

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import Sector, cosh_grid, reference_model, orbital_density_matrix
from qdjj_solver.common.algebra import create, annihilate
from qdjj_solver.common.problem import Hamiltonian
from qdjj_solver.qp_solver import CalculationCancelled, ExecutionContext, QPWorkspace, SolverOptions, solve
from qdjj_solver.qp_solver.solver import preflight
from qdjj_solver.qp_solver.workspace import hermitian_inertia


def junction(phi=0.8, delta=1., vacuum="isolated"):
    return reference_model(cosh_grid(1, delta, 5.), u=1.6, gamma=.3, phi=phi,
                           detuning=.11, field=.13, rho_wn=.02, symmetry=False,
                           compress=False, bath_reference=vacuum)


def one_particle(matrix):
    n = len(matrix)
    op = sum(matrix[i, j]*create(i)*annihilate(j) for i in range(n) for j in range(n))
    return Hamiltonian(op, 0, (1,)*n), Sector(1, 1, particle_number=1)


@pytest.mark.parametrize("method", ["dense", "sparse", "matrix-free", "block"])
def test_workspace_reuse_matches_fresh_and_preserves_snapshots(method):
    w = QPWorkspace()
    h = junction()
    sector = Sector(0, 0)
    options = SolverOptions(method=method, eigenpairs=2, maxiter=2000)
    fresh = solve(h, 2, sector, options)
    first = w.solve(h, 2, sector, options)
    old = first.vectors.copy()
    again = w.solve(h, 2, sector, options)
    changed = w.solve(junction(phi=.9), 2, sector, options)
    assert_allclose(first.energies, fresh.energies, atol=1e-9, rtol=0)
    assert_allclose(again.energies, fresh.energies, atol=1e-9, rtol=0)
    assert_allclose(changed.energies, solve(junction(phi=.9), 2, sector, options).energies, atol=1e-9, rtol=0)
    assert_allclose(old, first.vectors, atol=0, rtol=0)
    assert first.basis is again.basis is changed.basis
    assert not first.metadata["reuse"]["basis"]
    assert again.metadata["reuse"]["basis"]
    assert again.metadata["reuse"]["projected"]
    assert w.stats["basis_builds"] == 1
    if method == "dense":
        assert again.metadata["reuse"]["matrix"] and again.metadata["reuse"]["eigensystem"]
    else:
        assert w.stats["warm_starts"] == 2
    assert not first.vectors.flags.writeable
    assert w.live_bytes >= sum(r.vectors.nbytes for r in (first, changed))


@pytest.mark.parametrize("method", ["sparse", "matrix-free"])
@pytest.mark.parametrize("complex_operator", [False, True])
def test_workspace_shifted_arpack_keeps_hidden_zero_with_warm_start(method, complex_operator):
    # The isolated third mode has exactly zero energy between the +/-1 pair.
    # A small residual alone does not detect omission of this state.
    matrix = np.diag(np.r_[0., 0., 0., np.arange(1., 38.)]).astype(complex if complex_operator else float)
    matrix[0, 1] = 1j if complex_operator else 1.
    matrix[1, 0] = matrix[0, 1].conjugate()
    h, sector = one_particle(matrix)
    options = SolverOptions(method=method, eigenpairs=2)
    workspace = QPWorkspace()
    for _ in range(2):
        result = workspace.solve(h, 1, sector, options)
        assert_allclose(result.energies, [-1., 0.], rtol=0, atol=1e-12)
        assert_allclose(result.vectors.conj().T@result.vectors, np.eye(2), rtol=0, atol=1e-12)
        assert result.metadata["spectral_shift"] == pytest.approx(-2**-.5)
    assert workspace.stats["warm_starts"] == 1
    assert result.metadata["reuse"]["projected"]


def test_workspace_snapshots_support_transitions_and_physical_density():
    h = junction()
    options = SolverOptions(method="dense", eigenpairs=2)
    workspace = QPWorkspace()
    even = workspace.solve(h, 2, Sector(0, 0), options)
    odd = workspace.solve(h, 2, Sector(1, 1), options)
    reference_even = solve(h, 2, Sector(0, 0), options)
    reference_odd = solve(h, 2, Sector(1, 1), options)
    assert_allclose(abs(odd.matrix_elements("d_up_dagger", even))**2,
                    abs(reference_odd.matrix_elements("d_up_dagger", reference_even))**2, rtol=0, atol=1e-12)
    assert_allclose(even.overlaps(), np.eye(2), rtol=0, atol=1e-12)
    assert_allclose(orbital_density_matrix(even, state=None),
                    orbital_density_matrix(reference_even, state=None), rtol=0, atol=1e-12)
    assert not even.vectors.flags.writeable and not odd.vectors.flags.writeable


def test_block_solver_keeps_exact_complex_multiplicity():
    rng = np.random.default_rng(93)
    q, _ = np.linalg.qr(rng.normal(size=(24, 24))+1j*rng.normal(size=(24, 24)))
    expected = np.r_[[-2., -2., -2.], np.arange(21.)]
    matrix = (q*expected)@q.conj().T
    h, sector = one_particle(matrix)
    opts = SolverOptions(method="block", eigenpairs=3, maxiter=1000)
    result = solve(h, 1, sector, opts)
    assert_allclose(result.energies, [-2]*3, atol=1e-9, rtol=0)
    assert_allclose(result.vectors.conj().T@result.vectors, np.eye(3), atol=1e-10)
    assert_allclose(result.vectors@result.vectors.conj().T, q[:, :3]@q[:, :3].conj().T, atol=1e-8)
    reused = solve(h, 1, sector, opts, initial=result.vectors[:, :2])
    assert_allclose(reused.energies, result.energies, atol=1e-9)


@pytest.mark.parametrize("initial", [np.zeros((24, 2)), np.ones((23, 2)), np.full((24, 2), np.nan), np.ones((24, 1, 1))])
def test_invalid_starting_blocks(initial):
    h, sector = one_particle(np.diag(np.arange(24.)))
    with pytest.raises(ValueError, match="starting block"):
        solve(h, 1, sector, SolverOptions(method="block", eigenpairs=2), initial)


def test_block_preflight_prevents_implicit_large_dense_fallback():
    h, sector = one_particle(np.diag(np.arange(24.)))
    with pytest.raises(ValueError, match=r"5\*eigenpairs"):
        preflight(h, 1, sector, SolverOptions(method="block", eigenpairs=6))


def test_coordinate_change_prevents_warm_start_but_keeps_structural_basis():
    w = QPWorkspace()
    opts = SolverOptions(method="matrix-free", eigenpairs=1)
    first = w.solve(junction(delta=1), 2, options=opts)
    second = w.solve(junction(delta=1.2), 2, options=opts)
    assert first.basis is second.basis
    assert w.stats["warm_starts"] == 0
    w.solve(junction(phi=.9, delta=1.2), 2, options=opts)
    assert w.stats["warm_starts"] == 1
    w.solve(junction(vacuum="coupled"), 2, options=opts)
    w.solve(junction(phi=.9, vacuum="coupled"), 2, options=opts)
    assert w.stats["warm_starts"] == 1


def test_observable_selection_and_cached_dense_growth():
    w = QPWorkspace()
    h = junction()
    first = w.solve(h, 2, options=SolverOptions(method="dense"), observables=())
    assert not first.observables
    expanded = w.solve(h, 2, options=SolverOptions(method="dense", eigenpairs=4), observables=["impurity_charge"])
    assert expanded.metadata["reuse"]["eigensystem"]
    assert set(expanded.observables) == {"impurity_charge"}
    again = w.solve(h, 2, options=SolverOptions(method="dense", eigenpairs=4), observables=["impurity_charge"])
    assert_allclose(again.observables["impurity_charge"], expanded.observables["impurity_charge"])
    with pytest.raises(ValueError, match="unknown requested observable"):
        w.solve(h, 2, observables=["missing"])


@pytest.mark.parametrize("stage", ["workspace-preflight", "basis-build", "basis-built", "prepared", "operator",
                                    "hamiltonian", "diagonalization", "iteration", "residuals", "observables",
                                    "validated", "workspace-commit"])
def test_stage_cancellation_preserves_completed_results(stage):
    w = QPWorkspace()
    old = w.solve(junction(), 2, options=SolverOptions(method="dense"))
    previous = old.energies.copy()

    def stop(event):
        if event["stage"] == stage:
            raise CalculationCancelled()

    with pytest.raises(CalculationCancelled):
        w.solve(junction(phi=.9), 1, options=SolverOptions(method="matrix-free"),
                context=ExecutionContext(progress=stop, interval=0))
    assert_allclose(old.energies, previous)
    resumed = w.solve(junction(), 2, options=SolverOptions(method="dense"))
    assert_allclose(resumed.energies, previous)


def test_cancel_flag_is_checked_between_throttled_progress_events():
    cancelled = False
    events = []
    ctx = ExecutionContext(events.append, lambda: cancelled, interval=100)
    ctx.checkpoint("iteration")
    ctx.checkpoint("iteration")
    assert len(events) == 1
    cancelled = True
    with pytest.raises(CalculationCancelled):
        ctx.checkpoint("iteration")


def test_live_memory_counts_borrowed_results_after_eviction_and_clear():
    w = QPWorkspace(max_entries=1)
    a = w.solve(junction(), 2)
    b = w.solve(junction(phi=.9), 1)
    assert w.stats["evictions"] > 0
    w.clear()
    assert w.live_bytes >= a.vectors.nbytes+b.vectors.nbytes+a.basis.bytes+b.basis.bytes
    del a, b
    gc.collect()
    assert w.live_bytes == 0


def test_memory_preflight_and_capacity_eviction():
    w = QPWorkspace()
    h = junction()
    a = w.solve(h, 2, options=SolverOptions(method="dense"))
    w.budget_bytes = 1
    with pytest.raises(MemoryError, match="workspace live"):
        w.solve(h, 2)
    assert a.energies.size == 1
    with pytest.raises(MemoryError):
        QPWorkspace(max_memory_gib=1e-10).solve(h, 2)
    with pytest.raises(MemoryError, match="sector dimension"):
        w.prepare(h, 2, Sector(), SolverOptions(max_dimension=1))


@pytest.mark.parametrize("kwargs", [dict(max_memory_gib=0), dict(max_memory_gib=np.nan),
                                    dict(max_entries=0), dict(count_limit=-1)])
def test_workspace_invalid_limits(kwargs):
    with pytest.raises(ValueError):
        QPWorkspace(**kwargs)


def test_hermitian_inertia_counts_pivoted_complex_two_by_two_blocks():
    matrix = np.array([[0, 1j], [-1j, 0]])
    assert hermitian_inertia(matrix, 0)["count"] == 1
    assert hermitian_inertia(matrix, -2)["count"] == 0
    assert hermitian_inertia(matrix, 2)["count"] == 2
    assert hermitian_inertia(np.diag([0., 1.]), 0)["status"] == "unresolved"
    for bad in (np.array([1, 2]), np.array([[0, 1], [0, 1]]), np.array([[np.nan]])):
        with pytest.raises(ValueError):
            hermitian_inertia(bad, 0)


@pytest.mark.parametrize("method", ["dense", "matrix-free"])
def test_independent_window_counts_match_full_spectrum(method):
    h, sector = one_particle(np.diag(np.r_[[-3., -3., -.5], np.arange(21.)]))
    w = QPWorkspace(count_limit=24)
    opts = SolverOptions(method=method, eigenpairs=2)
    w.solve(h, 1, sector, opts)
    count = w.count_below(h, 1, sector, -1.)
    assert count["status"] == "checked" and count["below"] == count["through"] == 2
    boundary = w.count_below(h, 1, sector, -.5)
    assert boundary["below"] == 2 and boundary["through"] == 3
    assert w.count_below(h, 1, Sector(1, None), 0)["status"] == "unavailable"
    with pytest.raises(ValueError):
        w.count_below(h, 1, sector, 0, margin=0)


def test_count_budget_limits_and_singular_boundary_are_explicit():
    h, sector = one_particle(np.diag(np.arange(24.)))
    w = QPWorkspace(count_limit=0)
    opts = SolverOptions(method="matrix-free")
    key, _ = w.prepare(h, 1, sector, opts)
    w.count_limit = 24
    assert w.count_below(h, 1, sector, 1)["reason"] == "projected operator not prepared"
    w.solve(h, 1, sector, opts)
    w.count_limit = 0
    assert "dimension" in w.count_below(h, 1, sector, 1)["reason"]
    w.count_limit = 24
    # The lower bracket lands exactly on a repeated/singular pivot.
    assert w.count_below(h, 1, sector, 1., margin=1.)["status"] == "unresolved"
    w.budget_bytes = 1
    assert "memory" in w.count_below(h, 1, sector, 1)["reason"]
    assert key in w.problems


def test_real_matrix_cache_can_be_reused_with_complex_start():
    h, sector = one_particle(np.diag(np.arange(24.)))
    w = QPWorkspace()
    opts = SolverOptions(method="sparse", eigenpairs=2)
    w.solve(h, 1, sector, opts)
    complex_start = np.ones(24, complex)*(1+1j)
    result = w.solve(h, 1, sector, opts, complex_start)
    assert_allclose(result.energies, [0, 1], atol=1e-9)
    changed_options = replace(opts, eigenpairs=1)
    assert w.solve(h, 1, sector, changed_options).vectors.shape == (24, 1)


def test_projected_transition_reuse_and_admission():
    w = QPWorkspace()
    h = junction()
    sector = Sector(0, 0)
    operator = h.physical_operator(create(0))
    with pytest.raises(ValueError, match="prepare the source"):
        w.project_operator(h, 2, sector, operator)
    result = w.solve(h, 2, sector)
    projected = w.project_operator(h, 2, sector, operator)
    assert projected is w.project_operator(h, 2, sector, operator)
    assert projected.bytes > 0 and w.live_bytes > result.basis.bytes


def test_external_model_and_bath_cache_is_included_in_live_memory():
    bath = cosh_grid(1, 1., 5.)
    retained = [bath, junction()]
    w = QPWorkspace(retained=lambda: retained)
    assert w.live_bytes > bath.xi.nbytes+bath.weights.nbytes
    retained.clear()
    assert w.live_bytes == 0
    with pytest.raises(ValueError, match="retained"):
        QPWorkspace(retained=[])
    events = []
    w.reserve(1000, ExecutionContext(events.append))
    assert events[0]["estimated_bytes"] == 1000
    with pytest.raises(ValueError, match="reservation"):
        w.reserve(-1)
    with pytest.raises(MemoryError):
        w.reserve(w.budget_bytes+1)
