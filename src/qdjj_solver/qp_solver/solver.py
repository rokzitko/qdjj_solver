"""Ritz eigenstates of projected fermionic Hamiltonians."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter
from typing import ClassVar
import warnings

import numpy as np
from scipy.linalg import eigh, qr
from scipy.optimize import brentq
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, aslinearoperator, eigsh, lobpcg
from threadpoolctl import threadpool_limits

from . import _core
from ..common.algebra import FermionOperator
from ..common.problem import Hamiltonian, check_transition
from .basis import DEFAULT_SECTOR, FockBasis, dimension
from .execution import ExecutionContext


class ProjectedOperator:
    """A ``FermionOperator`` acting within a QP Fock basis, or from it into another compatible basis.

    Components outside the retained basis are removed by the QP projection.
    ``expectation`` and ``transition`` evaluate matrix elements; ``action_between``
    and ``sparse_matrix`` provide the operator action and matrix representation.
    """

    def __init__(self, operator, basis):
        self.basis = basis
        keys = list(operator.terms)
        if any(abs(op) > len(basis.spins) for key in keys for op in key):
            raise ValueError("operator mode outside the Fock space")
        self.core = _core.Operator(basis.core, [operator.terms[k] for k in keys], [list(k) for k in keys])
        self.dtype = np.dtype(float if self.core.real else complex)
        self.shape = (basis.dimension, basis.dimension)

    def action(self, vector):
        """The projected operator applied to a vector of the basis."""
        vector = np.asarray(vector)
        if self.core.real and not np.iscomplexobj(vector):
            return self.core.matvec_real(vector)
        return self.core.matvec_complex(vector)

    @property
    def bytes(self):
        """Storage owned by the native operator, in bytes, excluding its shared basis."""
        return self.core.bytes

    def linear_operator(self, dtype=None):
        """The action as a SciPy ``LinearOperator``, without storing a matrix."""
        return LinearOperator(self.shape, matvec=lambda x: self.action(np.asarray(x).reshape(-1)),
                              dtype=self.dtype if dtype is None else dtype)

    def sparse_matrix(self, max_nnz=100_000_000):
        """The operator as a sparse CSR matrix with at most ``max_nnz`` stored entries."""
        data, indices, indptr = self.core.sparse(int(max_nnz))
        if self.core.real:
            data = data.real.copy()
        return csc_matrix((data, indices, indptr), shape=self.shape).tocsr()

    def expectation(self, bra, ket=None):
        """``<bra|O|ket>`` within the basis; ``ket`` defaults to ``bra``."""
        return np.vdot(bra, self.action(bra if ket is None else ket))

    def diagonal(self):
        """Return diagonal matrix elements without constructing the sparse matrix or diagonalizing it."""
        return np.real_if_close(self.core.diagonal())

    def action_between(self, ket, destination_basis):
        """Apply the operator to ``ket`` and express the result in another basis (sector or QP cutoff)
        with the same one-particle coordinates."""
        ket = np.asarray(ket)
        if (self.basis.spins != destination_basis.spins or
                self.basis.nimp != destination_basis.nimp or ket.shape != (self.basis.dimension,)):
            raise ValueError("transition vectors or mode labels do not match their bases")
        return self.core.between(ket, destination_basis.core)

    def transition(self, bra, bra_basis, ket):
        """``<bra|O|ket>`` between sectors in the *same* one-particle coordinates.

        The bra and ket may have different parity, S_z, eta and QP cutoff.
        The supplied operator must already be expressed in these coordinates.
        """
        bra, ket = np.asarray(bra), np.asarray(ket)
        if bra.shape != (bra_basis.dimension,) or self.basis.spins != bra_basis.spins:
            raise ValueError("transition vectors or mode labels do not match their bases")
        return np.vdot(bra, self.action_between(ket, bra_basis))


@dataclass(frozen=True)
class SolverOptions:
    """Numerical settings of the QP solver; energies and tolerances share the Hamiltonian's units."""

    method: str = "auto"
    """Eigensolver: ``"auto"``, ``"dense"``, ``"sparse"``, ``"matrix-free"``, ``"schur"`` or ``"block"``. ``"auto"``
    diagonalizes densely up to dimension 128, uses sparse ARPACK up to 400,000 and matrix-free ARPACK
    above. ``"block"`` uses LOBPCG for several eigenstates, including degenerate multiplets."""
    eigenpairs: int = 1
    """Number of lowest eigenstates computed in the selected sector."""
    tolerance: float = 1e-11
    """Relative ARPACK eigenvalue tolerance, for the internally shifted operator H - sigma with sigma
    below the smallest diagonal element by an offset of order one; for ``"schur"`` the root-search
    tolerance in energy units. For ``"block"``, the residual tolerance in energy units, capped at
    ``residual_tolerance/10``. Unused by dense diagonalization."""
    residual_tolerance: float = 2e-9
    """Largest accepted projected residual norm, in energy units."""
    maxiter: int = 10000
    """Maximum number of ARPACK or block-LOBPCG iterations."""
    ncv: int = 40
    """Requested Krylov-space dimension; increased as needed for the number of eigenstates and limited by the sector
    dimension."""
    max_dimension: int = 5_000_000
    """Largest allowed number of Fock basis states, checked before the basis is enumerated."""
    max_nnz: int = 100_000_000
    """Largest number of stored sparse-matrix entries (dense, sparse and Schur methods); the remaining
    ``max_memory_gib`` budget can lower it."""
    max_memory_gib: float = 24.
    """Limit on the estimated basis and matrix storage, in GiB; the total memory use can be larger."""
    threads: int = 1
    """Maximum number of threads used by the numerical libraries."""
    seed: int = 1729
    """Seed of the random starting vector of the iterative eigensolvers."""

    def __post_init__(self):
        if self.method not in ("auto", "dense", "sparse", "matrix-free", "schur", "block"):
            raise ValueError("method must be auto, dense, sparse, matrix-free, schur or block")
        for name in ("eigenpairs", "maxiter", "ncv", "max_dimension", "max_nnz", "threads"):
            value = getattr(self, name)
            if not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("tolerance", "residual_tolerance", "max_memory_gib"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")


@dataclass
class Eigenstates:
    """Lowest eigenstates of a Hamiltonian in one sector, from the QP solver."""

    backend: ClassVar[str] = "qp"
    energies: np.ndarray
    """Eigenenergies, shape (k,), in the units of the Hamiltonian."""
    vectors: np.ndarray
    """Eigenvectors in the Fock basis, shape (D, k); state i is the column ``vectors[:, i]``."""
    residuals: np.ndarray
    """Projected residual norms, shape (k,), in energy units."""
    basis: FockBasis
    """The ``FockBasis`` that defines the row order of ``vectors`` and the canonical occupations."""
    qp_weights: np.ndarray
    """Shape (k, q+1); ``qp_weights[i, n]`` is the probability of n bath quasiparticles in state i."""
    observables: dict[str, np.ndarray]
    """Expectation values of the requested named model observables, each of shape (k,); all are measured by default."""
    timings: dict[str, float]
    """Elapsed times of the calculation stages, in seconds."""
    metadata: dict
    """Numerical settings, sector, Hilbert-space restrictions and convergence diagnostics."""
    hamiltonian: Hamiltonian | None = None
    """The ``Hamiltonian`` that was solved."""

    def matrix_elements(self, operator, ket=None, *, orbital=None):
        """All ``<self_i|operator|ket_j>`` in common canonical coordinates.

        ``operator`` is a FermionOperator in canonical coordinates or a probe name
        (``hamiltonian.probe``). The ket may belong to another sector or QP cutoff.
        Eigenvector phases are arbitrary: compare ``|elements|^2`` or complete multiplets.
        """
        ket = self if ket is None else ket
        _check_coordinates(self, ket)
        if isinstance(operator, str):
            operator = self.hamiltonian.probe(operator, orbital)
        h = self.hamiltonian
        check_transition(operator, h.spins, h.eta_labels, ket.basis.sector, self.basis.sector)
        projected = ProjectedOperator(operator, ket.basis)
        vectors = np.asarray(ket.vectors)
        if same_basis(self.basis, ket.basis):
            images = projected.sparse_matrix() @ vectors
        else:
            images = np.column_stack([projected.action_between(v, self.basis) for v in vectors.T])
        return np.asarray(self.vectors).conj().T @ images

    def overlaps(self, other=None):
        """All ``<self_i|other_j>``; the states may use different QP cutoffs of one sector."""
        return self.matrix_elements(FermionOperator.scalar(1.), other)


def same_basis(a, b):
    """Whether two Fock bases enumerate the same states in the same order."""
    return a is b or (a.nimp == b.nimp and a.spins == b.spins and a.cutoff == b.cutoff
                      and a.sector == b.sector and a.dimension == b.dimension
                      and np.array_equal(a.occupations(), b.occupations()))


def _check_coordinates(bra, ket):
    if not isinstance(ket, Eigenstates):
        raise ValueError("matrix elements need QP Eigenstates on both sides")
    if bra.hamiltonian is None or ket.hamiltonian is None:
        raise ValueError("matrix elements need results that keep their Hamiltonian")
    a, b = bra.hamiltonian, ket.hamiltonian
    if (a.nimp, a.spins, a.eta_labels) != (b.nimp, b.spins, b.eta_labels):
        raise ValueError("the states use different mode labels")
    if a is not b and a.coordinate_fingerprint() != b.coordinate_fingerprint():
        raise ValueError("the states use different canonical coordinates")


DEFAULT_OPTIONS = SolverOptions()


def _refine_ritz_subspace(projected, vectors):
    """Orthonormal Ritz states in ARPACK's retained subspace.

    Complex eigsh uses the nonsymmetric ARPACK driver, whose eigenvectors in
    degenerate eigenspaces need not be orthogonal. QR alone would mix nearby
    eigenstates; diagonalize H in the orthonormal span to recover Ritz pairs.
    Reject a deficient span rather than inventing missing eigenvectors.
    """
    orthogonal, triangular, _ = qr(vectors, mode="economic", pivoting=True)
    diagonal = abs(np.diag(triangular))
    threshold = np.finfo(float).eps * max(vectors.shape) * diagonal.max()
    if np.any(diagonal <= threshold):
        raise RuntimeError("Ritz eigenvectors are linearly dependent")
    # Contract one action at a time to avoid another dimension-by-roots array.
    reduced = np.column_stack([orthogonal.conj().T @ projected.action(v)
                               for v in orthogonal.T])
    energies, rotation = eigh((reduced + reduced.conj().T)/2)
    return energies, orthogonal @ rotation


def _arpack(operator, k, initial, options, ncv):
    """The k lowest Ritz pairs of a Hermitian operator, in ascending order."""
    try:
        energies, vectors = eigsh(operator, k=k, which="SA", v0=initial,
                                  tol=options.tolerance, maxiter=options.maxiter, ncv=ncv)
    except ArpackNoConvergence as error:
        raise RuntimeError(f"Ritz iteration did not converge ({len(error.eigenvalues)}/{k} eigenpairs)") from error
    order = np.argsort(energies)
    return energies[order], vectors[:, order]


def _within_diagonal_bound(energies, diagonal, tolerance):
    """Ky Fan: the sum of the j lowest eigenvalues cannot exceed the sum of the j smallest diagonal elements.

    A larger sum of Ritz values means that the iterative solver missed eigenvalues.
    """
    count = len(energies)
    bound = np.cumsum(np.sort(diagonal)[:count])
    slack = tolerance*np.arange(1, count + 1) + 1e-12*np.maximum(1., abs(bound))
    return bool(np.all(np.cumsum(np.sort(energies)) <= bound + slack))


# ARPACK accepts a Ritz value only when its residual bound is below tol*max(eps**(2/3), |value|).
# An exactly zero eigenvalue then never converges, and ARPACK silently returns converged unwanted
# ones instead (for example the decoupled doublet at Gamma=0, whose energy is exactly zero).
# The ARPACK solvers therefore work with H - sigma, sigma = (smallest diagonal element) - offset.
# Exactly structured eigenvalues are sums of model parameters; an irrational offset keeps them away
# from zero. If a result still violates the diagonal bound, a second, incommensurate offset is tried.
SHIFTS = (2**-.5, 3**-.5)


def _shifted(matrix, sigma):
    """H - sigma as a linear operator, without copying a sparse matrix."""
    base = aslinearoperator(matrix)
    return LinearOperator(base.shape, matvec=lambda x: base.matvec(x) - sigma*x, dtype=base.dtype)


def _shifted_arpack(projected, matrix, diagonal, k, initial, options, ncv):
    """Repeat the iterative eigensolution with the second shift; energies refer to H itself."""
    for offset in SHIFTS[1:]:
        sigma = float(diagonal.min()) - offset
        _, vectors = _arpack(_shifted(matrix, sigma), k, initial, options, ncv)
        if k > 1:
            energies, vectors = _refine_ritz_subspace(projected, vectors)
        else:
            energies = np.array([np.vdot(vectors[:, 0], projected.action(vectors[:, 0])).real])
        residuals = np.array([np.linalg.norm(projected.action(v) - e*v) for e, v in zip(energies, vectors.T, strict=True)])
        if np.all(np.isfinite(residuals)) and np.max(residuals) <= options.residual_tolerance \
                and _within_diagonal_bound(energies, diagonal, options.residual_tolerance):
            return energies, vectors, residuals, sigma
    raise RuntimeError("the iterative eigensolver missed eigenvalues below the diagonal bound; use method='dense'")


def _two_qp_schur(matrix, basis, tolerance):
    """Exact elimination of a diagonal two-QP block, below its lowest pole."""
    if basis.cutoff != 2:
        raise ValueError("the Schur eigensolver requires QP cutoff 2")
    counts = basis.qp_counts()
    low, high = np.flatnonzero(counts < 2), np.flatnonzero(counts == 2)
    if len(high) == 0:
        raise ValueError("the selected sector contains no two-QP configurations")
    if len(low) == 0:
        raise ValueError("Schur elimination requires configurations below two QPs")
    hh = matrix[high][:, high]
    diagonal = hh.diagonal().real
    off_diagonal = hh-diags(diagonal)
    if off_diagonal.nnz and np.max(abs(off_diagonal.data)) > 1e-13:
        raise ValueError("Schur elimination requires a diagonal two-QP block; use sparse diagonalization")
    ll = matrix[low][:, low].toarray()
    coupling = matrix[high][:, low]
    # Solve for the distance from the pole. Subtracting two absolute energies
    # close to a pole loses the precision needed to reconstruct the high block.
    pole = float(diagonal.min())
    diagonal = diagonal-pole
    ll = ll-pole*np.eye(len(low))
    calls = 0
    def schur(energy):
        return ll + (coupling.getH() @ coupling.multiply((1/(energy-diagonal))[:, None])).toarray()
    def equation(energy):
        nonlocal calls
        calls += 1
        return float(eigh(schur(energy), subset_by_index=(0, 0), eigvals_only=True)[0]-energy)
    # Keep the reciprocal finite at the endpoint without imposing an absolute
    # physical energy margin. Resolve the pole distance with relative accuracy.
    margin = np.sqrt(np.finfo(float).tiny)
    upper = min(float(np.linalg.eigvalsh(ll)[0]), -margin)
    lower = float(np.min(matrix.diagonal().real-np.asarray(abs(matrix).sum(axis=1)).ravel()))-1.-pole
    if equation(lower) < 0 or equation(upper) > 0:
        raise ValueError("the lowest state is not on a regular Schur branch below the two-QP threshold")
    energy = brentq(equation, lower, upper, xtol=min(tolerance/10, margin), rtol=4*np.finfo(float).eps)
    _, vector_low = eigh(schur(energy), subset_by_index=(0, 0))
    vector = np.zeros(basis.dimension, dtype=matrix.dtype)
    vector[low] = vector_low[:, 0]
    vector[high] = (coupling @ vector[low])/(energy-diagonal)
    vector /= np.linalg.norm(vector)
    return np.array([energy+pole]), vector[:, None], dict(schur_dimension=len(low), schur_evaluations=calls)


def preflight(hamiltonian, cutoff=None, sector=DEFAULT_SECTOR, options=DEFAULT_OPTIONS, initial=None):
    """Exact dimension and conservative working-storage estimate, without a basis."""
    hamiltonian.check_sector(sector)
    cutoff = hamiltonian.bath_modes if cutoff is None else cutoff
    dim = dimension(hamiltonian.nimp, hamiltonian.spins, cutoff, sector, hamiltonian.eta_labels)
    if dim == 0:
        raise ValueError("empty symmetry sector")
    if dim > options.max_dimension:
        raise MemoryError(f"sector dimension {dim:,} exceeds max_dimension={options.max_dimension:,}")
    k = options.eigenpairs
    if k > dim:
        raise ValueError("more eigenpairs requested than states in the sector")
    if options.method == "schur" and (cutoff != 2 or k != 1):
        raise ValueError("Schur elimination is implemented for the lowest state at QP cutoff 2")
    is_real = (all(v.imag == 0 for v in hamiltonian.operator.terms.values())
               and (initial is None or not np.iscomplexobj(initial)))
    itemsize = 8 if is_real else 16
    ncv = min(dim, max(options.ncv, 2*k + 2))
    estimate = dim * (8*((len(hamiltonian.spins)+63)//64) + 40 + itemsize*(2*ncv + k + 8))
    if options.method == "block":
        estimate += dim * itemsize * 12*k + 12*(3*k)**2*itemsize
    if options.method == "schur":
        low_dim = dimension(hamiltonian.nimp, hamiltonian.spins, 1, sector, hamiltonian.eta_labels)
        # Reserve the dense low-QP block and elimination/eigensolver workspaces.
        estimate += 4 * low_dim * low_dim * itemsize
    if estimate > options.max_memory_gib * 1024**3:
        raise MemoryError(f"the basis and eigensolver need {estimate/1024**3:.3g} GiB, "
                          f"more than max_memory_gib={options.max_memory_gib:g}")
    method = options.method
    if method == "auto":
        method = "dense" if dim <= 128 else ("sparse" if dim <= 400_000 else "matrix-free")
    if dim <= k + 1:
        method = "dense"
    if method == "block" and dim < 5*k:
        raise ValueError("block eigensolver requires dimension >= 5*eigenpairs; use dense explicitly")
    dense_bytes = dim*dim*itemsize*4 if method == "dense" else 0
    if dense_bytes > options.max_memory_gib*1024**3:
        raise MemoryError(f"dense diagonalization needs {dense_bytes/1024**3:.3g} GiB, more than "
                          f"max_memory_gib={options.max_memory_gib:g}; use method='sparse' or 'matrix-free'")
    return dict(dimension=dim, cutoff=cutoff, eigenpairs=k, ncv=ncv, real=is_real,
                itemsize=itemsize, method=method, estimated_solver_bytes=estimate,
                estimated_peak_bytes=estimate+dense_bytes)


def solve(hamiltonian: Hamiltonian, cutoff: int | None = None, sector=DEFAULT_SECTOR,
          options=DEFAULT_OPTIONS, initial=None, *, context=None, observables=None,
          _prepared=None) -> Eigenstates:
    """Lowest Ritz states in one sector. ``cutoff=None`` retains the complete finite bath.

    Tolerances and residuals use the same energy units as the Hamiltonian, and the
    returned energies and residuals refer to H as given. The ARPACK eigensolvers
    internally work with H - sigma (sigma is recorded as ``metadata["spectral_shift"]``)
    and check their results against the diagonal (Ky Fan) bound.

    Use ``QPWorkspace.solve`` for prepared execution and ``context`` for cooperative
    progress and cancellation. ``observables=None`` measures all named operators;
    an iterable of names selects a subset, and an empty iterable skips measurements.
    """
    total_start = perf_counter()
    context = context or ExecutionContext()
    context.checkpoint("preflight")
    plan = preflight(hamiltonian, cutoff, sector, options, initial)
    dim, cutoff, k = plan["dimension"], plan["cutoff"], plan["eigenpairs"]
    is_real, ncv, method = plan["real"], plan["ncv"], plan["method"]
    estimate = plan["estimated_solver_bytes"]
    names = list(hamiltonian.observables) if observables is None else list(observables)
    if any(name not in hamiltonian.observables for name in names):
        raise ValueError("unknown requested observable")
    reuse = dict(basis=False, projected=False, matrix=False, eigensystem=False)
    times = {}
    start = perf_counter()
    context.checkpoint("basis", dimension=dim)
    if _prepared is None:
        basis = FockBasis(hamiltonian.nimp, hamiltonian.spins, cutoff, sector,
                          hamiltonian.eta_labels, options.max_dimension)
    else:
        basis = _prepared.basis
        reuse["basis"] = True
    times["basis"] = perf_counter() - start
    start = perf_counter()
    context.checkpoint("operator")
    projected = None if _prepared is None else _prepared.projected
    if projected is None:
        projected = ProjectedOperator(hamiltonian.operator, basis)
        context.checkpoint("operator-ready")
        if _prepared is not None:
            _prepared.projected = projected
    else:
        reuse["projected"] = True
    matrix = None
    context.checkpoint("hamiltonian", method=method)
    if method in ("dense", "sparse", "schur"):
        matrix_key = (method, is_real)
        matrix = None if _prepared is None else _prepared.matrices.get(matrix_key)
        if matrix is None:
            # Native vector growth, exported arrays, and CSC/CSR conversion can
            # coexist. Reserve more than the final matrix's bytes per entry.
            allowed_nnz = min(options.max_nnz, int(max(0, options.max_memory_gib*1024**3 - estimate)//96))
            matrix = projected.sparse_matrix(allowed_nnz)
            if not is_real:
                matrix = matrix.astype(complex, copy=False)
            if method == "dense":
                matrix = matrix.toarray()
            context.checkpoint("hamiltonian-ready")
            if _prepared is not None:
                _prepared.matrices[matrix_key] = matrix
        else:
            reuse["matrix"] = True
    else:
        actions = 0

        def action(vector):
            nonlocal actions
            actions += 1
            context.checkpoint("iteration", matvecs=actions)
            return projected.action(np.asarray(vector).reshape(-1))

        matrix = LinearOperator((dim, dim), matvec=action, dtype=float if is_real else complex)
    times["hamiltonian"] = perf_counter() - start
    start = perf_counter()
    extra = {}
    full_eigenpairs = None
    context.checkpoint("diagonalization", method=method, eigenpairs=k)
    with threadpool_limits(limits=options.threads):
        if method == "dense":
            full_eigenpairs = None if _prepared is None else _prepared.dense_eigenpairs
            if full_eigenpairs is None:
                full_eigenpairs = np.linalg.eigh(matrix)
            else:
                reuse["eigensystem"] = True
            energies, vectors = full_eigenpairs
            energies, vectors = energies[:k], vectors[:, :k]
        elif method == "schur":
            energies, vectors, extra = _two_qp_schur(matrix, basis, options.tolerance)
        elif method == "block":
            rng = np.random.default_rng(options.seed)
            block = rng.normal(size=(dim, k))
            if not is_real:
                block = block + 1j*rng.normal(size=(dim, k))
            if initial is not None:
                supplied = np.asarray(initial)
                if supplied.ndim == 1:
                    supplied = supplied[:, None]
                if (supplied.ndim != 2 or supplied.shape[0] != dim or
                        not np.all(np.isfinite(supplied)) or np.linalg.norm(supplied) == 0):
                    raise ValueError("invalid starting block")
                columns = min(k, supplied.shape[1])
                block[:, :columns] = supplied[:, :columns]
            block, _ = qr(block, mode="economic")
            diagonal = np.maximum(abs(projected.diagonal()), 1.)
            preconditioner = LinearOperator((dim, dim), dtype=matrix.dtype,
                matvec=lambda x: np.asarray(x).reshape(-1)/diagonal,
                matmat=lambda x: x/diagonal[:, None])
            with warnings.catch_warnings(record=True) as caught:
                energies, vectors, history = lobpcg(matrix, block, M=preconditioner,
                    largest=False, tol=min(options.tolerance, options.residual_tolerance/10),
                    maxiter=options.maxiter, retResidualNormsHistory=True)
            extra.update(block_iterations=max(0, len(history)-2),
                         block_warnings=[str(w.message) for w in caught])
            order = np.argsort(energies)
            energies, vectors = energies[order], vectors[:, order]
        else:
            if initial is None:
                rng = np.random.default_rng(options.seed)
                initial = rng.normal(size=dim)
                if not is_real:
                    initial = initial + 1j*rng.normal(size=dim)
            initial = np.asarray(initial, dtype=float if is_real else complex)
            if initial.shape != (dim,) or not np.all(np.isfinite(initial)) or np.linalg.norm(initial) == 0:
                raise ValueError("invalid starting vector")
            diagonal = np.real(projected.diagonal())
            sigma = float(diagonal.min()) - SHIFTS[0]
            ritz, vectors = _arpack(_shifted(matrix, sigma), k, initial, options, ncv)
            energies = ritz + sigma
            if k == 1:
                energies = np.array([np.vdot(vectors[:, 0], projected.action(vectors[:, 0])).real])
            extra["spectral_shift"] = sigma
        if energies.shape != (k,) or vectors.shape != (dim, k):
            raise RuntimeError("eigensolution has invalid dimensions")
        if not np.all(np.isfinite(energies)) or not np.all(np.isfinite(vectors)):
            raise RuntimeError("eigensolution contains non-finite energies or amplitudes")
        if method in ("sparse", "matrix-free", "block") and k > 1:
            energies, vectors = _refine_ritz_subspace(projected, vectors)
            extra["ritz_refinement"] = "pivoted-qr-rayleigh-ritz"
    times["diagonalization"] = perf_counter() - start
    start = perf_counter()
    context.checkpoint("residuals")
    residuals = np.array([np.linalg.norm(projected.action(v) - e*v) for e, v in zip(energies, vectors.T, strict=True)])
    if not np.all(np.isfinite(residuals)) or np.max(residuals) > options.residual_tolerance:
        raise RuntimeError(f"Ritz residual {np.max(residuals):.3e} exceeds "
                           f"residual_tolerance={options.residual_tolerance:.3e}")
    if np.linalg.norm(vectors.conj().T @ vectors - np.eye(k)) > 1e-8:
        raise RuntimeError("eigenstates are not orthonormal")
    if method in ("sparse", "matrix-free") and not _within_diagonal_bound(energies, diagonal,
                                                                          options.residual_tolerance):
        with threadpool_limits(limits=options.threads):
            energies, vectors, residuals, sigma = _shifted_arpack(projected, matrix, diagonal, k, initial,
                                                                  options, ncv)
        extra["spectral_shift"] = sigma
    counts = basis.qp_counts() if _prepared is None else _prepared.counts
    qp_weights = np.array([np.bincount(counts, weights=abs(v)**2, minlength=cutoff+1) for v in vectors.T])
    measured = {}
    for index, name in enumerate(names):
        context.checkpoint("observables", completed=index, total=len(names))
        operator = hamiltonian.observables[name]
        obs = None if _prepared is None else _prepared.observables.get(name)
        if obs is None:
            obs = ProjectedOperator(operator, basis)
            context.checkpoint("observable-ready", name=name)
            if _prepared is not None:
                _prepared.observables[name] = obs
        values = np.array([obs.expectation(v) for v in vectors.T])
        measured[name] = np.real_if_close(values)
    times["observables_and_residuals"] = perf_counter() - start
    times["total"] = perf_counter() - total_start
    metadata = dict(hamiltonian.metadata)
    metadata.update(sector=asdict(sector), qp_cutoff=cutoff, dimension=dim,
                    bath_modes=hamiltonian.bath_modes, impurity_modes=hamiltonian.nimp,
                    method=method, options=asdict(options), basis_bytes=basis.bytes,
                    estimated_solver_bytes=estimate, real_arithmetic=is_real,
                    nnz=int(matrix.nnz) if hasattr(matrix, "nnz") else None)
    metadata.update(extra)
    metadata["reuse"] = reuse
    context.checkpoint("validated")
    if _prepared is not None and full_eigenpairs is not None:
        _prepared.dense_eigenpairs = full_eigenpairs
    return Eigenstates(energies, vectors, residuals, basis, qp_weights, measured, times, metadata, hamiltonian)
