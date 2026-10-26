"""Floquet quasienergies of periodically driven models in a fixed QP Fock space.

Sambe (extended-space) diagonalization, a one-period propagator and a few-level
rotating-wave model return results of the same form. With H(t)=sum_k H_k exp(ik omega t)
and psi(t)=exp(-i eps t) sum_m exp(im omega t) u_m, the Sambe matrix is
Q_{mm'}=H_{m-m'}+m omega delta_{mm'} for |m|,|m'|<=M, the Floquet cutoff (docs/floquet.md).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from time import perf_counter
from typing import ClassVar

import numpy as np
from scipy.linalg import eigh, qr, schur
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix, diags, eye, identity, kron
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigs, splu
from threadpoolctl import threadpool_limits

from ..common.algebra import is_integer
from ..common.floquet import PeriodicHamiltonian, fold_quasienergies, quasienergy_distance, without_scalar
from .basis import DEFAULT_SECTOR, FockBasis
from .solver import Eigenstates, ProjectedOperator, SolverOptions, same_basis, solve

_CHOICES = dict(basis=("fock", "eigen"), selection=("tracked", "all"),
                branch=("auto", "label", "centroid"), method=("auto", "dense", "shift-invert"))


@dataclass(frozen=True)
class FloquetOptions:
    """Floquet settings; energies and tolerances share the Hamiltonian's units.

    ``floquet_cutoff`` is M (blocks m=-M..M). ``basis='eigen'``
    truncates the many-body space to the ``basis_states`` lowest static
    eigenstates of the reference Hamiltonian. The lowest ``states`` static
    states are tracked through the drive and label the dressed states.
    """

    floquet_cutoff: int = 4
    """Floquet cutoff M: the Sambe matrix retains the harmonic blocks m = -M, ..., M."""
    basis: str = "fock"
    """``"fock"`` uses the full QP Fock space of the sector, ``"eigen"`` the ``basis_states`` lowest static
    eigenstates."""
    states: int = 4
    """Number of lowest static states tracked through the drive; they label the dressed states."""
    basis_states: int = 32
    """Size K of the static eigenbasis for ``basis="eigen"``; choose K at a gap of the static spectrum,
    retaining complete degenerate multiplets."""
    selection: str = "tracked"
    """``"tracked"`` returns the continuations of the tracked static states, ``"all"`` every distinct Floquet
    state of the truncated Sambe matrix (dense solver or explicit ``targets``)."""
    branch: str = "auto"
    """Replica that defines the reported quasienergy: ``"label"`` continues the assigned static state,
    ``"centroid"`` is the replica whose harmonic centroid lies in [-1/2, 1/2); ``"auto"`` uses ``"label"``
    for ``selection="tracked"`` and ``"centroid"`` for ``selection="all"``."""
    method: str = "auto"
    """Sambe eigensolver: ``"dense"``, ``"shift-invert"`` or ``"auto"`` (dense up to ``dense_dimension``)."""
    dense_dimension: int = 2000
    """Largest dimension diagonalized densely: the Sambe matrix for ``method="auto"``, the static solution
    and the propagator."""
    targets: tuple[float, ...] | None = None
    """Explicit shift-invert target quasienergies; by default the static energies of the tracked states."""
    eigenpairs_per_target: int = 8
    """Initial number of shift-invert eigenpairs per target; doubled for a target until its tracked state is
    assigned a distinct Floquet state with at least ``tracking_threshold`` weight."""
    max_eigenpairs_per_target: int = 64
    """Upper limit on ``eigenpairs_per_target`` during successive doubling."""
    target_shift: float = 1e-7
    """Offset added to the default shift-invert targets to separate them from the reference
    eigenvalues."""
    tracking_threshold: float = .25
    """Minimum time-averaged weight of a tracked static state in its assigned Floquet
    state."""
    steps: int = 256
    """Number of time steps per period of the one-period propagator; doubled for the reported result when
    ``step_check=True``."""
    step_check: bool = True
    """Repeat the propagation with twice as many steps and estimate the time-step error from the change of
    the quasienergies; the finer result is reported."""
    step_tolerance: float | None = None
    """If set, raise an error when the estimated propagator step error exceeds it, in energy units; requires
    ``step_check=True``."""
    zone_center: float = 0.
    """Centre c of the zone [c - omega/2, c + omega/2) of the ``folded`` quasienergies."""
    drop_scalar_harmonics: bool = False
    """Remove the c-number parts of the harmonics H_k with k != 0; the quasienergies are unchanged, but the
    harmonic weights of the modes are redistributed."""
    tolerance: float = 1e-11
    """Eigensolver tolerance of the static solution and of the shift-invert iterations."""
    residual_tolerance: float = 2e-9
    """Largest accepted window residual of the Sambe and static solutions, in energy units."""
    maxiter: int = 10000
    """Maximum number of iterations of the static and shift-invert eigensolvers."""
    max_dimension: int = 5_000_000
    """Largest allowed Fock-basis and Sambe dimension."""
    max_nnz: int = 100_000_000
    """Largest number of stored sparse-matrix entries."""
    max_memory_gib: float = 24.
    """Limit on the estimated storage of dense Sambe matrices, sparse LU factors, propagators and probe images,
    in GiB."""
    threads: int = 1
    """Maximum number of threads used by the numerical libraries."""
    seed: int = 1729
    """Seed of the random starting vectors of the iterative eigensolvers."""

    def __post_init__(self):
        for name, allowed in _CHOICES.items():
            if getattr(self, name) not in allowed:
                raise ValueError(f"{name} must be one of: {', '.join(allowed)}")
        if not is_integer(self.floquet_cutoff) or self.floquet_cutoff < 0:
            raise ValueError("floquet_cutoff must be a nonnegative integer")
        for name in ("states", "basis_states", "dense_dimension", "eigenpairs_per_target",
                     "max_eigenpairs_per_target", "steps", "maxiter", "max_dimension", "max_nnz", "threads"):
            if not is_integer(getattr(self, name)) or getattr(self, name) < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.max_eigenpairs_per_target < self.eigenpairs_per_target:
            raise ValueError("max_eigenpairs_per_target must be at least eigenpairs_per_target")
        for name in ("tolerance", "residual_tolerance", "max_memory_gib"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        for name in ("target_shift", "zone_center"):
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite")
        if not 0 <= self.tracking_threshold <= 1:
            raise ValueError("tracking_threshold must lie between zero and one")
        if self.step_tolerance is not None and (not np.isfinite(self.step_tolerance) or self.step_tolerance <= 0):
            raise ValueError("step_tolerance must be positive and finite")
        for name in ("step_check", "drop_scalar_harmonics"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a boolean")
        if self.step_tolerance is not None and not self.step_check:
            raise ValueError("step_tolerance requires step_check=True")
        if self.targets is not None:
            targets = tuple(float(x) for x in self.targets)
            if not targets or not np.all(np.isfinite(targets)):
                raise ValueError("targets must be a nonempty sequence of finite energies")
            object.__setattr__(self, "targets", targets)


@dataclass
class FloquetStates:
    """Floquet states psi(t)=exp(-i eps t) sum_m exp(i m omega t) u_m.

    ``modes[i, :, n]`` holds u_m of state n for m = ``harmonic_indices[i]``, in the Fock
    basis or in the static eigenbasis (``representation='static'``). The stored
    window ``|m| <= M'`` can exceed the Floquet cutoff M, because the branch choice
    shifts the harmonic index. Quasienergies belong to the chosen branch;
    ``mean_energies`` (eps-omega<m>) are branch invariant. ``static_weights[n, a]``
    is the time-averaged weight of static state a.
    """

    kind: ClassVar[str] = "floquet"
    method: str
    """Solver that produced the states: ``"sambe"``, ``"propagator"`` or ``"rwa"``."""
    representation: str
    """``"fock"`` or ``"static"``: the basis of ``modes``."""
    omega: float
    """Drive frequency."""
    quasienergies: np.ndarray
    """Quasienergies, shape (n,), on the chosen branch."""
    folded: np.ndarray
    """Quasienergies folded into the zone around ``zone_center``, shape (n,)."""
    mean_energies: np.ndarray
    """Replica-independent time-averaged energies eps - omega <m>, shape (n,)."""
    harmonic_indices: np.ndarray
    """Harmonic indices m of the stored modes, shape (2M'+1,)."""
    modes: np.ndarray
    """Mode harmonics, shape (2M'+1, d, n): u_m of state j is ``modes[i, :, j]`` with ``harmonic_indices[i] == m``, in
    the basis given by ``representation``."""
    harmonic_weights: np.ndarray
    """Harmonic weights, the squared norms of u_m, shape (n, 2M'+1)."""
    centroids: np.ndarray
    """Harmonic centroids <m>, shape (n,)."""
    boundary_weights: np.ndarray
    """Weight in the outermost harmonic blocks m = -M and M (Sambe) or outside the window -M <= m <= M
    (propagator); zero for the rotating-wave model; shape (n,)."""
    static_weights: np.ndarray
    """Time-averaged weights of the static eigenstates a, the sum over m of the squared overlaps of a with u_m;
    shape (n, s), with s = ``states`` for the Fock basis and s = K for the static basis."""
    labels: np.ndarray
    """Index of the static state assigned to each Floquet state, shape (n,)."""
    residuals: np.ndarray
    """Residual bound, shape (n,): the Sambe window residual combined with the leakage into the omitted
    harmonics, or the eigendecomposition residual of U(T) or of the rotating-wave matrix."""
    basis: FockBasis
    """The QP ``FockBasis`` of the sector."""
    static: Eigenstates
    """Static ``Eigenstates`` of the reference Hamiltonian ``periodic.reference``."""
    periodic: PeriodicHamiltonian
    """The ``PeriodicHamiltonian`` that was solved."""
    fock_residuals: np.ndarray | None = None
    """Fock-space residuals of static-basis results, shape (n,); ``None`` otherwise."""
    effective_matrix: np.ndarray | None = None
    """Rotating-wave matrix of ``rwa_model`` results; ``None`` otherwise."""
    timings: dict = field(default_factory=dict)
    """Elapsed times of the calculation stages, in seconds."""
    metadata: dict = field(default_factory=dict)
    """Settings and diagnostics, for example window residuals, leakage and step errors."""

    def fock_modes(self):
        """Mode harmonics u_m in the QP Fock basis, shape (2M'+1, Fock dimension, n)."""
        if self.representation == "fock":
            return self.modes
        vectors = np.asarray(self.static.vectors[:, :self.modes.shape[1]], complex)
        return np.einsum("Dd,mdn->mDn", vectors, self.modes)

    def mode_at(self, t):
        """Periodic part u(t)=sum_m exp(i m omega t) u_m in the stored representation."""
        return np.einsum("m,mdn->dn", np.exp(1j*self.harmonic_indices*self.omega*t), self.modes)

    def state_at(self, t):
        """Floquet state psi(t)=exp(-i eps t) u(t) in the stored representation."""
        return self.mode_at(t)*np.exp(-1j*self.quasienergies*t)

    def dressed_state(self, label=0):
        """Index of the Floquet state assigned to the static state ``label``.

        If no state, or more than one, carries this label, the state with the
        largest weight on the static state ``label`` is returned.
        """
        hits = np.flatnonzero(self.labels == label)
        if len(hits) == 1:
            return int(hits[0])
        if not 0 <= label < self.static_weights.shape[1]:
            raise ValueError(f"static state {label} is not tracked")
        return int(np.argmax(self.static_weights[:, label]))

    def matrix_elements(self, probe, ket=None, *, orbital=None, ket_states=None, bra_states=None,
                        max_memory_gib=None):
        """Probe matrix elements X^(q) from the ket states (initial) to these states (final).

        ``probe`` is a probe name (``PeriodicHamiltonian.probe``), a FermionOperator
        (treated as time independent) or a dict {k: FermionOperator} of harmonics;
        ``orbital`` selects the impurity orbital of a named probe. ``ket`` defaults to
        these states and may belong to another sector of the same drive.
        ``ket_states`` and ``bra_states`` restrict the state indices, and
        ``max_memory_gib`` limits the storage of the probe images X(t)|u(t)>.
        """
        from .transitions import floquet_matrix_elements
        return floquet_matrix_elements(self, probe, ket, orbital=orbital, ket_states=ket_states,
                                       bra_states=bra_states, max_memory_gib=max_memory_gib)


def _options(options):
    if options is None:
        return FloquetOptions()
    if isinstance(options, dict):
        return FloquetOptions(**options)
    if not isinstance(options, FloquetOptions):
        raise ValueError("options must be FloquetOptions, a dict or None")
    return options


def _static_options(options, eigenpairs, dense):
    return SolverOptions(method="dense" if dense else "auto", eigenpairs=eigenpairs,
                         tolerance=options.tolerance, residual_tolerance=options.residual_tolerance,
                         maxiter=options.maxiter, max_dimension=options.max_dimension,
                         max_nnz=options.max_nnz, max_memory_gib=options.max_memory_gib,
                         threads=options.threads, seed=options.seed)


@dataclass
class _Problem:
    periodic: PeriodicHamiltonian
    options: FloquetOptions
    cutoff: int
    sector: object
    basis: FockBasis
    static: Eigenstates
    fock: dict
    blocks: dict
    representation: str
    tracked: np.ndarray
    times: dict
    static_gap: float | None = None

    @property
    def dimension(self):
        return self.tracked.shape[0]


def _prepare(periodic, cutoff, sector, options, static, representation, needed, gap_check=True):
    """Fock basis, reference eigenstates and projected harmonics H_k for both signs of k."""
    if not isinstance(periodic, PeriodicHamiltonian):
        raise ValueError("provide a PeriodicHamiltonian")
    times = {}
    start = perf_counter()
    periodic.check_sector(sector)
    cutoff = periodic.reference.bath_modes if cutoff is None else cutoff
    basis = FockBasis(periodic.nimp, periodic.spins, cutoff, sector, periodic.eta_labels, options.max_dimension)
    times["basis"] = perf_counter() - start
    start = perf_counter()
    count = min(needed, basis.dimension)
    gap_check = gap_check and representation == "eigen" and count < basis.dimension
    request = count + 1 if gap_check else count
    dense = basis.dimension <= options.dense_dimension
    if static is None:
        # One extra eigenpair reveals whether the tracked set cuts a degenerate multiplet.
        static = solve(periodic.reference, cutoff, sector,
                       _static_options(options, min(count + 1, basis.dimension), dense))
    else:
        _check_static(static, periodic, basis, request)
    energies = static.energies
    static_gap = float(energies[count] - energies[count-1]) if len(energies) > count else None
    if gap_check and static_gap < 1e-8*max(1., abs(energies[count-1])):
        raise ValueError("basis_states splits a degenerate static multiplet; change basis_states")
    times["static"] = perf_counter() - start
    start = perf_counter()
    fock = {}
    for k, operator in periodic.harmonics.items():
        if k and options.drop_scalar_harmonics:
            operator = without_scalar(operator)
        if k and not operator.terms:
            continue
        matrix = ProjectedOperator(operator, basis).sparse_matrix(options.max_nnz).astype(complex)
        fock[k] = matrix
        if k:
            fock[-k] = matrix.conj().T.tocsr()
    vectors = np.asarray(static.vectors[:, :count], complex)
    if representation == "fock":
        blocks, tracked = fock, vectors
    else:
        blocks = {k: vectors.conj().T @ (matrix @ vectors) for k, matrix in fock.items()}
        blocks[0] = (blocks[0] + blocks[0].conj().T)/2
        tracked = np.eye(count, dtype=complex)
    times["harmonics"] = perf_counter() - start
    return _Problem(periodic, options, cutoff, sector, basis, static, fock, blocks,
                    "fock" if representation == "fock" else "static", tracked, times, static_gap)


def _check_static(static, periodic, basis, request):
    if not isinstance(static, Eigenstates):
        raise ValueError("static must be QP Eigenstates")
    fingerprints = {periodic.coordinate_fingerprint(), periodic.time_average().coordinate_fingerprint()}
    if static.hamiltonian is None or static.hamiltonian.coordinate_fingerprint() not in fingerprints:
        raise ValueError("the static states use different canonical coordinates")
    if not same_basis(static.basis, basis):
        raise ValueError("the static states belong to a different Fock basis")
    if len(static.energies) < request:
        raise ValueError(f"static states provide {len(static.energies)} eigenpairs; {request} are required")


def _sambe(blocks, harmonics, omega, dimension):
    size = 2*harmonics + 1
    q = kron(diags(omega*np.arange(-harmonics, harmonics+1, dtype=float)), identity(dimension), format="csr")
    q = q.astype(complex)
    for k, block in blocks.items():
        if abs(k) < size:
            q = q + kron(eye(size, k=-k), csr_matrix(block), format="csr")
    return q.tocsr()


def _sambe_residuals(blocks, modes, energies, omega, harmonics):
    """Window residual ||(Q-eps)u|| and leakage of H into blocks |m|>M."""
    top = max(abs(k) for k in blocks)
    inside, outside = np.zeros(len(energies)), np.zeros(len(energies))
    for m in range(-harmonics-top, harmonics+top+1):
        total = np.zeros(modes.shape[1:], complex)
        for mp in range(max(-harmonics, m-top), min(harmonics, m+top)+1):
            block = blocks.get(m-mp)
            if block is not None:
                total += block @ modes[mp+harmonics]
        if abs(m) <= harmonics:
            total += modes[m+harmonics]*(m*omega - energies)[None, :]
            inside += np.sum(abs(total)**2, axis=0)
        else:
            outside += np.sum(abs(total)**2, axis=0)
    return np.sqrt(inside), np.sqrt(outside)


def _shifted_overlap(a, b, shift):
    """|sum_m <a_m|b_{m-shift}>| for modes stored on the same harmonic window."""
    if shift >= 0:
        return abs(np.vdot(a[shift:], b[:len(b)-shift])) if shift < len(a) else 0.
    return _shifted_overlap(b, a, -shift)


def _distinct(energies, modes, centroids, omega):
    """Keep one replica (smallest |<m>|) of each Floquet state.

    A replica shifted by l harmonics has quasienergy eps+l*omega. It is recognized
    by its projection on the kept states of that energy, which also handles
    degenerate multiplets whose replicas are arbitrary rotations within the multiplet.
    Distinct Floquet states are orthogonal in the full Sambe space, so the overlap
    decides; the energy test only preselects pairs, because window truncation
    shifts replica energies (strongly near a resonant chain cut by the window edge).
    """
    order = np.argsort(abs(centroids), kind="stable")
    kept = []
    tolerance = .05*omega
    for i in order:
        if kept:
            difference = energies[i] - energies[kept]
            shifts = np.rint(difference/omega).astype(int)
            close = np.flatnonzero((shifts != 0) & (abs(difference - shifts*omega) < tolerance))
            projection = {}
            for j in close:
                overlap = _shifted_overlap(modes[:, :, i], modes[:, :, kept[j]], int(shifts[j]))
                projection[shifts[j]] = projection.get(shifts[j], 0.) + overlap**2
            if any(weight > .5 for weight in projection.values()):
                continue
        kept.append(i)
    return np.array(sorted(kept), dtype=int)


def _require_distinct(energies, modes, omega):
    """Refuse two tracked states that are replicas of one Floquet state."""
    for i in range(len(energies)):
        for j in range(i):
            shift = int(np.rint((energies[i]-energies[j])/omega))
            if shift and _shifted_overlap(modes[:, :, i], modes[:, :, j], shift)**2 > .5:
                raise RuntimeError("two tracked states are replicas of one Floquet state; increase floquet_cutoff")


def _static_weights(tracked, modes):
    """Time-averaged weights of the tracked static states, shape (Floquet states, tracked states)."""
    return np.sum(abs(np.einsum("da,mdn->man", tracked.conj(), modes))**2, axis=0).T


def _assign(static_weights, candidates, states):
    """One distinct candidate for each tracked state 0, ..., states-1, maximizing the total static weight.

    Returns the selected candidates and their labels, ordered by label.
    """
    rows, columns = linear_sum_assignment(-static_weights[candidates, :states])
    order = np.argsort(columns)
    return candidates[rows[order]], columns[order]


def _finish(problem, method, energies, modes, residuals, extra, per_state=None, boundary=None,
            fock_residuals=None, effective=None, central=None, mean=None):
    """Select states, apply the branch convention and assemble FloquetStates."""
    options, omega = problem.options, problem.periodic.omega
    harmonics = (modes.shape[0]-1)//2
    m = np.arange(-harmonics, harmonics+1)
    weights = np.sum(abs(modes)**2, axis=1).T
    centroids = weights @ m/np.maximum(weights.sum(axis=1), np.finfo(float).tiny)
    candidates = np.arange(len(energies)) if central is None else central
    tracked = problem.tracked
    projections = np.einsum("da,mdn->man", tracked.conj(), modes)
    static_weights = np.sum(abs(projections)**2, axis=0).T
    if options.selection == "tracked" and method != "rwa":
        states = min(options.states, tracked.shape[1])
        if len(candidates) < states:
            raise RuntimeError(f"{len(candidates)} Floquet candidates for {states} tracked states; "
                               "increase floquet_cutoff (or eigenpairs_per_target for shift-invert)")
        selected, labels = _assign(static_weights, candidates, states)
        _require_distinct(energies[selected], modes[:, :, selected], omega)
        weak = static_weights[selected, labels] < options.tracking_threshold
        if np.any(weak):
            a = int(labels[np.argmax(weak)])
            raise RuntimeError(f"tracked static state {a} has Floquet weight "
                               f"{static_weights[selected, labels][weak].min():.3f} below tracking_threshold; "
                               "increase floquet_cutoff, eigenpairs_per_target or basis_states")
    else:
        selected = candidates[np.argsort(energies[candidates], kind="stable")]
        labels = np.argmax(static_weights[selected], axis=1)
    branch = options.branch
    if branch == "auto":
        branch = "label" if options.selection == "tracked" else "centroid"
    weight_on_label = static_weights[selected, labels]
    if branch == "label":
        local = np.sum(abs(projections[:, labels, selected])**2*m[:, None], axis=0)
        local = np.where(weight_on_label > 0, local/np.maximum(weight_on_label, np.finfo(float).tiny),
                         centroids[selected])
    else:
        local = centroids[selected]
    shifts = -np.floor(local + .5).astype(int)
    mean = energies[selected] - omega*centroids[selected] if mean is None else np.asarray(mean)[selected]
    width = harmonics + (int(np.max(abs(shifts))) if len(shifts) else 0)
    stored = np.zeros((2*width+1, modes.shape[1], len(selected)), complex)
    for column, (state, shift) in enumerate(zip(selected, shifts, strict=True)):
        start = width - harmonics + shift
        stored[start:start+2*harmonics+1, :, column] = modes[:, :, state]
        label_vector = tracked[:, labels[column]]
        reference = np.vdot(label_vector, stored[width, :, column])
        if abs(reference) < 1e-12:
            flat = stored[:, :, column].ravel()
            reference = flat[np.argmax(abs(flat))]
        stored[:, :, column] *= np.conj(reference)/abs(reference)
    quasienergies = energies[selected] + omega*shifts
    window = np.arange(-width, width+1)
    stored_weights = np.sum(abs(stored)**2, axis=1).T
    if boundary is None:
        boundary = weights[:, 0] + (weights[:, -1] if harmonics else 0.)
    metadata = dict(problem.periodic.metadata.get("floquet", {}))
    metadata = dict(periodic=metadata, sector=asdict(problem.sector), qp_cutoff=problem.cutoff,
                    dimension=problem.basis.dimension, representation_dimension=problem.dimension,
                    method=method, options=asdict(options), branch=branch, branch_shifts=shifts.tolist(),
                    selection="all" if method == "rwa" else options.selection,
                    zone_center=options.zone_center, coordinates_sha256=problem.periodic.coordinate_fingerprint())
    # Complete: every Floquet state of the representation is present (needed for sum rules).
    complete = method != "rwa" and options.selection == "all" and len(selected) == problem.dimension
    gap = problem.static_gap
    scale = max(1., abs(float(problem.static.energies[0])))
    metadata.update(complete=bool(complete), static_gap=gap,
                    multiplet_split=bool(not complete and gap is not None and gap < 1e-8*scale))
    metadata.update(extra)
    metadata.update({name: np.asarray(values)[selected].tolist() for name, values in (per_state or {}).items()})
    norms = np.maximum(stored_weights.sum(axis=1), np.finfo(float).tiny)
    return FloquetStates(method, problem.representation, omega, quasienergies,
                         fold_quasienergies(quasienergies, omega, options.zone_center), mean, window, stored,
                         stored_weights, stored_weights @ window/norms, np.asarray(boundary)[selected],
                         static_weights[selected], labels, np.asarray(residuals)[selected], problem.basis,
                         problem.static,
                         problem.periodic, None if fock_residuals is None else fock_residuals[selected],
                         effective, problem.times, metadata)


def _needed(options):
    if options.basis == "eigen":
        if options.states > options.basis_states:
            raise ValueError("states, the number of tracked states, must not exceed basis_states")
        return options.basis_states
    return options.states


def solve_floquet(periodic, cutoff=None, sector=DEFAULT_SECTOR, options=None, static=None) -> FloquetStates:
    """Floquet states from the Sambe matrix truncated at the Floquet cutoff, ``|m| <= M``.

    ``residuals`` combine the window residual with the leakage of H into the
    omitted harmonics; since the infinite Sambe operator is self-adjoint, an
    exact quasienergy lies within that distance. ``basis='eigen'`` also reports
    ``fock_residuals``, the same bound in the full Fock space of the QP cutoff.
    """
    options = _options(options)
    problem = _prepare(periodic, cutoff, sector, options, static, options.basis, _needed(options))
    omega, harmonics, d = periodic.omega, options.floquet_cutoff, problem.dimension
    size = (2*harmonics+1)*d
    if size > options.max_dimension:
        raise MemoryError(f"Sambe dimension {size:,} exceeds max_dimension={options.max_dimension:,}")
    start = perf_counter()
    q = _sambe(problem.blocks, harmonics, omega, d)
    problem.times["sambe"] = perf_counter() - start
    method = options.method
    if method == "auto":
        method = "dense" if size <= options.dense_dimension else "shift-invert"
    if method == "shift-invert" and size < max(8, 2*options.eigenpairs_per_target + 2):
        method = "dense"
    if method == "shift-invert" and options.selection == "all" and options.targets is None:
        raise ValueError("selection='all' with shift-invert returns only the states near its targets; "
                         "use method='dense' (raise dense_dimension, or basis='eigen') or give explicit targets")
    start = perf_counter()
    extra = dict(sambe_dimension=size, sambe_nnz=int(q.nnz), floquet_cutoff=harmonics,
                 max_drive_harmonic=int(max(abs(k) for k in problem.blocks)))
    with threadpool_limits(limits=options.threads):
        if method == "dense":
            if 3*16*size*size > options.max_memory_gib*1024**3:
                raise MemoryError(f"dense Sambe diagonalization needs {3*16*size*size/1024**3:.3g} GiB, more than "
                                  f"max_memory_gib={options.max_memory_gib:g}; use method='shift-invert'")
            matrix = q.toarray()
            energies, vectors = eigh((matrix + matrix.conj().T)/2)
        else:
            energies, vectors, info = _shift_invert(q, problem)
            extra.update(info)
        modes = vectors.reshape(2*harmonics+1, d, -1)
        weights = np.sum(abs(modes)**2, axis=1).T
        centroids = weights @ np.arange(-harmonics, harmonics+1)
        # The dense spectrum contains every replica. Each Floquet state has one with |<m>|<=1;
        # a hard cut at |<m>|<1/2 can lose a state whose replicas are distorted by the window edges.
        pool = np.flatnonzero(abs(centroids) <= 1.) if method == "dense" else np.arange(len(energies))
        distinct = pool[_distinct(energies[pool], modes[:, :, pool], centroids[pool], omega)]
    problem.times["diagonalization"] = perf_counter() - start
    start = perf_counter()
    inside, outside = np.zeros(len(energies)), np.zeros(len(energies))
    window, leakage = _sambe_residuals(problem.blocks, modes[:, :, distinct], energies[distinct], omega, harmonics)
    inside[distinct], outside[distinct] = window, leakage
    if np.max(window, initial=0.) > options.residual_tolerance:
        raise RuntimeError(f"Sambe residual {np.max(window):.3e} exceeds "
                           f"residual_tolerance={options.residual_tolerance:.3e}")
    fock_residuals = None
    if problem.representation == "static":
        lifted = np.einsum("Dd,mdn->mDn", np.asarray(problem.static.vectors[:, :d], complex), modes[:, :, distinct])
        fock_window, fock_leakage = _sambe_residuals(problem.fock, lifted, energies[distinct], omega, harmonics)
        fock_residuals = np.zeros(len(energies))
        fock_residuals[distinct] = np.hypot(fock_window, fock_leakage)
    problem.times["residuals"] = perf_counter() - start
    extra.update(eigensolver=method)
    return _finish(problem, "sambe", energies, modes, np.hypot(inside, outside), extra,
                   per_state=dict(window_residuals=inside, leakage=outside),
                   fock_residuals=fock_residuals, central=distinct)


def _shift_invert(q, problem):
    """Eigenpairs near targets by sparse LU, merged by Rayleigh-Ritz.

    With the default targets, the number of eigenpairs of a target is doubled, up to
    ``max_eigenpairs_per_target``, until the selection of ``_finish`` finds its tracked
    state: a Floquet state with sufficient static weight that is not a replica of the
    state selected for another target. Near a resonance, the eigenpairs closest to one
    target can be replicas of the state found for another target, or states at the edge
    of the harmonic window.
    """
    options, static = problem.options, problem.static
    n = q.shape[0]
    tracking = options.targets is None
    if tracking:
        count = min(options.states, problem.tracked.shape[1])
        targets = [float(e) + options.target_shift for e in static.energies[:count]]
    else:
        targets = list(options.targets)
    limit = options.max_memory_gib*1024**3
    if 64*q.nnz > limit:
        raise MemoryError(f"the sparse Sambe matrix needs {64*q.nnz/1024**3:.3g} GiB, more than "
                          f"max_memory_gib={options.max_memory_gib:g}; use basis='eigen'")
    rng = np.random.default_rng(options.seed)
    largest = min(options.max_eigenpairs_per_target, n-2)
    used = [min(options.eigenpairs_per_target, largest)]*len(targets)
    found = [None]*len(targets)
    pending = range(len(targets))
    while True:
        for i in pending:
            found[i] = _eigenvectors_near(q, targets[i], used[i], rng, options)
        energies, vectors, rank = _rayleigh_ritz(q, found)
        if not tracking:
            break
        pending = [a for a in _missing(problem, energies, vectors) if used[a] < largest]
        if not pending:
            break
        for a in pending:
            used[a] = min(2*used[a], largest)
    return energies, vectors, dict(targets=targets, eigenpairs_per_target=used, ritz_dimension=rank)


def _eigenvectors_near(q, sigma, k, rng, options):
    """The k eigenvectors of q with eigenvalues closest to sigma (shift-invert Arnoldi)."""
    n = q.shape[0]
    # Q is structurally symmetric; minimum degree on A^T+A limits the fill-in.
    lu = splu((q - sigma*identity(n, format="csc", dtype=complex)).tocsc(), permc_spec="MMD_AT_PLUS_A")
    factor_bytes = 16*(lu.L.nnz + lu.U.nnz)
    if factor_bytes > options.max_memory_gib*1024**3:
        raise MemoryError(f"the sparse LU factors need {factor_bytes/1024**3:.3g} GiB, more than "
                          f"max_memory_gib={options.max_memory_gib:g}; use basis='eigen'")
    operator = LinearOperator((n, n), matvec=lu.solve, dtype=complex)
    start = rng.normal(size=n) + 1j*rng.normal(size=n)
    try:
        _, vectors = eigs(operator, k=k, which="LM", v0=start, tol=options.tolerance, maxiter=options.maxiter)
    except ArpackNoConvergence as error:
        raise RuntimeError(f"shift-invert iteration did not converge near {sigma:.6g}") from error
    return vectors


def _rayleigh_ritz(q, blocks):
    """Eigenpairs of q in the span of the column blocks: (energies, orthonormal vectors, rank)."""
    orthogonal, triangular, _ = qr(np.hstack(blocks), mode="economic", pivoting=True)
    diagonal = abs(np.diag(triangular))
    rank = int(np.sum(diagonal > 1e-10*diagonal.max()))
    orthogonal = orthogonal[:, :rank]
    reduced = orthogonal.conj().T @ (q @ orthogonal)
    energies, rotation = eigh((reduced + reduced.conj().T)/2)
    return energies, orthogonal @ rotation, rank


def _missing(problem, energies, vectors):
    """Tracked states that the selection of ``_finish`` would not find among these eigenpairs."""
    options = problem.options
    states = min(options.states, problem.tracked.shape[1])
    harmonics = options.floquet_cutoff
    modes = vectors.reshape(2*harmonics+1, problem.dimension, -1)
    centroids = np.sum(abs(modes)**2, axis=1).T @ np.arange(-harmonics, harmonics+1)
    candidates = _distinct(energies, modes, centroids, problem.periodic.omega)
    if len(candidates) < states:
        return list(range(states))
    weights = _static_weights(problem.tracked, modes)
    selected, labels = _assign(weights, candidates, states)
    return [int(a) for a in labels[weights[selected, labels] < options.tracking_threshold]]


def _dense_blocks(problem):
    return {k: (block.toarray() if hasattr(block, "toarray") else np.asarray(block, complex))
            for k, block in problem.blocks.items()}


def floquet_propagator(periodic, cutoff=None, sector=DEFAULT_SECTOR, options=None, static=None) -> FloquetStates:
    """Floquet states from the eigenphases of the one-period propagator U(T).

    The trigonometric polynomial H(t) is integrated with the fourth-order
    commutator-free Magnus scheme, without a Floquet cutoff. The harmonics of the
    modes are then truncated to ``|m| <= M``, as in the Sambe results, and
    ``boundary_weights`` report the weight outside that window.
    """
    options = _options(options)
    problem = _prepare(periodic, cutoff, sector, options, static, options.basis, _needed(options))
    d, omega = problem.dimension, periodic.omega
    if d > options.dense_dimension:
        raise MemoryError(f"propagator dimension {d:,} exceeds dense_dimension={options.dense_dimension:,}; "
                          "use basis='eigen'")
    harmonics = options.floquet_cutoff
    if options.steps <= 2*harmonics + 1:
        raise ValueError("steps must exceed 2*floquet_cutoff+1 to resolve the mode harmonics")
    blocks = _dense_blocks(problem)
    offset = float(np.trace(blocks[0]).real/d)
    blocks[0] = blocks[0] - offset*np.eye(d)
    period = 2*np.pi/omega

    def hamiltonian(t):
        return sum(block*np.exp(1j*k*omega*t) for k, block in blocks.items())

    start = perf_counter()
    with threadpool_limits(limits=options.threads):
        steps = options.steps
        coarse = None
        if options.step_check:
            coarse = _eigenphases(_propagate(hamiltonian, period, steps, False)[0], period, options)
            steps *= 2
        # Keep the step propagators, and later the sampled states, if each fits in a quarter of the memory limit.
        store = 16*steps*d*d < options.max_memory_gib*1024**3/4
        propagator, cache = _propagate(hamiltonian, period, steps, store)
        folded, vectors, decomposition = _eigenphases(propagator, period, options)
        errors = None
        if coarse is not None:
            distance = quasienergy_distance(folded[:, None], coarse[0][None, :], omega)
            rows, columns = linear_sum_assignment(distance)
            errors = np.zeros(d)
            errors[rows] = distance[rows, columns]/15
            if options.step_tolerance is not None and np.max(errors) > options.step_tolerance:
                raise RuntimeError(f"propagator step error {np.max(errors):.3e} exceeds "
                                   f"step_tolerance={options.step_tolerance:.3e}")
        # Mean energies fix the central replica; harmonics are then taken in that branch.
        h = period/steps
        mean = np.zeros(d)
        state = vectors.copy()
        samples = np.empty((steps, d, d), complex) if store else None
        for j in range(steps):
            if store:
                samples[j] = state
            mean += np.real(np.einsum("dn,dn->n", state.conj(), hamiltonian(j*h) @ state))/steps
            state = _step(hamiltonian, j*h, h, cache, j) @ state
        centroids = (folded - mean)/omega
        energies = folded - omega*np.floor(centroids + .5)
        m = np.arange(-harmonics, harmonics+1)
        if store:
            samples *= np.exp(1j*np.outer(np.arange(steps)*h, energies))[:, None, :]
            modes = (np.fft.fft(samples, axis=0)/steps)[m % steps]
        else:
            modes = np.zeros((2*harmonics+1, d, d), complex)
            state = vectors.copy()
            for j in range(steps):
                t = j*h
                sample = state*np.exp(1j*energies*t)[None, :]
                modes += np.exp(-1j*m*omega*t)[:, None, None]*sample[None]/steps
                state = _step(hamiltonian, t, h, cache, j) @ state
    problem.times["propagation"] = perf_counter() - start
    boundary = np.clip(1 - np.sum(abs(modes)**2, axis=(0, 1)), 0., 1.)
    extra = dict(steps=steps, energy_offset=offset, unitarity_error=decomposition["unitarity_error"],
                 schur_offdiagonal=decomposition["offdiagonal"], integrator="commutator-free Magnus, fourth order")
    per_state = {} if errors is None else dict(step_errors=errors)
    return _finish(problem, "propagator", energies + offset, modes, decomposition["residuals"], extra,
                   per_state=per_state, boundary=boundary, mean=mean + offset)


_CF4 = ((3-2*np.sqrt(3))/12, (3+2*np.sqrt(3))/12, .5-np.sqrt(3)/6, .5+np.sqrt(3)/6)


def _exponential(matrix, h):
    values, vectors = np.linalg.eigh((matrix + matrix.conj().T)/2)
    return (vectors*np.exp(-1j*h*values)) @ vectors.conj().T


def _step(hamiltonian, t, h, cache=None, index=None):
    if cache is not None:
        return cache[index]
    a1, a2, c1, c2 = _CF4
    h1, h2 = hamiltonian(t + c1*h), hamiltonian(t + c2*h)
    return _exponential(a1*h1 + a2*h2, h) @ _exponential(a2*h1 + a1*h2, h)


def _propagate(hamiltonian, period, steps, keep):
    h = period/steps
    cache = [] if keep else None
    result = None
    for j in range(steps):
        step = _step(hamiltonian, j*h, h)
        if keep:
            cache.append(step)
        result = step if result is None else step @ result
    return result, cache


def _eigenphases(propagator, period, options):
    d = propagator.shape[0]
    triangular, vectors = schur(propagator, output="complex")
    values = np.diag(triangular)
    offdiagonal = float(np.linalg.norm(np.triu(triangular, 1)))
    unitarity = float(np.linalg.norm(propagator.conj().T @ propagator - np.eye(d)))
    if unitarity > 1e-8 or offdiagonal > 1e-8:
        raise RuntimeError("the one-period propagator is not numerically unitary")
    residuals = np.linalg.norm(propagator @ vectors - vectors*values[None, :], axis=0)
    folded = -np.angle(values)/period
    return folded, vectors, dict(residuals=residuals, unitarity_error=unitarity, offdiagonal=offdiagonal)


def rwa_model(periodic, states, photons, cutoff=None, sector=DEFAULT_SECTOR, options=None,
              static=None) -> FloquetStates:
    """Few-level rotating-wave model: the Sambe matrix restricted to ``|a_i, n_i>``.

    Q_ij=(H_{n_i-n_j})_{a_i a_j}+n_i omega delta_ij, in static eigenstates of the
    reference Hamiltonian. Choose photon numbers so that E_a+n_a omega are nearly
    degenerate; at resonance the splitting is ``2|(H_{n_a-n_b})_{ab}|``.
    """
    options = _options(options)
    states = [int(a) for a in states]
    photons = [int(n) for n in photons]
    if not states or len(states) != len(photons) or min(states) < 0:
        raise ValueError("give one nonnegative static index and one photon number per level")
    if len(set(zip(states, photons, strict=True))) != len(states):
        raise ValueError("static states and photon numbers must form distinct pairs")
    problem = _prepare(periodic, cutoff, sector, options, static, "eigen", max(states)+1, gap_check=False)
    if max(states) >= problem.dimension:
        raise ValueError("a static state index exceeds the basis dimension")
    omega = periodic.omega
    size = len(states)
    matrix = np.zeros((size, size), complex)
    for i, (a, n) in enumerate(zip(states, photons, strict=True)):
        for j, (b, p) in enumerate(zip(states, photons, strict=True)):
            block = problem.blocks.get(n-p)
            if block is not None:
                matrix[i, j] = block[a, b]
        matrix[i, i] += n*omega
    matrix = (matrix + matrix.conj().T)/2
    energies, vectors = eigh(matrix)
    harmonics = max(abs(n) for n in photons)
    modes = np.zeros((2*harmonics+1, problem.dimension, size), complex)
    for i, (a, n) in enumerate(zip(states, photons, strict=True)):
        modes[n+harmonics, a, :] = vectors[i, :]
    residuals = np.linalg.norm(matrix @ vectors - vectors*energies[None, :], axis=0)
    extra = dict(states=states, photons=photons)
    return _finish(problem, "rwa", energies, modes, residuals, extra,
                   boundary=np.zeros(size), effective=matrix)


__all__ = ["FloquetOptions", "FloquetStates", "solve_floquet", "floquet_propagator", "rwa_model"]
