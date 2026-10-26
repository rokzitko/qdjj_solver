"""Bounded reusable QP preparation and independent finite-window state counts."""

from collections import OrderedDict
from dataclasses import dataclass, field, replace
import sys
import weakref

import numpy as np
from scipy.linalg import ldl
from scipy.sparse import issparse

from ..common.problem import Hamiltonian
from ..common.baths import DiscreteBath
from .basis import DEFAULT_SECTOR, FockBasis
from .execution import ExecutionContext
from .solver import DEFAULT_OPTIONS, ProjectedOperator, preflight, solve


def basis_key(h, cutoff, sector):
    cutoff = h.bath_modes if cutoff is None else cutoff
    return h.nimp, h.spins, cutoff, sector, h.eta_labels


@dataclass
class _Prepared:
    basis: FockBasis
    counts: np.ndarray
    projected: ProjectedOperator | None = None
    matrices: dict = field(default_factory=dict)
    dense_eigenpairs: tuple | None = None
    observables: dict = field(default_factory=dict)
    transitions: dict = field(default_factory=dict)


def hermitian_inertia(matrix, threshold, context=None):
    """Numerical inertia of H-threshold*I via pivoted Hermitian LDL*, not LU signs.

    A nearly singular pivot or poor reconstruction is unresolved. The caller
    brackets its boundary separately; this is not an interval-arithmetic proof.
    """
    context = context or ExecutionContext()
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.isfinite(threshold):
        raise ValueError("inertia needs a finite threshold and square matrix")
    shifted = matrix - threshold*np.eye(len(matrix))
    if not np.all(np.isfinite(shifted)) or not np.allclose(shifted, shifted.conj().T, rtol=0, atol=1e-12):
        raise ValueError("inertia needs a finite Hermitian matrix")
    context.checkpoint("count-factorization", dimension=len(matrix))
    lower, diagonal, _ = ldl(shifted, lower=True, hermitian=True, check_finite=False)
    context.checkpoint("count-validation")
    scale = max(1., float(np.linalg.norm(shifted, ord=np.inf)))
    error = float(np.linalg.norm(shifted-lower@diagonal@lower.conj().T, ord=np.inf))/scale
    tolerance = 64*np.finfo(float).eps*max(1, len(matrix))*scale
    values = []
    index = 0
    while index < len(diagonal):
        if index+1 < len(diagonal) and diagonal[index+1, index] != 0:
            values.extend(np.linalg.eigvalsh(diagonal[index:index+2, index:index+2]))
            index += 2
        else:
            values.append(diagonal[index, index].real)
            index += 1
    values = np.asarray(values)
    if error > 1e-10 or np.any(abs(values) <= tolerance):
        return dict(status="unresolved", reason="singular/uncertain LDL pivot", backward_error=error)
    return dict(status="checked", count=int(np.count_nonzero(values < 0)), backward_error=error)


class QPWorkspace:
    """Reuse QP preparation with least-recently-used eviction and bounded storage.

    Keys distinguish enumeration structure from physical canonical coordinates.
    Changing coefficients can reuse a basis, but only identical coordinates
    authorize a warm start. Returned arrays are read-only snapshots.

    ``max_memory_gib`` bounds accounted live data plus estimated working storage;
    ``max_entries`` limits each cache. ``count_limit`` is the largest dimension for
    an independent LDL state count. ``retained``, when supplied, is a callable
    returning additional numerical data to include in the memory accounting.
    """

    def __init__(self, max_memory_gib=0.5, max_entries=16, count_limit=512, retained=None):
        if not np.isfinite(max_memory_gib) or max_memory_gib <= 0:
            raise ValueError("workspace memory limit must be positive")
        if type(max_entries) is not int or max_entries < 1 or type(count_limit) is not int or count_limit < 0:
            raise ValueError("invalid workspace cache/count limits")
        if retained is not None and not callable(retained):
            raise ValueError("retained must be a callable returning externally owned numerical data")
        self.retained = retained
        self.budget_bytes = int(max_memory_gib*1024**3)
        self.max_entries = max_entries
        self.count_limit = count_limit
        self.bases = OrderedDict()
        self.problems = OrderedDict()
        self.warm = OrderedDict()
        self._live_bases = weakref.WeakValueDictionary()
        self._live_results = weakref.WeakValueDictionary()
        self.stats = dict(basis_builds=0, basis_hits=0, problem_hits=0, warm_starts=0, solves=0, evictions=0)

    @property
    def live_bytes(self):
        """Accounted bytes of caches and live results, including caller-retained arrays."""
        seen = set()

        def size(value):
            if value is None:
                return 0
            if isinstance(value, np.ndarray):
                while isinstance(value.base, np.ndarray):
                    value = value.base
            identity = id(value)
            if identity in seen:
                return 0
            seen.add(identity)
            if isinstance(value, np.ndarray):
                return int(value.nbytes)
            if issparse(value):
                return size(value.data)+size(value.indices)+size(value.indptr)
            if isinstance(value, FockBasis):
                return value.bytes+size(value.spins)
            if isinstance(value, ProjectedOperator):
                return value.bytes+size(value.basis)
            if isinstance(value, Hamiltonian):
                return size(value.operator.terms)+size(value.metadata)+sum(size(o.terms) for o in value.observables.values())
            if isinstance(value, DiscreteBath):
                return size(value.xi)+size(value.weights)+size(value.metadata)
            if isinstance(value, dict):
                return sys.getsizeof(value)+sum(size(k)+size(v) for k, v in value.items())
            if isinstance(value, (list, tuple)):
                return sys.getsizeof(value)+sum(size(v) for v in value)
            return sys.getsizeof(value)

        total = sum(size(b) for b in list(self._live_bases.values()))
        if self.retained is not None:
            total += sum(size(value) for value in self.retained())
        for prepared in self.problems.values():
            total += size(prepared.counts)+size(prepared.projected)+size(prepared.matrices)
            total += size(prepared.dense_eigenpairs)+size(prepared.observables)+size(prepared.transitions)
        for result in list(self._live_results.values()):
            total += size(result.vectors)+size(result.energies)+size(result.qp_weights)+size(result.residuals)
            total += size(result.observables)+size(result.hamiltonian)+size(result.basis)
        return total

    def _ensure(self, additional, keep=None):
        # Results still referenced by a caller stay in live_bytes after eviction.
        while self.live_bytes+additional > self.budget_bytes:
            candidate = next((key for key in self.problems if key != keep), None)
            if candidate is not None:
                del self.problems[candidate]
            elif self.warm:
                self.warm.popitem(last=False)
            elif self.bases:
                self.bases.popitem(last=False)
            else:
                raise MemoryError("workspace live data plus estimated working storage exceeds its memory limit")
            self.stats["evictions"] += 1

    def _limit(self, cache):
        while len(cache) > self.max_entries:
            cache.popitem(last=False)
            self.stats["evictions"] += 1

    def reserve(self, estimated_bytes, context=None):
        """Admit temporary model/operator construction alongside retained data."""
        if type(estimated_bytes) is not int or estimated_bytes < 0:
            raise ValueError("reservation must be a nonnegative byte count")
        self._ensure(estimated_bytes)
        (context or ExecutionContext()).checkpoint("workspace-reservation", estimated_bytes=estimated_bytes,
                                                  retained_bytes=self.live_bytes, budget_bytes=self.budget_bytes)

    def prepare(self, h, cutoff, sector, options=DEFAULT_OPTIONS, context=None):
        """Check resource limits and return the cache key and prepared basis entry for a problem."""
        context = context or ExecutionContext()
        plan = preflight(h, cutoff, sector, options)
        bkey = basis_key(h, cutoff, sector)
        key = (h.fingerprint(), bkey)
        self._ensure(plan["estimated_peak_bytes"], keep=key)
        context.checkpoint("workspace-preflight", **plan, retained_bytes=self.live_bytes,
                           budget_bytes=self.budget_bytes)
        if key in self.problems:
            self.problems.move_to_end(key)
            self.stats["problem_hits"] += 1
            self.stats["basis_hits"] += 1
            return key, self.problems[key]
        basis = self.bases.get(bkey)
        if basis is None:
            context.checkpoint("basis-build", dimension=plan["dimension"])
            basis = FockBasis(h.nimp, h.spins, plan["cutoff"], sector, h.eta_labels, options.max_dimension)
            context.checkpoint("basis-built")
            self._live_bases[id(basis)] = basis
            self.bases[bkey] = basis
            self.stats["basis_builds"] += 1
            self._limit(self.bases)
        else:
            self.bases.move_to_end(bkey)
            self.stats["basis_hits"] += 1
        counts = basis.qp_counts()
        counts.flags.writeable = False
        prepared = _Prepared(basis, counts)
        context.checkpoint("prepared")
        self.problems[key] = prepared
        self._limit(self.problems)
        return key, prepared

    def solve(self, h, cutoff=None, sector=DEFAULT_SECTOR, options=DEFAULT_OPTIONS,
              initial=None, *, context=None, observables=None):
        """Return validated QP eigenstates, reusing compatible completed preparation.

        Arguments follow ``qp_solver.solve``. A compatible prior result supplies
        the starting vector or block when ``initial`` is omitted for an explicitly
        iterative method. ``observables=None`` measures all named model operators.
        Returned numerical arrays are read-only; reuse diagnostics are in ``metadata``.
        """
        context = context or ExecutionContext()
        builds_before = self.stats["basis_builds"]
        key, prepared = self.prepare(h, cutoff, sector, options, context)
        warm_key = (h.coordinate_fingerprint(), basis_key(h, cutoff, sector))
        if initial is None and options.method in ("sparse", "matrix-free", "block") and warm_key in self.warm:
            previous = self.warm[warm_key]
            initial = previous.vectors if options.method == "block" else previous.vectors[:, 0]
            self.stats["warm_starts"] += 1
        plan = preflight(h, cutoff, sector, options, initial)
        self._ensure(plan["estimated_peak_bytes"], keep=key)
        available = max(1, self.budget_bytes-self.live_bytes)/1024**3
        bounded = replace(options, max_memory_gib=min(options.max_memory_gib, available))
        result = solve(h, cutoff, sector, bounded, initial, context=context,
                       observables=observables, _prepared=prepared)
        context.checkpoint("workspace-commit")
        for array in (result.vectors, result.energies, result.qp_weights, result.residuals):
            array.flags.writeable = False
        if prepared.dense_eigenpairs:
            for array in prepared.dense_eigenpairs:
                array.flags.writeable = False
        for array in result.observables.values():
            array.flags.writeable = False
        self._live_results[id(result)] = result
        self.warm[warm_key] = result
        self.warm.move_to_end(warm_key)
        self._limit(self.warm)
        self.stats["solves"] += 1
        result.metadata["reuse"]["basis"] = self.stats["basis_builds"] == builds_before
        result.metadata["workspace"] = dict(self.stats, live_bytes=self.live_bytes)
        return result

    def count_below(self, h, cutoff, sector, edge, margin=1e-8, *, context=None):
        """Count states below edge +/- margin; refuse oversized/uncertain certificates."""
        if not np.isfinite(edge) or not np.isfinite(margin) or margin <= 0:
            raise ValueError("count boundary and positive margin must be finite")
        context = context or ExecutionContext()
        key = h.fingerprint(), basis_key(h, cutoff, sector)
        prepared = self.problems.get(key)
        if prepared is None:
            return dict(status="unavailable", reason="problem not prepared")
        if prepared.dense_eigenpairs is not None:
            energies = prepared.dense_eigenpairs[0]
            return dict(status="checked", method="full-dense-spectrum",
                        below=int(np.count_nonzero(energies < edge-margin)),
                        through=int(np.count_nonzero(energies <= edge+margin)))
        dim = prepared.basis.dimension
        if dim > self.count_limit:
            return dict(status="unavailable", reason="independent count dimension limit")
        if prepared.projected is None:
            return dict(status="unavailable", reason="projected operator not prepared")
        try:
            self._ensure(8*dim*dim*16, keep=key)
        except MemoryError:
            return dict(status="unavailable", reason="independent count memory limit")
        context.checkpoint("count-matrix", dimension=dim)
        matrix = prepared.projected.sparse_matrix(max_nnz=dim*dim).toarray()
        lower = hermitian_inertia(matrix, edge-margin, context)
        upper = hermitian_inertia(matrix, edge+margin, context)
        if lower["status"] != "checked" or upper["status"] != "checked":
            return dict(status="unresolved", reason="boundary count could not be validated")
        return dict(status="checked", method="hermitian-ldl-inertia", below=lower["count"], through=upper["count"],
                    backward_error=max(lower["backward_error"], upper["backward_error"]))

    def project_operator(self, h, cutoff, sector, operator, *, context=None):
        """Cache a canonical-coordinate transition operator in an already prepared basis."""
        from ..common.io import fingerprint
        context = context or ExecutionContext()
        key = h.fingerprint(), basis_key(h, cutoff, sector)
        prepared = self.problems.get(key)
        if prepared is None:
            raise ValueError("prepare the source sector before projecting a transition operator")
        operator_key = fingerprint(operator.to_records())
        if operator_key not in prepared.transitions:
            self._ensure(32*prepared.basis.dimension+512*len(operator.terms), keep=key)
            context.checkpoint("transition-operator")
            projected = ProjectedOperator(operator, prepared.basis)
            context.checkpoint("transition-operator-ready")
            prepared.transitions[operator_key] = projected
        return prepared.transitions[operator_key]

    def clear(self):
        """Release cached preparation and warm starts; caller-owned results remain valid."""
        self.bases.clear()
        self.problems.clear()
        self.warm.clear()
