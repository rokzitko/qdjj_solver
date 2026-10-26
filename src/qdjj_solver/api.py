"""Lazy dispatch between numerical backends sharing the same Hamiltonian."""

from .common.problem import DEFAULT_SECTOR


def solve(hamiltonian, cutoff=None, sector=DEFAULT_SECTOR, options=None, initial=None, *, backend="qp"):
    if backend == "qp":
        from .qp_solver import SolverOptions, solve as backend_solve
    elif backend == "dmrg":
        from .dmrg_solver import SolverOptions, solve as backend_solve
    else:
        raise ValueError("backend must be qp or dmrg")
    if options is None:
        options = SolverOptions()
    elif isinstance(options, dict):
        options = SolverOptions(**options)
    if not isinstance(options, SolverOptions):
        raise ValueError("options do not belong to the selected backend")
    return backend_solve(hamiltonian, cutoff, sector, options, initial)


def save_result(result, path, save_states=False):
    """Write a common version-1 record and optionally a typed state artifact."""
    backend = getattr(result, "backend", None)
    if backend == "qp":
        from .qp_solver.io import save_result as save
        save(result, path, save_vectors=save_states, common_format=True)
    elif backend == "dmrg":
        from .dmrg_solver.io import save_result as save
        save(result, path, save_states=save_states)
    else:
        raise ValueError("save_result requires QP or DMRG Eigenstates")


def load_result(path):
    """Restore saved states, selecting the backend from the versioned record.

    A companion state artifact is required. Legacy QP records are also readable.
    Loading QP states does not import TeNPy.
    """
    import json
    from pathlib import Path
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(record, dict) or record.get("format_version") != 1:
        raise ValueError("unsupported result format")
    if record.get("format") == "bcs-qp-eigenstates":
        backend = "qp"
    elif record.get("format") == "qdjj-eigenstates":
        backend = record.get("backend")
    else:
        raise ValueError("unsupported result format")
    if backend == "qp":
        from .qp_solver.io import load_result as load
    elif backend == "dmrg":
        from .dmrg_solver.io import load_result as load
    else:
        raise ValueError("unsupported result backend")
    return load(path)


def orbital_density_matrix(result, orbital=0, *, state=0, weights=None, physical=True):
    """Physical impurity-orbital RDM for QP or DMRG states.

    The basis is [empty, up, down, up-down]. state is a zero-based eigenstate
    index; state=None returns all RDMs with shape (k,4,4). With state=None and
    nonnegative weights of length k, return their normalized incoherent mixture.
    Average complete degenerate multiplets for a basis-independent result.

    physical=True uses the declared electron-to-canonical map, including verified
    impurity-only eta compression. For a manually constructed canonical-mode
    Hamiltonian without an electron map use physical=False explicitly.
    Composite projectors are formed before QP projection. Local odd/even
    coherences vanish in the solvers' fixed total parity sectors.
    """
    import numpy as np
    from .common.density import orbital_density_operators

    backend = getattr(result, "backend", None)
    if backend not in ("qp", "dmrg"):
        raise ValueError("orbital_density_matrix requires QP or DMRG Eigenstates")
    if getattr(result, "hamiltonian", None) is None:
        raise ValueError("orbital density requires the result Hamiltonian")
    count = len(result.energies)
    if state is not None and (isinstance(state, (bool, np.bool_))
                              or not isinstance(state, (int, np.integer)) or not 0 <= state < count):
        raise ValueError("state must be None or a valid eigenstate index")
    if weights is not None:
        if state is not None:
            raise ValueError("weights require state=None")
        if np.iscomplexobj(weights):
            raise ValueError("weights must be real")
        weights = np.asarray(weights, dtype=float)
        if (weights.shape != (count,) or not np.all(np.isfinite(weights))
                or np.any(weights < 0) or weights.sum() <= 0 or not np.isfinite(weights.sum())):
            raise ValueError("weights must be finite nonnegative values with a positive sum, one per state")
        weights = weights/weights.sum()
    indices = range(count) if state is None else (int(state),)
    operators = orbital_density_operators(orbital, hamiltonian=result.hamiltonian, physical=physical)
    density = np.zeros((len(indices), 4, 4), complex)
    for (row, column), op in operators.items():
        if row > column:
            continue
        if backend == "qp":
            from .qp_solver import ProjectedOperator
            projected = ProjectedOperator(op, result.basis)
            values = np.array([projected.expectation(result.vectors[:, i]) for i in indices])
        else:
            values = np.array([result.prepared.matrix_element(op, result.states[i], result.states[i])
                               for i in indices])
        density[:, row, column] = values
        if row != column:
            density[:, column, row] = values.conj()
    if weights is not None:
        return np.einsum("k,kij->ij", weights, density)
    return density if state is None else density[0]
