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
