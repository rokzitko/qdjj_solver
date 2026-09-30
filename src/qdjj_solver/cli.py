"""Unified, versioned command-line interface for finite impurity calculations."""

import argparse
import json
from pathlib import Path
import sys

from . import __version__
from .api import save_result, solve
from .common.io import configured, model_from_config, validate_config
from .common.problem import Sector


def main(argv=None):
    parser = argparse.ArgumentParser(description="QP and DMRG solvers for superconducting quantum impurities")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("solve", help="solve a finite Hamiltonian")
    command.add_argument("configuration", type=Path)
    command.add_argument("--backend", choices=("qp", "dmrg"))
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--save-states", action="store_true")
    command.add_argument("--resume", type=Path, help="DMRG result JSON with an MPS checkpoint")
    args = parser.parse_args(argv)
    try:
        config = validate_config(json.loads(args.configuration.read_text(encoding="utf-8")))
        backend = args.backend or config.get("backend", "qp")
        if backend not in ("qp", "dmrg"):
            raise ValueError("backend must be qp or dmrg")
        cutoff = config.get("cutoff", 2 if backend == "qp" else None)
        if backend == "dmrg" and cutoff is not None:
            raise ValueError("DMRG requires cutoff=null or an omitted cutoff")
        if backend == "qp":
            from .qp_solver import SolverOptions
        else:
            from .dmrg_solver import SolverOptions
        sector = configured(Sector, config.get("sector", {}), "sector")
        options = configured(SolverOptions, config.get("solver", {}), "solver")
        # Backend selection changes the numerical approximation, never which
        # physical modes are present. Legacy qpsolver keeps its older default.
        h = model_from_config(config, default_compress=False)
        initial = None
        if backend == "dmrg":
            from .dmrg_solver.io import load_result
            if args.resume:
                initial = load_result(args.resume)
        elif args.resume:
            raise ValueError("--resume currently accepts DMRG checkpoints")
        result = solve(h, cutoff, sector, options, initial, backend=backend)
        save_result(result, args.output, args.save_states)
        print(f"Backend: {backend}; energy reference: {h.metadata.get('energy_reference', 'as specified')}")
        for i, energy in enumerate(result.energies):
            residual = "not measured" if result.residuals is None else f"{result.residuals[i]:.3e}"
            print(f"State {i}: E={energy:.15g}; residual={residual}")
        if backend == "dmrg":
            print(f"Finite-problem checks passed: {result.metadata['finite_problem_converged']}")
        print(f"Saved {args.output}")
        return 0
    except (ValueError, TypeError, KeyError, OSError, RuntimeError, MemoryError, ImportError) as error:
        print(f"Calculation failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
