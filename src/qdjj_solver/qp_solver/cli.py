"""Command-line calculations and bath fitting."""

import argparse
import json
from pathlib import Path
import sys

from . import __version__
from ..common.baths import fit_surrogate
from .basis import dimension
from .io import calculation_from_config, save_result, write_json
from .solver import solve


def main(argv=None):
    parser = argparse.ArgumentParser(description="Variational QP expansions for superconducting impurities")
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("solve", "dimension"):
        command = commands.add_parser(name, help="solve a Hamiltonian" if name == "solve" else "count the requested Hilbert space")
        command.add_argument("configuration", type=Path)
        if name == "solve":
            command.add_argument("--output", type=Path, required=True)
            command.add_argument("--save-vectors", action="store_true")
    command = commands.add_parser("fit", help="fit a particle-hole-symmetric BCS surrogate")
    command.add_argument("--levels", type=int, required=True)
    command.add_argument("--delta", type=float, default=1.)
    command.add_argument("--bandwidth", type=float, default=100.)
    command.add_argument("--frequency-cutoff", type=float, default=100.)
    command.add_argument("--starts", type=int, default=4)
    command.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "fit":
            bath = fit_surrogate(args.levels, args.delta, args.bandwidth, args.frequency_cutoff, starts=args.starts)
            write_json(args.output, bath.record())
            print(f"Saved {bath.levels}-level surrogate to {args.output}")
            print(f"Maximum relative g-function error on fitting mesh: {bath.metadata['max_relative_fit_error']:.6g}")
            return 0
        config = json.loads(args.configuration.read_text(encoding="utf-8"))
        h, cutoff, sector, options = calculation_from_config(config)
        if args.command == "dimension":
            q = h.bath_modes if cutoff is None else cutoff
            h.check_sector(sector)
            dim = dimension(h.nimp, h.spins, q, sector, h.eta_labels)
            print(f"Impurity spin orbitals: {h.nimp}; bath QP modes: {h.bath_modes}")
            print(f"QP cutoff: {q}; symmetry sector: {sector}")
            print(f"Hilbert-space dimension: {dim:,}")
            return 0
        result = solve(h, cutoff, sector, options)
        save_result(result, args.output, args.save_vectors)
        print(f"Hilbert-space dimension: {result.basis.dimension:,}; QP cutoff: {result.basis.cutoff}")
        print(f"Energy reference: {h.metadata.get('energy_reference', 'as specified in Hamiltonian')}")
        for i, (energy, residual) in enumerate(zip(result.energies, result.residuals, strict=True)):
            print(f"State {i}: E={energy:.15g}; Ritz residual={residual:.3e}")
            print("  QP probabilities: " + ", ".join(f"P({q})={w:.7g}" for q, w in enumerate(result.qp_weights[i])))
            for name, values in result.observables.items():
                print(f"  {name}={values[i]:.12g}")
        print(f"Calculation time: {result.timings['total']:.6g} s")
        print(f"Saved {args.output}")
        return 0
    except (ValueError, TypeError, KeyError, OSError, RuntimeError, MemoryError) as error:
        print(f"Calculation failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
