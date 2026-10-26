"""Tutorial 00: the main README's first calculation, with DMRG and exact diagonalization in the full profile.

    python -B examples/00_first_calculation/run.py                  # quick: QP cutoff two
    python -B examples/00_first_calculation/run.py --profile full   # also DMRG and exact (needs [dmrg])

Results go to runs/examples/00_first_calculation/<profile>/ unless --output is given;
results.json holds the key numbers that tools/check_tutorials.py compares with reference/.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from qdjj_solver import Sector, save_result, solve
from qdjj_solver.common.io import configured, model_from_config, validate_config, write_json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def calculate(name):
    config = validate_config(json.loads((HERE/"input"/f"{name}.json").read_text(encoding="utf-8")))
    h = model_from_config(config, default_compress=False)
    sector = configured(Sector, config["sector"], "sector")
    return solve(h, cutoff=config["cutoff"], sector=sector,
                 options=config["solver"], backend=config["backend"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=("quick", "full"), default="quick",
                        help="quick: QP only; full: also DMRG and exact diagonalization of the finite bath")
    parser.add_argument("--all", action="store_true", help="same as --profile full")
    parser.add_argument("--output", type=Path,
                        help="output directory (default: runs/examples/00_first_calculation/<profile>)")
    parser.add_argument("--save-states", action="store_true", help="also save the wavefunctions, for load_result")
    args = parser.parse_args(argv)
    args.all = args.all or args.profile == "full"
    profile = "full" if args.all else "quick"
    args.output = args.output or ROOT/"runs"/"examples"/HERE.name/profile

    qp = calculate("qp")
    assert abs(qp.energies[0] - (-0.2950595967187083)) < 1e-10
    assert max(qp.residuals) < 1e-9
    assert abs(qp.qp_weights[0].sum() - 1) < 1e-12
    results = {"qp": qp}

    if args.all:
        mps = calculate("dmrg")
        exact = calculate("exact")
        for result in (mps, exact):
            assert result.hamiltonian.fingerprint() == qp.hamiltonian.fingerprint()
            assert result.metadata["sector"] == qp.metadata["sector"]
        np.testing.assert_allclose(mps.energies, exact.energies, rtol=0, atol=1e-9)
        assert mps.metadata["finite_problem_converged"]
        spin = mps.matrix_elements(mps.hamiltonian.observables["impurity_spin_z"])
        results.update(dmrg=mps, exact=exact)
        comparison = dict(
            format_version=1,
            hamiltonian_sha256=qp.hamiltonian.fingerprint(),
            energy_unit="Delta=1",
            qp_minus_exact_lowest_energy=qp.energies[0]-exact.energies[0],
            dmrg_minus_exact_energies=mps.energies-exact.energies,
            dmrg_exact_energy_atol=1e-9,
            dmrg_exact_energy_rtol=0,
            dmrg_finite_problem_converged=mps.metadata["finite_problem_converged"],
            dmrg_impurity_spin_z_matrix=spin,
            impurity_spin_z_unit="hbar",
        )

    for name, result in results.items():
        save_result(result, args.output/f"{name}.json", save_states=args.save_states)
        print(f"{name}: energies={result.energies}; residuals={result.residuals}")
    print("QP weights:", qp.qp_weights)
    summary = dict(qp_energy=float(qp.energies[0]), qp_residual=float(max(qp.residuals)),
                   qp_weight_0=float(qp.qp_weights[0][0]))
    if args.all:
        summary.update(exact_energy_0=float(exact.energies[0]), exact_energy_1=float(exact.energies[1]),
                       dmrg_minus_exact=float(np.max(abs(mps.energies - exact.energies))),
                       qp_minus_exact=float(qp.energies[0] - exact.energies[0]))
    write_json(args.output/"results.json", dict(tutorial=HERE.name, profile=profile, values=summary))
    if args.all:
        write_json(args.output/"comparison.json", comparison)
        print("DMRG and exact energies agree within 1e-9.")
        print("DMRG impurity_spin_z matrix (hbar units):")
        print(spin)
    print(f"Saved results to {args.output}")


if __name__ == "__main__":
    main()
