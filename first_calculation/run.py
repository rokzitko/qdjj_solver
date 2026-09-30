"""Run the README's first calculation; --all adds its DMRG/ED comparison."""

import argparse
import json
from pathlib import Path

import numpy as np

from qdjj_solver import Sector, save_result, solve
from qdjj_solver.common.io import configured, model_from_config, validate_config, write_json


HERE = Path(__file__).resolve().parent


def calculate(name):
    config = validate_config(json.loads((HERE/"input"/f"{name}.json").read_text(encoding="utf-8")))
    h = model_from_config(config, default_compress=False)
    sector = configured(Sector, config["sector"], "sector")
    return solve(h, cutoff=config["cutoff"], sector=sector,
                 options=config["solver"], backend=config["backend"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="also run DMRG and unrestricted finite-bath ED")
    parser.add_argument("--output", type=Path, default=HERE.parent/"results"/"first_calculation",
                        help="output directory (default: results/first_calculation in the source tree)")
    parser.add_argument("--save-states", action="store_true", help="save companion wavefunctions for load_result")
    args = parser.parse_args(argv)

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
    if args.all:
        write_json(args.output/"comparison.json", comparison)
        print("DMRG and exact energies agree within 1e-9.")
        print("DMRG impurity_spin_z matrix (hbar units):")
        print(spin)
    print(f"Saved results to {args.output}")


if __name__ == "__main__":
    main()
