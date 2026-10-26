"""Reference-junction analysis; run from any directory with the installed package.

python examples/handbook/static.py --output runs/handbook/static
"""

import argparse
import json
from pathlib import Path

import numpy as np

from qdjj_solver import Sector, load_result, save_result, solve
from qdjj_solver.common.io import model_from_config, write_json

INPUT = Path(__file__).resolve().parents[1]/"01_subgap_states"/"input"/"reference.json"


def calculate(config, output):
    # BEGIN handbook: static
    h = model_from_config(config, default_compress=False)
    even = solve(h, cutoff=config["cutoff"], sector=Sector(0, 0), options={"eigenpairs": 2})
    odd = solve(h, cutoff=config["cutoff"], sector=Sector(1, 1))
    partner = solve(h, cutoff=config["cutoff"], sector=Sector(1, -1))
    signed_gap = float(odd.energies[0] - even.energies[0])
    pair_gap = float(even.energies[1] - even.energies[0])
    addition = odd.matrix_elements("d_up_dagger", even)
    pair = even.matrix_elements("current")
    current = float(np.real(even.observables["phase_derivative"][0]))
    save_result(even, output/"even.json", save_states=True)
    restored = load_result(output/"even.json")
    np.testing.assert_allclose(restored.energies, even.energies, rtol=0, atol=1e-13)
    # END handbook: static
    values = dict(signed_gap=signed_gap, pair_gap=pair_gap, ground_parity=int(signed_gap < 0),
                  current=current, addition_weight=float(abs(addition[0, 0])**2),
                  pair_weight=float(abs(pair[1, 0])**2),
                  doublet_splitting=float(abs(partner.energies[0]-odd.energies[0])),
                  max_residual=float(max(*even.residuals, *odd.residuals)))
    write_json(output/"summary.json", dict(input=config, values=values))
    for name, value in values.items():
        print(f"{name:22s} {value:.10g}")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=Path("results/handbook-static"))
    args = parser.parse_args()
    calculate(json.loads(args.input.read_text()), args.output)


if __name__ == "__main__":
    main()
