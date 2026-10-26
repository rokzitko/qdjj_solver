"""Match the explicit two-orbital Python Hamiltonian to custom.json."""

import argparse
import json
from pathlib import Path

import numpy as np

from qdjj_solver import Impurity, Sector, cosh_grid, make_model, number, save_result, solve
from qdjj_solver.common.io import model_from_config


def calculate(output):
    # BEGIN handbook: custom
    n = [number(i) for i in range(4)]
    impurity = Impurity(2, -.3*(n[0]+n[1]) + .2*(n[2]+n[3]) + n[0]*n[1])
    bath = cosh_grid(pairs=1, delta=1., bandwidth=5.)
    tunneling = np.array([[.3, 0., .1, 0.], [0., .3, 0., .1]])
    h = make_model(impurity, [bath], [tunneling], phases=[0.])
    result = solve(h, cutoff=2, sector=Sector(0, 0), options={"eigenpairs": 2})
    save_result(result, output, save_states=True)
    # END handbook: custom
    config = json.loads(Path(__file__).with_name("custom.json").read_text())
    from_json = model_from_config(config, default_compress=False)
    assert h.operator.terms == from_json.operator.terms
    check = solve(from_json, cutoff=2, sector=Sector(0, 0), options=config["solver"])
    np.testing.assert_allclose(result.energies, check.energies, rtol=0, atol=1e-12)
    print(result.energies)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/custom.json"))
    calculate(parser.parse_args().output)
