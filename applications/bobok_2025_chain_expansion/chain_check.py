"""Independent physical-electron-chain DMRG checks of the common-lead ChE model."""

import argparse
from pathlib import Path

import numpy as np

from qdjj_solver import Hamiltonian, Sector, annihilate, create, hermitian_pair, number
from applications._common import Run, eigenstate, read_json, scalar
from .models import double_dot, spin_operators, spin_square

CASE = Path(__file__).resolve().parent


def chain_model(bath, u, gamma_per_dot):
    """Use the paper's chain directly; physical-electron bath occupations are untruncated."""
    coefficients = np.array(bath.metadata["chain_h"])
    charges = [number(2*j)+number(2*j+1) for j in (0, 1)]
    op = sum(u/2*(n-1)*(n-1) for n in charges)
    op -= float(np.sum(bath.xi-bath.energies))
    for j in range(bath.levels):
        up, down = 4+2*j, 5+2*j
        op -= hermitian_pair(bath.delta*create(up)*create(down))
        for spin in (0, 1):
            if j:
                op -= hermitian_pair(bath.delta*np.sqrt(coefficients[j])*create(up+spin)*annihilate(up-2+spin))
            else:
                for dot in (0, 1):
                    op -= hermitian_pair(np.sqrt(gamma_per_dot*bath.delta*coefficients[0])
                                         *create(up+spin)*annihilate(2*dot+spin))
    z0, plus0 = spin_operators(0)
    z1, plus1 = spin_operators(1)
    observables = dict(total_spin_squared=spin_square(2+bath.levels),
                       spin_correlation=z0*z1+hermitian_pair(plus0*plus1.dagger())/2,
                       pairing=sum(hermitian_pair(create(2*j)*create(2*j+1))/4 for j in (0, 1)))
    return Hamiltonian(op, 4, (1, -1)*(2+bath.levels), observables=observables,
                        metadata=dict(kind="ChE-electron-chain-DQD", bath=bath.record(), u=u,
                                      gamma_per_dot=gamma_per_dot, restriction="none",
                                      energy_reference="isolated BCS reservoir subtracted"))


def calculate(config, output):
    run = Run(CASE, config, "dmrg", output)
    rows = []
    for settings in config["settings"]:
        bath = run.bath(settings)
        for u in config["u"]:
            star, chain = double_dot(bath, u, .5), chain_model(bath, u, .5)
            for sector in (Sector(0, 0), Sector(1, 1), Sector(0, 2)):
                exact = eigenstate(star, sector, {})
                mps = eigenstate(chain, sector, config["dmrg"])
                run.keep(dict(setting=settings["label"], u=u, coordinates="qp-star"), exact)
                run.keep(dict(setting=settings["label"], u=u, coordinates="electron-chain"), mps)
                rows.append(dict(setting=settings["label"], u=u, parity=sector.parity, twice_sz=sector.twice_sz,
                                 energy_error=float(mps.energies[0]-exact.energies[0]),
                                 pairing_error=scalar(mps, "pairing")-scalar(exact, "pairing"),
                                 spin_correlation_error=scalar(mps, "spin_correlation")-scalar(exact, "spin_correlation"),
                                 total_spin_squared_error=scalar(mps, "total_spin_squared")-scalar(exact, "total_spin_squared")))
            print(f"Completed independent DMRG {settings['label']}, U={u:g}", flush=True)
    summary = {f"maximum_{key}": max(abs(r[key]) for r in rows) for key in
               ("energy_error", "pairing_error", "spin_correlation_error", "total_spin_squared_error")}
    summary["maximum_residual"] = max(float(max(s["residuals"])) for s in run.states)
    run.finish({"comparison": rows}, summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=CASE/"input"/"dmrg_check.json")
    parser.add_argument("--output", type=Path, default=CASE/"runs"/"dmrg")
    args = parser.parse_args()
    config = read_json(args.input)
    if config.get("format_version") != 1:
        parser.error("unsupported application input version")
    calculate(config, args.output)


if __name__ == "__main__":
    main()
