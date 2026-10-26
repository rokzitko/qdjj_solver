"""A bounded DMRG bond ladder, compared with full same-bath diagonalization.

python examples/11_dmrg_convergence/run.py --profile quick --output runs/dmrg-ladder
"""

import argparse
from dataclasses import asdict
from pathlib import Path

from qdjj_solver import Sector, cosh_grid, load_result, reference_model, save_result
from qdjj_solver.common.io import numerical_environment, write_json
from qdjj_solver.dmrg_solver import SolverOptions, solve
from qdjj_solver.dmrg_solver.convergence import converge_bond_dimension
from qdjj_solver.qp_solver import solve as exact_solve

HERE = Path(__file__).resolve().parent


def compute(profile, output):
    h = reference_model(cosh_grid(pairs=1, delta=1., bandwidth=5.),
                        u=1.6, gamma=.3, phi=.8, symmetry=False, compress=False)
    sector = Sector(1, 1)
    exact = exact_solve(h, cutoff=None, sector=sector)
    rows = []

    def archive(result, entry):
        row = {**entry, "max_actual_bond": max(max(state.chi, default=1) for state in result.states),
               "charge": result.observables["impurity_charge"],
               "gram_error": result.metadata["gram_error"]}
        rows.append(row)
        write_json(output/"history.json", dict(rows=rows, exact_energy=exact.energies[0]))
        print(f"chi={entry['chi_max']:2d} E={result.energies[0]:.12f} "
              f"r={max(result.residuals):.3e} converged={entry['finite_problem_converged']}")

    # BEGIN handbook: dmrg-ladder
    options = SolverOptions(eigenpairs=1, seed_trials=1, seed=1729,
                            max_sweeps=40, min_sweeps=8, mixer_sweeps=4,
                            lanczos_probability_tolerance=1e-22)
    chis = [2, 4, 16, 32] if profile == "quick" else [2, 4, 8, 16, 32, 64]
    result, assessment = converge_bond_dimension(
        h, chis, sector=sector, options=options, stable_steps=1 if profile == "quick" else 2,
        energy_tolerance=1e-9, observable_tolerance=1e-8, callback=archive)
    if not assessment["empirical_convergence"]:
        raise RuntimeError("bond ladder did not converge; inspect history.json")
    save_result(result, output/"state.json", save_states=True)
    restored = load_result(output/"state.json")
    continued = solve(h, sector=sector, options=SolverOptions(
        eigenpairs=1, chi_max=chis[-1], seed_trials=1, require_convergence=True,
        lanczos_probability_tolerance=1e-22), initial=restored)
    # END handbook: dmrg-ladder
    values = dict(energy=float(result.energies[0]),
                  ed_error=float(abs(result.energies[0]-exact.energies[0])),
                  residual=float(max(result.residuals)),
                  continuation_error=float(abs(continued.energies[0]-result.energies[0])),
                  empirical_convergence=int(assessment["empirical_convergence"]),
                  coarse_error=float(abs(rows[0]["energies"][0]-exact.energies[0])))
    assert not rows[0]["finite_problem_converged"], "the coarse teaching step should expose truncation"
    write_json(output/"results.json", dict(tutorial=HERE.name, profile=profile, values=values,
               environment=numerical_environment(), sector=asdict(sector), options=asdict(options),
               hamiltonian=h.record(), assessment=assessment))
    plot(rows, float(exact.energies[0]), output/"convergence.png")
    return values


def plot(rows, exact, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(6.2, 2.6), layout="constrained")
    x = [r["chi_max"] for r in rows]
    axes[0].semilogy(x, [max(abs(r["energies"][0]-exact), 1e-14) for r in rows], "o-", color="#2a78d6")
    axes[1].semilogy(x, [max(max(r["residuals"]), 1e-14) for r in rows], "o-", color="#eb6834")
    axes[0].set_ylabel(r"$|E-E_{\mathrm{ED}}|/\Delta$")
    axes[1].set_ylabel(r"$\|(H-E)\psi\|/\Delta$")
    for axis in axes:
        axis.set_xlabel(r"bond limit $\chi_{\max}$")
        axis.grid(alpha=.25)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("quick", "full"), default="quick")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--figure", type=Path)
    args = parser.parse_args()
    output = args.output or HERE.parents[1]/"runs"/"examples"/HERE.name/args.profile
    if args.plot_only:
        import json
        record = json.loads((output/"history.json").read_text())
        plot(record["rows"], record["exact_energy"], args.figure or output/"convergence.png")
    else:
        compute(args.profile, output)


if __name__ == "__main__":
    main()
