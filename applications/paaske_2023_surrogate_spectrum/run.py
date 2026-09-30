"""Published surrogate-bath parity effects in the coupled subgap spectrum."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import reference_model
from qdjj_solver.common.baths import discrete_g, hybridization_g
from applications._common import Run, arguments, parity_states, read_csv, read_json, scalar

CASE = Path(__file__).resolve().parent


def published_settings(levels):
    return dict(levels=levels, frequency_cutoff=10.,
                bath_record=read_json(CASE/"reference"/"baths.json")[str(levels)])


def model(bath, u, gamma):
    # At phi=0 the even lead combination couples with TOTAL Gamma. The odd
    # combination is a spectator; removing it isolates the paper's Dg branch.
    return reference_model(bath, u=u, gamma=gamma, geometry="single", detuning=0.,
                           symmetry=False, compress=False)


def point(bath, u, gamma, settings):
    states = parity_states(model(bath, u, gamma), settings)
    even, odd = [float(s.energies[0]) for s in states]
    gap = odd-even
    return dict(gamma=gamma, even_energy=even, odd_energy=odd, signed_gap=gap,
                singlet_excitation=max(-gap, 0.), doublet_excitation=max(gap, 0.),
                doublet_dot_spin=scalar(states[1], "impurity_spin_z")), states


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    u = config["model"]["u"]
    references = read_csv(CASE/"reference"/"figure3.csv")
    if profile == "convergence":
        settings_list = spec["settings"]
        gammas = spec["gamma"]
    else:
        settings_list = [dict(label=f"L{n}", published=True, levels=n, frequency_cutoff=10.) for n in spec["levels"]]
        gammas = spec.get("gamma", np.geomspace(spec.get("gamma_min", .4), spec.get("gamma_max", 20),
                                                spec.get("gamma_points", 3)))
    rows, fits, crossings, published_comparisons = [], [], [], []
    for original in settings_list:
        settings = dict(original)
        if settings.get("published"):
            settings.update(published_settings(settings["levels"]))
        bath = run.bath(settings)
        label = settings["label"]
        omega = np.geomspace(.001, 10, 1000)
        fits.append(dict(setting=label, max_relative_fit_error=float(np.max(abs(discrete_g(bath, omega)/hybridization_g(omega, bandwidth=10)-1))),
                         min_bath_qp_energy=float(min(bath.energies))))
        cache = {}

        def evaluate(gamma, *, cache=cache, bath=bath, settings=settings, label=label):
            gamma = float(gamma)
            if gamma not in cache:
                row, states = point(bath, u, gamma, settings)
                cache[gamma] = row
                for state in states:
                    run.keep(dict(setting=label, gamma=gamma, channel="gerade"), state)
            return cache[gamma]

        scan = list(gammas)
        if profile == "paper" and settings["levels"] == 6:
            # Solve directly at selected original NRG vertices, not an interpolated
            # approximation to our own calculation near the transition.
            nrg = [r for r in references if r["method"] == "NRG"]
            scan += [float(r["gamma"]) for r in nrg[::8]]
        if profile == "paper":
            source_points = [r for r in references if r["method"] == label and float(r["excitation"]) > .001][::12]
            for reference in source_points:
                result = evaluate(float(reference["gamma"]))
                excitation = result["singlet_excitation" if reference["branch"] == "S" else "doublet_excitation"]
                published_comparisons.append(dict(setting=label, gamma=result["gamma"], branch=reference["branch"],
                                                  excitation=excitation, published_excitation=float(reference["excitation"]),
                                                  difference=excitation-float(reference["excitation"])))
        for gamma in sorted(set(scan)):
            rows.append(dict(setting=label, **evaluate(gamma)))
        if profile != "convergence":
            critical = brentq(lambda g: evaluate(g)["signed_gap"], 1., 5., xtol=1e-7)
            crossings.append(dict(setting=label, gamma_critical=critical))
        print(f"Completed {label}", flush=True)

    comparisons = []
    if profile == "paper":
        for row in rows:
            if row["setting"] != "L6":
                continue
            branch = "Dg" if row["signed_gap"] > 0 else "S"
            reference = [r for r in references if r["method"] == "NRG" and r["branch"] == branch]
            x = [float(r["gamma"]) for r in reference]
            if min(x) <= row["gamma"] <= max(x):
                value = float(np.interp(np.log(row["gamma"]), np.log(x), [float(r["excitation"]) for r in reference]))
                excitation = abs(row["signed_gap"])
                comparisons.append(dict(gamma=row["gamma"], branch=branch, excitation=excitation,
                                         nrg_excitation=value, difference=excitation-value))
    summary = dict(points=len(rows), crossings=crossings,
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    if comparisons:
        summary["L6_maximum_nrg_excitation_error"] = max(abs(r["difference"]) for r in comparisons)
    if published_comparisons:
        summary["maximum_published_surrogate_error"] = max(abs(r["difference"]) for r in published_comparisons)
    run.finish({"convergence" if profile == "convergence" else "spectrum": rows,
                "hybridization": fits, "crossings": crossings, "comparison": comparisons,
                "published_surrogates": published_comparisons}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
