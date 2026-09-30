"""Classic superconductivity/Kondo competition at the published bare parameters."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import reference_model
from applications._common import Run, arguments, parity_states, read_csv, scalar

CASE = Path(__file__).resolve().parent


def parameters(normal, delta_over_tk):
    u, gamma = normal["u_over_D"], normal["gamma_total_over_D"]
    # Eq. (4), evaluated at epsilon_d=-U/2. No empirical Kondo-scale rescaling.
    tk = np.sqrt(u*gamma/2)*np.exp(-np.pi*u/(8*gamma))
    delta = delta_over_tk*tk
    return dict(tk_over_D=float(tk), delta_over_D=float(delta), u=float(u/delta),
                gamma=float(gamma/delta), bandwidth=float(1/delta))


def model(bath, scales, phi):
    return reference_model(bath, u=scales["u"], gamma=scales["gamma"], phi=phi,
                           detuning=0., symmetry=True, compress=True)


def point(bath, scales, phi, settings):
    states = parity_states(model(bath, scales, phi), settings)
    energies = [float(s.energies[0]) for s in states]
    gap = energies[1]-energies[0]
    parity = int(gap < 0)
    current = 2*scalar(states[parity], "phase_derivative")
    even_current = 2*scalar(states[0], "phase_derivative") if abs(abs(phi)-np.pi) > 1e-12 else None
    return dict(phi_over_pi=phi/np.pi, even_energy=energies[0], odd_energy=energies[1],
                signed_gap=gap, ground_parity=parity, current=current,
                even_current=even_current, odd_current=2*scalar(states[1], "phase_derivative")), states


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    references = read_csv(CASE/"reference"/"figure3.csv")
    rows, crossings, checks, comparisons = [], [], [], []
    for original in settings_list:
        label = original.get("label", profile)
        for ratio in original.get("ratios", spec["ratios"]):
            scales = parameters(config["normal_state"], ratio)
            settings = original | {"frequency_cutoff": scales["bandwidth"]*original["frequency_factor"]}
            bath = run.bath(settings, bandwidth=scales["bandwidth"])
            phases = original.get("phase_over_pi", spec.get("phase_over_pi", np.linspace(0, 1, spec.get("phase_points", 3))))
            selected_ref = [r for r in references if float(r["delta_over_tk"]) == ratio and
                            0 <= float(r["phi_over_pi"]) <= 1]
            if spec.get("include_reference_phases"):
                phases = sorted(set(phases) | {float(r["phi_over_pi"]) for r in selected_ref})
            cache = {}

            def evaluate(p, *, cache=cache, bath=bath, scales=scales, settings=settings, label=label, ratio=ratio):
                p = float(p)
                if p not in cache:
                    result, states = point(bath, scales, p*np.pi, settings)
                    cache[p] = result
                    for state in states:
                        run.keep(dict(setting=label, delta_over_tk=ratio, phi_over_pi=p, **scales), state)
                return cache[p]

            for p in phases:
                rows.append(dict(setting=label, delta_over_tk=ratio, **evaluate(p)))
            if profile != "convergence":
                p = config["finite_difference_phase_over_pi"]
                step = config["finite_difference_step"]
                center, minus, plus = evaluate(p), evaluate(p-step/np.pi), evaluate(p+step/np.pi)
                for parity in ("even", "odd"):
                    fd = (plus[f"{parity}_energy"]-minus[f"{parity}_energy"])/step
                    checks.append(dict(delta_over_tk=ratio, parity=parity, current=center[f"{parity}_current"],
                                       energy_derivative_current=fd, error=fd-center[f"{parity}_current"]))
                if spec.get("locate_crossings"):
                    for a, b in zip(phases[:-1], phases[1:], strict=True):
                        if evaluate(a)["signed_gap"]*evaluate(b)["signed_gap"] < 0:
                            root = brentq(lambda p: evaluate(p)["signed_gap"], a, b,
                                          xtol=config["transition_phase_tolerance"])
                            crossings.append(dict(delta_over_tk=ratio, phi_critical_over_pi=root))
                for reference in selected_ref:
                    p = float(reference["phi_over_pi"])
                    if p in cache:
                        comparisons.append(dict(delta_over_tk=ratio, phi_over_pi=p,
                                                 current=cache[p]["current"], nrg_current=float(reference["current"]),
                                                 difference=cache[p]["current"]-float(reference["current"])))
            print(f"Completed {label}, Delta/T_K={ratio:g}", flush=True)
    summary = dict(tk_over_D=parameters(config["normal_state"], 1.)["tk_over_D"],
                   points=len(rows), crossings=crossings,
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    if checks:
        summary["maximum_current_derivative_error"] = max(abs(r["error"]) for r in checks)
    if comparisons:
        summary["maximum_published_current_difference"] = max(abs(r["difference"]) for r in comparisons)
    run.finish({"convergence" if profile == "convergence" else "current": rows,
                "crossings": crossings, "current_check": checks, "comparison": comparisons}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
