"""Unequal-gap Josephson currents and charge conservation, Zonda et al. Fig. 8."""

from pathlib import Path

import numpy as np

from qdjj_solver import Impurity, make_model
from applications._common import Run, arguments, parity_states, read_csv, scalar

CASE = Path(__file__).resolve().parent


def model(baths, parameters, detuning, phi):
    contacts = [np.sqrt(g/(np.pi*b.rho))*np.eye(2)
                for g, b in zip(parameters["gammas"], baths, strict=True)]
    impurity = Impurity.anderson(parameters["u"], detuning=detuning)
    h = make_model(impurity, baths, contacts, phases=[-phi/2, phi/2],
                   phase_velocities=[-.5, .5])
    for j in range(2):
        other = make_model(impurity, baths, contacts, phases=[-phi/2, phi/2],
                           phase_velocities=np.eye(2)[j])
        h.observables[f"lead_derivative_{j}"] = other.observables["phase_derivative"]
    return h


def point(baths, parameters, detuning, phi, settings):
    states = parity_states(model(baths, parameters, detuning, phi), settings)
    energies = [float(s.energies[0]) for s in states]
    parity = int(energies[1] < energies[0])
    state = states[parity]
    left, right = [2*scalar(state, f"lead_derivative_{j}") for j in range(2)]
    return dict(detuning=detuning, even_energy=energies[0], odd_energy=energies[1],
                signed_gap=energies[1]-energies[0], ground_parity=parity,
                current=2*scalar(state, "phase_derivative"),
                left_current=left, right_current=right, current_sum=left+right,
                charge=scalar(state, "impurity_charge")), states


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw, half_filling=False)
    spec, parameters = config["profiles"][profile], config["model"]
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    reference = read_csv(CASE/"reference"/"figure8.csv")
    rows, checks, comparison = [], [], []
    for settings in settings_list:
        label = settings.get("label", profile)
        bandwidth = settings.get("bandwidth", parameters["bandwidth"])
        for ratio in settings.get("gap_ratios", spec["gap_ratios"]):
            baths = [run.bath(settings, bandwidth=bandwidth, delta=d) for d in (1., ratio)]
            gates = spec["detunings"]
            cache = {}

            def evaluate(x, phi=parameters["phi"], *, baths=baths, cache=cache,
                         settings=settings, label=label, ratio=ratio):
                key = (float(x), float(phi))
                if key not in cache:
                    row, states = point(baths, parameters, x, phi, settings)
                    cache[key] = row
                    for state in states:
                        run.keep(dict(setting=label, gap_ratio=ratio, detuning=x, phi=phi), state)
                return cache[key]

            for x in gates:
                row = evaluate(x)
                rows.append(dict(setting=label, gap_ratio=ratio, **row))
                for method in ("NRG", "FDC_L", "DC_L", "DC_minus_R"):
                    selected = sorted((r for r in reference if float(r["gap_ratio"]) == ratio
                                       and r["method"] == method), key=lambda r: float(r["detuning"]))
                    if selected and float(selected[0]["detuning"])-.001 <= x <= float(selected[-1]["detuning"])+.001:
                        ref = float(np.interp(x, [float(r["detuning"]) for r in selected],
                                              [float(r["current"]) for r in selected]))
                        comparison.append(dict(setting=label, gap_ratio=ratio, detuning=x,
                                               method=method, current=row["current"],
                                               reference_current=ref, difference=row["current"]-ref))
            if profile != "convergence":
                x, step = config["check_detuning"], config["finite_difference_step"]
                center = evaluate(x)
                minus, plus = evaluate(x, parameters["phi"]-step), evaluate(x, parameters["phi"]+step)
                p = ("even", "odd")[center["ground_parity"]]
                fd = (plus[f"{p}_energy"]-minus[f"{p}_energy"])/step
                reflected = evaluate(-x)
                checks.append(dict(gap_ratio=ratio, current_derivative_error=fd-center["current"],
                                   charge_reflection_error=center["charge"]+reflected["charge"]-2,
                                   current_reflection_error=center["current"]-reflected["current"]))
            print(f"Completed {label}, Delta_R/Delta_L={ratio:g}", flush=True)
    summary = dict(points=len(rows), maximum_residual=max(max(s["residuals"]) for s in run.states),
                   maximum_current_conservation_error=max(abs(r["current_sum"]) for r in rows))
    for method in ("NRG", "FDC_L", "DC_L", "DC_minus_R"):
        data = [r for r in comparison if r["method"] == method]
        if data:
            summary[f"maximum_{method}_current_difference"] = max(abs(r["difference"]) for r in data)
    if checks:
        summary["maximum_current_derivative_error"] = max(abs(r["current_derivative_error"]) for r in checks)
    run.finish({"convergence" if profile == "convergence" else "current": rows,
                "checks": checks, "comparison": comparison}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
