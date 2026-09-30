"""Reproduce the two serial-double-dot CPRs of Zonda et al., Fig. 9(a)."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import (Hamiltonian, Impurity, annihilate, create, hermitian_pair,
                         make_model, number)
from applications._common import (Run, arguments, parity_states, read_csv, scalar)

CASE = Path(__file__).resolve().parent


def impurity_operator(u, hopping, detuning=0.):
    charges = [number(2*j)+number(2*j+1) for j in (0, 1)]
    return (sum(u/2*(n-1)*(n-1)+detuning*(n-1) for n in charges)
            -sum(hermitian_pair(hopping*create(s)*annihilate(2+s)) for s in (0, 1)))


def model(bath, u, curve, phi):
    left = np.zeros((2, 4))
    right = np.zeros((2, 4))
    left[:, :2] = np.sqrt(curve["gamma_left"]/(np.pi*bath.rho))*np.eye(2)
    right[:, 2:] = np.sqrt(curve["gamma_right"]/(np.pi*bath.rho))*np.eye(2)
    impurity = Impurity(2, impurity_operator(u, curve["hopping"]),
                       metadata=dict(kind="serial-double-dot", u=u, hopping=curve["hopping"]))
    return make_model(impurity, [bath, bath], [left, right], phases=[-phi/2, phi/2],
                      phase_velocities=[-.5, .5])


def gal_model(u, curve, phi):
    """Eqs. (9), (19)--(22), built as a four-mode effective Hamiltonian."""
    gammas = [curve["gamma_left"], curve["gamma_right"]]
    nu = 1/(1+np.array(gammas))
    charges = [number(2*j)+number(2*j+1) for j in (0, 1)]
    op = sum(u*nu[j]**2/2*(n-1)*(n-1) for j, n in enumerate(charges))
    op -= sum(hermitian_pair(curve["hopping"]*np.sqrt(np.prod(nu))*create(s)*annihilate(2+s))
              for s in (0, 1))
    derivative = 0*op
    for j, velocity in enumerate((-.5, .5)):
        pair = nu[j]*gammas[j]*np.exp(1j*velocity*phi)*create(2*j)*create(2*j+1)
        op += hermitian_pair(pair)
        derivative += hermitian_pair(1j*velocity*pair)
    return Hamiltonian(op, 4, (1, -1, 1, -1), observables={"phase_derivative": derivative},
                       metadata=dict(kind="GAL", u=u, curve=curve, phi=phi))


def point(h, settings):
    states = parity_states(h, settings)
    gap = float(states[1].energies[0]-states[0].energies[0])
    currents = [2*scalar(s, "phase_derivative") for s in states]
    parity = int(gap < 0)
    return dict(even_energy=float(states[0].energies[0]), odd_energy=float(states[1].energies[0]),
                signed_gap=gap, ground_parity=parity, current=currents[parity],
                even_current=currents[0], odd_current=currents[1]), states


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    u = config["model"]["u"]
    rows, gal_rows, crossings, fd_rows = [], [], [], []
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    references = read_csv(CASE/"reference"/"figure9a.csv")
    for settings in settings_list:
        bath = run.bath(settings)
        label = settings.get("label", profile)
        for curve in config["curves"]:
            if "curves" in settings and curve["label"] not in settings["curves"]:
                continue
            phases = settings.get("phase_over_pi", spec.get("phase_over_pi", np.linspace(0, 1, spec.get("phase_points", 3))))
            if spec.get("include_reference_phases"):
                phases = sorted(set(phases) | {float(r["phi_over_pi"]) for r in references
                                              if r["curve"] == curve["label"]})
            cache = {}

            def evaluate(p, *, cache=cache, bath=bath, curve=curve, settings=settings, label=label):
                p = float(p)
                if p not in cache:
                    result, states = point(model(bath, u, curve, p*np.pi), settings)
                    cache[p] = result
                    for state in states:
                        run.keep(dict(setting=label, curve=curve["label"], phi_over_pi=p), state)
                return cache[p]

            for p in phases:
                result = evaluate(p)
                rows.append(dict(setting=label, curve=curve["label"], phi_over_pi=float(p), **result))
            if profile != "convergence":
                for p in np.linspace(0, 1, 201):
                    result, _ = point(gal_model(u, curve, p*np.pi), {})
                    gal_rows.append(dict(curve=curve["label"], phi_over_pi=float(p), **result))
                if spec.get("locate_crossings"):
                    for a, b in zip(phases[:-1], phases[1:], strict=True):
                        if evaluate(a)["signed_gap"]*evaluate(b)["signed_gap"] < 0:
                            root = brentq(lambda p: evaluate(p)["signed_gap"], a, b, xtol=1e-6)
                            crossings.append(dict(curve=curve["label"], phi_over_pi=root,
                                                  signed_gap=evaluate(root)["signed_gap"]))
                step = config["finite_difference_step"]
                p = config["finite_difference_phi_over_pi"]
                center = evaluate(p)
                for parity, key in ((0, "even"), (1, "odd")):
                    minus, plus = evaluate(p-step/np.pi), evaluate(p+step/np.pi)
                    # J/J0 = 2 dE/dphi; the 2 cancels the centered denominator.
                    fd = (plus[f"{key}_energy"]-minus[f"{key}_energy"])/step
                    fd_rows.append(dict(curve=curve["label"], parity=parity,
                                        current=center[f"{key}_current"], finite_difference=fd,
                                        error=abs(fd-center[f"{key}_current"])))
            print(f"Completed {label}, {curve['label']}", flush=True)
    summary = dict(points=len(rows), maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    if profile != "convergence":
        comparisons = []
        for reference in references:
            p = float(reference["phi_over_pi"])
            matches = [r for r in rows if r["curve"] == reference["curve"] and abs(r["phi_over_pi"]-p) < 1e-9]
            if matches:
                r = matches[0]
                comparisons.append(dict(curve=r["curve"], phi_over_pi=p, current=r["current"],
                                         nrg_current=float(reference["current"]),
                                         difference=r["current"]-float(reference["current"])))
        summary.update(crossings=crossings,
                       maximum_current_finite_difference_error=max(r["error"] for r in fd_rows))
        for curve in config["curves"]:
            differences = [r["difference"] for r in comparisons if r["curve"] == curve["label"]]
            if differences:
                summary[f"{curve['label']}_maximum_nrg_current_error"] = max(abs(x) for x in differences)
        run.finish({"current": rows, "gal": gal_rows, "crossings": crossings,
                    "current_check": fd_rows, "comparison": comparisons}, summary)
    else:
        run.finish({"convergence": rows}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
