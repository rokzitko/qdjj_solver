"""Symmetric/asymmetric phase boundaries and exact current mapping."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import Impurity, make_model, reference_model
from applications._common import Run, arguments, parity_states, read_csv, scalar

CASE = Path(__file__).resolve().parent


def mapping(phi, asymmetry):
    """Eqs. (3b), (5a), (7), for 0 <= phi <= pi and positive a=Gamma_L/Gamma_R."""
    if not np.isfinite(asymmetry) or asymmetry <= 0 or not 0 <= phi <= np.pi:
        raise ValueError("positive finite asymmetry and phase in [0,pi] required")
    transparency = 4*asymmetry/(1+asymmetry)**2
    chi = max(0., 1-transparency*np.sin(phi/2)**2)
    phase = 2*np.arcsin(np.sqrt(transparency)*np.sin(phi/2))
    jacobian = 1. if asymmetry == 1 else np.sqrt(transparency)*np.cos(phi/2)/np.sqrt(chi)
    return float(chi), float(phase), float(jacobian)


def inverse_phase(phi_symmetric, asymmetry):
    argument = (asymmetry+1)/(2*np.sqrt(asymmetry))*np.sin(phi_symmetric/2)
    return float(2*np.arcsin(np.clip(argument, 0, 1))) if argument <= 1+1e-14 else None


def direct_model(bath, u, gamma, detuning, phi, asymmetry):
    gammas = gamma*np.array([asymmetry, 1.])/(1+asymmetry)
    contacts = [np.sqrt(g/(np.pi*bath.rho))*np.eye(2) for g in gammas]
    return make_model(Impurity.anderson(u, detuning=detuning), [bath, bath], contacts,
                      phases=[-phi/2, phi/2], phase_velocities=[-.5, .5])


def point(bath, u, gamma, x, phi, settings):
    h = reference_model(bath, u=u, gamma=gamma, detuning=x*u/2, phi=phi,
                        symmetry=True, compress=True)
    states = parity_states(h, settings)
    energies = [float(s.energies[0]) for s in states]
    return dict(tilde_epsilon=x, even_energy=energies[0], odd_energy=energies[1],
                signed_gap=energies[1]-energies[0]), states


def boundary(evaluate, tolerance):
    a, b = evaluate(0.)["signed_gap"], evaluate(1.)["signed_gap"]
    if a >= 0:
        return None
    if b <= 0:
        raise ValueError("transition not bracketed in tilde_epsilon=[0,1]")
    return float(brentq(lambda x: evaluate(x)["signed_gap"], 0., 1., xtol=tolerance))


def compare_mapping(bath, u, gamma, detuning, phi, asymmetry, settings):
    chi, effective, jacobian = mapping(phi, asymmetry)
    direct = parity_states(direct_model(bath, u, gamma, detuning, phi, asymmetry), settings)
    symmetric = parity_states(reference_model(bath, u=u, gamma=gamma, detuning=detuning,
                              phi=effective, symmetry=False, compress=False), settings)
    rows = []
    for p, (state, mapped) in enumerate(zip(direct, symmetric, strict=True)):
        current, predicted = 2*scalar(state, "phase_derivative"), 2*jacobian*scalar(mapped, "phase_derivative")
        rows.append(dict(asymmetry=asymmetry, phi_over_pi=phi/np.pi, chi=chi,
                         effective_phi_over_pi=effective/np.pi, parity=p,
                         energy_error=float(state.energies[0]-mapped.energies[0]),
                         charge_error=scalar(state, "impurity_charge")-scalar(mapped, "impurity_charge"),
                         current=current, predicted_current=predicted, current_error=current-predicted))
    return rows, direct, symmetric


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw, half_filling=False)
    spec, parameters = config["profiles"][profile], config["model"]
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    references = read_csv(CASE/"reference"/"figure1.csv")
    rows, comparisons, checks = [], [], []
    for settings in settings_list:
        label = settings.get("label", profile)
        bath = run.bath(settings, bandwidth=settings.get("bandwidth", parameters["bandwidth"]))
        for u_meV in spec["u_meV"]:
            u, gamma = u_meV/config["delta_meV"], config["gamma_meV"]/config["delta_meV"]
            for asymmetry in spec["asymmetries"]:
                for p in spec["phases_over_pi"]:
                    chi, effective, _ = mapping(p*np.pi, asymmetry)
                    cache = {}

                    def evaluate(x, *, cache=cache, bath=bath, u=u, gamma=gamma, effective=effective,
                                 settings=settings, label=label, u_meV=u_meV):
                        x = float(x)
                        if x not in cache:
                            row, states = point(bath, u, gamma, x, effective, settings)
                            cache[x] = row
                            for state in states:
                                run.keep(dict(setting=label, u_meV=u_meV, tilde_epsilon=x,
                                              effective_phi=effective), state)
                        return cache[x]

                    x = boundary(evaluate, config["root_tolerance"])
                    rows.append(dict(setting=label, u_meV=u_meV, asymmetry=asymmetry,
                                     phi_over_pi=p, chi=chi, tilde_epsilon=x,
                                     has_transition=int(x is not None)))
                    ref = sorted((r for r in references if float(r["u_meV"]) == u_meV),
                                 key=lambda r: float(r["phi_over_pi"]))
                    if x is not None and float(ref[0]["phi_over_pi"]) <= effective/np.pi <= float(ref[-1]["phi_over_pi"]):
                        xr = float(np.interp(effective/np.pi, [float(r["phi_over_pi"]) for r in ref],
                                             [float(r["tilde_epsilon"]) for r in ref]))
                        comparisons.append(dict(setting=label, u_meV=u_meV, asymmetry=asymmetry,
                                                phi_over_pi=p, tilde_epsilon=x, reference_tilde_epsilon=xr,
                                                difference=x-xr))
            print(f"Completed {label}, U={u_meV:g} meV", flush=True)
    if profile != "convergence":
        settings = config["mapping_settings"]
        bath = run.bath(settings)
        u, gamma = 3.2/config["delta_meV"], config["gamma_meV"]/config["delta_meV"]
        for a in (1., 4., 11., 1/11):
            for p in (.3, .8):
                data, direct, mapped = compare_mapping(bath, u, gamma, .2*u, p*np.pi, a, settings)
                checks.extend(data)
                for geometry, states in (("asymmetric", direct), ("symmetric", mapped)):
                    for state in states:
                        run.keep(dict(geometry=geometry, asymmetry=a, phi_over_pi=p), state)
    summary = dict(points=len(rows), maximum_residual=max(max(s["residuals"]) for s in run.states))
    if comparisons:
        summary["maximum_reference_gate_difference"] = max(abs(r["difference"]) for r in comparisons)
    for key in ("energy_error", "charge_error", "current_error"):
        if checks:
            summary[f"maximum_mapping_{key}"] = max(abs(r[key]) for r in checks)
    run.finish({"convergence" if profile == "convergence" else "boundary": rows,
                "comparison": comparisons, "mapping": checks}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
