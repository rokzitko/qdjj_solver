"""Three-terminal phase diagram and direct checks of the hidden-symmetry mapping."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import Impurity, Sector, make_model, reference_model
from applications._common import Run, arguments, eigenstate, parity_states, scalar

CASE = Path(__file__).resolve().parent


def geometry(relative, phases):
    relative, phases = np.asarray(relative, float), np.asarray(phases, float)
    if (relative.shape != phases.shape or relative.ndim != 1 or
            np.any(relative < 0) or not np.isclose(relative.sum(), 1.) or
            not np.all(np.isfinite(phases))):
        raise ValueError("normalized nonnegative couplings and matching finite phases required")
    terms = relative*np.exp(1j*phases)
    z = terms.sum()
    chi = float(np.clip(abs(z), 0, 1))
    gradient = -np.imag(terms*z.conjugate())/chi if chi > 1e-14 else None
    return chi, gradient


def direct_model(bath, u, gamma, relative, phases):
    """Independent reservoir modes; each lead-current operator is a partial derivative."""
    contacts = [np.sqrt(gamma*r/(np.pi*bath.rho))*np.eye(2) for r in relative]
    models = [make_model(Impurity.anderson(u), [bath]*len(relative), contacts,
                         phases=phases, phase_velocities=np.eye(len(relative))[j])
              for j in range(len(relative))]
    h = models[0]
    for j, other in enumerate(models):
        h.observables[f"lead_current_{j}"] = other.observables["phase_derivative"]
    return h


def effective_model(bath, u, gamma, chi, compress=True):
    if not 0 <= chi <= 1:
        raise ValueError("chi must lie in [0,1]")
    return reference_model(bath, u=u, gamma=gamma, phi=2*np.arccos(chi),
                           symmetry=True, compress=compress)


def compare_mapping(bath, parameters, relative, phases, settings):
    chi, gradient = geometry(relative, phases)
    if gradient is None or not 0 < chi < 1:
        raise ValueError("use nonsingular geometry for the current chain-rule test")
    direct = parity_states(direct_model(bath, parameters["u"], parameters["gamma"], relative, phases), settings)
    effective = parity_states(effective_model(bath, parameters["u"], parameters["gamma"], chi), settings)
    rows = []
    for p, (state, mapped) in enumerate(zip(direct, effective, strict=True)):
        currents = np.array([scalar(state, f"lead_current_{j}") for j in range(len(relative))])
        universal = -2*scalar(mapped, "phase_derivative")/np.sqrt(1-chi**2)
        prediction = universal*gradient
        rows.append(dict(chi=chi, parity=p, energy_direct=float(state.energies[0]),
                         energy_effective=float(mapped.energies[0]),
                         energy_error=float(state.energies[0]-mapped.energies[0]),
                         current_0=float(currents[0]), current_1=float(currents[1]), current_2=float(currents[2]),
                         current_mapping_error=float(np.max(abs(currents-prediction))),
                         current_sum=float(currents.sum()),
                         charge_error=scalar(state, "impurity_charge")-scalar(mapped, "impurity_charge")))
    return rows, direct, effective


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    parameters = config["model"]
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    curves, convergence = [], []
    for settings in settings_list:
        bath = run.bath(settings)
        label = settings.get("label", profile)
        cache = {}

        def evaluate(chi, *, cache=cache, bath=bath, settings=settings, label=label):
            chi = float(chi)
            if chi not in cache:
                h = effective_model(bath, parameters["u"], parameters["gamma"], chi)
                states = parity_states(h, settings)
                energies = [float(s.energies[0]) for s in states]
                cache[chi] = dict(chi=chi, even_energy=energies[0], odd_energy=energies[1],
                                  signed_gap=energies[1]-energies[0],
                                  # The two even levels are degenerate at chi=0;
                                  # an arbitrary eigenvector has no unique branch derivative.
                                  even_current=scalar(states[0], "phase_derivative") if chi > 1e-14 else None,
                                  odd_current=scalar(states[1], "phase_derivative"))
                for state in states:
                    run.keep(dict(setting=label, chi=chi, geometry="effective-two-lead"), state)
            return cache[chi]

        critical = brentq(lambda x: evaluate(x)["signed_gap"], .5, .95, xtol=2e-7)
        convergence.append(dict(setting=label, chi_critical=critical,
                                phase_critical_over_pi=float(2*np.arccos(critical)/np.pi),
                                gap_at_published_chi=evaluate(config["reference"]["chi_critical"])["signed_gap"]))
        if profile != "convergence":
            for chi in np.linspace(0, 1, spec["chi_points"]):
                curves.append(evaluate(chi))
        print(f"Completed {label}: chi_c={critical:.9f}", flush=True)

    if profile == "convergence":
        run.finish({"convergence": convergence}, dict(points=len(convergence),
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states)))
        return convergence

    relative = config["relative_couplings"]
    maps = []
    for p2 in np.linspace(-1, 1, spec["map_points"]):
        for p3 in np.linspace(-1, 1, spec["map_points"]):
            chi, _ = geometry(relative, [0, p2*np.pi, p3*np.pi])
            maps.append(dict(phi2_over_pi=float(p2), phi3_over_pi=float(p3), chi=chi,
                             ground_parity=int(chi < critical)))

    direct_settings = config["direct_settings"]
    direct_bath = run.bath(direct_settings)
    comparisons = []
    for p2, p3 in config["direct_phase_pairs_over_pi"]:
        phases = [0, p2*np.pi, p3*np.pi]
        rows, direct, effective = compare_mapping(direct_bath, parameters, relative, phases, direct_settings)
        for row in rows:
            comparisons.append(dict(phi2_over_pi=p2, phi3_over_pi=p3, **row))
        for name, states in (("direct-three-lead", direct), ("effective-two-lead", effective)):
            for state in states:
                run.keep(dict(geometry=name, phi2_over_pi=p2, phi3_over_pi=p3), state)

    # Equal second/third couplings: exact phasor closure without using the paper's misprinted Eq. (5).
    theta = float(np.arccos(-relative[0]/(2*relative[1])))
    phases = [0., theta, -theta]
    # Dense diagonalization on a small independent bath guarantees a complete
    # degenerate multiplet; single-vector Lanczos can miss an exactly repeated root.
    closure_settings = config["closure_settings"]
    closure_bath = run.bath(closure_settings)
    h = direct_model(closure_bath, parameters["u"], parameters["gamma"], relative, phases)
    even = eigenstate(h, Sector(0, 0), closure_settings, roots=2)
    odd = eigenstate(h, Sector(1, 1), closure_settings)
    for state in (even, odd):
        run.keep(dict(geometry="three-lead-closure", phases=phases), state)
    symmetry = dict(phi2_over_pi=theta/np.pi, phi3_over_pi=-theta/np.pi,
                    chi=geometry(relative, phases)[0],
                    even_level_splitting=float(even.energies[1]-even.energies[0]),
                    even_excitation=float(even.energies[0]-odd.energies[0]))
    summary = dict(critical=convergence[0], high_symmetry=symmetry,
                   maximum_energy_mapping_error=max(abs(r["energy_error"]) for r in comparisons),
                   maximum_current_mapping_error=max(r["current_mapping_error"] for r in comparisons),
                   maximum_current_conservation_error=max(abs(r["current_sum"]) for r in comparisons),
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    run.finish({"spectrum": curves, "phase_map": maps, "mapping": comparisons,
                "critical": convergence}, summary)
    return summary


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
