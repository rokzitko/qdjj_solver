"""Reproduce the junction-only gate/phase boundary of Bargerbos et al., Fig. S1(c)."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import reference_model
from applications._common import Run, arguments, parity_states, read_csv, scalar

CASE = Path(__file__).resolve().parent


def model(bath, u, gamma, gate_over_u, phi):
    return reference_model(bath, u=u, gamma=gamma, detuning=u*gate_over_u,
                           phi=phi, symmetry=False, compress=False)


def point(bath, parameters, gate_over_u, phi, settings):
    h = model(bath, parameters["u"], parameters["gamma"], gate_over_u, phi)
    states = parity_states(h, settings)
    even, odd = [float(s.energies[0]) for s in states]
    return dict(gate_over_u=gate_over_u, phi_over_pi=phi/np.pi,
                even_energy=even, odd_energy=odd, signed_gap=odd-even,
                ground_parity=int(odd < even), even_charge=scalar(states[0], "impurity_charge"),
                odd_charge=scalar(states[1], "impurity_charge")), states


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw, half_filling=False)
    spec = config["profiles"][profile]
    parameters = config["model"]
    reference = read_csv(CASE/"reference"/"boundary.csv")
    ref_phase = [float(r["phi_over_pi"]) for r in reference]
    ref_gate = [float(r["detuning_critical_over_u"]) for r in reference]
    settings_list = spec["settings"] if profile == "convergence" else [spec["settings"]]
    rows, gate_rows, checks = [], [], []
    for settings in settings_list:
        bath = run.bath(settings)
        label = settings.get("label", profile)
        phases = spec.get("phase_over_pi", [0., *ref_phase[::spec.get("reference_stride", 1)], 1.])
        cache = {}

        def evaluate(gate, p, *, cache=cache, bath=bath, settings=settings, label=label):
            key = float(gate), float(p)
            if key not in cache:
                row, states = point(bath, parameters, key[0], key[1]*np.pi, settings)
                cache[key] = row
                for state in states:
                    run.keep(dict(setting=label, gate_over_u=key[0], phi_over_pi=key[1]), state)
            return cache[key]

        for p in phases:
            root = brentq(lambda g, p=p: evaluate(g, p)["signed_gap"], *config["gate_bracket_over_u"],
                          xtol=config["gate_tolerance_over_u"])
            nrg = float(np.interp(p, ref_phase, ref_gate)) if ref_phase[0] <= p <= ref_phase[-1] else None
            rows.append(dict(setting=label, phi_over_pi=float(p), detuning_critical_over_u=root,
                             nrg_detuning_critical_over_u=nrg,
                             signed_gap_at_root=evaluate(root, p)["signed_gap"]))
        if profile != "convergence":
            for p in config["gate_scan_phases_over_pi"]:
                for gate in np.linspace(*config["gate_range_over_u"], spec["gate_points"]):
                    gate_rows.append(evaluate(gate, p))
            gate, p = config["charge_check_gate_over_u"], config["charge_check_phase_over_pi"]
            center, reflected = evaluate(gate, p), evaluate(-gate, p)
            step = config["charge_check_step"]
            minus, plus = evaluate(gate-step/parameters["u"], p), evaluate(gate+step/parameters["u"], p)
            for parity in ("even", "odd"):
                response = (plus[f"{parity}_energy"]-minus[f"{parity}_energy"])/(2*step)
                checks.append(dict(parity=parity, charge_minus_one=center[f"{parity}_charge"]-1,
                                   energy_derivative=response,
                                   derivative_error=response-(center[f"{parity}_charge"]-1),
                                   ph_energy_error=center[f"{parity}_energy"]-reflected[f"{parity}_energy"],
                                   ph_charge_sum_error=center[f"{parity}_charge"]+reflected[f"{parity}_charge"]-2))
        print(f"Completed {label}", flush=True)
    errors = [abs(r["detuning_critical_over_u"]-r["nrg_detuning_critical_over_u"]) for r in rows
              if r["nrg_detuning_critical_over_u"] is not None]
    summary = dict(boundary_points=len(rows), maximum_nrg_boundary_difference=max(errors),
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    if checks:
        summary["maximum_charge_derivative_error"] = max(abs(r["derivative_error"]) for r in checks)
        summary["maximum_particle_hole_energy_error"] = max(abs(r["ph_energy_error"]) for r in checks)
    run.finish({"convergence" if profile == "convergence" else "boundary": rows,
                "gate_scan": gate_rows, "charge_check": checks}, summary)
    return rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
