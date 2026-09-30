"""Bobok et al. (2025): GAL/ChE Josephson currents and the common-lead high-spin transition."""

from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from qdjj_solver import Sector
from qdjj_solver.qp_solver.basis import dimension
from qdjj_solver.qp_solver.solver import ProjectedOperator
from applications._common import (Run, arguments, eigenstate, parity_states, read_csv, scalar)
from .models import double_dot, gal_junction, junction

CASE = Path(__file__).resolve().parent


def current_point(bath, u, phi, settings, gamma_per_lead=1.):
    h = gal_junction(u, gamma_per_lead, phi) if bath is None else junction(bath, u, gamma_per_lead, phi)
    states = parity_states(h, settings)
    energies = [float(s.energies[0]) for s in states]
    parity = int(energies[1] < energies[0])
    return dict(even_energy=energies[0], odd_energy=energies[1], signed_gap=energies[1]-energies[0],
                ground_parity=parity, current=2*scalar(states[parity], "phase_derivative"),
                even_current=2*scalar(states[0], "phase_derivative") if phi < np.pi-1e-12 else None,
                odd_current=2*scalar(states[1], "phase_derivative")), states


def resolve_spin(state, tolerance=2e-9):
    """Resolve S^2 within degenerate Ritz subspaces before assigning multiplets."""
    energies = state.energies.copy()
    projected = ProjectedOperator(state.hamiltonian.observables["total_spin_squared"], state.basis)
    start = 0
    while start < len(energies):
        stop = start+1
        while stop < len(energies) and energies[stop]-energies[start] < tolerance:
            stop += 1
        if stop-start > 1:
            vectors = state.vectors[:, start:stop]
            actions = np.column_stack([projected.action(v) for v in vectors.T])
            matrix = vectors.conj().T@actions
            _, rotation = np.linalg.eigh((matrix+matrix.conj().T)/2)
            state.vectors[:, start:stop] = vectors@rotation
            h = ProjectedOperator(state.hamiltonian.operator, state.basis)
            for i in range(start, stop):
                v = state.vectors[:, i]
                action = h.action(v)
                state.energies[i] = np.vdot(v, action).real
                state.residuals[i] = np.linalg.norm(action-state.energies[i]*v)
                counts = state.basis.qp_counts()
                state.qp_weights[i] = np.bincount(counts, weights=abs(v)**2,
                                                 minlength=len(state.qp_weights[i]))
                for name, operator in state.hamiltonian.observables.items():
                    state.observables[name][i] = ProjectedOperator(operator, state.basis).expectation(v)
        start = stop
    casimir = np.real(state.observables["total_spin_squared"])
    spin = np.round(np.sqrt(np.maximum(0., 1+4*casimir))-1)/2
    error = float(max(abs(casimir-spin*(spin+1))))
    if error > 1e-6:
        raise RuntimeError("incomplete spin multiplet: increase roots or resolve further symmetries")
    if max(state.residuals) > 2e-8:
        raise RuntimeError("spin-resolved vectors failed the physical residual check")
    return spin, error


def shared_point(bath, u, settings, gamma_per_dot=.5):
    exchange_basis = settings.get("exchange_basis", False)
    h = double_dot(bath, u, gamma_per_dot, exchange_basis=exchange_basis)
    roots = settings.get("roots", 6)
    labels = (-1, 1) if exchange_basis else (None,)
    states = []
    for p in (0, 1):
        for eta in labels:
            sector = Sector(p, p, eta)
            cutoff = settings.get("cutoff", h.bath_modes)
            cutoff = h.bath_modes if cutoff is None else cutoff
            count = dimension(h.nimp, h.spins, cutoff, sector, h.eta_labels)
            states.append(eigenstate(h, sector, settings, roots=min(roots, count)))
    # Independent highest-weight representative checks the triplet assignment.
    triplet_states = [eigenstate(h, Sector(0, 2, eta), settings) for eta in labels]
    states.extend(triplet_states)
    records, spin_error = [], 0.
    for state in states:
        spins, error = resolve_spin(state)
        spin_error = max(spin_error, error)
        for i, spin in enumerate(spins):
            records.append(dict(energy=float(state.energies[i]), spin=float(spin),
                                twice_sz=state.basis.sector.twice_sz,
                                pairing=float(np.real(state.observables["pairing"][i])),
                                spin_correlation=float(np.real(state.observables["spin_correlation"][i]))))
    lowest = {s: min((r for r in records if r["spin"] == s), key=lambda r: r["energy"])
              for s in (0., .5, 1.)}
    ground = min(lowest.values(), key=lambda r: r["energy"])
    # Distinct multiplet energies, not the repeated magnetic projections.
    spectrum = []
    for record in sorted(records, key=lambda r: r["energy"]):
        if any(abs(record["energy"]-old["energy"]) < 2e-8 and record["spin"] == old["spin"] for old in spectrum):
            continue
        spectrum.append(record | {"excitation": record["energy"]-ground["energy"]})
    even_triplets = [r["energy"] for r in records if r["spin"] == 1 and r["twice_sz"] == 0]
    triplet_check = abs(min(even_triplets)-min(s.energies[0] for s in triplet_states)) if even_triplets else None
    row = dict(u=u, ground_spin=ground["spin"], singlet_energy=lowest[0.]["energy"],
               doublet_energy=lowest[.5]["energy"], triplet_energy=lowest[1.]["energy"],
               doublet_singlet_gap=lowest[.5]["energy"]-lowest[0.]["energy"],
               triplet_singlet_gap=lowest[1.]["energy"]-lowest[0.]["energy"],
               pairing=ground["pairing"], spin_correlation=ground["spin_correlation"],
               maximum_spin_quantization_error=spin_error, triplet_projection_error=triplet_check)
    return row, spectrum, states


def compare_current(rows, references):
    comparison = []
    for row in rows:
        for method in ("NRG", row["setting"]):
            reference = [r for r in references if r["method"] == method and float(r["u"]) == row["u"]]
            if method == "NRG":
                reference = [r for r in reference if abs(float(r["phi_over_pi"])-row["phi_over_pi"]) < 1e-7]
            else:
                reference = sorted(reference, key=lambda r: float(r["phi_over_pi"]))
            if not reference:
                continue
            ref = float(np.interp(row["phi_over_pi"], [float(r["phi_over_pi"]) for r in reference],
                                  [float(r["current"]) for r in reference]))
            comparison.append(dict(setting=row["setting"], u=row["u"], phi_over_pi=row["phi_over_pi"],
                                   reference=method, current=row["current"], reference_current=ref,
                                   difference=row["current"]-ref,
                                   same_sign=int(row["current"]*ref >= 0 or abs(ref) < 1e-5)))
    return comparison


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    gamma_lead = config["model"].get("junction_gamma_per_lead", 1.)
    gamma_dot = config["model"].get("shared_gamma_per_dot", .5)
    reference8 = read_csv(CASE/"reference"/"figure8.csv")
    reference17 = read_csv(CASE/"reference"/"figure17.csv")
    currents, checks, shared, spectra, crossings, comparisons17 = [], [], [], [], [], []
    for settings in spec["current_settings"]:
        label = settings["label"]
        bath = None if label == "GAL" else run.bath(settings)
        for u in spec["current_u"]:
            cache = {}

            def evaluate(p, *, cache=cache, bath=bath, u=u, settings=settings, label=label):
                p = float(p)
                if p not in cache:
                    row, states = current_point(bath, u, p*np.pi, settings, gamma_lead)
                    cache[p] = row
                    for state in states:
                        run.keep(dict(figure=8, setting=label, u=u, phi_over_pi=p), state)
                return cache[p]

            phases = spec["phases_over_pi"]
            if spec.get("include_reference"):
                phases = sorted(set(phases) | {float(r["phi_over_pi"]) for r in reference8
                                               if r["method"] == "NRG" and float(r["u"]) == u})
            for phase in phases:
                currents.append(dict(setting=label, u=u, phi_over_pi=phase, **evaluate(phase)))
            if spec.get("locate_crossings"):
                root = brentq(lambda p: evaluate(p)["signed_gap"], 0., 1., xtol=2e-7)
                crossings.append(dict(figure=8, setting=label, u=u, transition="singlet-doublet",
                                      critical_parameter=root, parameter="phi_over_pi"))
            if profile != "convergence":
                p, step = .3, 1e-4
                center, minus, plus = evaluate(p), evaluate(p-step/np.pi), evaluate(p+step/np.pi)
                for branch in ("even", "odd"):
                    fd = (plus[f"{branch}_energy"]-minus[f"{branch}_energy"])/step
                    checks.append(dict(setting=label, u=u, branch=branch,
                                       derivative_error=fd-center[f"{branch}_current"]))
        print(f"Completed Fig. 8 {label}", flush=True)

    for settings in spec["shared_settings"]:
        label = settings["label"]
        bath = run.bath(settings)
        cache = {}

        def evaluate(u, *, cache=cache, bath=bath, settings=settings, label=label):
            u = float(u)
            if u not in cache:
                row, spectrum, states = shared_point(bath, u, settings, gamma_dot)
                cache[u] = (row, spectrum)
                for state in states:
                    run.keep(dict(figure=17, setting=label, u=u), state)
            return cache[u]

        interactions = settings.get("u", spec["shared_u"])
        if spec.get("include_reference"):
            interactions = sorted(set(interactions) | {round(float(r["u"]), 4) for r in reference17
                                                      if r["method"] == "NRG"})
        for u in interactions:
            row, spectrum = evaluate(u)
            shared.append(dict(setting=label, **row))
            spectra.extend(dict(setting=label, u=u, **r) for r in spectrum if r["excitation"] <= 1.02)
            for ref in reference17:
                if ref["method"] != "NRG" or abs(float(ref["u"])-u) > 3e-5:
                    continue
                observable = ref["observable"]
                expected = float(ref["value"])
                if observable == "excitation":
                    value = min((r["excitation"] for r in spectrum), key=lambda e: abs(e-expected))
                else:
                    value = row[observable]
                comparisons17.append(dict(setting=label, u=u, observable=observable, value=value,
                                          reference_value=expected, difference=value-expected,
                                          ground_spin=row["ground_spin"]))
        if spec.get("locate_crossings"):
            for name, key, bracket in (("doublet-singlet", "doublet_singlet_gap", (.1, 3.)),
                                       ("singlet-triplet", "triplet_singlet_gap", (3., 10.))):
                left, right = bracket
                if evaluate(left)[0][key]*evaluate(right)[0][key] < 0:
                    root = brentq(lambda u, key=key: evaluate(u)[0][key], left, right, xtol=2e-5)
                    crossings.append(dict(figure=17, setting=label, u=None, transition=name,
                                          critical_parameter=root, parameter="u"))
        print(f"Completed Fig. 17 {label}", flush=True)
    comparisons8 = compare_current(currents, reference8)
    summary = dict(current_points=len(currents), shared_points=len(shared), crossings=crossings,
                   maximum_residual=max(float(max(s["residuals"])) for s in run.states),
                   maximum_spin_quantization_error=max((r["maximum_spin_quantization_error"] for r in shared), default=0.),
                   maximum_triplet_projection_error=max((r["triplet_projection_error"] for r in shared
                                                         if r["triplet_projection_error"] is not None), default=0.))
    if checks:
        summary["maximum_current_derivative_error"] = max(abs(r["derivative_error"]) for r in checks)
    run.finish(dict(current=currents, shared=shared, spectrum=spectra, crossings=crossings,
                    current_checks=checks, current_comparison=comparisons8, shared_comparison=comparisons17), summary)
    return summary


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
