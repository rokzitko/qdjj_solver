"""Reproduce the weak/moderate-coupling Knight shift and its phase dependence."""

from pathlib import Path

import numpy as np
from scipy.integrate import quad

from qdjj_solver import Sector, reference_model
from applications._common import (Run, arguments, eigenstate, read_csv, scalar)

CASE = Path(__file__).resolve().parent


def model(bath, u, gamma, phi, field=0., compress=True):
    return reference_model(bath, u=u, gamma=gamma, phi=phi, field=field,
                           detuning=0., symmetry=True, compress=compress)


def leading_coefficient(u, bandwidth):
    """Eq. (26): kappa/Gamma at Gamma -> 0 for the continuous finite band."""
    return 2/np.pi*quad(lambda x: 1/(u/2+np.hypot(1., x))**2,
                        0, bandwidth, epsabs=1e-12, epsrel=1e-12)[0]


def knight(bath, u, gamma, phi, settings, field=0., twice_sz=1):
    h = model(bath, u, gamma, phi, field, settings.get("compress", True))
    result = eigenstate(h, Sector(1, twice_sz), settings)
    return 1-2*scalar(result, "impurity_spin_z"), result


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    specification = config["profiles"][profile]
    u, bandwidth = config["model"]["u"], config["model"]["bandwidth"]
    slope = leading_coefficient(u, bandwidth)
    if profile == "convergence":
        rows = []
        for settings in specification["settings"]:
            bath = run.bath(settings)
            for g in specification["gamma_over_u"]:
                for p in specification["phase_over_pi"]:
                    k, state = knight(bath, u, g*u, p*np.pi, settings)
                    run.keep(dict(settings=settings["label"], gamma_over_u=g, phi_over_pi=p), state)
                    rows.append(dict(setting=settings["label"], gamma_over_u=g,
                                     phi_over_pi=p, kappa=k, residual=float(state.residuals[0])))
            print(f"Completed {settings['label']}", flush=True)
        run.finish({"convergence": rows}, dict(points=len(rows), leading_kappa_over_gamma=slope))
        return rows

    settings = specification["settings"]
    bath = run.bath(settings)
    ref9 = read_csv(CASE/"reference"/"figure9.csv")
    ref10 = read_csv(CASE/"reference"/"figure10.csv")
    if profile == "paper":
        available = [r for r in ref9 if float(r["gamma_over_u"]) <= specification["gamma_over_u_max"]]
        selected = available[::specification["reference_stride"]]
        if selected[-1] != available[-1]:
            selected.append(available[-1])
        gammas = [float(r["gamma_over_u"]) for r in selected]
        phases = [float(r["phi_over_pi"]) for r in ref10[::specification["phase_stride"]]]
        if phases[-1] != float(ref10[-1]["phi_over_pi"]):
            phases.append(float(ref10[-1]["phi_over_pi"]))
    else:
        gammas = specification["gamma_over_u"]
        phases = specification["phase_over_pi"]

    coupling, phase_rows, field_rows = [], [], []
    for g in gammas:
        values = []
        for p in (0., 1.):
            k, state = knight(bath, u, g*u, p*np.pi, settings)
            run.keep(dict(scan="coupling", gamma_over_u=g, phi_over_pi=p), state)
            values.append(k)
        coupling.append(dict(gamma_over_u=g, kappa_phi0=values[0], kappa_phipi=values[1],
                             leading_kappa=g*u*slope, modulation=values[1]-values[0]))
    g = config["phase_scan_gamma_over_u"]
    for p in phases:
        k, state = knight(bath, u, g*u, p*np.pi, settings)
        run.keep(dict(scan="phase", gamma_over_u=g, phi_over_pi=p), state)
        phase_rows.append(dict(phi_over_pi=p, gamma_over_u=g, kappa=k))
    p = config["finite_field_phi_over_pi"]
    zero, state = knight(bath, u, g*u, p*np.pi, settings)
    run.keep(dict(scan="field", field=0., phi_over_pi=p), state)
    for field in config["finite_fields"]:
        states = [knight(bath, u, g*u, p*np.pi, settings, field, sz)[1] for sz in (1, -1)]
        for sz, state in zip((1, -1), states, strict=True):
            run.keep(dict(scan="field", field=field, twice_sz=sz, phi_over_pi=p), state)
        finite = 1-(states[0].energies[0]-states[1].energies[0])/field
        field_rows.append(dict(field=field, kappa_energy=float(finite), kappa_spin=zero,
                               difference=float(finite-zero)))

    ref_x = np.array([float(r["gamma_over_u"]) for r in ref9])
    ref_p = np.array([float(r["phi_over_pi"]) for r in ref10])
    errors = [abs(r[key]-np.interp(r["gamma_over_u"], ref_x,
              [float(t[key]) for t in ref9])) for r in coupling for key in ("kappa_phi0", "kappa_phipi")]
    phase_errors = [abs(r["kappa"]-np.interp(r["phi_over_pi"], ref_p,
                    [float(t["kappa"]) for t in ref10])) for r in phase_rows]
    summary = dict(leading_kappa_over_gamma=slope,
                   maximum_coupling_scan_nrg_absolute_error=max(errors),
                   maximum_phase_scan_nrg_absolute_error=max(phase_errors),
                   smallest_field_identity_error=abs(field_rows[-1]["difference"]),
                   maximum_residual=max(float(np.max(s["residuals"])) for s in run.states))
    run.finish({"coupling": coupling, "phase": phase_rows, "field_check": field_rows}, summary)
    return coupling, phase_rows, field_rows


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
