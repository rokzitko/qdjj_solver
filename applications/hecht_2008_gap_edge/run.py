"""Near-gap spectra: resolved quadratic limit and interacting finite-bath diagnostics."""

from pathlib import Path

import numpy as np

from applications._common import Run, arguments, read_csv
from .spectral import (analytic_continuum, broaden_continuum, fock_measure,
                       logarithmic_bath, quadratic_green)

CASE = Path(__file__).resolve().parent


def calculate(config, profile, output, raw=False):
    run = Run(CASE, config, profile, output, raw)
    spec = config["profiles"][profile]
    references = read_csv(CASE/"reference"/"spectra.csv")
    curves, convergence, comparison, poles_table, krylov = [], [], [], [], []
    for settings in spec["quadratic_settings"]:
        for delta, eps_ratio in spec["quadratic_cases"]:
            gamma, epsilon = .008/delta, eps_ratio*.008/delta
            bath = logarithmic_bath(settings["cells"], 1/delta)
            run.bath({"bath_record": bath.record()}, bandwidth=1/delta)
            # Relative edge distances resolve even the epsilon/Gamma=8 peak.
            offsets = np.geomspace(1e-10, 1e2, spec["frequency_points"])
            z = 1+offsets+1j*settings["eta_ratio"]*offsets
            discrete = -np.imag(quadratic_green(z, gamma, epsilon, bath=bath))*gamma
            continuum = -np.imag(quadratic_green(z, gamma, epsilon, bandwidth=1/delta))*gamma
            analytic = analytic_continuum(offsets, gamma, epsilon)*np.pi*gamma
            scale = max(analytic)
            convergence.append(dict(setting=settings["label"], kind="quadratic", delta_over_D=delta,
                                    epsilon_over_gamma=eps_ratio, cells=settings["cells"], eta_ratio=settings["eta_ratio"],
                                    max_quadrature_error_over_peak=float(max(abs(discrete-continuum))/scale),
                                    max_broadening_band_error_over_peak=float(max(abs(continuum-analytic))/scale),
                                    peak_offset_over_D=float(offsets[np.argmax(discrete)]*delta),
                                    peak_pi_gamma_A=float(max(discrete))))
            for offset, value, exact, limit in zip(offsets, discrete, continuum, analytic, strict=True):
                curves.append(dict(setting=settings["label"], delta_over_D=delta, epsilon_over_gamma=eps_ratio,
                                   omega_prime_over_D=offset*delta, pi_gamma_A=float(value),
                                   finite_band_continuum=float(exact), wide_band_analytic=float(limit)))
            if profile != "convergence":
                for ref in references:
                    if (int(ref["figure"]) != 6 or float(ref["delta_over_D"]) != delta or
                            not np.isclose(float(ref["epsilon_over_D"])/.008, eps_ratio)):
                        continue
                    offset = float(ref["omega_prime_over_D"])/delta
                    # Near-edge comparison, excluding the finite-band high-frequency range.
                    if not 1e-10 <= offset <= .1:
                        continue
                    exact = float(analytic_continuum(offset, gamma, epsilon)*np.pi*gamma)
                    value = float(-quadratic_green(1+offset+1j*settings["eta_ratio"]*offset,
                                                   gamma, epsilon, bath=bath).imag*gamma)
                    published = float(ref["pi_gamma_A"])
                    comparison.append(dict(kind="quadratic", method=ref["method"], delta_over_D=delta,
                                           epsilon_over_gamma=eps_ratio, omega_prime_over_D=offset*delta,
                                           pi_gamma_A=value, reference_pi_gamma_A=published,
                                           relative_difference=value/published-1,
                                           wide_band_pi_gamma_A=exact,
                                           wide_band_relative_difference=exact/published-1))
            print(f"Quadratic {settings['label']}: Delta/D={delta:g}, epsilon/Gamma={eps_ratio:g}", flush=True)

    # Connect the many-body transition calculation with the exact same finite
    # quadratic bath resolvent, rather than only comparing two continuum formulas.
    check_bath = run.bath(config["fock_check_settings"], bandwidth=100.)
    check_poles, states, info = fock_measure(check_bath, 0., 2., .3, config["fock_check_settings"], 256)
    for state in states:
        run.keep(dict(kind="quadratic-Fock-check"), state)
    z = np.array([.2, 1.2, 3.])+0.07j
    lehmann = sum(p["weight"]/(z-p["omega"]) for p in check_poles)
    check_error = float(max(abs(lehmann-quadratic_green(z, 2., .3, bath=check_bath))))

    for settings in spec["interacting_settings"]:
        for delta in spec["interacting_gaps_over_D"]:
            gamma, u = .049/delta, .6/delta
            bath = run.bath(settings | {"frequency_cutoff": 1/delta}, bandwidth=1/delta)
            poles, states, diagnostics = fock_measure(bath, u, gamma, 0., settings, settings["steps"])
            label = settings["label"]
            for state in states:
                run.keep(dict(kind="interacting", setting=label, delta_over_D=delta), state)
            for p in poles:
                poles_table.append(dict(setting=label, delta_over_D=delta, **p))
            krylov.extend(dict(setting=label, delta_over_D=delta, **d) for d in diagnostics)
            positive = [p for p in poles if p["omega"] > 1 and p["weight"] > 1e-10]
            first_edge = min(p["omega"]-1 for p in positive)
            convergence.append(dict(setting=label, kind="interacting", delta_over_D=delta,
                                    epsilon_over_gamma=-.3/.049, cells=0, eta_ratio=0.,
                                    max_quadrature_error_over_peak=None, max_broadening_band_error_over_peak=None,
                                    peak_offset_over_D=None, peak_pi_gamma_A=None))
            selected = [r for r in references if int(r["figure"]) == 7 and float(r["delta_over_D"]) == delta]
            offsets = np.array([float(r["omega_prime_over_D"])/delta for r in selected])
            for width in spec["log_widths"]:
                values = broaden_continuum(offsets, poles, width)*np.pi*gamma
                for ref, value in zip(selected, values, strict=True):
                    comparison.append(dict(kind=f"interacting-{label}-b{width}", method="NRG", delta_over_D=delta,
                                           epsilon_over_gamma=-.3/.049,
                                           omega_prime_over_D=float(ref["omega_prime_over_D"]),
                                           pi_gamma_A=float(value), reference_pi_gamma_A=float(ref["pi_gamma_A"]),
                                           relative_difference=float(value)/float(ref["pi_gamma_A"])-1,
                                           wide_band_pi_gamma_A=None, wide_band_relative_difference=None))
            krylov.append(dict(setting=label, delta_over_D=delta, spin="sum", addition="both", dimension=0,
                               steps=settings["steps"], weight=sum(p["weight"] for p in poles), terminal_beta=0.,
                               first_continuum_offset_over_delta=first_edge,
                               signed_gap=float(states[1].energies[0]-states[0].energies[0])))
            print(f"Interacting {label}: Delta/D={delta:g}, first continuum offset/Delta={first_edge:.5g}", flush=True)
    # Uniform schemas keep the compact CSVs usable without custom readers.
    for row in krylov:
        row.setdefault("first_continuum_offset_over_delta", None)
        row.setdefault("signed_gap", None)
    summary = dict(quadratic_fock_resolvent_error=check_error,
                   quadratic_fock_sum_rule_error=abs(sum(p["weight"] for p in check_poles)-1),
                   maximum_residual=max(max(s["residuals"]) for s in run.states),
                   maximum_interacting_sum_rule_error=max((abs(r["weight"]-1) for r in krylov if r["spin"] == "sum"), default=0.),
                   interacting_continuum_converged=False)
    for method in ("analytic", "NRG"):
        data = [r for r in comparison if r["kind"] == "quadratic" and r["method"] == method]
        if data:
            summary[f"maximum_quadratic_{method}_relative_difference"] = max(abs(r["relative_difference"]) for r in data)
            summary[f"maximum_wide_band_{method}_relative_difference"] = max(abs(r["wide_band_relative_difference"]) for r in data)
    run.finish({"spectra": curves, "convergence": convergence, "comparison": comparison,
                "poles": poles_table, "krylov": krylov}, summary)
    return summary


def main():
    args, config = arguments(CASE, __doc__)
    calculate(config, args.profile, args.output, args.raw)


if __name__ == "__main__":
    main()
