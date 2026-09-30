"""Convention-matched comparison with externally supplied NRG observables.

This module does not run or validate an NRG discretization. The supplied
reference must document its own convergence. QP error estimates must include
the QP cutoff and quadrature, not only the Ritz residual.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from .io import write_json


def compare(reference, result, qp_uncertainty, state=0):
    if reference.get("format") != "bcs-qp-nrg-reference" or reference.get("format_version") != 1:
        raise ValueError("unsupported NRG reference format")
    if reference.get("status") != "converged-reference":
        raise ValueError("the NRG reference is a template or has not been marked converged")
    if result.get("format") != "bcs-qp-eigenstates" or result.get("format_version") != 1:
        raise ValueError("unsupported QP result format")
    if not np.isfinite(qp_uncertainty) or qp_uncertainty < 0:
        raise ValueError("QP uncertainty must be nonnegative and finite")
    conventions = reference["conventions"]
    if conventions.get("hybridization") != "total" or conventions.get("phase_bias") != "phi_R-phi_L":
        raise ValueError("NRG hybridization or phase convention is not matched")
    if conventions.get("impurity_convention") != "charge-square":
        raise ValueError("match the impurity gate and additive constant to the charge-square convention")
    settings = reference["nrg"]
    if not isinstance(settings.get("version"), str) or not settings["version"]:
        raise ValueError("NRG software version is required")
    if not isinstance(settings.get("kept_states"), int) or settings["kept_states"] <= 0:
        raise ValueError("NRG retained-state count is required")
    if not isinstance(settings.get("Lambda"), (int, float)) or not np.isfinite(settings["Lambda"]) or settings["Lambda"] <= 1:
        raise ValueError("NRG Lambda must exceed one")
    temperature = settings.get("temperature_over_delta")
    if not isinstance(temperature, (int, float)) or not np.isfinite(temperature) or temperature < 0:
        raise ValueError("NRG temperature must be specified and nonnegative")
    if not isinstance(settings.get("symmetry"), str) or not settings["symmetry"]:
        raise ValueError("the NRG symmetry sector must be specified")
    zs = np.asarray(settings.get("z_values", []), dtype=float)
    if zs.ndim != 1 or zs.size == 0 or not np.all(np.isfinite(zs)) or np.any((zs <= 0) | (zs > 1)):
        raise ValueError("NRG z values must lie in (0,1]")
    if not settings.get("convergence_notes"):
        raise ValueError("describe the convergence of the NRG reference")
    calculation = result["calculation"]
    params = reference["parameters"]
    if not all(name in calculation for name in ("u", "gamma", "phi", "rho_ws", "rho_wn", "detuning", "field")):
        raise ValueError("this NRG helper expects a built-in single-dot reference_model calculation")
    actual = {name: calculation[name] for name in ("u", "gamma", "phi", "rho_ws", "rho_wn", "detuning", "field")}
    actual.update(delta=calculation["bath"][0]["delta"], bandwidth=calculation["bath"][0]["bandwidth"],
                  geometry=calculation["geometry"], parity=calculation["sector"]["parity"],
                  twice_sz=calculation["sector"]["twice_sz"])
    for name, value in actual.items():
        other = params.get(name)
        if isinstance(value, str) or value is None:
            matched = value == other
        else:
            matched = isinstance(other, (int, float)) and np.isfinite(other) and np.isclose(value, other, rtol=1e-12, atol=1e-14)
        if not matched:
            raise ValueError(f"parameter mismatch for {name}: QP={value}, NRG={other}")
    if not 0 <= state < len(result["energies"]):
        raise ValueError("state index outside the QP result")
    comparisons = []
    for name, entry in reference["observables"].items():
        value, uncertainty = entry["value"], entry["uncertainty"]
        if value is None or uncertainty is None or not np.all(np.isfinite([value, uncertainty])) or uncertainty < 0:
            raise ValueError(f"missing finite value or uncertainty for {name}")
        if name == "energy":
            if conventions.get("energy_reference") != calculation.get("energy_reference"):
                raise ValueError("absolute-energy references do not match")
            calculated = result["energies"][state]
        else:
            calculated = result["observables"][name][state]
        difference = float(calculated-value)
        combined = float(uncertainty+qp_uncertainty)
        comparisons.append(dict(observable=name, qp_value=float(calculated), nrg_value=float(value),
                                difference=difference, combined_uncertainty=combined,
                                consistent=abs(difference) <= combined))
    if not comparisons:
        raise ValueError("no NRG observables were supplied")
    return dict(comparisons=comparisons, consistent=all(c["consistent"] for c in comparisons),
                parameters=actual, nrg=settings, qp_uncertainty=qp_uncertainty)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nrg_reference", type=Path)
    parser.add_argument("qp_result", type=Path)
    parser.add_argument("--qp-uncertainty", type=float, required=True)
    parser.add_argument("--state", type=int, default=0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    comparison = compare(json.loads(args.nrg_reference.read_text()), json.loads(args.qp_result.read_text()),
                         args.qp_uncertainty, args.state)
    for row in comparison["comparisons"]:
        print(f"{row['observable']}: QP={row['qp_value']:.12g}, NRG={row['nrg_value']:.12g}, "
              f"difference={row['difference']:.3e}, combined uncertainty={row['combined_uncertainty']:.3e}")
    if args.output:
        write_json(args.output, comparison)
    return 0 if comparison["consistent"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
