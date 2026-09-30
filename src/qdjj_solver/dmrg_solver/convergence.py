"""Explicit finite-bond ladders and state/subspace overlap diagnostics."""

from dataclasses import replace

import numpy as np
from scipy.optimize import linear_sum_assignment

from ..common.problem import DEFAULT_SECTOR
from .solver import DEFAULT_OPTIONS, solve


def compare_states(previous, current):
    """Compare roots in identical coordinates, reporting both root and subspace overlap.

    The assignment maximizes individual overlaps; singular values are invariant
    under rotations within a degenerate retained subspace. No phase-dependent
    amplitude is used as an individual-state convergence criterion.
    """
    overlap = previous.overlaps(current)
    rows, columns = linear_sum_assignment(-abs(overlap)**2)
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    return dict(assignment=list(zip(rows.tolist(), columns.tolist(), strict=True)),
                overlap_squared=abs(overlap)**2, subspace_singular_values=singular_values,
                energy_changes=current.energies[columns]-previous.energies[rows],
                observable_changes={name: current.observables[name][columns]-values[rows]
                                    for name, values in previous.observables.items()
                                    if name in current.observables})


def converge_bond_dimension(hamiltonian, chi_values, *, sector=DEFAULT_SECTOR, options=DEFAULT_OPTIONS,
                            energy_tolerance=1e-8, observable_tolerance=1e-7,
                            stable_steps=2, callback=None):
    """Run the supplied increasing chi ladder, retaining explicit evidence.

    Requires finite-problem residual/sweep convergence and stable successive
    energies/observables. The returned assessment is empirical; it neither
    establishes bath convergence nor proves that no lower root was missed.
    ``callback(result, entry)`` can archive each completed step immediately.
    """
    chis = tuple(chi_values)
    if (not chis or any(not isinstance(c, int) or isinstance(c, bool) or c < 1 for c in chis)
            or any(b <= a for a, b in zip(chis, chis[1:], strict=False))):
        raise ValueError("chi_values must be strictly increasing positive integers")
    if (not isinstance(stable_steps, int) or isinstance(stable_steps, bool) or stable_steps < 1
            or not np.isfinite(energy_tolerance) or energy_tolerance <= 0
            or not np.isfinite(observable_tolerance) or observable_tolerance <= 0):
        raise ValueError("invalid convergence targets")
    previous, history, stable = None, [], 0
    for chi in chis:
        result = solve(hamiltonian, sector=sector,
                       options=replace(options, chi_max=chi, require_convergence=False), initial=previous)
        entry = dict(chi_max=chi, energies=result.energies, residuals=result.residuals,
                     finite_problem_converged=result.metadata["finite_problem_converged"])
        if previous is not None:
            comparison = compare_states(previous, result)
            entry["comparison"] = comparison
            ok = (entry["finite_problem_converged"]
                  and np.max(abs(comparison["energy_changes"])) <= energy_tolerance
                  and all(np.max(abs(delta)) <= observable_tolerance
                          for delta in comparison["observable_changes"].values()))
            stable = stable+1 if ok else 0
        entry["stable_steps"] = stable
        history.append(entry)
        if callback is not None:
            callback(result, entry)
        previous = result
    assessment = dict(format="qdjj-bond-convergence", format_version=1, history=history,
                      empirical_convergence=stable >= stable_steps,
                      required_stable_steps=stable_steps, energy_tolerance=energy_tolerance,
                      observable_tolerance=observable_tolerance, bath_convergence_checked=False)
    return previous, assessment
