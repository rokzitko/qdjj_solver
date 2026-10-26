"""Tutorial 05: numerical convergence of a driven calculation.

    python -B examples/05_convergence/run.py                  # quick profile, about 7 s
    python -B examples/05_convergence/run.py --profile full   # the physics manual's junction, about 3 minutes

The phase-driven junction phi(t) = phi0 + A cos(Omega t) is checked in four independent
approximations: (1) the Floquet cutoff M, with the rigorous residual bound, for a strong
drive; (2) a cross-check of (1) with the one-period propagator, which has no Floquet cutoff;
(3) the size K of the static basis against the full QP Fock space, with the Fock-space
residual, for a weak drive; (4) the bath discretization, by following the dressed singlet
across the one-photon threshold of the continuum for several bath sizes. Results go to
runs/examples/05_convergence/<profile>/ (results.json, and convergence.png if matplotlib is
installed).
"""

import argparse
import dataclasses
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, floquet_propagator, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The reference junction of the physics manual, energies in units of Delta = 1.
JUNCTION = dict(u=1., gamma=.5, phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
BANDWIDTH = 10.
CUTOFF = 2               # QP cutoff of every calculation
SECTOR = Sector(0, 0)    # the singlets S and S*
OMEGA = .45              # drive frequency of sections 1-3
STRONG = 1.5             # phase-drive amplitude of sections 1 and 2
WEAK = .3                # phase-drive amplitude of sections 3 and 4
BOUND = 1e-6             # accuracy target for the residual bound in section 1
WINDOW = 2.4             # section 4 keeps the static states up to E_S + WINDOW
REFERENCE_PAIRS = 4      # bath of the threshold quoted in the summary; present in both profiles
PROFILES = {
    # bath level pairs per lead of sections 1-3; Floquet cutoffs and the reference cutoff;
    # smallest static basis of section 1; static-basis sizes of section 3; baths and frequencies of section 4
    "quick": dict(pairs=2, harmonics=range(1, 11), reference=14, basis_states=20,
                  basis_sizes=(4, 8, 16, 32, 64, 96), baths=(2, 4), frequencies=np.linspace(.95, 1.75, 33)),
    "full": dict(pairs=4, harmonics=range(1, 11), reference=16, basis_states=60,
                 basis_sizes=(4, 8, 16, 32, 64, 128, 256, 400), baths=(2, 4, 8),
                 frequencies=np.linspace(.95, 1.75, 321)),
}


def junction(pairs):
    return dict(JUNCTION, bath=cosh_grid(pairs=pairs, bandwidth=BANDWIDTH))


def phase_drive(model, amplitude, omega):
    """phi(t) = phi0 + A cos(Omega t). 24 samples per period resolve the harmonics for A <= 1.5;
    from_parameters checks this at the midpoints between the samples."""
    return PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("phi", amplitude),
                                               omega, samples=24)


def basis_size(energies, minimum):
    """Smallest K >= minimum that does not split a degenerate static multiplet."""
    gaps = np.flatnonzero(np.diff(energies) > 1e-7) + 1
    gaps = gaps[gaps >= minimum]
    return int(gaps[0]) if len(gaps) else len(energies)


def all_static_states(model):
    """Every static eigenstate of the sector: needed for any static basis size up to the Fock dimension."""
    h = reference_model(**model)
    dimension = solve(h, CUTOFF, SECTOR, SolverOptions(method="dense", eigenpairs=1)).basis.dimension
    return solve(h, CUTOFF, SECTOR, SolverOptions(method="dense", eigenpairs=dimension))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/05_convergence/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    model = junction(settings["pairs"])
    static = all_static_states(model)
    dimension = static.basis.dimension
    print(f"{settings['pairs']} bath level pairs per lead: Fock dimension {dimension}, "
          f"E_S = {static.energies[0]:.6f}, E_S* = {static.energies[1]:.6f}")

    # 1. Floquet cutoff M for a strong drive, in a static basis of K states. A much larger cutoff
    #    is the reference; the residual bound must lie above the actual error.
    strong = phase_drive(model, STRONG, OMEGA)
    kept = basis_size(static.energies, settings["basis_states"])
    options = dict(states=2, basis="eigen", basis_states=kept, dense_dimension=6000)
    reference = solve_floquet(strong, CUTOFF, SECTOR, dict(options, floquet_cutoff=settings["reference"]), static=static)
    cutoff = []
    for harmonics in settings["harmonics"]:
        result = solve_floquet(strong, CUTOFF, SECTOR, dict(options, floquet_cutoff=harmonics), static=static)
        cutoff.append(dict(harmonics=harmonics, error=abs(result.quasienergies - reference.quasienergies).tolist(),
                           bound=result.residuals.tolist(), boundary=result.boundary_weights.tolist()))
    error = np.array([row["error"] for row in cutoff])
    bound = np.array([row["bound"] for row in cutoff])
    m = np.array(settings["harmonics"])
    converged = int(m[np.flatnonzero(bound.max(axis=1) < BOUND)[0]])
    decay = -np.polyfit(m[1:7], np.log10(error[1:7].max(axis=1)), 1)[0]    # M = 2..7, above rounding
    print(f"K = {kept}: residual bound below {BOUND:g} from M = {converged}; error falls by {decay:.2f} decades "
          f"per harmonic; largest error/bound = {np.max(error/bound):.2e}")

    # 2. The one-period propagator: no Floquet cutoff, fourth-order Magnus steps checked by step doubling.
    propagator = floquet_propagator(strong, CUTOFF, SECTOR, dict(options, floquet_cutoff=8), static=static)
    deviation = abs(np.angle(np.exp(2j*np.pi*(propagator.quasienergies - reference.quasienergies)/OMEGA)))*OMEGA/(2*np.pi)
    step_error = float(np.max(propagator.metadata["step_errors"]))
    print(f"propagator - Sambe: {deviation.max():.2e} (estimated step error {step_error:.2e})")

    # 3. Static basis K against the full Fock space at a weaker drive. The Fock-space result, from
    #    sparse shift-invert iterations, is the reference; fock_residuals bound the static-basis errors.
    weak = phase_drive(model, WEAK, OMEGA)
    fock = solve_floquet(weak, CUTOFF, SECTOR, dict(floquet_cutoff=4, states=2, method="shift-invert"), static=static)
    basis = []
    for minimum in settings["basis_sizes"]:
        size = basis_size(static.energies, minimum)
        result = solve_floquet(weak, CUTOFF, SECTOR, dict(floquet_cutoff=4, states=2, basis="eigen", basis_states=size,
                                                          dense_dimension=6000), static=static)
        basis.append(dict(basis_states=size, error=abs(result.quasienergies - fock.quasienergies).tolist(),
                          fock_residuals=result.fock_residuals.tolist()))
    small = min(basis, key=lambda row: abs(np.log(row["basis_states"]/10)))       # K near 10
    large = max((row for row in basis if row["basis_states"] < dimension), key=lambda row: row["basis_states"])
    shift = float(fock.quasienergies[0] - static.energies[0])
    print(f"AC Stark shift of S (Fock space) {shift:.6f}; static-basis error {max(small['error']):.2e} at "
          f"K = {small['basis_states']}, {max(large['error']):.2e} at K = {large['basis_states']}")

    # 4. The dressed singlet across the one-photon threshold E_D + Delta - E_S for several bath sizes.
    continuum = {}
    for pairs in settings["baths"]:
        bath_model = junction(pairs)
        h = reference_model(**bath_model)
        states = solve(h, CUTOFF, SECTOR, SolverOptions(method="dense", eigenpairs=120))
        doublet = solve(h, CUTOFF, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1)).energies[0]
        size = basis_size(states.energies, int(np.sum(states.energies < states.energies[0] + WINDOW)))
        periodic = phase_drive(bath_model, WEAK, 1.)
        rows = []
        for omega in settings["frequencies"]:
            # The harmonics of phi0 + A cos(Omega t) do not depend on Omega: only the frequency changes.
            result = solve_floquet(dataclasses.replace(periodic, omega=float(omega)), CUTOFF, SECTOR,
                                   dict(floquet_cutoff=3, states=1, basis="eigen", basis_states=size,
                                        dense_dimension=6000, tracking_threshold=0.), static=states)
            rows.append(dict(omega=float(omega), shift=float(result.quasienergies[0] - states.energies[0]),
                             weight=float(result.static_weights[0, 0])))
        threshold = float(doublet + 1. - states.energies[0])
        weights = np.array([row["weight"] for row in rows])
        above = np.array([row["omega"] for row in rows]) > threshold
        below = solve_floquet(dataclasses.replace(periodic, omega=1.), CUTOFF, SECTOR,
                              dict(floquet_cutoff=3, states=1, basis="eigen", basis_states=size), static=states)
        continuum[str(pairs)] = dict(rows=rows, basis_states=size, dimension=states.basis.dimension,
                                     threshold=threshold, shift_below=float(below.quasienergies[0] - states.energies[0]),
                                     min_weight_above=float(weights[above].min()))
        print(f"{pairs} pairs (Fock dimension {states.basis.dimension}, K = {size}): threshold {threshold:.6f}, "
              f"shift at Omega = 1: {continuum[str(pairs)]['shift_below']:.6f}, "
              f"smallest weight above threshold {weights[above].min():.3f}")

    values = dict(
        fock_dimension=dimension, energy_S=static.energies[0], energy_Sstar=static.energies[1],
        basis_states=kept, strong_S=reference.quasienergies[0], strong_Sstar=reference.quasienergies[1],
        converged_harmonics=converged, decades_per_harmonic=decay, max_error_over_bound=np.max(error/bound),
        propagator_deviation=deviation.max(), propagator_step_error=step_error,
        fock_S=fock.quasienergies[0], fock_Sstar=fock.quasienergies[1], stark_shift=shift,
        basis_error_small=max(small["error"]), basis_error_large=max(large["error"]),
        max_error_over_fock_residual=max(e/r for row in basis for e, r in zip(row["error"], row["fock_residuals"],
                                                                              strict=True)),
        threshold=continuum[str(REFERENCE_PAIRS)]["threshold"],
    )
    for pairs, data in continuum.items():
        values[f"shift_below_{pairs}"] = data["shift_below"]
        values[f"min_weight_above_{pairs}"] = data["min_weight_above"]
    values = {key: float(value) for key, value in values.items()}
    for key, value in values.items():
        print(f"{key:30s} {value:.6g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values, cutoff=cutoff,
                  propagator=propagator.quasienergies.tolist(), basis=basis, continuum=continuum)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"convergence.png")
    print(f"Saved {output}")


def plot(record, path):
    """(a) Floquet cutoff, (b) static basis, (c) AC Stark shift and (d) static weight across the threshold."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "0.55"
    values = record["values"]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6), layout="constrained")
    axis = axes[0, 0]
    m = [row["harmonics"] for row in record["cutoff"]]
    for state, color, name in ((0, blue, "S"), (1, orange, "S^*")):
        axis.semilogy(m, [max(row["error"][state], 1e-16) for row in record["cutoff"]], "o-", ms=4, mfc="none",
                      color=color, label=f"error, ${name}$")
        axis.semilogy(m, [row["bound"][state] for row in record["cutoff"]], "--", lw=1, color=color,
                      label=f"bound, ${name}$")
    axis.axhline(max(values["propagator_deviation"], 1e-16), color=gray, lw=.8, ls=":")
    axis.text(m[-1], 2*max(values["propagator_deviation"], 1e-16), "propagator", ha="right", fontsize=7, color=gray)
    axis.set(xlabel="Floquet cutoff $M$", ylabel="quasienergy error")
    axis.legend(frameon=False, fontsize=7, ncols=2, loc="upper right")
    axis = axes[0, 1]
    k = [row["basis_states"] for row in record["basis"]]
    for state, color, name in ((0, blue, "S"), (1, orange, "S^*")):
        axis.loglog(k, [row["error"][state] for row in record["basis"]], "o-", ms=4, mfc="none", color=color,
                    label=f"error, ${name}$")
        axis.loglog(k, [row["fock_residuals"][state] for row in record["basis"]], "--", lw=1, color=color,
                    label=f"Fock residual, ${name}$")
    axis.axhline(abs(values["stark_shift"]), color=gray, lw=.8, ls=":")
    axis.text(k[-1], 1.3*abs(values["stark_shift"]), "AC Stark shift of $S$", ha="right", fontsize=7, color=gray)
    axis.set(xlabel="static basis size $K$", ylabel="error against the Fock space")
    axis.legend(frameon=False, fontsize=7, loc="lower left")
    for (pairs, data), color in zip(record["continuum"].items(), (blue, orange, aqua, gray, "k"), strict=False):
        omega = [row["omega"] for row in data["rows"]]
        axes[1, 0].plot(omega, [1e3*row["shift"] for row in data["rows"]], color=color, lw=1,
                        label=f"{pairs} pairs per lead (dimension {data['dimension']})")
        axes[1, 1].plot(omega, [row["weight"] for row in data["rows"]], color=color, lw=1)
    for axis in axes[1]:
        axis.axvline(values["threshold"], color=gray, lw=.8, ls=":")
        axis.set(xlabel=r"drive frequency $\Omega$", xlim=(min(omega), max(omega)))
    axes[1, 0].text(values["threshold"] + .01, 5, r"$E_D+\Delta-E_S$", fontsize=7, color=gray)
    axes[1, 0].set(ylabel=r"dressed-singlet shift $(\varepsilon_S-E_S)\times10^3$", ylim=(-6, 6))
    axes[1, 0].legend(frameon=False, fontsize=7, loc="lower left")
    axes[1, 1].set(ylabel="static weight on $S$", ylim=(.4, 1.02))
    for axis, label in zip(axes.flat, "abcd", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=100)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
