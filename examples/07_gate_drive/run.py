"""Tutorial 07: the gate-driven junction.

    python -B examples/07_gate_drive/run.py                  # quick profile, about 15 s
    python -B examples/07_gate_drive/run.py --profile full   # the complete calculation, under 2 minutes

A microwave gate voltage modulates the dot level, epsilon(t) = epsilon0 + A cos(Omega t), and
couples the singlets S and S* through the dot charge n_d. The script computes (1) the S-S*
resonance of the gate drive next to that of a phase drive, (2) the resonant splittings of both
drives versus the operating point epsilon0, where the eta selection rule makes them
complementary, (3) the AC Stark shift of the dressed singlet below the resonance, compared with
second-order perturbation theory, and (4) the one- and two-photon splittings versus the
amplitude. Results go to runs/examples/07_gate_drive/<profile>/ (results.json, and
gate_drive.png if matplotlib is installed).
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
from qdjj_solver.qp_solver import SolverOptions, rwa_model, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The junction of the physics manual: energies in units of the gap, Delta = 1, phases in radians.
JUNCTION = dict(u=1., gamma=.5, phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
SINGLETS = Sector(0, 0)     # even parity, S_z = 0: S, S* and their continuum
HARMONICS = 4               # Floquet cutoff M
# Weak drives (linear response): gate amplitude A_eps and phase amplitude A_phi (rad).
DRIVES = dict(gate=("detuning", .01), phase=("phi", .1))
STARK_OMEGA = .6            # drive frequency below the resonance for the AC Stark shift
STARK_AMPLITUDE = .05       # gate amplitude at which the AC Stark shift is reported
STRONG_AMPLITUDE = .1       # gate amplitude at which the one- and two-photon splittings are reported
PROFILES = {
    # Omega - omega_SS* of the sweep, operating points epsilon0, amplitudes of the AC Stark and photon scans
    "quick": dict(sweep=np.linspace(-.02, .02, 7), positions=(0., .2), stark=(.05,), photons=(.1,)),
    "full": dict(sweep=np.linspace(-.02, .02, 41), positions=np.arange(13)/20, stark=np.arange(1, 16)/100,
                 photons=(.02, .03, .05, .07, .1, .14, .2)),
}


def junction(**changes):
    """Keyword arguments of reference_model: the manual's junction with four bath level pairs per lead."""
    return dict(JUNCTION, **changes, bath=cosh_grid(pairs=4, delta=1., bandwidth=10.))


def static_singlets(model, count=60):
    """The ``count`` lowest even-parity eigenstates of the undriven junction (dense diagonalization)."""
    return solve(reference_model(**model), 2, SINGLETS, SolverOptions(method="dense", eigenpairs=count))


def drive(model, kind, amplitude=None, omega=1.):
    """H(t) for a gate (``"gate"``) or phase (``"phase"``) drive p(t) = p0 + A cos(Omega t).

    The gate drive is linear in epsilon, so H(t) has only the harmonics H_0 and H_{+-1}
    and a few samples per period reproduce it exactly.
    """
    parameter, default = DRIVES[kind]
    return PeriodicHamiltonian.from_parameters(
        reference_model, model, ParameterDrive.cosine(parameter, default if amplitude is None else amplitude),
        omega, samples=8 if kind == "gate" else 16)


def at(periodic, omega):
    """The same drive at another frequency.

    The harmonics H_k of H(p0 + A cos(Omega t)) depend on A but not on Omega, so a
    frequency scan builds them once and only replaces the frequency.
    """
    return dataclasses.replace(periodic, omega=float(omega))


def wrap(value, omega):
    """Map an energy difference into [-Omega/2, Omega/2): quasienergies are defined modulo Omega."""
    return (value + omega/2) % omega - omega/2


def dressed_pair(periodic, static):
    """Quasienergies of the dressed S and S*, and the S* quasienergy moved next to eps_S by multiples of Omega."""
    result = solve_floquet(periodic, 2, SINGLETS, dict(floquet_cutoff=HARMONICS, states=2, method="shift-invert"),
                           static=static)
    eps_s, eps_star = result.quasienergies[np.argsort(result.labels)]   # the states labelled S and S*
    return eps_s, eps_s + wrap(eps_star - eps_s, periodic.omega)


def splitting(periodic, static):
    """Quasienergy distance of the dressed pair modulo Omega: the splitting at an n-photon resonance."""
    eps_s, eps_star = dressed_pair(periodic, static)
    return abs(eps_star - eps_s)


def rwa(periodic, static):
    """Two-level rotating-wave model |S, 0>, |S*, -1>: (coupling 2|(H_1)_SS*|, resonance, eigenvalues)."""
    matrix = rwa_model(periodic, [0, 1], [0, -1], 2, SINGLETS, static=static).effective_matrix
    resonance = periodic.omega + (matrix[1, 1] - matrix[0, 0]).real   # where the two diagonal entries meet
    return 2*abs(matrix[0, 1]), float(resonance), np.linalg.eigvalsh(matrix)


def resonance(periodic, static, center, width, tolerance=1e-6):
    """Resonance frequency and splitting from a fit of the dressed-pair hyperbola.

    Near an n-photon resonance the dressed states are split by
    w(Omega) = sqrt(n^2 (Omega - Omega_r)^2 + delta^2), so w^2 is a parabola in Omega. Three
    points around ``center`` fix it; the fit is repeated around the new vertex until it no
    longer moves. AC Stark and Bloch-Siegert shifts move Omega_r away from omega_SS*/n.
    """
    offsets = np.array([-width, 0., width])
    for _ in range(8):
        squares = [splitting(at(periodic, center + x), static)**2 for x in offsets]
        a, b, c = np.polyfit(offsets, squares, 2)
        shift = -b/(2*a)
        center += shift
        if abs(shift) < tolerance:
            break
    return float(center), float(np.sqrt(max(c - b*b/(4*a), 0.)))


def operating_point(position):
    """Static matrix elements and the resonant splittings of both weak drives at epsilon0 = ``position``."""
    model = junction(detuning=float(position))
    static = static_singlets(model)
    row = dict(position=float(position), omega0=float(static.energies[1] - static.energies[0]),
               charge=float(abs(static.matrix_elements("charge")[0, 1])),
               current=float(abs(static.matrix_elements("current")[0, 1])))
    for kind, (_, amplitude) in DRIVES.items():
        periodic = drive(model, kind, omega=row["omega0"])
        coupling, guess, _ = rwa(periodic, static)
        omega_r, delta = resonance(periodic, static, guess, max(coupling, 1e-3))
        row[kind] = dict(resonance=omega_r, splitting=delta, rwa=coupling, per_amplitude=delta/amplitude)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/07_gate_drive/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Static junction: the pair transition and the matrix elements of both drive operators.
    model = junction()
    static = static_singlets(model)
    ground = float(static.energies[0])
    omega0 = float(static.energies[1] - static.energies[0])
    charge = float(abs(static.matrix_elements("charge")[0, 1]))
    current = float(abs(static.matrix_elements("current")[0, 1]))
    print(f"omega_SS* = {omega0:.6f}, |<S|n_d|S*>| = {charge:.6f}, |<S|dH/dphi|S*>| = {current:.6f}")

    # 2. Rotating-frame quasienergies through the resonance, gate and phase drive (Sambe and RWA).
    sweep = []
    drives = {kind: drive(model, kind) for kind in DRIVES}
    for omega in omega0 + settings["sweep"]:
        row = dict(omega=float(omega))
        for kind, periodic in drives.items():
            row[kind] = sorted(float(e) for e in dressed_pair(at(periodic, omega), static))
            row[f"{kind}_rwa"] = rwa(at(periodic, omega), static)[2].tolist()
        sweep.append(row)

    # 3. Resonant splittings versus the operating point: the eta selection rule.
    positions = [operating_point(position) for position in settings["positions"]]
    reference = min(positions, key=lambda row: abs(row["position"] - JUNCTION["detuning"]))
    symmetric = min(positions, key=lambda row: abs(row["position"]))

    # 4. AC Stark shift at Omega = 0.6 < omega_SS*. The gate drive leaves the time average H_0 = H(epsilon0)
    #    unchanged, so the shift is purely second order: sum over all static states b of the sector,
    #    delta eps_S = -(A^2/4) sum_b |<S|n_d|b>|^2 [1/(E_b - E_S - Omega) + 1/(E_b - E_S + Omega)].
    complete = static_singlets(model, count=static.basis.dimension)          # every state of the sector
    excitation = complete.energies[1:] - complete.energies[0]
    weight = np.abs(complete.matrix_elements("charge")[0, 1:])**2
    response = weight*(1/(excitation - STARK_OMEGA) + 1/(excitation + STARK_OMEGA))
    stark = []
    for amplitude in settings["stark"]:
        periodic = drive(model, "gate", float(amplitude), STARK_OMEGA)
        dressed = solve_floquet(periodic, 2, SINGLETS, dict(floquet_cutoff=HARMONICS, states=1, method="shift-invert"),
                                static=static)
        stark.append(dict(amplitude=float(amplitude), sambe=float(dressed.quasienergies[0]) - ground,
                          second_order=float(-amplitude**2/4*response.sum()),
                          pair_only=float(-amplitude**2/4*response[0])))
    reported = min(stark, key=lambda row: abs(row["amplitude"] - STARK_AMPLITUDE))
    averaged = solve(drive(model, "gate", STARK_AMPLITUDE, STARK_OMEGA).time_average(), 2, SINGLETS,
                     SolverOptions(method="dense", eigenpairs=1))

    # 5. One- and two-photon resonances versus the gate amplitude.
    photons = []
    for amplitude in settings["photons"]:
        periodic = drive(model, "gate", float(amplitude), omega0)
        one = resonance(periodic, static, omega0, amplitude*charge)
        two = resonance(periodic, static, omega0/2, 1e-3)
        photons.append(dict(amplitude=float(amplitude), one=one[1], one_resonance=one[0],
                            two=two[1], two_resonance=two[0]))
    strong = min(photons, key=lambda row: abs(row["amplitude"] - STRONG_AMPLITUDE))

    values = dict(
        omega0=omega0, charge_element=charge, current_element=current, element_ratio=charge/current,
        gate_resonance=reference["gate"]["resonance"], gate_splitting=reference["gate"]["splitting"],
        gate_ratio=reference["gate"]["splitting"]/(DRIVES["gate"][1]*charge),
        phase_resonance=reference["phase"]["resonance"], phase_splitting=reference["phase"]["splitting"],
        phase_ratio=reference["phase"]["splitting"]/(DRIVES["phase"][1]*current),
        symmetric_charge_element=symmetric["charge"], symmetric_current_element=symmetric["current"],
        symmetric_gate_per_amplitude=symmetric["gate"]["per_amplitude"],
        symmetric_phase_per_amplitude=symmetric["phase"]["per_amplitude"],
        averaged_shift=abs(float(averaged.energies[0]) - ground),
        stark_sambe=reported["sambe"], stark_second_order=reported["second_order"],
        stark_pair_only=reported["pair_only"],
        one_photon_strong=strong["one"], two_photon_strong=strong["two"],
        two_photon_resonance=strong["two_resonance"],
    )
    if len(positions) > 2:
        values["gate_max_deviation"] = max(abs(row["gate"]["per_amplitude"] - row["charge"]) for row in positions)
        for position, key in ((.5, "phase_ratio_05"), (.6, "phase_ratio_06")):
            row = min(positions, key=lambda entry: abs(entry["position"] - position))
            values[key] = row["phase"]["per_amplitude"]/row["current"]
    if len(photons) > 2:
        # Power laws from the weak-drive part of the scan: delta_1 ~ A, delta_2 ~ A^2.
        weak = [row for row in photons if row["amplitude"] <= STRONG_AMPLITUDE]
        logs = np.log([row["amplitude"] for row in weak])
        values["one_photon_exponent"] = float(np.polyfit(logs, np.log([row["one"] for row in weak]), 1)[0])
        values["two_photon_exponent"] = float(np.polyfit(logs, np.log([row["two"] for row in weak]), 1)[0])
    for key, value in values.items():
        print(f"{key:30s} {value:.6g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values, ground=ground,
                  amplitudes={kind: amplitude for kind, (_, amplitude) in DRIVES.items()},
                  sweep=sweep, positions=positions, stark=dict(omega=STARK_OMEGA, rows=stark), photons=photons)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"gate_drive.png")
    print(f"Saved {output}")


def plot(record, path):
    """Four panels: the two resonances, the selection rule, the AC Stark shift, one- and two-photon splittings."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    values, amplitudes = record["values"], record["amplitudes"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.6), layout="constrained")

    axis = axes[0, 0]
    rows = record["sweep"]
    omega = np.array([r["omega"] for r in rows])
    for kind, color, label in (("gate", "#2a78d6", fr"gate, $A_\epsilon={amplitudes['gate']}$"),
                               ("phase", "#eb6834", fr"phase, $A_\phi={amplitudes['phase']}$")):
        axis.plot(omega, 1e3*(np.array([r[kind] for r in rows]) - record["ground"]), color=color, label=[label, None])
        axis.plot(omega, 1e3*(np.array([r[f"{kind}_rwa"] for r in rows]) - record["ground"]), ":", color=color,
                  lw=1)
    axis.axvline(values["omega0"], color="gray", lw=.6, ls=":")
    axis.set(xlabel=r"drive frequency $\Omega$", ylabel=r"$(\varepsilon-E_S)\times10^3$, rotating frame",
             ylim=(-12, 12))
    axis.legend(frameon=False, fontsize=8, loc="upper left", title="solid: Sambe, dotted: RWA", title_fontsize=8)

    axis = axes[0, 1]
    rows = record["positions"]
    position = np.array([r["position"] for r in rows])
    axis.plot(position, [r["charge"] for r in rows], color="#2a78d6", lw=1.1,
              label=r"gate: $|\langle S|n_d|S^*\rangle|$")
    axis.plot(position, [r["gate"]["per_amplitude"] for r in rows], "o", mfc="none", color="#2a78d6")
    axis.plot(position, [10*r["current"] for r in rows], color="#eb6834", lw=1.1,
              label=r"phase: $10\,|\langle S|\partial_\phi H|S^*\rangle|$")
    axis.plot(position, [10*r["phase"]["per_amplitude"] for r in rows], "s", mfc="none", color="#eb6834")
    axis.set(xlabel=r"detuning $\epsilon_0$", ylabel="resonant splitting / amplitude", ylim=(0, 1.))
    axis.legend(frameon=False, fontsize=8, loc="lower right", title="lines: static, markers: Sambe",
                title_fontsize=8)

    axis = axes[1, 0]
    rows = record["stark"]["rows"]
    amplitude = np.array([r["amplitude"] for r in rows])
    line = np.linspace(0, amplitude.max(), 100)
    scale = rows[-1]["second_order"]/rows[-1]["amplitude"]**2
    pair = rows[-1]["pair_only"]/rows[-1]["amplitude"]**2
    axis.plot(amplitude, 1e3*np.array([r["sambe"] for r in rows]), "o", mfc="none", color="#2a78d6", label="Sambe")
    axis.plot(line, 1e3*scale*line**2, color="#1baf7a", lw=1.1, label="second order, all states")
    axis.plot(line, 1e3*pair*line**2, "--", color="gray", lw=1, label=r"second order, $S^*$ only")
    axis.set(xlabel=r"gate amplitude $A_\epsilon$", ylabel=r"$(\varepsilon_S-E_S)\times10^3$",
             title=fr"AC Stark shift at $\Omega={record['stark']['omega']}$")
    axis.legend(frameon=False, fontsize=8, loc="lower left")

    axis = axes[1, 1]
    rows = record["photons"]
    amplitude = np.array([r["amplitude"] for r in rows])
    for key, order, color, label in (("one", 1, "#2a78d6", r"one photon, $\Omega\approx\omega_{SS^*}$"),
                                     ("two", 2, "#eb6834", r"two photons, $\Omega\approx\omega_{SS^*}/2$")):
        splittings = np.array([r[key] for r in rows])
        axis.loglog(amplitude, splittings, "o" if order == 1 else "s", mfc="none", color=color, label=label)
        axis.loglog(amplitude, splittings[0]*(amplitude/amplitude[0])**order, ":", color=color, lw=1)
    axis.set(xlabel=r"gate amplitude $A_\epsilon$", ylabel="splitting at resonance")
    axis.set_xticks([.02, .05, .1, .2], ["0.02", "0.05", "0.1", "0.2"])
    axis.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axis.legend(frameon=False, fontsize=8, loc="upper left", title=r"dotted: $\propto A_\epsilon$, $A_\epsilon^2$",
                title_fontsize=8)
    for axis, label in zip(axes.ravel(), "abcd", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=100)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
