"""Tutorial 08: two-tone spectroscopy of the pair transition of a gate-driven junction.

    python -B examples/08_two_tone/run.py                  # quick profile, a few seconds
    python -B examples/08_two_tone/run.py --profile full   # the complete calculation, about a minute

A strong gate pump drives the S-S* pair transition; a weak probe measures the dressed
states. The script computes (1) the probe lines versus pump detuning, compared with the
rotating-wave prediction, (2) the Mollow pattern at resonance and (3) the Autler-Townes
doublet seen by a tunnelling probe from the doublet. Results go to
runs/examples/08_two_tone/<profile>/ (results.json, and two_tone.png if matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The junction of the physics manual: energies in units of the gap, Delta = 1.
JUNCTION = dict(u=1., gamma=.5, phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
PUMP = .03          # gate amplitude A_eps (amplitude of the detuning)
WIDTH = .0015       # half width of the plotted lines
PROFILES = {
    # bath level pairs per lead, Floquet cutoff M, minimal static basis K, pump detunings
    "quick": dict(pairs=2, harmonics=5, basis_states=20, detunings=np.linspace(-.04, .04, 9)),
    "full": dict(pairs=4, harmonics=6, basis_states=43, detunings=np.linspace(-.06, .06, 61)),
}


def basis_size(energies, minimum):
    """Smallest K >= minimum that does not split a degenerate static multiplet."""
    gaps = np.flatnonzero(np.diff(energies) > 1e-7) + 1
    gaps = gaps[gaps >= minimum]
    return int(gaps[0]) if len(gaps) else len(energies)


def dressed(model, omega, settings, static_even, size, static_odd=None):
    """All even Floquet states in the static basis and, if ``static_odd`` is given, the dressed doublet."""
    drive = ParameterDrive.cosine("detuning", PUMP)
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model, drive, omega, samples=16)
    even = solve_floquet(periodic, 2, Sector(0, 0),
                         dict(floquet_cutoff=settings["harmonics"], selection="all", method="dense", basis="eigen",
                              basis_states=size, dense_dimension=20000),
                         static=static_even)
    if static_odd is None:
        return even, None
    odd = solve_floquet(periodic, 2, Sector(1, 1),
                        dict(floquet_cutoff=settings["harmonics"], states=1, basis="eigen",
                             basis_states=basis_size(static_odd.energies, 16)), static=static_odd)
    return even, odd


def strongest_near(frequencies, strengths, center, window):
    """The strongest line within ``window`` of ``center``: (frequency, strength)."""
    near = np.flatnonzero(abs(frequencies - center) < window)
    best = near[np.argmax(strengths[near])]
    return float(frequencies[best]), float(strengths[best])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/08_two_tone/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Static junction: the pair transition and its gate matrix element.
    model = dict(JUNCTION, bath=cosh_grid(pairs=settings["pairs"], delta=1., bandwidth=10.))
    h = reference_model(**model)
    static_even = solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=settings["basis_states"] + 40))
    static_odd = solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=40))
    size = basis_size(static_even.energies, settings["basis_states"])
    omega0 = float(static_even.energies[1] - static_even.energies[0])
    element = float(abs(static_even.matrix_elements("charge")[0, 1]))
    rabi = PUMP*element                     # rotating-wave splitting delta = A_eps |<S|n_d|S*>|
    print(f"omega_SS* = {omega0:.6f}, |<S|n_d|S*>| = {element:.6f}, delta = {rabi:.6f}, K = {size}")

    # 2. Probe lines versus pump detuning: the two dressed lines follow Omega -+ w,
    #    w = sqrt((Omega - omega_SS*)^2 + delta^2).
    sweep, deviation = [], 0.
    for detuning in settings["detunings"]:
        omega = omega0 + float(detuning)
        even, _ = dressed(model, omega, settings, static_even, size)
        lines = even.matrix_elements("charge").lines(initial=even.dressed_state(0))
        w = np.hypot(detuning, rabi)
        keep = np.flatnonzero((abs(abs(lines.frequencies) - omega) < 3*rabi + abs(detuning)) & ~lines.elastic)
        strongest = keep[np.argmax(lines.strengths[keep])]                # the dominant dressed line
        predicted = np.array([omega - w, omega + w])
        deviation = max(deviation, float(np.min(abs(abs(lines.frequencies[strongest]) - predicted))))
        sweep.append(dict(detuning=float(detuning), frequencies=lines.frequencies.tolist(),
                          strengths=lines.strengths.tolist()))

    # 3. At resonance: the Mollow pattern of the pumped transition, out of the dressed ground state.
    even, odd = dressed(model, omega0, settings, static_even, size, static_odd)
    mollow = even.matrix_elements("charge").lines(initial=even.dressed_state(0))
    quarter = element**2/4
    f, s = mollow.frequencies, mollow.strengths
    lower_f, lower = strongest_near(-f, s, omega0 - rabi, rabi/2)      # stimulated emission at -(Omega - delta)
    upper_f, upper = strongest_near(f, s, omega0 + rabi, rabi/2)       # absorption at Omega + delta
    _, central = strongest_near(f, s, omega0, rabi/2)                  # elastic line at +Omega

    # 4. A spectator: removing a spin-up electron from the doublet probes the dressed singlets.
    removal = even.matrix_elements("d_up", odd).lines(initial=odd.dressed_state(0))
    center = float(static_odd.energies[0] - static_even.energies[0])
    static_weight = float(abs(static_even.matrix_elements("d_up", static_odd)[0, 0])**2)
    near = np.flatnonzero(abs(-removal.frequencies - center) < 2*rabi)
    pair = near[np.argsort(removal.strengths[near])[-2:]]
    doublet = np.sort(-removal.frequencies[pair])

    values = dict(
        omega0=omega0, element=element, rabi=rabi, basis_states=size, removal_energy=center,
        rwa_deviation=float(deviation),
        mollow_splitting=float((upper_f - lower_f)/2),
        mollow_lower=lower/quarter, mollow_central=central/quarter, mollow_upper=upper/quarter,
        autler_townes_splitting=float(doublet[1] - doublet[0]),
        autler_townes_weight=float(removal.strengths[pair].sum()/static_weight),
    )
    for key, value in values.items():
        print(f"{key:24s} {value:.6g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values,
                  sweep=sweep, mollow=mollow.record(), removal=removal.record())
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"two_tone.png")
    print(f"Saved {output}")


def plot(record, path):
    """Three panels: lines versus detuning, the Mollow rate and net spectra, the Autler-Townes doublet."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    from qdjj_solver import line_shape
    values = record["values"]
    omega0, rabi = values["omega0"], values["rabi"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), layout="constrained")
    x, y, s = [], [], []
    for row in record["sweep"]:
        x += [row["detuning"]]*len(row["frequencies"])
        y += [abs(f) for f in row["frequencies"]]
        s += row["strengths"]
    s = np.array(s)
    order = np.argsort(s)
    axes[0].scatter(1e3*np.array(x)[order], np.array(y)[order], c=s[order]/s.max(), cmap="Blues", vmin=0, vmax=1,
                    s=6, lw=0)
    detuning = np.linspace(min(x), max(x), 200)
    for sign in (-1, 1):
        axes[0].plot(1e3*detuning, omega0 + detuning + sign*np.hypot(detuning, rabi), "--", color="#eb6834", lw=.8)
    axes[0].set(xlabel=r"pump detuning $(\Omega-\omega_{SS^*})\times10^3$", ylabel=r"probe frequency $|\omega|$",
                ylim=(omega0 - .08, omega0 + .08))
    grid = np.linspace(omega0 - .05, omega0 + .05, 1001)
    lines = record["mollow"]
    def signed(x, freqs, strengths):
        return sum(st*line_shape(x - f, WIDTH) for f, st in zip(freqs, strengths, strict=True))
    rate = signed(grid, lines["frequencies"], lines["strengths"]) + signed(-grid, lines["frequencies"], lines["strengths"])
    net = signed(grid, lines["frequencies"], lines["strengths"]) - signed(-grid, lines["frequencies"], lines["strengths"])
    axes[1].plot(grid, rate, color="#2a78d6", label="rate")
    axes[1].plot(grid, net, color="#eb6834", label="absorption $-$ emission")
    axes[1].set(xlabel=r"probe frequency $\omega$", ylabel="gate-probe spectrum")
    axes[1].legend(frameon=False, fontsize=8)
    removal = record["removal"]
    center = values["removal_energy"]
    energies = np.linspace(center - .04, center + .04, 801)
    axes[2].plot(energies, signed(-energies, removal["frequencies"], removal["strengths"]), color="#2a78d6")
    axes[2].set(xlabel=r"removal energy $-\omega$", ylabel=r"$d_\uparrow$ spectrum from $D$")
    for axis, label in zip(axes, "abc", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=110)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
