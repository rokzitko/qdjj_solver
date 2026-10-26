"""Tutorial 01: subgap states of the reference junction, the singlet and the doublet.

    python -B examples/01_subgap_states/run.py                  # quick profile, about 3 s
    python -B examples/01_subgap_states/run.py --profile full   # dense scans, about 10 s

A quantum dot between two superconducting leads binds a singlet S, an excited singlet S*
(even parity) and a spin doublet D (odd parity) inside the gap. The script computes
(1) their energies at the reference point of the physics manual, (2) the spectrum versus
the phase difference and the singlet-doublet (0-pi) transition, (3) the spectrum versus the
detuning, and (4) the supercurrent <dH/dphi> of S and D, checked against a finite difference
of the energy. Results go to runs/examples/01_subgap_states/<profile>/ (results.json, and
subgap_states.png if matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from qdjj_solver import Sector, cosh_grid, reference_model  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, solve  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The reference junction of the physics manual: energies in units of the gap, Delta = 1.
# gamma is the total hybridization (0.25 per lead), detuning = epsilon_d + U/2, phases in radians.
CONFIG = json.loads((HERE/"input"/"reference.json").read_text())
JUNCTION = {k: v for k, v in CONFIG["model"].items() if k != "kind"}
BATH = {k: v for k, v in CONFIG["bath"].items() if k != "kind"}
CUTOFF = CONFIG["cutoff"]                        # at most two bath quasiparticles (QP cutoff)
STEP = 1e-4                                      # phase step of the finite-difference check
PROFILES = {
    # points of the phase scan over [0, 2 pi], of the detuning scan over [-1.2, 1.2], phases of the current check
    "quick": dict(phases=25, detunings=25, checks=np.array([.4, .8, 1.6, 2.4])),
    "full": dict(phases=121, detunings=97, checks=np.linspace(.1, 2*np.pi - .1, 24)),
}


def model(**changes):
    """The reference junction with some parameters changed."""
    return reference_model(cosh_grid(**BATH), **{**JUNCTION, **changes})


def subgap(method="auto", **changes):
    """Energies and currents of S and S* (even sector) and of D (odd sector, S_z = +1/2)."""
    h = model(**changes)
    even = solve(h, CUTOFF, Sector(0, 0), SolverOptions(method=method, eigenpairs=2))
    odd = solve(h, CUTOFF, Sector(1, 1), SolverOptions(method=method, eigenpairs=1))
    current = np.real(even.observables["phase_derivative"])     # <dH/dphi> in units of Delta
    return dict(singlet=float(even.energies[0]), excited=float(even.energies[1]), doublet=float(odd.energies[0]),
                current_singlet=float(current[0]), current_excited=float(current[1]),
                current_doublet=float(np.real(odd.observables["phase_derivative"][0])))


def scan(name, values):
    """``subgap`` along one parameter, as a dict of arrays."""
    rows = [subgap(**{name: float(value)}) for value in values]
    return {key: np.array([row[key] for row in rows]) for key in rows[0]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/01_subgap_states/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. The reference point. The global ground state is the lower of the two sector minima.
    point = subgap()
    ground_parity = int(point["doublet"] < point["singlet"])
    print(f"E_S = {point['singlet']:.6f}, E_S* - E_S = {point['excited'] - point['singlet']:.6f}, "
          f"E_D - E_S = {point['doublet'] - point['singlet']:.6f}, ground-state parity {ground_parity}")
    # At zero field the doublet is a spin-1/2 multiplet: the S_z = -1/2 partner is degenerate.
    partner = solve(model(), CUTOFF, Sector(1, -1), SolverOptions(eigenpairs=1)).energies[0]

    # 2. Phase scan: the doublet becomes the ground state around phi = pi, where |Gamma_phi| is small.
    phases = np.linspace(0., 2*np.pi, settings["phases"])
    by_phase = scan("phi", phases)
    gap = by_phase["doublet"] - by_phase["singlet"]
    # Bracket the first sign change of E_D - E_S on the scan, then refine it with Brent's method.
    first = int(np.flatnonzero(np.diff(np.sign(gap)) != 0)[0])
    critical = brentq(lambda phi: (lambda r: r["doublet"] - r["singlet"])(subgap(phi=phi)),
                      phases[first], phases[first + 1], xtol=1e-13)
    at_pi = subgap(phi=np.pi)

    # 3. Detuning scan at the reference phase; the reference junction is particle-hole symmetric in epsilon.
    detunings = np.linspace(-1.2, 1.2, settings["detunings"])
    by_detuning = scan("detuning", detunings)
    asymmetry = max(float(np.max(abs(by_detuning[key] - by_detuning[key][::-1])))
                    for key in ("singlet", "excited", "doublet"))
    half_filling = subgap(detuning=0.)

    # 4. Supercurrent: Hellmann-Feynman <dH/dphi> against the central difference of the energy.
    deviation, checks = 0., []
    for phi in settings["checks"]:
        plus, minus, here = subgap(phi=phi + STEP), subgap(phi=phi - STEP), subgap(phi=phi)
        for state in ("singlet", "doublet"):
            slope = (plus[state] - minus[state])/(2*STEP)
            deviation = max(deviation, abs(slope - here["current_" + state]))
            checks.append(dict(phi=float(phi), state=state, finite_difference=slope))

    # 5. Atomic limit Gamma = 0: E_|0> = U/2 - epsilon, E_|2> = U/2 + epsilon, E_|sigma> = 0. The decoupled
    # Hamiltonian is diagonal with highly degenerate levels and an exactly zero doublet energy, a case that
    # plain ARPACK gets wrong; the solver's shifted iteration handles it (see docs/numerics.md).
    atomic = subgap(gamma=0.)
    u, epsilon = JUNCTION["u"], JUNCTION["detuning"]
    atomic_deviation = max(abs(atomic["singlet"] - (u/2 - epsilon)), abs(atomic["excited"] - (u/2 + epsilon)),
                           abs(atomic["doublet"]))

    values = dict(
        singlet_energy=point["singlet"],
        excited_energy=point["excited"] - point["singlet"],
        doublet_energy=point["doublet"] - point["singlet"],
        ground_parity=ground_parity,
        doublet_sz_splitting=float(abs(partner - point["doublet"])),
        current_singlet=point["current_singlet"],
        current_doublet=point["current_doublet"],
        current_deviation=float(deviation),
        critical_phase=float(critical),
        ground_parity_pi=int(at_pi["doublet"] < at_pi["singlet"]),
        doublet_energy_pi=at_pi["doublet"] - at_pi["singlet"],
        doublet_energy_half_filling=half_filling["doublet"] - half_filling["singlet"],
        particle_hole_asymmetry=asymmetry,
        atomic_deviation=float(atomic_deviation),
    )
    for key, value in values.items():
        print(f"{key:28s} {value:.10g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values,
                  parameters=dict(JUNCTION, bath=BATH, cutoff=CUTOFF),
                  phase_scan=dict(phi=phases.tolist(), **{k: v.tolist() for k, v in by_phase.items()}),
                  detuning_scan=dict(detuning=detunings.tolist(), **{k: v.tolist() for k, v in by_detuning.items()}),
                  current_checks=checks)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"subgap_states.png")
    print(f"Saved {output}")


def plot(record, path):
    """Spectrum versus phase and detuning, measured from the ground state, and the supercurrents."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    critical = record["values"]["critical_phase"]/np.pi
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), layout="constrained")
    for axis, (name, data) in zip(axes[:2], (("phi", record["phase_scan"]), ("detuning", record["detuning_scan"])),
                                  strict=True):
        x = np.array(data[name])/(np.pi if name == "phi" else 1.)
        singlet, excited, doublet = (np.array(data[k]) for k in ("singlet", "excited", "doublet"))
        ground = np.minimum(singlet, doublet)
        axis.axhspan(1., 2.2, color="#f0efec", lw=0, zorder=0)       # ground state + one quasiparticle (Delta)
        axis.text(.03, .95, "continuum", transform=axis.transAxes, color=gray, fontsize=8, va="top")
        axis.plot(x, singlet - ground, color=blue, label=r"singlet $S$")
        axis.plot(x, excited - ground, color=blue, ls="--", label=r"excited singlet $S^*$")
        axis.plot(x, doublet - ground, color=orange, label=r"doublet $D$")
        axis.set(xlim=(x[0], x[-1]), ylim=(-.05, 1.6), ylabel=r"$E-E_{\mathrm{ground}}$ ($\Delta$)")
    for value in (critical, 2 - critical):
        axes[0].axvline(value, color=gray, lw=.8, ls=":")
    axes[0].set_xlabel(r"phase $\phi/\pi$")
    axes[0].legend(frameon=False, fontsize=8, loc="center left")
    axes[1].set_xlabel(r"detuning $\epsilon/\Delta$  ($\phi=0.8$)")
    data = record["phase_scan"]
    x = np.array(data["phi"])/np.pi
    axis = axes[2]
    axis.axhline(0., color=gray, lw=.6)
    axis.plot(x, data["current_singlet"], color=blue, label=r"$S$")
    axis.plot(x, data["current_doublet"], color=orange, label=r"$D$")
    odd = np.array(data["doublet"]) < np.array(data["singlet"])
    ground = np.where(odd, data["current_doublet"], data["current_singlet"])
    for part in np.split(np.arange(len(x)), np.flatnonzero(np.diff(odd)) + 1):   # the current jumps at phi_c
        axis.plot(x[part], ground[part], color=aqua, lw=2.4, alpha=.45, label="ground state" if part[0] == 0 else None)
    for state in ("singlet", "doublet"):
        rows = [r for r in record["current_checks"] if r["state"] == state]
        axis.plot([r["phi"]/np.pi for r in rows], [r["finite_difference"] for r in rows], "o", ms=3.5, mfc="none",
                  mew=.9, color=gray, label=r"$\partial E/\partial\phi$" if state == "singlet" else None)
    for value in (critical, 2 - critical):
        axis.axvline(value, color=gray, lw=.8, ls=":")
    axis.set(xlim=(x[0], x[-1]), xlabel=r"phase $\phi/\pi$",
             ylabel=r"current $I$ ($2e\Delta/\hbar$)")
    axis.legend(frameon=False, fontsize=8, loc="lower left")
    for axis, label in zip(axes, "abc", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=110)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
