"""Plot the computed parity boundary and gate-dependent subgap excitation."""

from pathlib import Path

import numpy as np

from applications._common import plot_arguments, read_csv, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    rows = read_csv(args.output/"boundary.csv")
    phase = np.array([float(r["phi_over_pi"]) for r in rows])
    gate = np.array([float(r["detuning_critical_over_u"]) for r in rows])
    phases = np.r_[-phase[:0:-1], phase]
    gates = np.r_[gate[:0:-1], gate]
    axes[0].fill_betweenx(phases, -gates, gates, color="#c0b1db", label="Doublet")
    axes[0].set_facecolor("#ffe3bf")
    for sign in (-1, 1):
        axes[0].plot(sign*gates, phases, color="C0", lw=1.5, label="QP boundary" if sign == 1 else None)
        ref = read_csv(CASE/"reference"/"boundary.csv")
        for side in (-1, 1):
            axes[0].plot([sign*float(r["detuning_critical_over_u"]) for r in ref],
                         [side*float(r["phi_over_pi"]) for r in ref], "k.", ms=3,
                         label="Deposited NRG" if sign == 1 and side == 1 else None)
    axes[0].text(.49, 0, "Singlet", rotation=90, ha="center", va="center")
    axes[0].set(xlabel=r"$\xi/U=(\epsilon_d+U/2)/U$", ylabel=r"$\phi/\pi$", xlim=(-.6, .6),
                title="Fig. S1(c): gate--phase parity diagram")
    axes[0].legend(fontsize=7, loc="upper left")
    rows = read_csv(args.output/"gate_scan.csv")
    for p in sorted({float(r["phi_over_pi"]) for r in rows}):
        data = [r for r in rows if float(r["phi_over_pi"]) == p]
        axes[1].plot([float(r["gate_over_u"]) for r in data],
                     [float(r["signed_gap"]) for r in data], label=rf"$\phi/\pi={p:g}$")
    axes[1].axhline(0, color="0.5", lw=.5)
    axes[1].set(xlabel=r"$\xi/U$", ylabel=r"$(E_D-E_S)/\Delta$", title="Parity-changing excitation branches")
    axes[1].legend(fontsize=8)
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
