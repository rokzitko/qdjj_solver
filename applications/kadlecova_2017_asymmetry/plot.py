"""Plot the symmetric/asymmetric phase-boundary comparison."""

import numpy as np

from applications._common import plot_arguments, read_csv, save_figure
from .run import CASE, inverse_phase


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    rows = read_csv(args.output/"boundary.csv")
    reference = read_csv(CASE/"reference"/"figure1.csv")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    for ax, a in zip(axes, (1., 11.), strict=True):
        for u, color in zip((2., 3.2, 5.), ("C0", "C1", "C2"), strict=True):
            data = [r for r in rows if float(r["asymmetry"]) == a and float(r["u_meV"]) == u and r["tilde_epsilon"]]
            ax.plot([float(r["tilde_epsilon"]) for r in data], [float(r["phi_over_pi"]) for r in data],
                    "o-", color=color, ms=3, label=rf"QP, $U={u:g}$ meV")
            ref = [(float(r["tilde_epsilon"]), inverse_phase(float(r["phi_over_pi"])*np.pi, a))
                   for r in reference if float(r["u_meV"]) == u]
            ref = [(x, p/np.pi) for x, p in ref if p is not None]
            ax.plot([r[0] for r in ref], [r[1] for r in ref], ".--", color=color, lw=1,
                    label="published NRG / mapping" if u == 2 else None)
        ax.set(xlabel=r"$\tilde\epsilon=2\,\mathrm{detuning}/U$", title=rf"$a={a:g}$", ylim=(0, 1))
        ax.legend(fontsize=8)
    axes[0].set_ylabel(r"$\phi_c/\pi$")
    fig.suptitle("Kadlecová et al. (2017), selected Fig. 1 curves")
    fig.tight_layout()
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
