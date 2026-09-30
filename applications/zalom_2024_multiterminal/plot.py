"""Plot the solver's phase boundary and three-terminal phase diagram."""

from pathlib import Path

import numpy as np

from applications._common import plot_arguments, read_csv, read_json, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    data = read_csv(args.output/"spectrum.csv")
    summary = read_json(args.output/"summary.json")
    critical = summary["critical"]["chi_critical"]
    axes[0].plot([float(r["chi"]) for r in data], [float(r["signed_gap"]) for r in data], "o-", ms=3)
    axes[0].axhline(0, color="0.6", lw=.5)
    axes[0].axvline(.721, color="k", ls="--", label="Published NRG")
    axes[0].axvline(critical, color="C1", label="QP crossing")
    axes[0].set(xlabel=r"$\chi$", ylabel=r"$(E_D-E_S)/\Delta$", title="Computed parity crossing")
    axes[0].legend(fontsize=8)
    data = read_csv(args.output/"phase_map.csv")
    size = int(round(np.sqrt(len(data))))
    x = np.array([float(r["phi2_over_pi"]) for r in data]).reshape(size, size)
    y = np.array([float(r["phi3_over_pi"]) for r in data]).reshape(size, size)
    chi = np.array([float(r["chi"]) for r in data]).reshape(size, size)
    parity = np.array([int(r["ground_parity"]) for r in data]).reshape(size, size)
    from matplotlib.colors import ListedColormap
    axes[1].pcolormesh(x, y, parity, cmap=ListedColormap(["#b4d3ef", "#efa9a4"]), shading="nearest")
    axes[1].contour(x, y, chi, levels=[critical], colors="C0", linewidths=1.5)
    axes[1].contour(x, y, chi, levels=[.721], colors="k", linestyles="--", linewidths=.8)
    theta = summary["high_symmetry"]["phi2_over_pi"]
    axes[1].plot([theta, -theta], [-theta, theta], "*", color="green", ms=9)
    axes[1].text(0, 0, "Singlet", ha="center")
    axes[1].text(.65, -.25, "Doublet", ha="center")
    axes[1].set(xlabel=r"$\phi_2/\pi$", ylabel=r"$\phi_3/\pi$", aspect="equal",
                title=r"Fig. 2(d): $\gamma=(0.30,0.35,0.35)$")
    axes[1].legend(handles=[Line2D([], [], color="C0", label="QP boundary"),
                             Line2D([], [], color="k", ls="--", label="Published NRG boundary")],
                    fontsize=7, loc="upper right")
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
