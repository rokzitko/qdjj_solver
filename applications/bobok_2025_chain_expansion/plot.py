"""Plot the archived GAL/ChE current and common-lead spin benchmarks."""

import numpy as np

from applications._common import plot_arguments, read_csv, read_json, save_figure
from .run import CASE


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    current = read_csv(args.output/"current.csv")
    crossings = read_csv(args.output/"crossings.csv") if (args.output/"crossings.csv").exists() else []
    reference8 = read_csv(CASE/"reference"/"figure8.csv")
    fig, axes = plt.subplots(3, 2, figsize=(9, 9), sharex=True, sharey=True)
    for ax, setting in zip(axes.flat, ("GAL", "ChE-F1", "ChE-W2", "ChE-F2", "ChE-W4", "ChE-F4"), strict=True):
        for u, color in ((2., "black"), (4., "C3"), (8., "C0")):
            rows = sorted((r for r in current if r["setting"] == setting and float(r["u"]) == u),
                          key=lambda r: float(r["phi_over_pi"]))
            crossing = next((float(r["critical_parameter"]) for r in crossings if int(r["figure"]) == 8
                             and r["setting"] == setting and float(r["u"]) == u), None)
            if crossing is None:
                x, y = [float(r["phi_over_pi"]) for r in rows], [float(r["current"]) for r in rows]
            else:
                # Use the refined transition, rather than drawing a sloped
                # interpolation across a first-order current jump.
                limits = []
                for branch in ("even_current", "odd_current"):
                    data = [r for r in rows if r[branch]]
                    limits.append(float(np.interp(crossing, [float(r["phi_over_pi"]) for r in data],
                                                  [float(r[branch]) for r in data])))
                left = [(float(r["phi_over_pi"]), float(r["current"])) for r in rows if float(r["phi_over_pi"]) < crossing]
                right = [(float(r["phi_over_pi"]), float(r["current"])) for r in rows if float(r["phi_over_pi"]) > crossing]
                points = left+[(crossing, value) for value in limits]+right
                x, y = zip(*points, strict=True)
            ax.plot(x, y, color=color, lw=1.4, label=rf"$U/\Delta={u:g}$")
            ref = [r for r in reference8 if r["method"] == "NRG" and float(r["u"]) == u]
            ax.scatter([float(r["phi_over_pi"]) for r in ref], [float(r["current"]) for r in ref],
                       facecolors="none", edgecolors=color, s=12)
        ax.set(title=setting, xlim=(0, 1), ylim=(-.2, .7))
    axes[0, 0].legend(fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$I/(e\Delta/\hbar)$")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\phi/\pi$")
    fig.suptitle("Bobok et al. (2025), Fig. 8: lines QP/GAL; markers published NRG")
    fig.tight_layout()
    save_figure(fig, args.output/"comparison.svg")

    rows = read_csv(args.output/"shared.csv")
    states = read_json(args.output/"eigenstates.json")
    for crossing in (r for r in crossings if int(r["figure"]) == 17):
        u, setting = float(crossing["critical_parameter"]), crossing["setting"]
        candidates = [s for s in states if s["label"]["figure"] == 17 and s["label"]["setting"] == setting]
        nearest = min(candidates, key=lambda s: abs(s["label"]["u"]-u))["label"]["u"]
        candidates = [s for s in candidates if s["label"]["u"] == nearest]
        # Scalar branch limits come from the archived spin-resolved states
        # at the root. A ground-state expectation is nonunique at coexistence.
        for spin in ((.5, 0.) if crossing["transition"] == "doublet-singlet" else (0., 1.)):
            options = [(s["energies"][i], s, i) for s in candidates for i, value in enumerate(s["observables"]["total_spin_squared"])
                       if abs(value-spin*(spin+1)) < 1e-6]
            _, state, index = min(options, key=lambda item: item[0])
            rows.append(dict(setting=setting, u=u, pairing=state["observables"]["pairing"][index],
                             spin_correlation=state["observables"]["spin_correlation"][index]))
    spectrum = read_csv(args.output/"spectrum.csv")
    reference17 = read_csv(CASE/"reference"/"figure17.csv")
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    for setting, color in (("ChE-W2", "C1"), ("ChE-W4", "C0"), ("ChE-W6", "green"), ("ChE-W8", "C3"), ("SM8", "C4")):
        selected = sorted((r for r in rows if r["setting"] == setting), key=lambda r: float(r["u"]))
        for rank in range(1, 6):
            x, y = [], []
            for u in sorted({float(r["u"]) for r in spectrum if r["setting"] == setting}):
                energies = sorted({round(float(r["excitation"]), 7) for r in spectrum if r["setting"] == setting and float(r["u"]) == u})
                if len(energies) > rank:
                    x.append(u)
                    y.append(energies[rank])
                else:
                    x.append(u)
                    y.append(np.nan)
            axes[0].plot(x, y, color=color, lw=1.1, ls="--" if setting == "SM8" else "-", label=setting if rank == 1 else None)
        for ax, observable in zip(axes[1:], ("pairing", "spin_correlation"), strict=True):
            ax.plot([float(r["u"]) for r in selected], [float(r[observable]) for r in selected], color=color,
                    lw=1.2, ls="--" if setting == "SM8" else "-")
    for ax, observable, label in zip(axes, ("excitation", "pairing", "spin_correlation"),
                                     (r"$(E-E_g)/\Delta$", r"$\nu$", r"$\langle\mathbf{S}_1\cdot\mathbf{S}_2\rangle$"), strict=True):
        ref = [r for r in reference17 if r["method"] == "NRG" and r["observable"] == observable]
        ax.scatter([float(r["u"]) for r in ref], [float(r["value"]) for r in ref], color="black", s=8, label="NRG")
        ax.set_ylabel(label)
        ax.set_xlim(0, 10)
    axes[0].set_ylim(-.01, 1.02)
    axes[0].legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel(r"$U/\Delta$")
    fig.suptitle("Bobok et al. (2025), Fig. 17: coherent common-lead double dot")
    fig.tight_layout()
    save_figure(fig, args.output/"shared_lead.svg")


if __name__ == "__main__":
    main()
