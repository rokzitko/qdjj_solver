"""Plot the unequal-gap current comparison from archived CSVs."""

import numpy as np

from applications._common import plot_arguments, read_csv, save_figure
from .run import CASE


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    rows = read_csv(args.output/"current.csv")
    reference = read_csv(CASE/"reference"/"figure8.csv")
    fig, ax = plt.subplots(figsize=(7, 4))
    for ratio, color in zip((1., .5, .25), ("C0", "C1", "C2"), strict=True):
        data = [r for r in rows if float(r["gap_ratio"]) == ratio]
        ax.plot([float(r["detuning"]) for r in data], [float(r["current"]) for r in data],
                "o-", color=color, ms=3, label=rf"QP, $\Delta_R/\Delta_L={ratio:g}$")
        for method in ("FDC_L", "DC_L"):
            ref = sorted((r for r in reference if float(r["gap_ratio"]) == ratio and
                          r["method"] == method), key=lambda r: float(r["detuning"]))
            ax.plot([float(r["detuning"]) for r in ref], [float(r["current"]) for r in ref],
                    "--", color=color, lw=1, label="published (F)DC" if ratio == 1 and method == "FDC_L" else None)
    nrg = [r for r in reference if r["method"] == "NRG"]
    ax.scatter([float(r["detuning"]) for r in nrg], [float(r["current"]) for r in nrg],
               facecolors="none", edgecolors="k", s=20, label="published NRG, ratio 0.5")
    ax.set(xlabel=r"$(\epsilon_d+U/2)/\Delta_L$", ylabel=r"$I/(e\Delta_L/\hbar)$",
           xticks=np.arange(-4, 5), title="Žonda et al. (2016), Fig. 8")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
