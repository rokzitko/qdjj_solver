"""Plot the reproduced coupled spectrum, dot spin and independent NRG curve."""

from pathlib import Path

from applications._common import plot_arguments, read_csv, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    rows = read_csv(args.output/"spectrum.csv")
    labels = list(dict.fromkeys(r["setting"] for r in rows))
    for i, label in enumerate(labels):
        data = [r for r in rows if r["setting"] == label]
        gamma = [float(r["gamma"]) for r in data]
        for key, style in (("singlet_excitation", "-"), ("doublet_excitation", "--")):
            axes[0].plot(gamma, [float(r[key]) for r in data], style, color=f"C{i}",
                         label=label if style == "-" else None, lw=1)
        axes[1].plot(gamma, [float(r["doublet_dot_spin"]) for r in data], color=f"C{i}", label=label)
    ref = read_csv(CASE/"reference"/"figure3.csv")
    for branch in ("S", "Dg"):
        data = [r for r in ref if r["method"] == "NRG" and r["branch"] == branch]
        axes[0].plot([float(r["gamma"]) for r in data], [float(r["excitation"]) for r in data],
                     "k", lw=1.5, ls="-" if branch == "S" else "--", label="NRG" if branch == "S" else None)
    axes[0].axhline(1, color="0.5", ls=":", lw=.8)
    axes[0].set(xscale="log", xlim=(.4, 20), ylim=(0, 2), xlabel=r"$\Gamma/\Delta$",
                ylabel=r"$(E-E_{\rm ground})/\Delta$", title="Fig. 3: S (solid), coupled Dg (dashed)")
    axes[1].set(xscale="log", xlim=(.4, 20), xlabel=r"$\Gamma/\Delta$",
                ylabel=r"$\langle S^z_{\rm dot}\rangle_{D_g}$", title="Screening in the coupled doublet")
    for ax in axes:
        ax.legend(fontsize=7, ncol=2)
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
