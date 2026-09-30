"""Plot full finite-bath, GAL and published NRG double-dot currents."""

from pathlib import Path

from applications._common import plot_arguments, read_csv, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    for ax, label in zip(axes, ("weak", "strong"), strict=True):
        for name, path, style in (("Full finite bath", args.output/"current.csv", "o"),
                                   ("GAL", args.output/"gal.csv", "-"),
                                   ("Published NRG", CASE/"reference"/"figure9a.csv", ".")):
            data = [r for r in read_csv(path) if r["curve"] == label]
            # Draw each parity branch independently; do not interpolate across a QPT.
            if name == "GAL":
                for parity in (0, 1):
                    branch = [r for r in data if int(r["ground_parity"]) == parity]
                    ax.plot([float(r["phi_over_pi"]) for r in branch],
                            [float(r["current"]) for r in branch], style,
                            color="C1", label=name if parity == 0 else None)
            else:
                ax.plot([float(r["phi_over_pi"]) for r in data],
                        [float(r["current"]) for r in data], style, ms=4,
                        color="C0" if name == "Full finite bath" else "black", label=name)
        ax.axhline(0, color="0.7", lw=.5)
        ax.set(xlabel=r"$\phi/\pi$", ylabel=r"$J/(e\Delta/\hbar)$",
               title=(r"$\Gamma_L=\Gamma_R=\Delta,\ t_d=0.1\Delta$" if label == "weak" else
                      r"$\Gamma_L=\Gamma_R=1.4\Delta,\ t_d=\Delta$"))
        ax.legend(fontsize=8)
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
