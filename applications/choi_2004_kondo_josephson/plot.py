"""Show the reproduction and the original currents on identical, untuned axes."""

from pathlib import Path

from applications._common import plot_arguments, read_csv, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
    rows = read_csv(args.output/"current.csv")
    references = read_csv(CASE/"reference"/"figure3.csv")
    for ax, ratios in zip(axes, ([10.], [1.6, 1.8, 2., 2.2], [.1]), strict=True):
        for i, ratio in enumerate(ratios):
            color = f"C{i}"
            data = [r for r in rows if float(r["delta_over_tk"]) == ratio]
            if not data:
                continue
            ref = [r for r in references if float(r["delta_over_tk"]) == ratio and float(r["phi_over_pi"]) >= 0]
            ax.plot([float(r["phi_over_pi"]) for r in ref], [float(r["current"]) for r in ref],
                    "--", color=color, lw=1, label=rf"NRG $\Delta/T_K={ratio:g}$")
            for parity in (0, 1):
                branch = [r for r in data if int(r["ground_parity"]) == parity]
                ax.plot([float(r["phi_over_pi"]) for r in branch], [float(r["current"]) for r in branch],
                        "-", color=color, lw=1.5, label=rf"QP $\Delta/T_K={ratio:g}$" if parity == int(data[0]["ground_parity"]) else None)
        ax.axhline(0, color="0.6", lw=.5)
        ax.set(xlabel=r"$\phi/\pi$", ylabel=r"$I/(e\Delta/\hbar)$", xlim=(0, 1))
        ax.legend(fontsize=6)
    for ax, title in zip(axes, (r"Weak coupling: $\pi$ junction", "Intermediate: phase-driven crossing", "Strong coupling: positive supercurrent"), strict=True):
        ax.set_title(title, fontsize=10)
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
