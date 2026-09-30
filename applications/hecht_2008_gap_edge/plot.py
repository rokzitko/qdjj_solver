"""Plot resolved quadratic spectra and the interacting resolution diagnostic."""

from applications._common import plot_arguments, read_csv, save_figure
from .run import CASE


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    rows = read_csv(args.output/"spectra.csv")
    ref = read_csv(CASE/"reference"/"spectra.csv")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for delta, color in ((.01, "C0"), (.0001, "C1"), (.000001, "C2")):
        data = [r for r in rows if float(r["delta_over_D"]) == delta and float(r["epsilon_over_gamma"]) == 0.]
        axes[0].loglog([float(r["omega_prime_over_D"]) for r in data], [float(r["pi_gamma_A"]) for r in data],
                       color=color, label=rf"$\Delta/D={delta:g}$")
    axes[0].set(xlim=(1e-16, .1), ylim=(.01, 1e4), title="Fig. 6(a): quadratic gap variation")
    for epsilon, color in zip((0., 2., 4., 8.), ("C0", "C1", "C2", "C3"), strict=True):
        data = [r for r in rows if float(r["delta_over_D"]) == .0001 and float(r["epsilon_over_gamma"]) == epsilon]
        axes[1].loglog([float(r["omega_prime_over_D"]) for r in data], [float(r["pi_gamma_A"]) for r in data],
                       color=color, label=rf"$\epsilon_d/\Gamma={epsilon:g}$")
        if epsilon in (0., 4.):
            nr = [r for r in ref if int(r["figure"]) == 6 and r["method"] == "NRG" and
                  abs(float(r["epsilon_over_D"])/.008-epsilon) < 1e-8]
            axes[1].loglog([float(r["omega_prime_over_D"]) for r in nr], [float(r["pi_gamma_A"]) for r in nr],
                           "--", color=color, alpha=.7)
    axes[1].set(xlim=(1e-13, 1e-3), ylim=(.01, 100), title="Fig. 6(b): quadratic gate variation\n(dashed: published NRG)")
    comparisons = read_csv(args.output/"comparison.csv")
    for delta, color in ((.03, "C0"), (.00003, "C1")):
        for width, style in ((.2, "-"), (.4, "--")):
            data = [r for r in comparisons if r["kind"] == f"interacting-L8_k240-b{width}" and float(r["delta_over_D"]) == delta]
            axes[2].loglog([float(r["omega_prime_over_D"]) for r in data],
                           [max(float(r["pi_gamma_A"]), 1e-20) for r in data], style, color=color,
                           label=rf"$\Delta/D={delta:g}$, $b={width}$")
        nr = [r for r in ref if int(r["figure"]) == 7 and float(r["delta_over_D"]) == delta]
        axes[2].loglog([float(r["omega_prime_over_D"]) for r in nr], [float(r["pi_gamma_A"]) for r in nr],
                       ".", color=color, ms=2)
    axes[2].set(xlim=(1e-9, 1.), ylim=(1e-3, 1e3), title="Fig. 7 diagnostic: L=8 finite bath\n(dots: published NRG; unconverged continuum)")
    for ax in axes:
        ax.set(xlabel=r"$(\omega-\Delta)/D$", ylabel=r"$\pi\Gamma A(\omega)$")
        ax.legend(fontsize=7)
    fig.tight_layout()
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
