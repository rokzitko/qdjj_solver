"""Plot Fig. 9(b)/10(a) reproductions together with the deposited NRG data."""

from pathlib import Path

import numpy as np

from applications._common import plot_arguments, read_csv, save_figure

CASE = Path(__file__).resolve().parent


def main():
    args = plot_arguments(CASE)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), layout="constrained")
    rows = read_csv(args.output/"coupling.csv")
    ref = read_csv(CASE/"reference"/"figure9.csv")
    x = np.array([float(r["gamma_over_u"]) for r in rows])
    rx = np.array([float(r["gamma_over_u"]) for r in ref])
    mask = rx <= x.max()*1.01
    for key, label, color in (("kappa_phi0", r"$\phi=0$", "C0"),
                              ("kappa_phipi", r"$\phi=\pi$", "C1")):
        axes[0].plot(rx[mask], np.array([float(r[key]) for r in ref])[mask]/rx[mask],
                     color=color, alpha=.6, label=label+" NRG")
        axes[0].plot(x, [float(r[key])/float(r["gamma_over_u"]) for r in rows],
                     "o", ms=4, color=color, label=label+" QP")
    axes[0].plot(x, [float(r["leading_kappa"])/float(r["gamma_over_u"]) for r in rows],
                 "k--", label="Leading perturbation theory")
    axes[0].set(xscale="log", xlabel=r"$\Gamma/U$", ylabel=r"$\kappa/(\Gamma/U)$",
                title="Fig. 9(b): leading shift and phase correction")
    axes[0].legend(fontsize=7)
    for name, path, style in (("NRG", CASE/"reference"/"figure10.csv", "-"),
                               ("QP", args.output/"phase.csv", "o")):
        data = read_csv(path)
        p = [float(r["phi_over_pi"]) for r in data]
        k = np.array([float(r["kappa"]) for r in data])
        axes[1].plot(p, (k-k.min())/(k.max()-k.min()), style, ms=4, label=name)
    p = np.linspace(0, 1, 201)
    axes[1].plot(p, (1-np.cos(np.pi*p))/2, "k--", label="Leading cosine")
    axes[1].set(xlabel=r"$\phi/\pi$", ylabel="Normalized Knight shift",
                title=r"Fig. 10(a): $\Gamma/U=0.1$")
    axes[1].legend(fontsize=8)
    save_figure(fig, args.output/"comparison.svg")


if __name__ == "__main__":
    main()
