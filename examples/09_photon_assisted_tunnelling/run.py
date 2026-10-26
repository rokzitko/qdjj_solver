"""Tutorial 09: photon-assisted tunnelling into a gate-driven junction.

    python -B examples/09_photon_assisted_tunnelling/run.py                  # quick profile, a few seconds
    python -B examples/09_photon_assisted_tunnelling/run.py --profile full   # the complete calculation, about 10 s

A slow gate drive eps(t) = eps0 + A cos(Omega t) dresses the singlet ground state S and the
doublet D. A weak tunnelling probe d_up^dagger adds an electron (S -> D); its line splits into
photon sidebands at E_D - E_S + q Omega. The script computes (1) the static charges and the
tunnelling matrix element, (2) the Floquet sideband weights versus A/Omega, (3) the plain
Tien-Gordon weights J_q^2(dn A/Omega) and (4) their adiabatic generalization from static
calculations along the drive, including the asymmetry between q = +1 and q = -1. Results go to
runs/examples/09_photon_assisted_tunnelling/<profile>/ (results.json, and
photon_assisted_tunnelling.png if matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402
from numpy.polynomial import Chebyshev  # noqa: E402

from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The junction of the physics manual's photon-assisted tunnelling chapter: Delta = 1, two bath level pairs per lead.
JUNCTION = dict(u=1., gamma=.5, phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
PAIRS = 2
OMEGA = .05         # drive frequency, far below the gaps inside each parity sector
REPORT = .4         # drive amplitude A at which the sideband weights are reported (A/Omega = 8)
ORDERS = np.arange(-3, 4)
WIDTH = .004        # half width of the plotted lines
PROFILES = {
    # drive amplitudes A, Floquet cutoff M = ceil(A/Omega) + extra, Chebyshev nodes of the static scan,
    # time samples per period in the adiabatic formula
    "quick": dict(amplitudes=np.round(np.linspace(0., .4, 5), 12), extra=4, nodes=24, samples=128),
    "full": dict(amplitudes=np.round(np.linspace(0., .6, 25), 12), extra=10, nodes=40, samples=256),
}


def floquet_sidebands(model, amplitude, harmonics, statics):
    """Lines of d_up^dagger from the dressed singlet into the dressed doublet, and the largest residual."""
    drive = ParameterDrive.cosine("detuning", amplitude)
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model, drive, OMEGA, samples=16)
    # One tracked state per sector in the full QP Fock space (shift-invert): no static-basis truncation.
    options = dict(floquet_cutoff=harmonics, states=1, method="shift-invert")
    singlet = solve_floquet(periodic, 2, Sector(0, 0), options, static=statics[0])
    doublet = solve_floquet(periodic, 2, Sector(1, 1), options, static=statics[1])
    lines = doublet.matrix_elements("d_up_dagger", singlet).lines(initial=singlet.dressed_state(0))
    return lines, float(max(singlet.residuals.max(), doublet.residuals.max()))


def static_scan(model, low, high, nodes):
    """Chebyshev interpolants of Delta E = E_D - E_S and |M| = |<D|d_up^dagger|S>| on detunings [low, high].

    Also returns the largest phase change of M along the scan in the parallel-transport gauge, where the
    overlap of the ground states at neighbouring detunings is real and positive. The adiabatic formula
    needs M in this gauge; when its phase is constant, |M| can be used.
    """
    detunings = (low + high)/2 + (high - low)/2*np.polynomial.chebyshev.chebpts1(nodes)
    gaps, elements, previous = [], [], None
    gauge_s = gauge_d = 1.
    for value in detunings:
        h = reference_model(**{**model, "detuning": float(value)})
        singlet = solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=1))
        doublet = solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1))
        if previous is not None:      # rotate the new eigenvectors so that overlaps with the previous ones are real
            overlap_s = complex(previous[0].overlaps(singlet)[0, 0])
            overlap_d = complex(previous[1].overlaps(doublet)[0, 0])
            gauge_s *= np.conj(overlap_s)/abs(overlap_s)
            gauge_d *= np.conj(overlap_d)/abs(overlap_d)
        previous = (singlet, doublet)
        gaps.append(doublet.energies[0] - singlet.energies[0])
        elements.append(np.conj(gauge_d)*gauge_s*complex(doublet.matrix_elements("d_up_dagger", singlet)[0, 0]))
    elements = np.array(elements)
    phase = float(np.max(np.abs(np.angle(elements/elements[0]))))
    gap = Chebyshev.fit(detunings, gaps, nodes - 1, domain=[low, high])
    element = Chebyshev.fit(detunings, np.abs(elements), nodes - 1, domain=[low, high])
    return gap, element, phase


def adiabatic_sidebands(gap, element, amplitude, samples, detuning=JUNCTION["detuning"]):
    """Weights |F_q|^2 of F(t) = M(eps(t)) exp[i int_0^t (Delta E - mean Delta E) dt'] for q in ORDERS.

    ``gap`` and ``element`` are functions of the detuning. Also returns the mean of Delta E over the period.
    """
    theta = 2*np.pi*np.arange(samples)/samples                    # theta = Omega t over one period
    eps = detuning + amplitude*np.cos(theta)
    splitting, matrix = gap(eps), element(eps)
    # Integrate Delta E - mean in Fourier space: exp(ik theta) integrates to exp(ik theta)/(ik Omega).
    coefficients = np.fft.fft(splitting)/samples
    k = np.fft.fftfreq(samples, 1./samples)
    integral = np.zeros(samples, complex)
    integral[k != 0] = coefficients[k != 0]/(1j*k[k != 0]*OMEGA)
    phase = np.real(np.fft.ifft(integral)*samples)
    transform = np.fft.fft(matrix*np.exp(1j*phase))/samples       # F_q multiplies exp(+i q Omega t)
    return np.abs(transform[ORDERS % samples])**2, float(np.real(coefficients[0]))


def bessel(q, x, samples=256):
    """J_q(x) from its integral representation (1/2pi) int exp[i(x sin theta - q theta)] d theta."""
    theta = 2*np.pi*np.arange(samples)/samples
    return float(np.real(np.mean(np.exp(1j*(x*np.sin(theta) - q*theta)))))


def extremes(series, low, high):
    """Minimum and maximum of a Chebyshev series on [low, high] (end points and stationary points)."""
    roots = series.deriv().roots()
    points = [low, high] + [r.real for r in roots if abs(r.imag) < 1e-9 and low < r.real < high]
    values = series(np.array(points))
    return float(values.min()), float(values.max())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path,
                        help="output directory (default: runs/examples/09_photon_assisted_tunnelling/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)
    epsilon = JUNCTION["detuning"]

    # 1. Static junction: charges, the tunnelling line S -> D and the gaps inside each parity sector.
    model = dict(JUNCTION, bath=cosh_grid(pairs=PAIRS, delta=1., bandwidth=10.))
    h = reference_model(**model)
    statics = (solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2)),
               solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=2)))
    n_s = float(np.real(statics[0].observables["impurity_charge"][0]))
    n_d = float(np.real(statics[1].observables["impurity_charge"][0]))
    dn = n_d - n_s
    excitation = float(statics[1].energies[0] - statics[0].energies[0])
    element0 = float(abs(statics[1].matrix_elements("d_up_dagger", statics[0])[0, 0]))
    strength = element0**2
    gap_singlets = float(statics[0].energies[1] - statics[0].energies[0])
    gap_doublet = float(statics[1].energies[1] - statics[1].energies[0])
    print(f"n_S = {n_s:.6f}, n_D = {n_d:.6f}, E_D - E_S = {excitation:.6f}, |<D|d_up^+|S>|^2 = {strength:.6f}")
    print(f"excitation gaps: singlets {gap_singlets:.4f}, doublet {gap_doublet:.4f}; Omega = {OMEGA}")

    # 2. Static calculations along the drive, for the adiabatic (generalized Tien-Gordon) formula.
    reach = settings["amplitudes"].max()
    low, high = epsilon - reach, epsilon + reach
    gap, element, element_phase = static_scan(model, low, high, settings["nodes"])
    interpolation = max(abs(gap(epsilon) - excitation), abs(element(epsilon) - element0))
    # Hellmann-Feynman: dE/d eps = <n_d> - 1, so the slope of Delta E is the charge difference dn.
    slope = float(gap.deriv()(epsilon))
    linear_gap = lambda eps: excitation + dn*(eps - epsilon)        # noqa: E731
    constant_element = lambda eps: np.full_like(eps, element0)      # noqa: E731

    # 3. Floquet sidebands versus amplitude, compared with the adiabatic formula and plain Tien-Gordon.
    rows, deviation, residual = [], 0., 0.
    for amplitude in settings["amplitudes"]:
        ratio = amplitude/OMEGA
        harmonics = int(np.ceil(ratio - 1e-9)) + settings["extra"]
        lines, largest = floquet_sidebands(model, amplitude, harmonics, statics)
        residual = max(residual, largest)
        floquet = np.array([lines.strengths[lines.photons == q].sum() for q in ORDERS])/strength
        adiabatic, mean_gap = adiabatic_sidebands(gap, element, amplitude, settings["samples"])
        adiabatic /= strength
        tien_gordon = np.array([bessel(q, dn*ratio)**2 for q in ORDERS])
        deviation = max(deviation, float(np.max(abs(floquet - adiabatic))))
        rows.append(dict(amplitude=float(amplitude), ratio=float(ratio), harmonics=harmonics, floquet=floquet.tolist(),
                         adiabatic=adiabatic.tolist(), tien_gordon=tien_gordon.tolist(),
                         center=float(lines.frequencies[lines.photons == 0][0]), mean_gap=mean_gap,
                         total=float(lines.strengths.sum()/strength), lines=lines.record()))
        print(f"A/Omega = {ratio:4.1f}  M = {harmonics:2d}  q=-1,0,+1: Floquet "
              + " ".join(f"{w:.4f}" for w in floquet[2:5]) + "  adiabatic "
              + " ".join(f"{w:.4f}" for w in adiabatic[2:5]))

    # 4. At A = REPORT: which dependence produces the +-1 asymmetry?
    (row,) = [r for r in rows if r["amplitude"] == REPORT]
    amplitude = REPORT
    plus, zero, minus = (list(ORDERS).index(q) for q in (1, 0, -1))
    def asymmetry(weights):
        return float(weights[plus] - weights[minus])
    variants = {name: adiabatic_sidebands(g, m, amplitude, settings["samples"])[0]/strength
                for name, g, m in (("constant_element", gap, constant_element), ("linear_gap", linear_gap, element),
                                   ("both", linear_gap, constant_element))}
    gap_min, gap_max = extremes(gap, epsilon - amplitude, epsilon + amplitude)
    element_min, element_max = extremes(element, epsilon - amplitude, epsilon + amplitude)
    theta = 2*np.pi*np.arange(settings["samples"])/settings["samples"]
    mean_square = float(np.mean(element(epsilon + amplitude*np.cos(theta))**2)/strength)

    values = dict(
        n_S=n_s, n_D=n_d, dn=dn, excitation=excitation, static_strength=strength,
        gap_singlets=gap_singlets, gap_doublet=gap_doublet,
        hellmann_feynman=abs(slope - dn), interpolation=float(interpolation), element_phase=element_phase,
        gap_min=gap_min, gap_max=gap_max, element_min=element_min, element_max=element_max,
        weight_minus=row["floquet"][minus], weight_zero=row["floquet"][zero], weight_plus=row["floquet"][plus],
        asymmetry=asymmetry(row["floquet"]), asymmetry_adiabatic=asymmetry(row["adiabatic"]),
        asymmetry_constant_element=asymmetry(variants["constant_element"]),
        asymmetry_linear_gap=asymmetry(variants["linear_gap"]),
        tien_gordon_zero=row["tien_gordon"][zero], tien_gordon_one=row["tien_gordon"][plus],
        tien_gordon_check=float(np.max(abs(variants["both"] - np.array(row["tien_gordon"])))),
        line_center=row["center"], line_center_adiabatic=row["mean_gap"],
        total_weight=row["total"], total_weight_adiabatic=mean_square,
        max_deviation_adiabatic=deviation, residual=residual,
    )
    for key, value in values.items():
        print(f"{key:28s} {value:.6g}")

    # Smooth adiabatic and Tien-Gordon curves for the figure (cheap: they use the interpolants).
    fine = np.linspace(0., settings["amplitudes"].max()/OMEGA, 121)
    curves = dict(ratios=fine.tolist(),
                  adiabatic=[(adiabatic_sidebands(gap, element, r*OMEGA, settings["samples"])[0]/strength).tolist()
                             for r in fine],
                  constant_element=[asymmetry(adiabatic_sidebands(gap, constant_element, r*OMEGA,
                                                                  settings["samples"])[0]/strength) for r in fine],
                  linear_gap=[asymmetry(adiabatic_sidebands(linear_gap, element, r*OMEGA,
                                                            settings["samples"])[0]/strength) for r in fine],
                  tien_gordon=[[bessel(q, dn*r)**2 for q in ORDERS] for r in fine])
    record = dict(tutorial=HERE.name, profile=args.profile, values=values, omega=OMEGA, orders=ORDERS.tolist(),
                  rows=rows, curves=curves)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"photon_assisted_tunnelling.png")
    print(f"Saved {output}")


def plot(record, path):
    """Three panels: sideband spectra, sideband weights versus A/Omega, and the +-1 asymmetry."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    from qdjj_solver import line_shape
    values, rows, curves = record["values"], record["rows"], record["curves"]
    orders = list(record["orders"])
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.3), layout="constrained")

    center = values["excitation"]
    grid = np.linspace(center - .2, center + .2, 1601)
    picks = [r for r in rows if r["amplitude"] in (0., .1, .2, .4, .6)]
    for offset, row in enumerate(picks):
        lines = row["lines"]
        spectrum = sum(s*line_shape(grid - f, WIDTH) for f, s in zip(lines["frequencies"], lines["strengths"],
                                                                     strict=True))
        axes[0].plot(grid, np.pi*WIDTH*spectrum + .55*offset, color="#2a78d6", lw=1.)
        axes[0].text(grid[0], .55*offset + .12, rf"$A/\Omega={row['ratio']:.0f}$", fontsize=8, color="gray")
    axes[0].set(xlabel=r"addition energy $\omega$", ylabel=r"$d^\dagger_\uparrow$ spectrum from $S$ (offset)",
                yticks=[])

    ratios = [r["ratio"] for r in rows]
    for q, color, marker, style in ((0, "#2a78d6", "o", "-"), (1, "#eb6834", "s", "-"),
                                    (-1, "#eb6834", "^", "--"), (2, "#1baf7a", "D", "-")):
        index = orders.index(q)
        axes[1].plot(ratios, [r["floquet"][index] for r in rows], marker, ms=4, mfc="none", color=color,
                     label=f"$q={q}$")
        axes[1].plot(curves["ratios"], [w[index] for w in curves["adiabatic"]], style, color=color, lw=.9)
        if q >= 0:
            axes[1].plot(curves["ratios"], [w[index] for w in curves["tien_gordon"]], ":", color="gray", lw=.9)
    axes[1].set(xlabel=r"drive amplitude $A/\Omega$",
                ylabel=r"sideband weight / $\vert\langle D\vert d^\dagger_\uparrow\vert S\rangle\vert^2$")
    axes[1].legend(frameon=False, fontsize=8, title="markers: Floquet\nlines: adiabatic\ndotted: Tien–Gordon",
                   title_fontsize=7)

    plus, minus = orders.index(1), orders.index(-1)
    axes[2].axhline(0., color="gray", ls=":", lw=.9, label="plain Tien–Gordon")
    axes[2].plot(curves["ratios"], [w[plus] - w[minus] for w in curves["adiabatic"]], color="#2a78d6", lw=1.,
                 label="adiabatic")
    axes[2].plot(curves["ratios"], curves["constant_element"], "--", color="#eb6834", lw=1.,
                 label=r"adiabatic, constant $\mathcal{M}$")
    axes[2].plot(curves["ratios"], curves["linear_gap"], "-.", color="#1baf7a", lw=1.,
                 label=r"adiabatic, linear $\Delta E$")
    axes[2].plot(ratios, [r["floquet"][plus] - r["floquet"][minus] for r in rows], "o", ms=4, mfc="none",
                 color="#2a78d6", label="Floquet")
    axes[2].set(xlabel=r"drive amplitude $A/\Omega$", ylabel=r"asymmetry $(w_{+1}-w_{-1})/\vert M\vert^2$")
    axes[2].legend(frameon=False, fontsize=7.5)
    for axis, label in zip(axes, "abc", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=105)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
