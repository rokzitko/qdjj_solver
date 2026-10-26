"""Tutorial 03: the proximitized dot as a driven two-level system.

    python -B examples/03_two_level_system/run.py                  # quick profile, about 4 s
    python -B examples/03_two_level_system/run.py --profile full   # the complete calculation, about 30 s

In the large-gap limit the even sector of the dot is a two-level system: the empty and the doubly
occupied dot, coupled by the induced pairing Gamma_phi. The script drives its gate,
eps(t) = eps0 + A cos(Omega t), and compares the exact Sambe quasienergies with simple theories:
(1) the Autler-Townes avoided crossing versus Omega against the rotating-wave model (rwa_model),
(2) the Rabi splitting and the Bloch-Siegert shift of the resonance, and (3) the k-photon
resonances of a strongly driven, weakly coupled dot, with splittings 2 Gamma_phi |J_k(2A/Omega)|
and their collapse at a zero of J_k. Results go to runs/examples/03_two_level_system/<profile>/
(results.json, and two_level_system.png if matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402
from scipy.optimize import brentq, minimize_scalar  # noqa: E402
from scipy.special import jn_zeros, jv  # noqa: E402

from qdjj_solver import Hamiltonian, Impurity, PeriodicHamiltonian, Sector, create, hermitian_pair  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, rwa_model, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

U = 1.
WEAK = dict(epsilon=0., pairing=.5, amplitude=.1)     # Rabi regime: omega_0 = 2 Gamma_phi = 1
STRONG = dict(pairing=.05, omega=1.)                  # multiphoton regime: Gamma_phi << Omega <~ A
STRONGER = .5                                         # an amplitude comparable with omega_0
TWO_PHOTON = 1.                                       # amplitude of the reported two-photon resonance (2A/Omega = 2)
J2_ZERO = float(jn_zeros(2, 1)[0])/2                  # A/Omega at the first zero of J_2(2A/Omega)
# Floquet cutoff M: the modes extend over the modulation index A/Omega (up to 6) around a centre that the
# detuning can displace by a few photons; M = 32 keeps the residual bound near machine precision throughout.
HARMONICS = 32
PROFILES = {
    # Omega of the Autler-Townes sweep, amplitudes of the resonance search, amplitudes A/Omega of the
    # multiphoton curves, detunings eps0/Omega of the scan at the zero of J_2
    "quick": dict(sweep=np.linspace(.85, 1.15, 31), amplitudes=np.array([.1, .3, .5]),
                  strong=np.array([.5, 1., 2., 3.]), biases=np.linspace(0., 2.2, 45)),
    "full": dict(sweep=np.linspace(.85, 1.15, 121), amplitudes=np.linspace(.05, .5, 10),
                 strong=np.linspace(.1, 6., 60), biases=np.linspace(0., 2.2, 441)),
}


def proximitized(detuning, pairing):
    """Large-gap dot: U/2 (n-1)^2 + eps (n-1) - Gamma_phi (d_up^+ d_dn^+ + h.c.)."""
    impurity = Impurity.anderson(U, detuning)
    operator = impurity.operator + hermitian_pair(-pairing*create(0)*create(1))
    return Hamiltonian(operator, 2, (1, -1), observables=impurity.observables)


def driven(eps0, pairing, amplitude, omega):
    """Gate drive eps(t) = eps0 + amplitude cos(omega t); eight samples are exact for harmonics 0, +-1."""
    return PeriodicHamiltonian.from_callable(lambda t: proximitized(eps0 + amplitude*np.cos(omega*t), pairing),
                                             omega, samples=8, assume_fixed_coordinates=True)


RESIDUALS = []      # residual bounds of every Sambe calculation, the evidence for the Floquet cutoff


def floquet(eps0, pairing, amplitude, omega, **options):
    """Both Floquet states of the even sector (|0> and |2>)."""
    result = solve_floquet(driven(eps0, pairing, amplitude, omega), 0, Sector(0, 0),
                           dict(floquet_cutoff=HARMONICS, states=2, **options))
    RESIDUALS.append(float(np.max(result.residuals)))
    return result


def zone_distance(a, b, omega):
    """Distance between two quasienergies modulo omega: the splitting of an avoided crossing."""
    d = np.mod(a - b, omega)
    return min(d, omega - d)


def splitting(eps0, pairing, amplitude, omega):
    first, second = floquet(eps0, pairing, amplitude, omega).quasienergies
    return float(zone_distance(first, second, omega))


def minimum(function, center, half_width):
    """Location and value of the smallest splitting within center +- half_width."""
    found = minimize_scalar(function, bounds=(center - half_width, center + half_width), method="bounded",
                            options=dict(xatol=1e-10))
    return float(found.x), float(found.fun)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path,
                        help="output directory (default: runs/examples/03_two_level_system/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)
    epsilon, pairing, amplitude = WEAK["epsilon"], WEAK["pairing"], WEAK["amplitude"]

    # 1. The undriven two-level system: E_{S,S*} = U/2 -+ sqrt(eps^2 + Gamma_phi^2).
    static = solve(proximitized(epsilon, pairing), 0, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2))
    singlet, excited = static.energies
    omega0 = float(excited - singlet)
    element = float(abs(static.matrix_elements("charge")[0, 1]))    # |<S|n_d|S*>| = sin(theta) = 1 here
    print(f"E_S = {singlet:.6f}, omega_0 = {omega0:.6f}, |<S|n_d|S*>| = {element:.6f}")

    # 2. Autler-Townes avoided crossing in the rotating frame: eps_- and eps_+ - Omega, folded around E_S.
    #    The rotating-wave model keeps only |S, 0> and |S*, -1> (one photon absorbed).
    sweep, rwa_deviation = [], 0.
    for omega in settings["sweep"]:
        exact = floquet(epsilon, pairing, amplitude, float(omega), zone_center=singlet)
        rwa = rwa_model(driven(epsilon, pairing, amplitude, float(omega)), states=[0, 1], photons=[0, -1], cutoff=0,
                        sector=Sector(0, 0), options=dict(zone_center=singlet))
        pair, model = np.sort(exact.folded), np.sort(rwa.folded)
        rwa_deviation = max(rwa_deviation, float(np.max(abs(pair - model))))
        sweep.append(dict(omega=float(omega), exact=pair.tolist(), rwa=model.tolist()))
    at_resonance = rwa_model(driven(epsilon, pairing, amplitude, omega0), [0, 1], [0, -1], 0, Sector(0, 0))
    rabi_rwa = float(np.ptp(at_resonance.folded))
    bare = splitting(epsilon, pairing, amplitude, omega0)

    # 3. The true resonance moves to omega_0 + A^2/(4 omega_0): the Bloch-Siegert shift. Two definitions agree to
    #    this order: the smallest quasienergy splitting, and equal mixing of S and S* in the time-averaged dressed
    #    states, where the time-averaged probability of a transition from S to S* equals 1/2.
    def mixing(omega, a):
        result = floquet(epsilon, pairing, a, omega, zone_center=singlet)
        lower = int(np.argmin(result.folded))            # mostly S below the resonance, mostly S* above it
        return result.static_weights[lower, 0] - .5

    resonances = []
    for value in sorted(set(settings["amplitudes"].tolist()) | {amplitude, STRONGER}):
        location, gap = minimum(lambda omega, a=value: splitting(epsilon, pairing, a, omega), omega0 + .045, .075)
        mixed = brentq(mixing, omega0 - .03, omega0 + .12, args=(value,), xtol=1e-14)
        resonances.append(dict(amplitude=value, resonance=location, splitting=gap, mixing=mixed))
    weak = next(r for r in resonances if r["amplitude"] == amplitude)
    stronger = next(r for r in resonances if r["amplitude"] == STRONGER)
    bloch_siegert = amplitude**2/(4*omega0)

    # 4. Strong drive, weak pairing: the k-photon resonance near 2 eps0 = k Omega has the splitting
    #    2 Gamma_phi |J_k(2A/Omega)|. The search over eps0 absorbs the small O(Gamma_phi^2/Omega) level shifts.
    weak_pairing, omega = STRONG["pairing"], STRONG["omega"]

    def photon_resonance(k, a):
        return minimum(lambda eps0: splitting(eps0, weak_pairing, a, omega), k*omega/2, .02)

    curves, multiphoton = {}, 0.
    for k in (1, 2, 3):
        rows = []
        for ratio in sorted(set(settings["strong"].tolist()) | {TWO_PHOTON, J2_ZERO}):
            bias, gap = photon_resonance(k, ratio*omega)
            theory = 2*weak_pairing*abs(jv(k, 2*ratio))
            multiphoton = max(multiphoton, abs(gap - theory))
            rows.append(dict(ratio=ratio, bias=bias, splitting=gap, theory=float(theory)))
        curves[str(k)] = rows
    two_photon = next(r for r in curves["2"] if r["ratio"] == TWO_PHOTON)
    collapsed = next(r for r in curves["2"] if r["ratio"] == J2_ZERO)

    # 5. Quasienergy splitting versus detuning at the zero of J_2: the two-photon resonance has disappeared.
    biases = []
    for bias in settings["biases"]:
        biases.append(dict(bias=float(bias), splitting=splitting(float(bias)*omega, weak_pairing, J2_ZERO*omega, omega)))

    values = dict(
        omega0=omega0,
        element=element,
        rabi_rwa=rabi_rwa,
        splitting_bare_resonance=bare,
        rwa_deviation=rwa_deviation,
        bloch_siegert_shift=weak["resonance"] - omega0,
        bloch_siegert_ratio=(weak["resonance"] - omega0)/bloch_siegert,
        resonance_splitting=weak["splitting"],
        mixing_shift_ratio=(weak["mixing"] - omega0)/bloch_siegert,
        stronger_shift_ratio=(stronger["resonance"] - omega0)/(STRONGER**2/(4*omega0)),
        stronger_mixing_ratio=(stronger["mixing"] - omega0)/(STRONGER**2/(4*omega0)),
        stronger_splitting_ratio=stronger["splitting"]/STRONGER,
        two_photon_splitting=two_photon["splitting"],
        two_photon_ratio=two_photon["splitting"]/two_photon["theory"],
        collapsed_splitting=collapsed["splitting"],
        multiphoton_deviation=float(multiphoton),
        max_residual=max(RESIDUALS),
    )
    for key, value in values.items():
        print(f"{key:26s} {value:.10g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values,
                  parameters=dict(u=U, weak=WEAK, strong=STRONG, stronger=STRONGER, j2_zero=J2_ZERO),
                  sweep=sweep, resonances=resonances, multiphoton=curves, biases=biases)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"two_level_system.png")
    print(f"Saved {output}")


def plot(record, path):
    """Autler-Townes crossing, Bloch-Siegert shift, multiphoton splittings and the detuning scan at a zero of J_2."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    values = record["values"]
    omega0 = values["omega0"]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.), layout="constrained")
    axis = axes[0, 0]
    sweep = record["sweep"]
    omega = np.array([r["omega"] for r in sweep])
    axis.plot(omega, 0*omega, color=gray, lw=.8, ls=":")
    axis.plot(omega, omega0 - omega, color=gray, lw=.8, ls=":")
    exact, rwa = np.array([r["exact"] for r in sweep]), np.array([r["rwa"] for r in sweep])
    for branch in range(2):
        axis.plot(omega, exact[:, branch], color=blue, label="Sambe" if branch == 0 else None)
        axis.plot(omega, rwa[:, branch], color=orange, ls="--", lw=1.1, label="RWA" if branch == 0 else None)
    axis.set(xlim=(omega[0], omega[-1]), xlabel=r"drive frequency $\Omega$",
             ylabel=r"$\varepsilon_S$ and $\varepsilon_{S^*}-\Omega$")
    axis.legend(frameon=False, fontsize=8, loc="lower left")
    axis = axes[0, 1]
    rows = record["resonances"]
    amplitude = np.array([r["amplitude"] for r in rows])
    fine = np.linspace(0., amplitude.max(), 200)
    axis.plot(fine, fine**2/(4*omega0), color=orange, lw=1.1, label=r"Bloch–Siegert $A^2/4\omega_0$")
    axis.plot(amplitude, [r["resonance"] - omega0 for r in rows], "o", ms=4, mfc="none", mew=.9, color=blue,
              label="smallest splitting")
    axis.plot(amplitude, [r["mixing"] - omega0 for r in rows], "s", ms=3.5, mfc="none", mew=.9, color=aqua,
              label="equal mixing")
    axis.set(xlabel=r"drive amplitude $A$", ylabel=r"resonance shift $\Omega_{\mathrm{res}}-\omega_0$")
    axis.legend(frameon=False, fontsize=8, loc="upper left")
    axis = axes[1, 0]
    fine = np.linspace(0., 6., 400)
    pairing = record["parameters"]["strong"]["pairing"]
    for k, color in zip(("1", "2", "3"), (blue, orange, aqua), strict=True):
        rows = record["multiphoton"][k]
        axis.plot(fine, 2*pairing*abs(jv(int(k), 2*fine)), color=color, lw=1.1, label=fr"$k={k}$")
        axis.plot([r["ratio"] for r in rows], [r["splitting"] for r in rows], "o", ms=3.5, mfc="none", mew=.9,
                  color=color)
    axis.set(xlim=(0, 6), xlabel=r"drive amplitude $A/\Omega$", ylabel=r"splitting at the $k$-photon resonance")
    axis.legend(frameon=False, fontsize=8, loc="upper right", title=r"lines: $2\Gamma_\phi|J_k(2A/\Omega)|$",
                title_fontsize=7.5)
    axis = axes[1, 1]
    rows = record["biases"]
    axis.semilogy([r["bias"] for r in rows], [r["splitting"] for r in rows], color=blue, lw=1.1, label="Sambe")
    ratio = record["parameters"]["j2_zero"]
    for k in range(5):
        theory = 2*pairing*abs(jv(k, 2*ratio))
        axis.plot([k/2], [max(theory, 2e-5)], "v" if theory < 1e-4 else "o", ms=4.5, color=orange, ls="none",
                  label=r"$2\Gamma_\phi|J_k(2A/\Omega)|$" if k == 0 else None)
    axis.set(xlim=(0, 2.2), ylim=(1e-5, .6), xlabel=r"detuning $\epsilon_0/\Omega$",
             ylabel="quasienergy splitting")
    axis.text(1., 3e-5, r"$J_2=0$", ha="center", fontsize=8, color=gray)
    axis.legend(frameon=False, fontsize=8, loc="lower right")
    for axis, label in zip(axes.flat, "abcd", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=105)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
