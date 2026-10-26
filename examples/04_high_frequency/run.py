"""Tutorial 04: the high-frequency limit of a phase-driven junction.

    python -B examples/04_high_frequency/run.py                  # quick profile, about 13 s
    python -B examples/04_high_frequency/run.py --profile full   # the complete calculation, about 30 s

A phase drive phi(t) = phi0 + A cos(Omega t), much faster than every internal process,
renormalizes the hybridization to Gamma J0(A/4)^2. The script computes (1) the Floquet
quasienergies of the singlets S, S* and the doublet D versus A, compared with the static
junction with the renormalized coupling: the drive-induced singlet-doublet crossing and the
decoupling of the dot at the first zero of J0; (2) the van Vleck effective Hamiltonians of
orders 0, 1 and 2 against exact Sambe quasienergies; (3) the scaling of the van Vleck error
with Omega for a drive of the phase and the level with a relative phase. Results go to
runs/examples/04_high_frequency/<profile>/ (results.json, and high_frequency.png if
matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

from qdjj_solver import (  # noqa: E402
    Impurity, ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, make_model, reference_model)
from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The junction of the physics manual's high-frequency chapter, energies in units of Delta = 1:
# one bath level pair per lead and a small bandwidth, so that Omega = 40 exceeds the whole
# many-body spectrum of the QP-cutoff-2 space.
BANDWIDTH = 2.
JUNCTION = dict(u=1., gamma=.8, phi=.8, detuning=0., field=0., symmetry=False, compress=False)
OMEGA = 40.          # drive frequency
CUTOFF = 2           # QP cutoff of the junction calculations
HARMONICS = 4        # Floquet cutoff M of the amplitude scan
SAMPLES = 32         # samples per period; resolve the Bessel harmonics of e^{i A cos(Omega t)/4} up to A = 12
VV_AMPLITUDE = 6.    # amplitude of the van Vleck comparison
VV_HARMONICS = 6     # harmonics |k| <= 6 enter the van Vleck commutators; the rest are below J_7(A/4)
# Section 4: a dot on one lead, driven through the lead phase and the level with a relative phase.
SCALING_DRIVES = (ParameterDrive.cosine("phase", 1.), ParameterDrive.cosine("detuning", .5, phase=np.pi/2))
PROFILES = {
    "quick": dict(amplitudes=np.linspace(0., 12., 9), omegas=(20., 40., 80.)),
    "full": dict(amplitudes=np.linspace(0., 12., 49), omegas=tuple(10*2**(np.arange(9)/2))),
}
FIRST_ZERO = 2.404825557695773    # first zero of the Bessel function J0


def bessel_j0(z):
    """J0(z) as the time average of cos(z cos(Omega t)); exact to rounding for |z| < 20."""
    theta = 2*np.pi*np.arange(64)/64
    return float(np.mean(np.cos(z*np.cos(theta))))


def secant(function, a, b, tolerance=1e-10):
    """Root of a smooth function from two nearby starting points."""
    fa, fb = function(a), function(b)
    while abs(b - a) > tolerance and fb != fa:
        a, fa, b = b, fb, b - fb*(b - a)/(fb - fa)
        fb = function(b)
    return b


def gap(levels):
    """eps_D - eps_S; negative when the doublet is the lowest state."""
    return levels[2] - levels[0]


def junction(gamma=JUNCTION["gamma"]):
    return dict(JUNCTION, gamma=gamma, bath=cosh_grid(pairs=1, bandwidth=BANDWIDTH))


def driven(amplitude):
    """The phase-driven junction phi(t) = phi0 + A cos(Omega t)."""
    return PeriodicHamiltonian.from_parameters(reference_model, junction(), ParameterDrive.cosine("phi", amplitude),
                                               OMEGA, samples=SAMPLES)


def static_levels(hamiltonian, cutoff=CUTOFF):
    """Lowest two even (S, S*) and lowest odd (D) energies of a static Hamiltonian."""
    even = solve(hamiltonian, cutoff, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2)).energies
    odd = solve(hamiltonian, cutoff, Sector(1, 1), SolverOptions(method="dense", eigenpairs=1)).energies
    return np.array([even[0], even[1], odd[0]])


def floquet_levels(periodic, harmonics=HARMONICS):
    """Quasienergies of S, S* and D: the lowest Floquet states of each sector.

    Omega is far above the many-body spectrum, so the central replicas are not folded.
    selection="all" avoids tracking where S and S* become degenerate (at the zero of J0).
    """
    options = dict(floquet_cutoff=harmonics, selection="all", method="dense")
    even = np.sort(solve_floquet(periodic, CUTOFF, Sector(0, 0), options).quasienergies)
    odd = np.sort(solve_floquet(periodic, CUTOFF, Sector(1, 1), options).quasienergies)
    return np.array([even[0], even[1], odd[0]])


def renormalized(amplitude):
    """Static levels of the junction with Gamma -> Gamma J0(A/4)^2."""
    return static_levels(reference_model(**junction(JUNCTION["gamma"]*bessel_j0(amplitude/4)**2)))


def one_lead(phase, detuning):
    """A dot with the junction's U, Gamma and bath, coupled to a single lead of phase ``phase``."""
    bath = cosh_grid(pairs=1, bandwidth=BANDWIDTH)
    hopping = np.sqrt(JUNCTION["gamma"]/(np.pi*bath.rho))*np.eye(2)
    return make_model(Impurity.anderson(JUNCTION["u"], detuning), [bath], [hopping], phases=[phase])


def van_vleck_errors(omega):
    """Largest error of the van Vleck orders 0, 1, 2 for the ground states of both parities."""
    periodic = PeriodicHamiltonian.from_parameters(one_lead, dict(phase=0., detuning=JUNCTION["detuning"]),
                                                   SCALING_DRIVES, omega, samples=16)
    errors = np.zeros(3)
    for sector in (Sector(0, 0), Sector(1, 1)):
        # The full Fock space (cutoff None), so that the expansion and the Sambe matrix agree as omega -> oo.
        exact = solve_floquet(periodic, None, sector, dict(floquet_cutoff=6, states=1, method="dense")).quasienergies[0]
        for order in range(3):
            energy = solve(periodic.van_vleck(order), None, sector,
                           SolverOptions(method="dense", eigenpairs=1)).energies[0]
            errors[order] = max(errors[order], abs(energy - exact))
    return errors


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path,
                        help="output directory (default: runs/examples/04_high_frequency/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Quasienergies versus drive amplitude, against the junction with Gamma J0(A/4)^2.
    static = static_levels(reference_model(**junction()))
    print(f"undriven: E_S = {static[0]:.6f}, E_S* = {static[1]:.6f}, E_D = {static[2]:.6f}")
    amplitudes = settings["amplitudes"]
    floquet = np.array([floquet_levels(driven(a)) for a in amplitudes])
    effective = np.array([renormalized(a) for a in amplitudes])
    deviation = abs(floquet - effective)
    print(f"largest |quasienergy - renormalized junction| over the scan: {deviation.max():.3e}")

    # 2. The drive-induced singlet-doublet crossing: first sign change of eps_D - eps_S.
    j = int(np.flatnonzero(gap(floquet.T) < 0)[0])
    crossing = secant(lambda a: gap(floquet_levels(driven(a))), amplitudes[j-1], amplitudes[j])
    # In the renormalized junction it happens when Gamma J0(A/4)^2 reaches the static critical coupling.
    bracket = amplitudes[j-1]/4, amplitudes[j]/4
    critical = secant(lambda g: gap(static_levels(reference_model(**junction(g)))),
                      *(JUNCTION["gamma"]*bessel_j0(x)**2 for x in bracket))
    expected = 4*secant(lambda x: JUNCTION["gamma"]*bessel_j0(x)**2 - critical, *bracket)
    print(f"crossing at A = {crossing:.6f}; renormalized junction: Gamma_c = {critical:.6f}, A = {expected:.6f}")

    # 3. Van Vleck Hamiltonians: static Hamiltonians in the model's coordinates, solved like any other.
    #    At the first zero of J0 the dot decouples, and the atomic levels U/2 (S, S*) and 0 (D) return.
    decoupling = 4*FIRST_ZERO
    vv = {}
    for name, amplitude in (("A6", VV_AMPLITUDE), ("decoupled", decoupling)):
        periodic = driven(amplitude)
        exact = floquet_levels(periodic, harmonics=6)
        orders = [static_levels(periodic.van_vleck(order, max_harmonic=VV_HARMONICS)) for order in range(3)]
        vv[name] = dict(amplitude=float(amplitude), exact=exact.tolist(), orders=[o.tolist() for o in orders],
                        errors=[float(abs(o - exact).max()) for o in orders],
                        first_order=float(abs(orders[1] - orders[0]).max()),
                        renormalized=float(abs(orders[0] - renormalized(amplitude)).max()))
        print(f"A = {amplitude:.3f}: van Vleck errors (orders 0, 1, 2) = {np.round(vv[name]['errors'], 9)}")
    decoupled = np.array(vv["decoupled"]["exact"])
    print(f"at A = {decoupling:.4f}: eps_S, eps_S*, eps_D = {np.round(decoupled, 6)}")

    # 4. Error scaling with Omega, for a drive with H_{-k} != H_k (finite first-order term).
    omegas = np.array(settings["omegas"])
    scaling = np.array([van_vleck_errors(omega) for omega in omegas])
    fit = omegas >= 20.     # exponents of power-law fits; below Omega = 20 higher orders still matter
    slopes = [float(np.polyfit(np.log(omegas[fit]), np.log(scaling[fit, order]), 1)[0]) for order in range(3)]
    at40 = scaling[np.flatnonzero(omegas == 40.)[0]]
    print(f"van Vleck error exponents (expected -1, -2, -3): {np.round(slopes, 3)}")

    values = dict(
        energy_S=static[0], energy_Sstar=static[1], energy_D=static[2],
        max_renormalized_deviation=float(deviation.max()),
        crossing_amplitude=float(crossing), critical_gamma=float(critical), crossing_renormalized=float(expected),
        decoupled_S=decoupled[0], decoupled_Sstar=decoupled[1], decoupled_D=decoupled[2],
        vv0_error=vv["A6"]["errors"][0], vv2_error=vv["A6"]["errors"][2],
        vv2_error_decoupled=vv["decoupled"]["errors"][2],
        vv_first_order=vv["A6"]["first_order"], vv0_renormalized=vv["A6"]["renormalized"],
        scaling_error0=at40[0], scaling_error1=at40[1], scaling_error2=at40[2],
        scaling_slope0=slopes[0], scaling_slope1=slopes[1], scaling_slope2=slopes[2],
    )
    values = {key: float(value) for key, value in values.items()}
    for key, value in values.items():
        print(f"{key:28s} {value:.6g}")
    fine = np.linspace(0., 12., 121)
    record = dict(tutorial=HERE.name, profile=args.profile, values=values,
                  scan=dict(amplitudes=amplitudes.tolist(), floquet=floquet.tolist(), effective=effective.tolist()),
                  renormalized=dict(amplitudes=fine.tolist(), levels=[renormalized(a).tolist() for a in fine]),
                  van_vleck=vv, scaling=dict(omegas=omegas.tolist(), errors=scaling.tolist()))
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"high_frequency.png")
    print(f"Saved {output}")


def plot(record, path):
    """Quasienergies versus A, deviation from the renormalized junction, van Vleck error versus Omega."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "0.55"
    values = record["values"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), layout="constrained")
    axis = axes[0]
    fine = np.array(record["renormalized"]["amplitudes"])
    levels = np.array(record["renormalized"]["levels"])
    a = np.array(record["scan"]["amplitudes"])
    floquet = np.array(record["scan"]["floquet"])
    for column, color, style, name in ((0, blue, "-", "$S$"), (1, blue, "--", "$S^*$"), (2, orange, "-", "$D$")):
        axis.plot(fine, levels[:, column], color=color, ls=style, lw=1, label=name)
        axis.plot(a, floquet[:, column], "o" if column < 2 else "s", ms=3.5, mfc="none", mew=.9, color=color)
    axis.axvline(4*FIRST_ZERO, color=gray, lw=.8, ls=":")
    axis.axvline(values["crossing_amplitude"], color=gray, lw=.8)
    axis.text(4*FIRST_ZERO - .2, .1, "$J_0(A/4)=0$", ha="right", fontsize=8, color=gray)
    axis.text(values["crossing_amplitude"] + .2, -.6, "$D$ lowest", fontsize=8, color=gray)
    axis.set(xlabel="phase-drive amplitude $A$", ylabel="quasienergy", xlim=(0, 12))
    axis.legend(frameon=False, fontsize=8, loc="lower right", title="lines: $\\Gamma J_0(A/4)^2$\nmarkers: Floquet",
                title_fontsize=7)
    axis = axes[1]
    deviation = abs(floquet - np.array(record["scan"]["effective"]))
    for column, color, marker, name in ((0, blue, "o", "$S$"), (1, aqua, "^", "$S^*$"), (2, orange, "s", "$D$")):
        axis.semilogy(a[1:], deviation[1:, column], marker=marker, ms=3, mfc="none", mew=.8, lw=.8,
                      color=color, label=f"order 0, {name}")
    for row in record["van_vleck"].values():
        axis.semilogy([row["amplitude"]], [row["errors"][2]], "*", ms=8, color="k")
    axis.semilogy([], [], "*", ms=8, color="k", label="order 2, max")
    axis.set(xlabel="phase-drive amplitude $A$", ylabel="|van Vleck $-$ Sambe|", xlim=(0, 12), ylim=(1e-8, 3e-3))
    axis.legend(frameon=False, fontsize=7, loc="lower right", ncols=2)
    axis = axes[2]
    omegas = np.array(record["scaling"]["omegas"])
    errors = np.array(record["scaling"]["errors"])
    for order, color, marker in ((0, blue, "o"), (1, orange, "s"), (2, aqua, "D")):
        axis.loglog(omegas, errors[:, order], marker=marker, ms=4, mfc="none", mew=.9, color=color,
                    label=f"order {order}")
        axis.loglog(omegas, errors[-1, order]*(omegas/omegas[-1])**(-(order+1)), color=gray, lw=.8, ls=":")
    axis.set(xlabel=r"drive frequency $\Omega$", ylabel="largest error of the van Vleck energy")
    axis.legend(frameon=False, fontsize=8, title=r"dotted: $\Omega^{-(n+1)}$", title_fontsize=7)
    for axis, label in zip(axes, "abc", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=105)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
