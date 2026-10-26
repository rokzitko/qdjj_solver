"""Tutorial 06: the phase-driven junction.

    python -B examples/06_phase_drive/run.py                  # quick profile, about 15 s
    python -B examples/06_phase_drive/run.py --profile full   # the complete calculation, about 3.5 minutes

A microwave phase drive phi(t) = phi0 + A cos(Omega t) of the reference junction couples the
singlet ground state S to the excited singlet S* through dH/dphi. The script computes
(1) the S-S* avoided crossing in the rotating frame, compared with the rotating-wave model,
(2) the splitting at the resonance versus the drive amplitude, (3) the eta selection rule at
the particle-hole-symmetric point, where the two replicas cross instead of repelling, and
(4) the Floquet spectrum versus the drive amplitude at a fixed frequency. Results go to
runs/examples/06_phase_drive/<profile>/ (results.json, and phase_drive.png if matplotlib is
installed).
"""

import argparse
import dataclasses
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model  # noqa: E402
from qdjj_solver.qp_solver import SolverOptions, rwa_model, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The junction of the physics manual: energies in units of the gap, Delta = 1, phases in radians.
JUNCTION = dict(u=1., gamma=.5, phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
SINGLETS = Sector(0, 0)     # even parity, S_z = 0: S, S* and their continuum
HARMONICS = 4               # Floquet cutoff M
SWEEP_AMPLITUDE = .4        # phase amplitude A (rad) of the frequency sweeps
FAN_OMEGA = .6              # drive frequency of the spectrum versus amplitude
WINDOW = .22                # plotted quasienergies within E_S +- WINDOW
PROFILES = {
    # Omega - omega_SS* of the sweeps, amplitudes of the resonance scan and of the Floquet spectrum
    "quick": dict(sweep=np.linspace(-.03, .03, 7), parity_sweep=np.linspace(-.03, .03, 5),
                  amplitudes=(1.,), fan=(1.,)),
    "full": dict(sweep=np.linspace(-.03, .03, 41), parity_sweep=np.linspace(-.03, .03, 25),
                 amplitudes=np.arange(1, 13)/10, fan=np.linspace(0., 3., 61)),
}


def junction(pairs, **changes):
    """Keyword arguments of reference_model: the manual's junction with ``pairs`` bath level pairs per lead."""
    return dict(JUNCTION, **changes, bath=cosh_grid(pairs=pairs, delta=1., bandwidth=10.))


def static_singlets(model, count):
    """The ``count`` lowest even-parity eigenstates of the undriven junction (dense diagonalization)."""
    return solve(reference_model(**model), 2, SINGLETS, SolverOptions(method="dense", eigenpairs=count))


def phase_drive(model, amplitude, omega=1., samples=16):
    """H(t) for phi(t) = phi0 + A cos(Omega t), sampled at ``samples`` times per period."""
    return PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("phi", amplitude),
                                               omega, samples=samples)


def at(periodic, omega):
    """The same drive at another frequency.

    The harmonics H_k of H(phi0 + A cos(Omega t)) depend on A but not on Omega, so a
    frequency scan builds them once and only replaces the frequency.
    """
    return dataclasses.replace(periodic, omega=float(omega))


def wrap(value, omega):
    """Map an energy difference into [-Omega/2, Omega/2): quasienergies are defined modulo Omega."""
    return (value + omega/2) % omega - omega/2


def rotating_pair(periodic, static):
    """Dressed S and S* in the rotating frame: (eps_S, eps_S* - Omega, static weights of both states)."""
    result = solve_floquet(periodic, 2, SINGLETS, dict(floquet_cutoff=HARMONICS, states=2, method="shift-invert"),
                           static=static)
    order = np.argsort(result.labels)                   # rows: the states labelled S and S*
    eps_s, eps_star = result.quasienergies[order]
    omega = periodic.omega
    return eps_s, eps_s + wrap(eps_star - omega - eps_s, omega), result.static_weights[order][:, :2]


def splitting(periodic, static):
    """Distance between the two rotating-frame quasienergies."""
    eps_s, eps_star, _ = rotating_pair(periodic, static)
    return abs(eps_star - eps_s)


def rwa(periodic, static):
    """Two-level rotating-wave model |S, 0>, |S*, -1>: (coupling 2|(H_1)_SS*|, resonance, eigenvalues)."""
    matrix = rwa_model(periodic, [0, 1], [0, -1], 2, SINGLETS, static=static).effective_matrix
    resonance = periodic.omega + (matrix[1, 1] - matrix[0, 0]).real   # where the two diagonal entries meet
    return 2*abs(matrix[0, 1]), float(resonance), np.linalg.eigvalsh(matrix)


def resonance(periodic, static, center, width, tolerance=1e-6):
    """Resonance frequency and splitting from a fit of the dressed-pair hyperbola.

    Near a resonance the two dressed states are split by w(Omega) = sqrt((Omega - Omega_r)^2 + delta^2),
    so w^2 is a parabola in Omega. Three points around ``center`` fix it; the fit is repeated around
    the new vertex until it no longer moves. The AC Stark shift moves Omega_r away from omega_SS*.
    """
    offsets = np.array([-width, 0., width])
    for _ in range(8):
        squares = [splitting(at(periodic, center + x), static)**2 for x in offsets]
        a, b, c = np.polyfit(offsets, squares, 2)
        shift = -b/(2*a)
        center += shift
        if abs(shift) < tolerance:
            break
    return float(center), float(np.sqrt(max(c - b*b/(4*a), 0.)))


def sweep(periodic, static, frequencies):
    """Rotating-frame quasienergies of the Sambe solution and of the RWA, frequency by frequency."""
    rows = []
    for omega in frequencies:
        eps_s, eps_star, weights = rotating_pair(at(periodic, omega), static)
        rows.append(dict(omega=float(omega), sambe=sorted([float(eps_s), float(eps_star)]),
                         rwa=rwa(at(periodic, omega), static)[2].tolist(),
                         mixing=float(max(weights[0, 1], weights[1, 0]))))
    return rows


def secant_crossing(periodic, static, a, b, steps=4):
    """Frequency at which eps_S* - Omega crosses eps_S (a true crossing: the difference changes sign)."""
    def signed(omega):
        eps_s, eps_star, _ = rotating_pair(at(periodic, omega), static)
        return eps_star - eps_s
    fa, fb = signed(a), signed(b)
    for _ in range(steps):
        a, b, fa = b, b - fb*(b - a)/(fb - fa), fb
        fb = signed(b)
    return float(b)


def adiabatic_average(model, amplitude, samples=32):
    """Time average of the instantaneous ground-state energy E_S(phi0 + A cos(Omega t)): the slow-drive limit."""
    phases = model["phi"] + amplitude*np.cos(2*np.pi*np.arange(samples)/samples)
    return float(np.mean([static_singlets(dict(model, phi=float(phi)), 1).energies[0] for phi in phases]))


def bessel_j0(x):
    """J_0(x) as the time average of cos(x sin t) (the trapezoidal rule is exact to rounding here)."""
    t = np.linspace(0., 2*np.pi, 64, endpoint=False)
    return float(np.mean(np.cos(x*np.sin(t))))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/06_phase_drive/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Static junction: the pair transition and the phase matrix element <S|dH/dphi|S*>.
    model = junction(pairs=4)
    static = static_singlets(model, 60)
    ground = float(static.energies[0])
    omega0 = float(static.energies[1] - static.energies[0])
    element = float(abs(static.matrix_elements("current")[0, 1]))
    print(f"omega_SS* = {omega0:.6f}, |<S|dH/dphi|S*>| = {element:.6f}")

    # 2. The avoided crossing at A = 0.4: Sambe in the full cutoff-2 Fock space and the two-level RWA.
    periodic = phase_drive(model, SWEEP_AMPLITUDE)
    coupling, rwa_resonance, _ = rwa(at(periodic, omega0), static)
    crossing = sweep(periodic, static, omega0 + settings["sweep"])
    resonance_omega, split = resonance(periodic, static, rwa_resonance, coupling)

    # 3. Splitting at the resonance versus amplitude: Sambe, RWA and the linear law A |<S|dH/dphi|S*>|.
    amplitudes = []
    for amplitude in settings["amplitudes"]:
        drive = phase_drive(model, float(amplitude))
        rwa_coupling, rwa_omega, _ = rwa(at(drive, omega0), static)
        omega_r, delta = resonance(drive, static, rwa_omega, rwa_coupling)
        amplitudes.append(dict(amplitude=float(amplitude), resonance=omega_r, splitting=delta, rwa=rwa_coupling,
                               linear=float(amplitude)*element))
    strong = min(amplitudes, key=lambda row: abs(row["amplitude"] - 1.))       # the row at A = 1

    # 4. Eta selection rule: at epsilon = 0 the phase drive cannot connect S and S*.
    symmetric = junction(pairs=4, detuning=0.)
    static0 = static_singlets(symmetric, 60)
    omega0_symmetric = float(static0.energies[1] - static0.energies[0])
    element_symmetric = float(abs(static0.matrix_elements("current")[0, 1]))
    charge_symmetric = float(abs(static0.matrix_elements("charge")[0, 1]))
    periodic0 = phase_drive(symmetric, SWEEP_AMPLITUDE)
    coupling_symmetric, guess, _ = rwa(at(periodic0, omega0_symmetric), static0)
    parity = sweep(periodic0, static0, omega0_symmetric + settings["parity_sweep"])
    cross = secant_crossing(periodic0, static0, guess - .005, guess + .005)

    # 5. Floquet spectrum versus amplitude at Omega = 0.6, in the complete Fock space of a smaller bath
    #    (two level pairs per lead, 146 states), so that all AC Stark shifts are included.
    small = junction(pairs=2)
    complete = static_singlets(small, 146)
    ground_small = float(complete.energies[0])
    fan = []
    for amplitude in settings["fan"]:
        drive = phase_drive(small, float(amplitude), FAN_OMEGA, samples=32)
        states = solve_floquet(drive, 2, SINGLETS, dict(floquet_cutoff=HARMONICS, selection="all", method="dense"),
                               static=complete)
        points = []
        for eps, weight in zip(states.quasienergies, states.static_weights[:, 0], strict=True):
            for shift in range(-6, 7):
                value = eps + shift*FAN_OMEGA - ground_small
                if abs(value) < WINDOW:
                    points.append([round(float(value), 8), round(float(weight), 6)])
        dressed = int(np.argmax(states.static_weights[:, 0]))      # the dressed singlet: most weight on S
        # The time-averaged Hamiltonian of a phase drive is the junction with Gamma -> Gamma J_0(A/4)^2.
        factor = bessel_j0(amplitude/4)**2
        averaged = static_singlets(junction(pairs=2, gamma=JUNCTION["gamma"]*factor), 1)
        time_average = solve(drive.time_average(), 2, SINGLETS, SolverOptions(method="dense", eigenpairs=1))
        fan.append(dict(amplitude=float(amplitude), points=points,
                        dressed=float(wrap(states.quasienergies[dressed] - ground_small, FAN_OMEGA)),
                        weight=float(states.static_weights[dressed, 0]), factor=factor,
                        averaged=float(averaged.energies[0] - ground_small),
                        time_average=float(time_average.energies[0] - ground_small),
                        adiabatic=adiabatic_average(small, float(amplitude)) - ground_small))
    reference = min(fan, key=lambda row: abs(row["amplitude"] - 1.))           # the spectrum at A = 1

    values = dict(
        omega0=omega0, element=element,
        rwa_splitting=coupling, linear_splitting=SWEEP_AMPLITUDE*element, rwa_resonance=rwa_resonance,
        resonance_omega=resonance_omega, resonance_shift=resonance_omega - omega0, splitting=split,
        strong_linear_ratio=strong["splitting"]/strong["linear"],
        strong_rwa_ratio=strong["rwa"]/strong["linear"],
        strong_resonance_shift=strong["resonance"] - omega0,
        symmetric_omega0=omega0_symmetric, symmetric_element=element_symmetric,
        symmetric_charge_element=charge_symmetric, symmetric_rwa_coupling=coupling_symmetric,
        symmetric_mixing=max(row["mixing"] for row in parity), symmetric_crossing=cross,
        fan_omega0=float(complete.energies[1] - complete.energies[0]), hybridization_factor=reference["factor"],
        dressed_shift=reference["dressed"], averaged_shift=reference["averaged"],
        averaged_check=abs(reference["averaged"] - reference["time_average"]),
        ac_stark_shift=reference["dressed"] - reference["averaged"], adiabatic_shift=reference["adiabatic"],
    )
    if len(fan) > 2:
        # The first avoided crossing of the dressed singlet versus amplitude: its weight on S dips below 0.9.
        weights = [row["weight"] for row in fan]
        dip = next(i for i in range(1, len(fan) - 1)
                   if weights[i - 1] > weights[i] <= weights[i + 1] and weights[i] < .9)
        values["fan_crossing"] = fan[dip]["amplitude"]
    for key, value in values.items():
        print(f"{key:26s} {value:.6g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values, ground=ground,
                  sweep_amplitude=SWEEP_AMPLITUDE, crossing=crossing, amplitudes=amplitudes,
                  parity=dict(omega0=omega0_symmetric, ground=float(static0.energies[0]), sweep=parity),
                  fan=dict(omega=FAN_OMEGA, window=WINDOW, rows=fan))
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"phase_drive.png")
    print(f"Saved {output}")


def plot(record, path):
    """Four panels: avoided crossing, splitting versus amplitude, eta selection rule, spectrum versus amplitude."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    values = record["values"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.6), layout="constrained")

    axis = axes[0, 0]
    rows = record["crossing"]
    omega = np.array([r["omega"] for r in rows])
    sambe = 1e3*(np.array([r["sambe"] for r in rows]) - record["ground"])
    approximation = 1e3*(np.array([r["rwa"] for r in rows]) - record["ground"])
    axis.plot(omega, sambe, color="#2a78d6", label=["Sambe", None])
    axis.plot(omega, approximation, "--", color="#eb6834", lw=1.1, label=["RWA", None])
    axis.axvline(values["omega0"], color="gray", lw=.6, ls=":")
    axis.set(xlabel=r"drive frequency $\Omega$", ylabel=r"$(\varepsilon-E_S)\times10^3$, rotating frame",
             title=fr"$\epsilon=0.2$, $A={record['sweep_amplitude']}$")
    axis.legend(frameon=False, fontsize=8, loc="lower left")

    axis = axes[0, 1]
    rows = record["amplitudes"]
    amplitude = np.array([r["amplitude"] for r in rows])
    axis.plot(amplitude, 1e3*np.array([r["splitting"] for r in rows]), "o", mfc="none", color="#2a78d6",
              label="Sambe, at the shifted resonance")
    axis.plot(amplitude, 1e3*np.array([r["rwa"] for r in rows]), "--", color="#eb6834", lw=1.1,
              label=r"RWA $2|(H_1)_{SS^*}|$")
    line = np.linspace(0, amplitude.max(), 50)
    axis.plot(line, 1e3*line*values["element"], color="#1baf7a", lw=1.1,
              label=r"$A\,|\langle S|\partial_\phi H|S^*\rangle|$")
    axis.set(xlabel=r"phase amplitude $A$ (rad)", ylabel=r"splitting at resonance $\times10^3$")
    axis.legend(frameon=False, fontsize=8, loc="upper left")

    axis = axes[1, 0]
    parity = record["parity"]
    omega = np.array([r["omega"] for r in parity["sweep"]])
    sambe = 1e3*(np.array([r["sambe"] for r in parity["sweep"]]) - parity["ground"])
    axis.plot(omega, sambe, color="#2a78d6")
    axis.axvline(values["symmetric_crossing"], color="gray", lw=.6, ls=":")
    axis.set(xlabel=r"drive frequency $\Omega$", ylabel=r"$(\varepsilon-E_S)\times10^3$, rotating frame",
             title=fr"$\epsilon=0$, $A={record['sweep_amplitude']}$: the replicas cross")

    axis = axes[1, 1]
    fan = record["fan"]
    x, y, weight = [], [], []
    for row in fan["rows"]:
        for value, w in row["points"]:
            x.append(row["amplitude"])
            y.append(value)
            weight.append(w)
    order = np.argsort(weight)
    points = axis.scatter(np.array(x)[order], 1e3*np.array(y)[order], c=np.array(weight)[order], cmap="Blues",
                          vmin=0, vmax=1, s=4, lw=0)
    amplitude = [row["amplitude"] for row in fan["rows"]]
    axis.plot(amplitude, [1e3*row["averaged"] for row in fan["rows"]], "--", color="gray", lw=1,
              label=r"$E_S$ of $H_0$: $\Gamma\to\Gamma J_0(A/4)^2$")
    axis.plot(amplitude, [1e3*row["adiabatic"] for row in fan["rows"]], ":", color="gray", lw=1.2,
              label=r"adiabatic: time average of $E_S(\phi(t))$")
    axis.set(xlabel=r"phase amplitude $A$ (rad)", ylabel=r"$(\varepsilon+\ell\Omega-E_S)\times10^3$",
             title=fr"$\Omega={fan['omega']}$, two bath level pairs", ylim=(-1e3*fan["window"], 1e3*fan["window"]))
    axis.legend(frameon=False, fontsize=8, loc="lower left")
    fig.colorbar(points, ax=axis, label=r"weight on $S$")
    for axis, label in zip(axes.ravel(), "abcd", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=100)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
