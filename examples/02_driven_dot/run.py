"""Tutorial 02: the gate-driven isolated dot, Bessel sidebands and Floquet zones.

    python -B examples/02_driven_dot/run.py                  # quick profile, about 1.5 s
    python -B examples/02_driven_dot/run.py --profile full   # the complete calculation, about 3 s

A quantum dot without leads is driven by its gate, eps(t) = eps0 + A cos(Omega t). The problem is
solvable exactly, and the script compares the Sambe calculation with the Bessel-function solution:
(1) the quasienergies equal the static energies for every drive, (2) the harmonic weights of a state
of charge n are J_m((n-1)A/Omega)^2, or J_m(nA/Omega)^2 when the c-number harmonics are dropped,
(3) replicas eps + m Omega and their folding into the first Floquet zone, and (4) the Tien-Gordon
sidebands J_q(A/Omega)^2 of a tunnelling probe. Results go to runs/examples/02_driven_dot/<profile>/
(results.json, and driven_dot.png if matplotlib is installed).
"""

import argparse
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402
from scipy.special import jn_zeros, jv  # noqa: E402

from qdjj_solver import Hamiltonian, Impurity, PeriodicHamiltonian, Sector, create  # noqa: E402
from qdjj_solver.qp_solver import solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The dot of the physics manual: H_d = U/2 (n_d-1)^2 + eps (n_d-1), energies in arbitrary units.
U, EPSILON = 1., .2
FAN_AMPLITUDE = .6          # drive amplitude A of the replica fan
PROBE_OMEGA = .25           # drive frequency of the tunnelling-probe example
HARMONICS = 24              # Floquet cutoff M, well above the largest modulation index (8)
J0_ZERO = float(jn_zeros(0, 1)[0])      # 2.405: the central Tien-Gordon line vanishes
PROFILES = {
    # modulation indices A/Omega of the Bessel scan (Omega = 1), frequencies of the fan, probe indices
    "quick": dict(ratios=np.linspace(0., 4., 9), frequencies=np.linspace(.15, 1.6, 30),
                  probes=(1., J0_ZERO)),
    "full": dict(ratios=np.linspace(0., 4., 41), frequencies=np.linspace(.15, 1.6, 146),
                 probes=(.5, 1., 1.5, 2., J0_ZERO, 3.)),
}


def dot(detuning):
    """The isolated dot as a Hamiltonian: two spin orbitals (up, down), no reservoir."""
    impurity = Impurity.anderson(U, detuning)
    return Hamiltonian(impurity.operator, 2, (1, -1), observables=impurity.observables)


def driven(amplitude, omega):
    """eps(t) = EPSILON + amplitude cos(omega t); eight samples reproduce H(t) exactly (harmonics 0, +-1)."""
    # The bare dot declares no canonical coordinates, so they are accepted explicitly.
    return PeriodicHamiltonian.from_callable(lambda t: dot(EPSILON + amplitude*np.cos(omega*t)), omega, samples=8,
                                             assume_fixed_coordinates=True)


def floquet(periodic, sector, states, drop=False):
    return solve_floquet(periodic, 0, sector, dict(floquet_cutoff=HARMONICS, states=states, drop_scalar_harmonics=drop))


def charges(result):
    """Dot charge n of each Floquet state, from its static label."""
    return np.rint(np.real(result.static.observables["impurity_charge"][result.labels]))


def zone_distance(a, b, omega):
    """Distance between two quasienergies modulo omega."""
    d = np.mod(np.asarray(a) - np.asarray(b), omega)
    return np.minimum(d, omega - d)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/02_driven_dot/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Static atomic energies: E_|0> = U/2 - eps, E_|2> = U/2 + eps (even), E_|sigma> = 0 (odd).
    static = {0: U/2 - EPSILON, 1: 0., 2: U/2 + EPSILON}

    # 2. Harmonic weights versus the modulation index, with the c-number harmonics kept and dropped.
    #    Kept: q_n = n - 1 (charge from half filling); dropped: q_n = n.
    quasienergy_deviation, deviation, bessel = 0., dict(kept=0., dropped=0.), dict(kept=[], dropped=[])
    for ratio in settings["ratios"]:
        periodic = driven(float(ratio), 1.)
        for gauge, drop in (("kept", False), ("dropped", True)):
            result = floquet(periodic, Sector(0, 0), 2, drop)
            for state, n in enumerate(charges(result)):
                quasienergy_deviation = max(quasienergy_deviation, abs(result.quasienergies[state] - static[n]))
                z = ((n - 1) if gauge == "kept" else n)*ratio
                exact = jv(result.harmonic_indices, z)**2
                deviation[gauge] = max(deviation[gauge], float(np.max(abs(result.harmonic_weights[state] - exact))))
                if n == 2:
                    center = int(np.flatnonzero(result.harmonic_indices == 0)[0])
                    bessel[gauge].append(dict(ratio=float(ratio), index=float(z),
                                              weights=result.harmonic_weights[state, center:center + 3].tolist()))

    # The example of the physics manual: A = 0.6, Omega = 1.
    listing = floquet(driven(FAN_AMPLITUDE, 1.), Sector(0, 0), 2)
    m0 = int(np.flatnonzero(listing.harmonic_indices == 0)[0])                  # column of the harmonic m = 0
    dropped = floquet(driven(FAN_AMPLITUDE, 1.), Sector(0, 0), 2, drop=True)
    doubly = int(np.flatnonzero(charges(dropped) == 2)[0])
    print("quasienergies", listing.quasienergies, "m=0 weights", listing.harmonic_weights[:, m0])

    # 3. Replicas versus the drive frequency at fixed amplitude; folded quasienergies lie in |eps| < Omega/2.
    fan, folding = [], 0.
    for omega in settings["frequencies"]:
        result = floquet(driven(FAN_AMPLITUDE, float(omega)), Sector(0, 0), 2)
        n = charges(result)
        for state in range(2):
            quasienergy_deviation = max(quasienergy_deviation, abs(result.quasienergies[state] - static[n[state]]))
            folding = max(folding, float(zone_distance(result.folded[state], static[n[state]], omega)),
                          float(abs(result.folded[state]) - omega/2))
        center = int(np.flatnonzero(result.harmonic_indices == 0)[0])
        fan.append(dict(omega=float(omega), charge=n.tolist(), quasienergies=result.quasienergies.tolist(),
                        folded=result.folded.tolist(),
                        weights=result.harmonic_weights[:, center - 8:center + 9].tolist()))

    # 4. Tien-Gordon: add a spin-up electron to the doublet |down> (the ground state), |down> -> |2>.
    #    Lines at E_2 - E_down + q Omega with strengths J_q(A/Omega)^2, in both gauges.
    probe = create(0)            # d_up^dagger; mode 0 is spin up, mode 1 spin down
    addition = static[2] - static[1]
    strength_deviation, frequency_deviation, sum_rule, spectra = 0., 0., 0., []
    for ratio in settings["probes"]:
        periodic = driven(ratio*PROBE_OMEGA, PROBE_OMEGA)
        for drop in (False, True):
            even = floquet(periodic, Sector(0, 0), 2, drop)
            odd = floquet(periodic, Sector(1, -1), 1, drop)
            table = even.matrix_elements(probe, odd)                    # initial: odd, final: even
            lines = table.lines(initial=odd.dressed_state(0))
            photons = np.rint((lines.frequencies - addition)/PROBE_OMEGA)
            frequency_deviation = max(frequency_deviation, float(np.max(abs(
                lines.frequencies - addition - photons*PROBE_OMEGA))))
            found = dict(zip(photons.astype(int).tolist(), lines.strengths.tolist(), strict=True))
            orders = np.arange(-HARMONICS, HARMONICS + 1)
            strengths = np.array([found.get(int(q), 0.) for q in orders])
            strength_deviation = max(strength_deviation, float(np.max(abs(strengths - jv(orders, ratio)**2))))
            sum_rule = max(sum_rule, abs(float(lines.strengths.sum()) - 1.))     # |<2|d_up^dagger|down>|^2 = 1
            if not drop:
                spectra.append(dict(ratio=float(ratio), central=found.get(0, 0.), lines=lines.record()))
    unit = next(row for row in spectra if row["ratio"] == 1.)
    zero = next(row for row in spectra if row["ratio"] == J0_ZERO)

    values = dict(
        quasienergy_deviation=float(quasienergy_deviation),
        weight_m0=float(listing.harmonic_weights[0, m0]),
        weight_m1=float(listing.harmonic_weights[0, m0 + 1]),
        weight_m0_dropped=float(dropped.harmonic_weights[doubly, int(np.flatnonzero(dropped.harmonic_indices == 0)[0])]),
        bessel_deviation_kept=deviation["kept"],
        bessel_deviation_dropped=deviation["dropped"],
        folding_deviation=folding,
        sideband_0=unit["central"],
        sideband_1=next(s for f, s in zip(unit["lines"]["frequencies"], unit["lines"]["strengths"], strict=True)
                        if abs(f - addition - PROBE_OMEGA) < 1e-9),
        central_line_at_zero=zero["central"],
        sum_rule_deviation=sum_rule,
        sideband_deviation=strength_deviation,
        sideband_frequency_deviation=frequency_deviation,
    )
    for key, value in values.items():
        print(f"{key:28s} {value:.10g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values,
                  parameters=dict(u=U, epsilon=EPSILON, fan_amplitude=FAN_AMPLITUDE, probe_omega=PROBE_OMEGA,
                                  harmonics=HARMONICS),
                  bessel=bessel, fan=fan, spectra=spectra)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"driven_dot.png")
    print(f"Saved {output}")


def plot(record, path):
    """Bessel weights versus the modulation index, the replica fan and the Tien-Gordon sidebands."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    from matplotlib.collections import LineCollection
    from matplotlib.colors import to_rgba
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), layout="constrained")
    axis = axes[0]
    fine = np.linspace(0., 8., 400)
    for m, color in enumerate((blue, orange, aqua)):
        axis.plot(fine, jv(m, fine)**2, color=color, lw=1.2, label=fr"$J_{m}^2(z)$")
        for gauge, marker in (("kept", "o"), ("dropped", "x")):
            rows = record["bessel"][gauge]
            axis.plot([r["index"] for r in rows], [r["weights"][m] for r in rows], marker, ms=3.5, mfc="none",
                      mew=.9, color=color)
    axis.plot([], [], "o", mfc="none", color=gray, label=r"kept, $z=A/\Omega$")
    axis.plot([], [], "x", color=gray, label=r"dropped, $z=2A/\Omega$")
    axis.set(xlabel=r"modulation index $z$", ylabel=r"harmonic weight $\Vert u_m\Vert^2$ of $|2\rangle$",
             xlim=(0, 8), ylim=(0, 1.03))
    axis.legend(frameon=False, fontsize=7.5, loc="upper right")
    axis = axes[1]
    fan = record["fan"]
    omega = np.array([r["omega"] for r in fan])
    axis.fill_between(omega, -omega/2, omega/2, color="#f0efec", lw=0, zorder=0)
    for state in range(2):
        color = blue if fan[0]["charge"][state] == 2 else orange
        energy = np.array([r["quasienergies"][state] for r in fan])
        for m in range(-8, 9):
            weight = np.array([r["weights"][state][m + 8] for r in fan])
            points = np.column_stack([omega, energy + m*omega])
            segments = np.stack([points[:-1], points[1:]], axis=1)
            middle = np.clip((weight[:-1] + weight[1:])/2, 0., 1.)
            keep = middle > .02                     # opacity and width follow the harmonic weight
            colors = [to_rgba(color, .15 + .85*np.sqrt(w)) for w in middle[keep]]
            axis.add_collection(LineCollection(segments[keep], colors=colors, linewidths=.4 + 1.8*np.sqrt(middle[keep]),
                                               capstyle="round"))
    axis.plot([], [], color=blue, label=r"$|2\rangle$ and replicas")
    axis.plot([], [], color=orange, label=r"$|0\rangle$ and replicas")
    axis.text(1.58, -.7, r"first zone $|\varepsilon|<\Omega/2$", color=gray, fontsize=7.5, ha="right")
    axis.set(xlim=(omega[0], omega[-1]), ylim=(-.9, 1.1), xlabel=r"drive frequency $\Omega$",
             ylabel=r"$\varepsilon+m\Omega$")
    axis.legend(frameon=False, fontsize=7.5, loc="upper right")
    axis = axes[2]
    from qdjj_solver import line_shape
    width, addition = .01, U/2 + EPSILON          # Lorentzian half width; the undriven line E_2 - E_down
    energies = np.linspace(addition - 1.1, addition + 1.1, 2201)
    axis.axvline(addition, color=gray, lw=.8, ls=":")
    for row in record["spectra"]:
        color = {1.: blue, J0_ZERO: orange}.get(row["ratio"])
        if color is None:
            continue
        lines = row["lines"]
        spectrum = sum(s*line_shape(energies - f, width) for f, s in zip(lines["frequencies"], lines["strengths"],
                                                                        strict=True))
        axis.plot(energies, np.pi*width*spectrum, color=color, lw=1.1, label=fr"$A/\Omega={row['ratio']:.4g}$")
        q = np.arange(-4, 5)
        axis.plot(addition + q*PROBE_OMEGA, jv(q, row["ratio"])**2, "o", ms=3.5, mfc="none", mew=.9, color=gray)
    axis.plot([], [], "o", mfc="none", color=gray, label=r"$J_q^2(A/\Omega)$")
    axis.set(xlabel=r"addition energy $\omega$", ylabel=r"$d^\dagger_\uparrow$ line strength from $|\!\downarrow\rangle$",
             xlim=(energies[0], energies[-1]), ylim=(0, 1.05))
    axis.legend(frameon=False, fontsize=7.5, loc="upper right")
    for axis, label in zip(axes, "abc", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=110)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
