"""Gate-driven junction: complementary selection rules and the two-photon resonance."""

import numpy as np

import manual_common as mc

# Drive amplitudes per drive type for the selection-rule scan (small: linear response).
AMPLITUDES = dict(gate=("detuning", .01), phase=("phi", .1))
PHASE_SCALE = 10  # panel (a) magnifies the phase-drive splitting and matrix element


def splitting(periodic, omega, static, harmonics=2):
    from qdjj_solver import Sector
    from qdjj_solver.common.floquet import quasienergy_distance
    from qdjj_solver.qp_solver import solve_floquet
    result = solve_floquet(mc.with_omega(periodic, omega), 2, Sector(0, 0), dict(floquet_cutoff=harmonics, states=2),
                           static=static)
    return float(quasienergy_distance(*result.quasienergies, omega))


def compute():
    from scipy.optimize import minimize_scalar
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import ProjectedOperator, SolverOptions, solve
    sector = Sector(0, 0)
    # (a) Splitting per unit amplitude at the one-photon resonance, versus the static detuning.
    rows = []
    for position in np.linspace(-.6, .6, 25):
        model = mc.junction(detuning=float(position))
        static = solve(reference_model(**model), 2, sector, SolverOptions(method="dense", eigenpairs=4))
        h = static.hamiltonian
        bra, ket = static.vectors[:, 0], static.vectors[:, 1]
        omega = float(static.energies[1] - static.energies[0])
        row = dict(position=float(position), omega=omega,
                   charge=float(abs(ProjectedOperator(h.observables["impurity_charge"], static.basis)
                                    .expectation(bra, ket))),
                   current=float(abs(ProjectedOperator(h.observables["phase_derivative"], static.basis)
                                     .expectation(bra, ket))))
        for label, (parameter, amplitude) in AMPLITUDES.items():
            periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                           ParameterDrive.cosine(parameter, amplitude), 1., samples=32)
            row[label] = splitting(periodic, omega, static)/amplitude
        rows.append(row)
    # (b) One- and two-photon splittings versus the gate amplitude, at the shifted resonances.
    model = mc.junction()
    static = solve(reference_model(**model), 2, sector, SolverOptions(method="dense", eigenpairs=4))
    omega0 = float(static.energies[1] - static.energies[0])
    photons = []
    for amplitude in np.geomspace(.01, .3, 12):
        periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                       ParameterDrive.cosine("detuning", float(amplitude)), 1.,
                                                       samples=32)
        row = dict(amplitude=float(amplitude))
        for order in (1, 2):
            guess = omega0/order
            found = minimize_scalar(lambda omega, periodic=periodic: splitting(periodic, omega, static, harmonics=4),
                                    bounds=(guess-.03/order, guess+.03/order), method="bounded",
                                    options=dict(xatol=1e-10))
            row[f"gap_{order}"] = float(found.fun)
            row[f"resonance_{order}"] = float(found.x)
        photons.append(row)
    mc.save("gate", dict(parameters={k: v for k, v in model.items() if k != "bath"}, bath=mc.BATH,
                         amplitudes={k: v[1] for k, v in AMPLITUDES.items()}, selection=rows, omega0=omega0,
                         photons=photons), __file__)


def plot():
    record = mc.load("gate")
    fig, axes = mc.figure(1, 2, height=2.5)
    axis = axes[0]
    rows = record["selection"]
    position = np.array([r["position"] for r in rows])
    axis.plot(position, [r["charge"] for r in rows], color=mc.BLUE, lw=1.1, label=r"gate: $|\langle S|n_d|S^*\rangle|$")
    axis.plot(position[::2], [r["gate"] for r in rows][::2], "o", ms=3.2, mfc="none", mew=.9, color=mc.BLUE)
    axis.plot(position, [PHASE_SCALE*r["current"] for r in rows], color=mc.ORANGE, lw=1.1,
              label=fr"phase: ${PHASE_SCALE}\,|\langle S|\partial_\phi H|S^*\rangle|$")
    axis.plot(position[::2], [PHASE_SCALE*r["phase"] for r in rows][::2], "s", ms=3.2, mfc="none", mew=.9,
              color=mc.ORANGE)
    axis.set_xlabel(r"detuning $ϵ_0/\Delta$")
    axis.set_ylabel("resonant splitting / amplitude")
    axis.set_ylim(0, 1.08)
    axis.legend(loc="upper center", ncols=1, title="lines: static matrix elements; markers: Sambe",
                title_fontsize=7)
    mc.panel_label(axis, "a")
    axis = axes[1]
    rows = record["photons"]
    amplitude = np.array([r["amplitude"] for r in rows])
    for order, color, marker in ((1, mc.BLUE, "o"), (2, mc.ORANGE, "s")):
        gap = np.array([r[f"gap_{order}"] for r in rows])
        axis.loglog(amplitude, gap, marker, ms=3.5, mfc="none", mew=.9, color=color,
                    label=r"one-photon, $\Omega\approx\omega_{SS^*}$" if order == 1
                    else r"two-photon, $\Omega\approx\omega_{SS^*}/2$")
        axis.loglog(amplitude, gap[0]*(amplitude/amplitude[0])**order, color=color, lw=.8, ls=":")
    axis.set_xlabel(r"gate amplitude $A_ϵ$")
    axis.set_ylabel("splitting at resonance")
    axis.legend(loc="lower right", title=r"dotted: $\propto A_ϵ$ and $\propto A_ϵ^2$", title_fontsize=7)
    mc.panel_label(axis, "b")
    mc.finish(fig, "gate")


def _scientific(value):
    """One significant digit in LaTeX, for example 8\\times10^{-5}."""
    exponent = int(np.floor(np.log10(value)))
    mantissa = round(value/10**exponent)
    if mantissa == 10:
        mantissa, exponent = 1, exponent + 1
    return f"10^{{{exponent}}}" if mantissa == 1 else f"{mantissa}\\times10^{{{exponent}}}"


def quoted(record):
    """Numbers of the gate-drive chapter, from the record and the plot scale."""
    rows = record["selection"]
    position = np.array([r["position"] for r in rows])
    gate = max(abs(r["gate"] - r["charge"]) for r in rows)
    phase = np.array([abs(r["phase"] - r["current"]) for r in rows])
    # Apart from the symmetric point, where the element vanishes, one pair of positions stands out.
    outlier = abs(position[np.argmax(np.where(position != 0, phase, 0))])
    regular = phase[(position != 0) & (abs(abs(position) - outlier) > 1e-9)]
    if phase[abs(abs(position) - outlier) < 1e-9].min() < 3*regular.max():
        raise ValueError("the text singles out one pair of gate positions for the phase drive")
    return dict(gate_amplitude=(record["amplitudes"]["gate"], "g"),
                phase_amplitude=(record["amplitudes"]["phase"], "g"), phase_scale=(PHASE_SCALE, "d"),
                gate_deviation=(_scientific(gate), None), phase_deviation=(_scientific(regular.max()), None),
                outlier_position=(outlier, "g"),
                outlier_deviation=(_scientific(phase[abs(abs(position) - outlier) < 1e-9].max()), None))


if __name__ == "__main__":
    mc.main(compute, plot)
