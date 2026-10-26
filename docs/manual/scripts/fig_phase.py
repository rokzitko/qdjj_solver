"""Phase-driven junction: S-S* resonance, Rabi splitting and the adiabatic limit."""

import numpy as np

import manual_common as mc

CUTOFF = 2                   # QP cutoff
ADIABATIC_OMEGAS = (.1, .3)  # drive frequencies of the adiabatic limit, panel (c)
FADED_WEIGHT = .5            # panel (c) fades the dressed state when its weight on S is lower


def static_states(model, sector, count):
    from qdjj_solver import reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    return solve(reference_model(**model), CUTOFF, sector, SolverOptions(method="dense", eigenpairs=count))


def phase_periodic(model, amplitude, samples=32):
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, reference_model
    return PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("phi", amplitude), 1.,
                                               samples=samples)


def rotating_pair(periodic, omega, static, harmonics=2):
    """Sorted (eps_S, eps_S* - omega) of the tracked pair in the full Fock space, and their static weights."""
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    result = solve_floquet(mc.with_omega(periodic, omega), CUTOFF, Sector(0, 0), dict(floquet_cutoff=harmonics, states=2),
                           static=static)
    order = np.argsort(result.labels)
    eps = result.quasienergies[order]
    return np.sort([eps[0], eps[1]-omega]), result.static_weights[order][:, :2]


def resonance(periodic, static, guess, width=.02):
    """Frequency of the smallest S-S* splitting near ``guess``; the bracket widens until the minimum is inside."""
    from scipy.optimize import minimize_scalar

    def gap(omega):
        # Away from the pair resonance the tracked states may hybridize with other replicas; such
        # frequencies are not the resonance, so they count as a large splitting.
        try:
            lower, upper = rotating_pair(periodic, omega, static)[0]
        except RuntimeError:
            return 1.
        return upper - lower
    low, high = guess - width, guess + width
    for _ in range(10):
        found = minimize_scalar(gap, bounds=(low, high), method="bounded", options=dict(xatol=1e-8))
        if found.x - low < 1e-4:
            low -= width
        elif high - found.x < 1e-4:
            high += width
        else:
            return float(found.x), float(found.fun)
    raise RuntimeError(f"no S-S* resonance found near {guess}")


def compute():
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import ProjectedOperator, rwa_model
    model = mc.junction()
    sector = Sector(0, 0)
    static = static_states(model, sector, 300)
    h = static.hamiltonian
    current = ProjectedOperator(h.observables["phase_derivative"], static.basis)
    element = float(abs(current.expectation(static.vectors[:, 0], static.vectors[:, 1])))
    omega0 = float(static.energies[1] - static.energies[0])
    ground = float(static.energies[0])
    # (a) Avoided crossing in the rotating frame. The RWA matrix is already in the rotating frame.
    sweep_amplitude = .4
    periodic = phase_periodic(model, sweep_amplitude)
    sweep = []
    for omega in np.linspace(omega0-.03, omega0+.03, 97):
        (lower, upper), weights = rotating_pair(periodic, omega, static)
        rwa = rwa_model(mc.with_omega(periodic, omega), [0, 1], [0, -1], CUTOFF, sector, static=static)
        sweep.append(dict(omega=float(omega), branches=[float(lower), float(upper)],
                          rwa=np.linalg.eigvalsh(rwa.effective_matrix).tolist(), weight=float(weights[0, 0])))
    # (b) Splitting at resonance versus amplitude.
    rows, guess = [], omega0
    # Up to A = 1 the driven pair stays a clean two-level system; beyond it the pair hybridizes with
    # replicas of other states and the splitting of the pair is no longer defined.
    for amplitude in np.linspace(.05, 1., 20):
        periodic = phase_periodic(model, float(amplitude))
        omega, gap = resonance(periodic, static, guess)
        guess = omega              # follow the AC-Stark-shifted resonance to larger amplitudes
        rwa = rwa_model(mc.with_omega(periodic, omega), [0, 1], [0, -1], CUTOFF, sector, static=static)
        # The RWA splitting at its own resonance is twice the coupling 2|(H_1)_{S S*}|.
        rows.append(dict(amplitude=float(amplitude), resonance=omega, gap=gap,
                         rwa=float(2*abs(rwa.effective_matrix[0, 1]))))
    # (c) Low frequencies approach the adiabatic average of the instantaneous energy. The AC Stark
    # shift collects virtual high-energy contributions, so use the complete Fock space of a
    # two-pair bath rather than a truncated static basis.
    model = mc.junction(pairs=2)
    static = static_states(model, sector, 146)
    adiabatic = []
    phases = np.linspace(-np.pi, np.pi, 721)
    from qdjj_solver import reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet
    curve = []
    for phi in np.linspace(model["phi"]-3.2, model["phi"]+3.2, 161):
        energy = solve(reference_model(**{**model, "phi": float(phi)}), CUTOFF, sector, SolverOptions(eigenpairs=1))
        curve.append(dict(phi=float(phi), energy=float(energy.energies[0])))
    grid = np.array([c["phi"] for c in curve])
    energies = np.array([c["energy"] for c in curve])
    amplitudes = np.linspace(0., 3., 25)
    for amplitude in amplitudes:
        average = float(np.mean(np.interp(model["phi"] + amplitude*np.cos(phases[:-1]), grid, energies)))
        row = dict(amplitude=float(amplitude), adiabatic=average)
        for omega in ADIABATIC_OMEGAS:
            periodic = phase_periodic(model, float(amplitude), samples=48)
            # At low frequency the zone is densely packed with replicas of all states: search widely,
            # and record no point where no Floquet state keeps 30% of S (the dressed singlet is gone).
            try:
                result = solve_floquet(mc.with_omega(periodic, omega), CUTOFF, sector,
                                       dict(floquet_cutoff=int((10 + 6*amplitude) if omega < .2 else (6 + 4*amplitude)),
                                            states=1, method="shift-invert", eigenpairs_per_target=32,
                                            max_eigenpairs_per_target=512, tracking_threshold=.3), static=static)
            except RuntimeError:
                row[str(omega)] = row[f"residual_{omega}"] = row[f"weight_{omega}"] = None
                continue
            row[str(omega)] = float(result.quasienergies[0])
            row[f"residual_{omega}"] = float(result.residuals[0])
            row[f"weight_{omega}"] = float(result.static_weights[0, 0])
        adiabatic.append(row)
    mc.save("phase", dict(parameters={k: v for k, v in model.items() if k != "bath"}, bath=mc.BATH,
                          omega0=omega0, element=element, sweep_amplitude=sweep_amplitude, sweep=sweep,
                          resonances=rows, ground=ground, adiabatic=adiabatic, energy_curve=curve,
                          adiabatic_pairs=2, adiabatic_ground=float(static.energies[0])), __file__)


def plot():
    record = mc.load("phase")
    fig, axes = mc.figure(1, 3, height=2.35)
    axis = axes[0]
    sweep = record["sweep"]
    omega = np.array([r["omega"] for r in sweep])
    branches = np.array([r["branches"] for r in sweep])
    rwa = np.array([r["rwa"] for r in sweep])
    ground = record["ground"]
    axis.plot(omega, 1e3*(branches[:, 0]-ground), color=mc.BLUE, label="Sambe")
    axis.plot(omega, 1e3*(branches[:, 1]-ground), color=mc.BLUE)
    axis.plot(omega, 1e3*(rwa[:, 0]-ground), color=mc.ORANGE, ls="--", lw=1.1, label="RWA")
    axis.plot(omega, 1e3*(rwa[:, 1]-ground), color=mc.ORANGE, ls="--", lw=1.1)
    axis.set_xlabel(r"drive frequency $\Omega$")
    axis.set_ylabel(r"$(\varepsilon-E_S)\times10^3$, rotating frame")
    axis.set_xlim(omega[0], omega[-1])
    axis.legend(loc="lower left")
    mc.panel_label(axis, "a")
    axis = axes[1]
    rows = record["resonances"]
    amplitude = np.array([r["amplitude"] for r in rows])
    axis.plot(amplitude, [1e3*r["gap"] for r in rows], "o", ms=3, mfc="none", mew=.9, color=mc.BLUE, label="Sambe")
    axis.plot(amplitude, [1e3*r["rwa"] for r in rows], ls="--", lw=1.1, color=mc.ORANGE, label="RWA")
    axis.plot(amplitude, 1e3*amplitude*record["element"], color=mc.AQUA, lw=1.1,
              label=r"$A|\langle S|\partial_\phi H|S^*\rangle|$")
    axis.set_xlabel(r"phase amplitude $A$ (rad)")
    axis.set_ylabel(r"splitting $\times10^3$")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "b")
    axis = axes[2]
    rows = record["adiabatic"]
    amplitude = np.array([r["amplitude"] for r in rows])
    ground = record["adiabatic_ground"]
    axis.plot(amplitude, [r["adiabatic"]-ground for r in rows], color=mc.INK, lw=1.1,
              label=r"$\langle E_S(\phi(t))\rangle$")
    # Faded markers: the dressed state keeps less than half of its weight on S.
    for omega, color, marker in zip(map(str, ADIABATIC_OMEGAS), (mc.BLUE, mc.ORANGE), "os", strict=True):
        found = np.array([r[omega] is not None for r in rows])
        values = np.array([r[omega]-ground if r[omega] is not None else np.nan for r in rows])
        kept = found & np.array([(r[f"weight_{omega}"] or 0.) >= FADED_WEIGHT for r in rows])
        faded = found & ~kept
        axis.plot(amplitude[kept], values[kept], marker, ms=3, mfc="none", mew=.9, color=color,
                  label=fr"$\Omega={omega}$")
        axis.plot(amplitude[faded], values[faded], marker, ms=3, mfc="none", mew=.6, color=color, alpha=.3)
    axis.set_xlabel(r"phase amplitude $A$ (rad)")
    axis.set_ylabel(r"$\varepsilon_S-E_S$")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "c")
    mc.finish(fig, "phase")


def quoted(record):
    """Numbers of the phase-drive chapter, from the record and the parameters above."""
    slow, fast = ADIABATIC_OMEGAS
    # Centre of the first run of faded markers at the faster drive (the resonance with S*).
    faded = []
    for row in record["adiabatic"]:
        weight = row[f"weight_{fast}"]
        if weight is not None and weight < FADED_WEIGHT:
            faded.append(row["amplitude"])
        elif faded:
            break
    pairs = "zero one two three four five six seven eight nine ten".split()[record["adiabatic_pairs"]]
    # The splitting of the pair falls smoothly with the amplitude until the pair starts to hybridize with
    # replicas of other states: the first decrement of the ratio that more than doubles, or a rise.
    rows = record["resonances"]
    ratios = np.array([row["gap"]/(row["amplitude"]*record["element"]) for row in rows])
    steps = -np.diff(ratios)
    onset = next(i + 1 for i in range(1, len(steps)) if steps[i] > 2*steps[i-1] or steps[i] < 0)
    smooth = rows[onset - 1]
    linear = smooth["amplitude"]*record["element"]
    if not smooth["gap"] < smooth["rwa"] < linear or smooth["amplitude"] < .5:
        raise ValueError("the text says that the RWA and the full splitting fall below the linear law")
    return dict(
        ratio_amplitude=(smooth["amplitude"], "g"), sambe_ratio=(smooth["gap"]/linear, ".3f"),
        rwa_ratio=(smooth["rwa"]/linear, ".3f"), resonance_shift=(record["omega0"] - smooth["resonance"], ".3f"),
        hybridization_amplitude=(rows[onset]["amplitude"], "g"),
        cutoff=(CUTOFF, "d"), detuning=(record["parameters"]["detuning"], "g"),
        omega_star=(record["omega0"], ".3f"), element=(record["element"], ".3g"),
        sweep_amplitude=(record["sweep_amplitude"], "g"), adiabatic_pairs=(rf"\text{{{pairs}}}", None),
        slow_omega=(slow, "g"), fast_omega=(fast, "g"), mixing_amplitude=((faded[0] + faded[-1])/2, ".1f"))


if __name__ == "__main__":
    mc.main(compute, plot)
