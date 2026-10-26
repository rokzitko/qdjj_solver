"""Floquet fan diagrams of the junction: replicas near the ground state versus frequency and amplitude."""

import numpy as np

import manual_common as mc

CUTOFF = 2            # QP cutoff
WINDOW = .22          # quasienergies within E_S +- WINDOW
HARMONICS = 4
DRIVES = dict(phase=dict(parameter="phi", sweep_amplitude=1., amplitudes=(0., 3.), omega=.6),
              gate=dict(parameter="detuning", sweep_amplitude=.05, amplitudes=(0., .3), omega=.6))


def replicas(result, ground, omega):
    """All replicas eps + l*omega inside the window, with their weight on S and on S*."""
    points = []
    for state in range(len(result.quasienergies)):
        weight_s, weight_star = result.static_weights[state, :2]
        for shift in range(-6, 7):
            value = result.quasienergies[state] + shift*omega - ground
            if abs(value) < WINDOW:
                # Plot precision: quasienergies to 1e-7, weights to 1e-4.
                points.append((round(float(value), 7), round(float(weight_s), 4), round(float(weight_star), 4)))
    return points


def solve_all(periodic, static):
    """All Floquet states in the complete cutoff-2 Fock space (AC Stark shifts need every state)."""
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    return solve_floquet(periodic, CUTOFF, Sector(0, 0),
                         dict(floquet_cutoff=HARMONICS, selection="all", method="dense", dense_dimension=8000),
                         static=static)


def compute():
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    # Two level pairs per lead keep the full Fock space (dimension 146) affordable for every point.
    model = mc.junction(pairs=2)
    static = solve(reference_model(**model), CUTOFF, Sector(0, 0), SolverOptions(method="dense", eigenpairs=146))
    ground = float(static.energies[0])
    record = dict(window=WINDOW, harmonics=HARMONICS, pairs=2, ground=ground,
                  levels=(static.energies - ground).tolist())
    for name, drive in DRIVES.items():
        periodic = PeriodicHamiltonian.from_parameters(
            reference_model, model, ParameterDrive.cosine(drive["parameter"], drive["sweep_amplitude"]), 1.,
            samples=48)
        frequency = []
        for omega in np.linspace(.25, 1.45, 321):
            result = solve_all(mc.with_omega(periodic, omega), static)
            frequency.append(dict(omega=float(omega), points=replicas(result, ground, omega)))
        amplitude = []
        for value in np.linspace(*drive["amplitudes"], 161):
            periodic = PeriodicHamiltonian.from_parameters(
                reference_model, model, ParameterDrive.cosine(drive["parameter"], float(value)), drive["omega"],
                samples=48)
            result = solve_all(periodic, static)
            amplitude.append(dict(amplitude=float(value), points=replicas(result, ground, drive["omega"])))
        record[name] = dict(drive=drive, frequency=frequency, amplitude=amplitude)
    record["decomposition"] = decomposition(model, static, ground)
    mc.save("spectrum", record, __file__, compact=True)


def decomposition(model, static, ground):
    """Shift of the dressed singlet under the phase drive, split into the time-averaged H_0 and the rest."""
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    drive = DRIVES["phase"]
    amplitude, omega = drive["sweep_amplitude"], drive["omega"]
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("phi", amplitude),
                                                   omega, samples=48)
    average = solve(periodic.time_average(), CUTOFF, Sector(0, 0), SolverOptions(method="dense", eigenpairs=1))
    result = solve_all(periodic, static)
    dressed = int(np.argmax(result.static_weights[:, 0]))
    shift = float(result.quasienergies[dressed] - ground)
    shift -= omega*round(shift/omega)
    # Adiabatic limit: the time average of the instantaneous singlet energy.
    phases = model["phi"] + amplitude*np.cos(2*np.pi*np.arange(64)/64)
    adiabatic = np.mean([solve(reference_model(**{**model, "phi": float(phi)}), CUTOFF, Sector(0, 0),
                               SolverOptions(method="dense", eigenpairs=1)).energies[0] for phi in phases])
    return dict(amplitude=amplitude, omega=omega, h0=float(average.energies[0] - ground), dressed=shift,
                weight=float(result.static_weights[dressed, 0]), adiabatic=float(adiabatic - ground))


def scatter(axis, rows, key):
    x, y, weight = [], [], []
    for row in rows:
        for value, weight_s, _ in row["points"]:
            x.append(row[key])
            y.append(value)
            weight.append(weight_s)
    order = np.argsort(weight)
    x, y, weight = np.array(x)[order], np.array(y)[order], np.array(weight)[order]
    return axis.scatter(x, 1e3*y, c=weight, s=1.2, cmap=mc.weight_colormap(), vmin=0, vmax=1, lw=0,
                        rasterized=True)


def plot():
    record = mc.load("spectrum")
    names = (("phase", r"phase drive, $A_\phi=1$", r"phase amplitude $A_\phi$ (rad)"),
             ("gate", r"gate drive, $A_ϵ=0.05$", r"gate amplitude $A_ϵ$"))
    fig, axes = mc.figure(2, 2, height=4.9)
    for row, (name, title, label) in enumerate(names):
        data = record[name]
        points = scatter(axes[row, 0], data["frequency"], "omega")
        axes[row, 0].set_title(title, loc="left", fontsize=8.5)
        axes[row, 0].set_xlabel(r"drive frequency $\Omega$")
        axes[row, 0].set_ylabel(r"$(\varepsilon+\ell\Omega-E_S)\times10^3$")
        scatter(axes[row, 1], data["amplitude"], "amplitude")
        axes[row, 1].set_title(fr"$\Omega={data['drive']['omega']}$", loc="left", fontsize=8.5)
        axes[row, 1].set_xlabel(label)
        for axis in axes[row]:
            axis.set_ylim(-1e3*record["window"], 1e3*record["window"])
            axis.grid(False)
        omega = [r["omega"] for r in data["frequency"]]
        axes[row, 0].set_xlim(omega[0], omega[-1])
        span = [r["amplitude"] for r in data["amplitude"]]
        axes[row, 1].set_xlim(span[0], span[-1])
        star = record["levels"][1]
        for photons in (1, 2):
            axes[row, 0].axvline(star/photons, color=mc.MUTED, lw=.6, ls=":", zorder=0)
    axes[0, 0].text(record["levels"][1], 1e3*record["window"]*.92, r" $S^*$, 1 photon", fontsize=6.5,
                    color=mc.SECONDARY, va="top")
    axes[0, 0].text(record["levels"][1]/2, 1e3*record["window"]*.92, r" 2 photons", fontsize=6.5,
                    color=mc.SECONDARY, va="top")
    bar = fig.colorbar(points, ax=axes, shrink=.5, aspect=25, pad=.01)
    bar.set_label(r"weight on the undriven singlet $S$")
    bar.outline.set_visible(False)
    for axis, label in zip(axes.ravel(), "abcd", strict=True):
        mc.panel_label(axis, label)
    mc.finish(fig, "spectrum")


def dressed_singlet(row):
    """(value, weight on S, weight on S*) of the replica with the largest weight on S."""
    return max(row["points"], key=lambda point: point[1])


def quoted(record):
    """Numbers of the fan-diagram discussion, from the record and the parameters above."""
    from scipy.special import j0
    phase, gate = record["phase"]["drive"], record["gate"]["drive"]
    star = record["levels"][1]
    # Away from resonances: the frequencies at which the dressed singlet keeps 95% of its weight on S.
    offset = np.median([value for value, weight, _ in map(dressed_singlet, record["phase"]["frequency"])
                        if weight >= .95])
    # The first avoided crossing versus the amplitude: the first dip of that weight below 0.9.
    rows = record["phase"]["amplitude"]
    weights = [dressed_singlet(row)[1] for row in rows]
    dip = next(i for i in range(1, len(rows)-1) if weights[i-1] > weights[i] <= weights[i+1] and weights[i] < .9)
    # The continuum threshold E_D+Delta-E_S of the same junction is computed in fig_numerics.
    threshold = mc.load("numerics")["continuum"][str(record["pairs"])]["threshold"]
    pairs = "zero one two three four five six seven eight nine ten".split()[record["pairs"]]
    parts = record["decomposition"]
    stark = parts["dressed"] - parts["h0"]
    if not (parts["h0"] > 0 > stark and abs(stark) > parts["h0"]/4):
        raise ValueError("the text says that H_0 raises the singlet and the AC Stark shifts lower it comparably")
    return dict(
        decomposition_omega=(parts["omega"], "g"),
        h0_shift=(1e3*parts["h0"], ".0f"), stark_shift=(-1e3*stark, ".0f"), dressed_shift=(1e3*parts["dressed"], ".0f"),
        adiabatic_shift=(1e3*parts["adiabatic"], ".0f"),
        cutoff=(CUTOFF, "d"), harmonics=(record["harmonics"], "d"), window=(record["window"], ".2f"),
        pairs=(rf"\text{{{pairs}}}", None), dimension=(len(record["levels"]), "d"),
        omega_star=(star, ".3f"), one_photon=(star, ".2f"), two_photon=(star/2, ".2f"), threshold=(threshold, ".1f"),
        phase_amplitude=(phase["sweep_amplitude"], "g"), phase_omega=(phase["omega"], "g"),
        gate_amplitude=(gate["sweep_amplitude"], "g"), gate_omega=(gate["omega"], "g"),
        # Gamma -> Gamma J_0^2(A/4) in the time-averaged Hamiltonian (eq:phase-h0 of the manual).
        hybridization_factor=(j0(phase["sweep_amplitude"]/4)**2, ".3f"),
        offset=(1e3*offset, ".0f"), crossing=(rows[dip]["amplitude"], ".1f"))


if __name__ == "__main__":
    mc.main(compute, plot)
