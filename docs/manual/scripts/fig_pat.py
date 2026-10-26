"""Photon-assisted tunnelling into a gate-driven junction: generalized Tien-Gordon sidebands."""

import numpy as np

import manual_common as mc

OMEGA, WIDTH, PAIRS = .05, .004, 2
AMPLITUDES = np.linspace(0., .6, 13)
RANGE_AMPLITUDE = .4    # amplitude at which the text quotes the ranges of Delta E and |M|


def floquet_lines(model, amplitude, statics, size):
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import solve_floquet
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                   ParameterDrive.cosine("detuning", float(amplitude)), OMEGA,
                                                   samples=16)
    harmonics = int(np.ceil(amplitude/OMEGA)) + 10
    even = solve_floquet(periodic, 2, Sector(0, 0), dict(floquet_cutoff=harmonics, states=1), static=statics[0])
    odd = solve_floquet(periodic, 2, Sector(1, 1),
                        dict(floquet_cutoff=harmonics, selection="all", method="dense", basis="eigen",
                             basis_states=size, dense_dimension=8000), static=statics[1])
    addition = odd.matrix_elements("d_up_dagger", even).lines(initial=even.dressed_state(0))
    removal = odd.matrix_elements("d_dn", even).lines(initial=even.dressed_state(0))
    ground = odd.dressed_state(0)
    def doublet(lines):
        keep = np.array([ground in final for final in lines.final])
        return dict(frequencies=lines.frequencies[keep].tolist(), strengths=lines.strengths[keep].tolist(),
                    orders=lines.photons[keep].tolist())
    return doublet(addition), doublet(removal), harmonics


def adiabatic_sidebands(model, amplitude, orders):
    """Generalized Tien-Gordon: Fourier series of M(eps(t)) exp[i int (dE - mean dE)]."""
    from qdjj_solver import Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    samples = 256
    theta = 2*np.pi*np.arange(samples)/samples
    grid = np.linspace(model["detuning"] - amplitude - 1e-9, model["detuning"] + amplitude + 1e-9, 41)
    gap, element = [], []
    for value in grid:
        h = reference_model(**{**model, "detuning": float(value)})
        even = solve(h, 2, Sector(0, 0), SolverOptions(eigenpairs=1))
        odd = solve(h, 2, Sector(1, 1), SolverOptions(eigenpairs=1))
        gap.append(odd.energies[0] - even.energies[0])
        element.append(abs(odd.matrix_elements("d_up_dagger", even)[0, 0]))
    detuning = model["detuning"] + amplitude*np.cos(theta)
    splitting = np.interp(detuning, grid, gap)
    matrix = np.interp(detuning, grid, element)

    def sidebands(splitting, matrix):
        # Phase accumulated relative to the mean, integrated over one period (theta = omega t).
        phase = np.cumsum(splitting - splitting.mean())*(2*np.pi/samples)/OMEGA
        phase -= phase[0]
        transform = np.fft.fft(matrix*np.exp(1j*phase))/samples
        return {int(q): float(abs(transform[q % samples])**2) for q in orders}

    # Which dependence makes the sidebands asymmetric: a constant M, or a linearized Delta E.
    linear = np.polyval(np.polyfit(grid, gap, 1), detuning)
    variants = dict(constant_element=sidebands(splitting, np.full_like(matrix, matrix.mean())),
                    linear_gap=sidebands(linear, matrix))
    ranges = dict(gap=[float(min(gap)), float(max(gap))], element=[float(min(element)), float(max(element))])
    return sidebands(splitting, matrix), variants, ranges


def compute():
    from qdjj_solver import Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    model = mc.junction(pairs=PAIRS)
    h = reference_model(**model)
    statics = (solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=2)),
               solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=109)))
    odd_gaps = np.diff(statics[1].energies)
    gaps = dict(singlets=float(statics[0].energies[1] - statics[0].energies[0]),
                doublet=float(statics[1].energies[1:][odd_gaps > 1e-7][0] - statics[1].energies[0]))
    size = mc.basis_size(statics[1].energies, 60)
    charge = (float(np.real(statics[0].observables["impurity_charge"][0])),
              float(np.real(statics[1].observables["impurity_charge"][0])))
    rows = []
    for amplitude in AMPLITUDES:
        addition, removal, harmonics = floquet_lines(model, amplitude, statics, size)
        adiabatic, variants, ranges = adiabatic_sidebands(model, float(amplitude), range(-3, 4))
        rows.append(dict(amplitude=float(amplitude), harmonics=harmonics, addition=addition, removal=removal,
                         adiabatic=adiabatic, variants=variants, ranges=ranges))
    mc.save("pat", dict(omega=OMEGA, basis_states=size, charges=dict(S=charge[0], D=charge[1]), gaps=gaps,
                        static_element=float(abs(statics[1].matrix_elements("d_up_dagger", statics[0])[0, 0])**2),
                        excitation=float(statics[1].energies[0] - statics[0].energies[0]), rows=rows), __file__)


def plot():
    from scipy.special import jv
    from qdjj_solver import line_shape
    record = mc.load("pat")
    rows = record["rows"]
    excitation = record["excitation"]
    fig, axes = mc.figure(1, 2, height=2.6)
    axis = axes[0]
    grid = np.linspace(excitation - .2, excitation + .2, 1601)
    picks = [0, 2, 4, 8, 12]
    for offset, index in enumerate(picks):
        row = rows[index]
        spectrum = sum(s*line_shape(grid - f, WIDTH) for f, s in zip(row["addition"]["frequencies"],
                                                                     row["addition"]["strengths"], strict=True))
        color = mc.RAMP[min(len(mc.RAMP)-1, 2 + offset)]
        axis.plot(grid, spectrum/ (1/(np.pi*WIDTH)) + 1.1*offset, color=color, lw=1.)
        axis.text(grid[-1], 1.1*offset + .15, fr"$A/\Omega={row['amplitude']/record['omega']:.0f}$",
                  ha="right", fontsize=7, color=mc.SECONDARY)
    axis.set_xlabel(r"addition energy $\omega$")
    axis.set_ylabel(r"$d_\uparrow^\dagger$ spectrum from $S$ (offset)")
    axis.set_yticks([])
    mc.panel_label(axis, "a")
    axis = axes[1]
    amplitude = np.array([r["amplitude"] for r in rows])/record["omega"]
    total = record["static_element"]
    shift = record["charges"]["D"] - record["charges"]["S"]
    fine = np.linspace(0, amplitude[-1], 300)
    for q, color, marker in ((0, mc.BLUE, "o"), (1, mc.ORANGE, "s"), (-1, mc.ORANGE, "^"), (2, mc.AQUA, "D")):
        weights = []
        for row in rows:
            match = [s for s, o in zip(row["addition"]["strengths"], row["addition"]["orders"], strict=True) if o == q]
            weights.append(sum(match))
        axis.plot(amplitude, np.array(weights)/total, marker, ms=3.2, mfc="none", mew=.9, color=color,
                  label=fr"$q={q}$")
        axis.plot(amplitude, [row["adiabatic"][str(q)]/total for row in rows], color=color, lw=.9,
                  ls="-" if q >= 0 else "--")
    for q in (0, 1, 2):
        axis.plot(fine, jv(q, shift*fine)**2, color=mc.MUTED, lw=.7, ls=":")
    axis.set_xlabel(r"drive amplitude $A/\Omega$")
    axis.set_ylabel(r"sideband weight / $|\langle D|d^\dagger_\uparrow|S\rangle|^2$")
    axis.set_xlim(0, amplitude[-1])
    axis.legend(loc="upper right", title="markers: Floquet; lines: adiabatic\ndotted: plain Tien–Gordon",
                title_fontsize=6.5)
    mc.panel_label(axis, "b")
    mc.finish(fig, "pat")


def quoted(record):
    charges = record["charges"]
    row = next(r for r in record["rows"] if abs(r["amplitude"] - RANGE_AMPLITUDE) < 1e-9)
    asymmetry = {name: weights["1"] - weights["-1"] for name, weights in row["variants"].items()}
    if np.sign(asymmetry["constant_element"]) == np.sign(asymmetry["linear_gap"]):
        raise ValueError("the text says that M(eps) and Delta E(eps) alone give opposite asymmetries")
    words = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
    return {
        "omega": (record["omega"], ".2f"),
        "pairs": (rf"\text{{{words[PAIRS]}}}", None),
        "u": (mc.JUNCTION["u"], ".0f"),
        "gamma": (mc.JUNCTION["gamma"], ".1f"),
        "epsilon": (mc.JUNCTION["detuning"], ".1f"),
        "phi": (mc.JUNCTION["phi"], ".1f"),
        "width": (WIDTH, ".3f"),
        "n_S": (charges["S"], ".3f"),
        "n_D": (charges["D"], ".3f"),
        "dn": (charges["D"] - charges["S"], ".2f"),
        "gap_singlets": (record["gaps"]["singlets"], ".2f"),
        "gap_doublet": (record["gaps"]["doublet"], ".2f"),
        "range_amplitude": (RANGE_AMPLITUDE, "g"),
        "gap_min": (row["ranges"]["gap"][0], ".2f"),
        "gap_max": (row["ranges"]["gap"][1], ".2f"),
        "element_min": (row["ranges"]["element"][0], ".2f"),
        "element_max": (row["ranges"]["element"][1], ".2f"),
    }


if __name__ == "__main__":
    mc.main(compute, plot)
