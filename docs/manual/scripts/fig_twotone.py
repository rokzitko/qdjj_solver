"""Two-tone spectroscopy: Mollow sidebands of the pumped pair and the Autler-Townes doublet of a spectator."""

import numpy as np

import manual_common as mc

PUMP, HARMONICS, WINDOW = .03, 6, (.55, .9)


def solve_pair(model, omega, static_even, static_odd, size):
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import solve_floquet
    periodic = PeriodicHamiltonian.from_parameters(reference_model, model, ParameterDrive.cosine("detuning", PUMP),
                                                   omega, samples=16)
    options = dict(floquet_cutoff=HARMONICS, selection="all", method="dense", basis="eigen", basis_states=size,
                   dense_dimension=6000)
    even = solve_floquet(periodic, 2, Sector(0, 0), options, static=static_even)
    odd = solve_floquet(periodic, 2, Sector(1, 1), dict(floquet_cutoff=HARMONICS, states=1), static=static_odd)
    return even, odd


def gate_lines(even):
    lines = even.matrix_elements("charge").lines(initial=even.dressed_state(0))
    keep = (abs(lines.frequencies) > WINDOW[0]) & (abs(lines.frequencies) < WINDOW[1])
    return [(float(f), float(s), bool(e)) for f, s, e in zip(lines.frequencies[keep], lines.strengths[keep],
                                                             lines.elastic[keep], strict=True)]


def compute():
    from qdjj_solver import Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    model = mc.junction()
    h = reference_model(**model)
    static_even = solve(h, 2, Sector(0, 0), SolverOptions(method="dense", eigenpairs=120))
    static_odd = solve(h, 2, Sector(1, 1), SolverOptions(method="dense", eigenpairs=2))
    size = mc.basis_size(static_even.energies, 43)
    omega0 = float(static_even.energies[1] - static_even.energies[0])
    element = float(abs(static_even.matrix_elements("charge")[0, 1]))
    sweep = []
    for omega in np.linspace(omega0 - .06, omega0 + .06, 121):
        even, _ = solve_pair(model, float(omega), static_even, static_odd, size)
        sweep.append(dict(omega=float(omega), lines=gate_lines(even)))
    even, odd = solve_pair(model, omega0, static_even, static_odd, size)
    resonant = even.matrix_elements("charge").lines(initial=even.dressed_state(0))
    spectator = even.matrix_elements("d_up", odd).lines(initial=odd.dressed_state(0))
    undriven_d = abs(static_even.matrix_elements("d_up", static_odd)[:2, 0])**2
    current = even.matrix_elements("current").lines(initial=even.dressed_state(0))
    def pack(lines):
        return dict(frequencies=lines.frequencies.tolist(), strengths=lines.strengths.tolist(),
                    elastic=lines.elastic.tolist())
    mc.save("twotone", dict(pump=PUMP, harmonics=HARMONICS, basis_states=size, omega0=omega0, element=element,
                            rabi=PUMP*element, energies=dict(S=float(static_even.energies[0]),
                                                             Sstar=float(static_even.energies[1]),
                                                             D=float(static_odd.energies[0])),
                            sweep=sweep, resonant=pack(resonant), spectator=pack(spectator),
                            undriven_spectator=undriven_d.tolist(), current=pack(current)), __file__)


def plot():
    from qdjj_solver import SpectralLines
    record = mc.load("twotone")
    omega0, rabi = record["omega0"], record["rabi"]
    fig, axes = mc.figure(1, 3, height=2.45)
    axis = axes[0]
    x, y, weight = [], [], []
    for row in record["sweep"]:
        for frequency, strength, _ in row["lines"]:
            x.append(row["omega"] - omega0)
            y.append(abs(frequency))
            weight.append(strength)
    order = np.argsort(weight)
    weight = np.array(weight)[order]
    points = axis.scatter(1e3*np.array(x)[order], np.array(y)[order], c=weight/weight.max(), s=2.,
                          cmap=mc.weight_colormap(), vmin=0, vmax=1, lw=0, rasterized=True)
    detuning = np.linspace(-.06, .06, 200)
    width = np.hypot(detuning, rabi)
    for sign in (-1, 1):
        axis.plot(1e3*detuning, omega0 + detuning + sign*width, color=mc.ORANGE, lw=.7, ls="--")
    axis.plot(1e3*detuning, omega0 + detuning, color=mc.MUTED, lw=.7, ls=":")
    axis.set_xlabel(r"pump detuning $(\Omega-\omega_{SS^*})\times10^3$")
    axis.set_ylabel(r"probe frequency $|\omega|$")
    axis.set_ylim(.64, .8)
    axis.text(.03, .95, r"RWA: $\Omega\pm w$", transform=axis.transAxes, fontsize=7, color=mc.ORANGE, va="top")
    mc.panel_label(axis, "a")
    axis = axes[1]
    resonant = record["resonant"]
    lines = SpectralLines(np.array(resonant["frequencies"]), np.array(resonant["strengths"]), [], [],
                          np.zeros(len(resonant["frequencies"]), int), np.array(resonant["elastic"]))
    grid = np.linspace(omega0 - .05, omega0 + .05, 1001)
    rate = lines.spectrum(grid, .0015, kind="rate")
    net = lines.spectrum(grid, .0015, kind="net")
    axis.plot(grid, rate, color=mc.BLUE, label="rate")
    axis.plot(grid, net, color=mc.ORANGE, lw=1., label="absorption $-$ emission")
    for value in (omega0 - rabi, omega0, omega0 + rabi):
        axis.axvline(value, color=mc.GRID, lw=.7, zorder=0)
    axis.set_xlabel(r"probe frequency $\omega$")
    axis.set_ylabel("charge-probe spectrum")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "b")
    axis = axes[2]
    spectator = record["spectator"]
    lines = SpectralLines(np.array(spectator["frequencies"]), np.array(spectator["strengths"]), [], [],
                          np.zeros(len(spectator["frequencies"]), int), np.zeros(len(spectator["frequencies"]), bool))
    energies = record["energies"]
    center = energies["D"] - energies["S"]
    grid = np.linspace(-center - .04, -center + .04, 1001)
    axis.plot(-grid, lines.spectrum(grid, .0015), color=mc.BLUE, label="pumped")
    reference = record["undriven_spectator"][0]*mc_lorentz(-grid, center, .0015)
    axis.plot(-grid, reference, color=mc.MUTED, lw=.9, ls="--", label="no pump")
    axis.set_xlabel(r"removal energy $-\omega$")
    axis.set_ylabel(r"$d_\uparrow$ spectrum from $D$")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "c")
    del points
    mc.finish(fig, "twotone")


def mc_lorentz(x, center, width):
    from qdjj_solver import line_shape
    return line_shape(np.asarray(x) - center, width)


def lines_near(lines, frequency, tolerance):
    """Strengths and elastic flags of the lines with |omega| within tolerance of frequency."""
    frequencies = np.abs(np.array(lines["frequencies"]))
    keep = np.abs(frequencies - frequency) < tolerance
    return np.array(lines["strengths"])[keep], np.array(lines["elastic"], bool)[keep]


def quoted(record):
    omega0, rabi = record["omega0"], record["rabi"]
    strengths, elastic = lines_near(record["resonant"], omega0, rabi/2)
    lower, _ = lines_near(record["resonant"], omega0 - rabi, rabi/2)
    upper, _ = lines_near(record["resonant"], omega0 + rabi, rabi/2)
    # Current probe: the four lines near the pump (Omega, Omega -+ delta) share |<S|dH/dphi|S*>|^2 in the RWA;
    # the line at delta between the hybrids has strength [(I_S - I_S*)/2]^2.
    pump, _ = lines_near(record["current"], omega0, 1.5*rabi)
    current = record["current"]
    frequencies = np.array(current["frequencies"])
    hybrids = np.array(current["strengths"])[np.abs(frequencies - rabi) < rabi/2].sum()
    return {
        "pump": (record["pump"], ".2f"),
        "omega0": (omega0, ".3f"),
        "rabi": (rabi, ".4f"),
        "harmonics": (record["harmonics"], "d"),
        "basis": (record["basis_states"], "d"),
        "epsilon": (mc.JUNCTION["detuning"], ".1f"),
        "mollow_rwa": (record["element"]**2/4, ".3f"),
        "mollow_lower": (float(lower.sum()), ".3f"),
        "mollow_elastic": (float(strengths[elastic].mean()), ".3f"),
        "mollow_upper": (float(upper.sum()), ".3f"),
        "current_element": (float(np.sqrt(pump.sum())), ".3f"),
        "current_rabi": (float(hybrids), ".4f"),
    }


if __name__ == "__main__":
    mc.main(compute, plot)
