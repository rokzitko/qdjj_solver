"""ESR probe of the Rabi-driven doublet: Mollow triplet in S_x and the Rabi line in S_z."""

import numpy as np

import fig_esr
import manual_common as mc

GAMMA, AMPLITUDE, HARMONICS = .3, .004, 8
WIDTH = .006  # line half width of the plotted spectra, in units of the Zeeman splitting


def driven_doublet(omega, static, size):
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    periodic = mc.with_omega(fig_esr.periodic(GAMMA, AMPLITUDE), omega)
    return solve_floquet(periodic, 2, Sector(1, None),
                         dict(floquet_cutoff=HARMONICS, selection="all", method="dense", basis="eigen",
                              basis_states=size, dense_dimension=8000), static=static)


def compute():
    static = fig_esr.doublet(GAMMA)
    size = mc.basis_size(static.energies, 24)
    zeeman = float(static.energies[1] - static.energies[0])
    spin = float(abs(np.real(static.observables["impurity_spin_z"][0])))
    # Within the doublet: the polarizations of both Zeeman levels and the transverse element that drives them.
    polarizations = [float(np.real(value)) for value in static.observables["impurity_spin_z"][:2]]
    transverse = float(abs(static.matrix_elements("spin_x")[0, 1]))
    sweep = []
    for omega in zeeman*np.linspace(.8, 1.2, 81):
        states = driven_doublet(float(omega), static, size)
        lines = states.matrix_elements("spin_x").lines(initial=states.dressed_state(0))
        keep = abs(abs(lines.frequencies) - omega) < .5*zeeman
        sweep.append(dict(omega=float(omega), frequencies=lines.frequencies[keep].tolist(),
                          strengths=lines.strengths[keep].tolist()))
    states = driven_doublet(zeeman, static, size)
    def pack(name, window):
        lines = states.matrix_elements(name).lines(initial=states.dressed_state(0))
        keep = abs(lines.frequencies) < window
        return dict(frequencies=lines.frequencies[keep].tolist(), strengths=lines.strengths[keep].tolist())
    mc.save("esr_probe", dict(gamma=GAMMA, amplitude=AMPLITUDE, harmonics=HARMONICS, basis_states=size,
                              zeeman=zeeman, spin=spin, polarizations=polarizations, transverse=transverse,
                              rabi=AMPLITUDE*transverse, sweep=sweep,
                              spin_x=pack("spin_x", 2*zeeman), spin_z=pack("spin_z", .5*zeeman)), __file__)


def plot():
    from qdjj_solver import SpectralLines
    record = mc.load("esr_probe")
    zeeman, rabi = record["zeeman"], record["rabi"]
    fig, axes = mc.figure(1, 3, height=2.45)
    axis = axes[0]
    x, y, weight = [], [], []
    for row in record["sweep"]:
        for frequency, strength in zip(row["frequencies"], row["strengths"], strict=True):
            x.append(row["omega"]/zeeman)
            y.append(abs(frequency)/zeeman)
            weight.append(strength)
    order = np.argsort(weight)
    weight = np.array(weight)[order]
    axis.scatter(np.array(x)[order], np.array(y)[order], c=weight/weight.max(), s=3., cmap=mc.weight_colormap(),
                 vmin=0, vmax=1, lw=0, rasterized=True)
    detuning = np.linspace(-.2, .2, 200)*zeeman
    width = np.hypot(detuning, rabi)
    for sign in (-1, 1):
        axis.plot(1 + detuning/zeeman, 1 + (detuning + sign*width)/zeeman, color=mc.ORANGE, lw=.7, ls="--")
    axis.set_xlabel(r"drive frequency $\Omega/Z$")
    axis.set_ylabel(r"probe frequency $|\omega|/Z$")
    axis.set_ylim(.7, 1.3)
    mc.panel_label(axis, "a")
    def lines(name):
        data = record[name]
        count = len(data["frequencies"])
        return SpectralLines(np.array(data["frequencies"]), np.array(data["strengths"]), [], [],
                             np.zeros(count, int), np.zeros(count, bool))
    axis = axes[1]
    grid = np.linspace(.8, 1.2, 801)*zeeman
    axis.plot(grid/zeeman, lines("spin_x").spectrum(grid, WIDTH*zeeman, kind="rate"), color=mc.BLUE)
    for value in (zeeman - rabi, zeeman, zeeman + rabi):
        axis.axvline(value/zeeman, color=mc.GRID, lw=.7, zorder=0)
    axis.set_xlabel(r"probe frequency $\omega/Z$")
    axis.set_ylabel(r"$S^x_d$ probe (rate)")
    mc.panel_label(axis, "b")
    axis = axes[2]
    grid = np.linspace(0., .3, 601)*zeeman
    axis.plot(grid/zeeman, lines("spin_z").spectrum(grid, WIDTH*zeeman, kind="rate"), color=mc.BLUE)
    axis.axvline(rabi/zeeman, color=mc.GRID, lw=.7, zorder=0)
    axis.set_xlabel(r"probe frequency $\omega/Z$")
    axis.set_ylabel(r"$S^z_d$ probe (rate)")
    mc.panel_label(axis, "c")
    mc.finish(fig, "esr_probe")


def quoted(record):
    zeeman, rabi = record["zeeman"], record["rabi"]
    def strength(name, frequency):
        """Total strength of the lines with |omega| within half a Rabi splitting of frequency."""
        lines = record[name]
        keep = np.abs(np.abs(np.array(lines["frequencies"])) - frequency) < rabi/2
        return np.array(lines["strengths"])[keep]
    # The resonant drive is at Omega = Z: elastic lines at +-Z, Mollow sidebands at Z -+ delta.
    return {
        "gamma": (record["gamma"], ".1f"),
        "amplitude": (record["amplitude"], ".3f"),
        "harmonics": (record["harmonics"], "d"),
        "basis": (record["basis_states"], "d"),
        "zeeman": (zeeman, ".4f"),
        "rabi": (rabi, ".5f"),
        "ratio": (rabi/zeeman, ".1f"),
        "width": (WIDTH, ".3f"),
        "sx_rwa": (record["transverse"]**2/4, ".4f"),
        "sx_lower": (float(strength("spin_x", zeeman - rabi).sum()), ".3f"),
        "sx_elastic": (float(strength("spin_x", zeeman).mean()), ".3f"),
        "sx_upper": (float(strength("spin_x", zeeman + rabi).sum()), ".3f"),
        "sz_rwa": (((record["polarizations"][1] - record["polarizations"][0])/2)**2, ".4f"),
        "sz_line": (float(strength("spin_z", rabi).sum()), ".4f"),
    }


if __name__ == "__main__":
    mc.main(compute, plot)
