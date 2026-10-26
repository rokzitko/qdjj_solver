"""Gate-driven isolated dot: Bessel sidebands and Floquet zones (Tien-Gordon physics)."""

import numpy as np

import manual_common as mc

U, EPSILON = 1., .2


def model(detuning):
    from qdjj_solver import Hamiltonian, Impurity
    impurity = Impurity.anderson(U, detuning)
    return Hamiltonian(impurity.operator, 2, (1, -1), observables=impurity.observables)


def periodic(amplitude, omega, samples=16):
    from qdjj_solver import PeriodicHamiltonian
    # An isolated impurity declares no canonical coordinates: accept them explicitly.
    return PeriodicHamiltonian.from_callable(lambda t: model(EPSILON + amplitude*np.cos(omega*t)), omega,
                                             samples=samples, assume_fixed_coordinates=True)


def compute():
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    omega = 1.
    ratios = np.linspace(0., 4., 41)
    weights = {}
    for keep in (False, True):
        rows = []
        for ratio in ratios:
            result = solve_floquet(periodic(ratio*omega, omega), 0, Sector(0, 0),
                                   dict(floquet_cutoff=14, states=2, drop_scalar_harmonics=not keep))
            charge = np.rint(np.real(result.static.observables["impurity_charge"][result.labels]))
            index = int(np.flatnonzero(charge == 2)[0])
            center = int(np.flatnonzero(result.harmonic_indices == 0)[0])
            rows.append(dict(ratio=float(ratio), weights=result.harmonic_weights[index, center:center+4].tolist(),
                             quasienergy=float(result.quasienergies[index])))
        weights["kept" if keep else "dropped"] = rows
    # Zone folding: quasienergies do not depend on the drive, harmonic weights do.
    amplitude = .6
    frequencies = np.linspace(.15, 1.6, 146)
    fan = []
    for omega in frequencies:
        result = solve_floquet(periodic(amplitude, omega), 0, Sector(0, 0),
                               dict(floquet_cutoff=12, states=2, drop_scalar_harmonics=False))
        center = int(np.flatnonzero(result.harmonic_indices == 0)[0])
        fan.append(dict(omega=float(omega), quasienergies=result.quasienergies.tolist(),
                        weights=result.harmonic_weights.tolist(), center=center,
                        charge=np.rint(np.real(result.static.observables["impurity_charge"][result.labels])).tolist()))
    mc.save("dot", dict(parameters=dict(u=U, detuning=EPSILON), bessel=weights, fan_amplitude=amplitude, fan=fan),
            __file__)


def plot():
    from scipy.special import jv
    record = mc.load("dot")
    fig, axes = mc.figure(1, 2, height=2.5)
    axis = axes[0]
    rows = record["bessel"]["kept"]
    ratio = np.array([r["ratio"] for r in rows])
    fine = np.linspace(0., ratio[-1], 400)
    for m, color in zip(range(3), (mc.BLUE, mc.ORANGE, mc.AQUA), strict=True):
        axis.plot(fine, jv(m, fine)**2, color=color, lw=1.2, label=fr"$J_{m}^2(A/\Omega)$")
        axis.plot(ratio[::2], [r["weights"][m] for r in rows][::2], "o", ms=3.2, mfc="none", mew=.9, color=color)
    axis.set_xlabel(r"modulation index $A/\Omega$")
    axis.set_ylabel(r"harmonic weight $\Vert u_m\Vert^2$ of $|2\rangle$")
    axis.set_xlim(0, ratio[-1])
    axis.set_ylim(0, 1.03)
    axis.legend(loc="upper right", title="lines: Bessel; circles: Sambe", title_fontsize=7)
    mc.panel_label(axis, "a")
    axis = axes[1]
    fan = record["fan"]
    omega = np.array([r["omega"] for r in fan])
    axis.fill_between(omega, -omega/2, omega/2, color="#f0efec", lw=0, zorder=0)
    for state in range(2):
        energy = np.array([r["quasienergies"][state] for r in fan])
        color = mc.BLUE if fan[0]["charge"][state] == 2 else mc.ORANGE
        for m in range(-8, 9):
            weight = np.array([r["weights"][state][r["center"]+m] for r in fan])
            mc.weighted_line(axis, omega, energy + m*omega, weight, color)
    axis.set_xlim(omega[0], omega[-1])
    axis.set_ylim(-.9, 1.1)
    axis.set_xlabel(r"drive frequency $\Omega$")
    axis.set_ylabel(r"$\varepsilon+m\Omega$")
    axis.text(1.58, -.62, r"first zone $|\varepsilon|<\Omega/2$", color=mc.SECONDARY, fontsize=7, ha="right")
    axis.plot([], [], color=mc.BLUE, label=r"$|2\rangle$ and replicas")
    axis.plot([], [], color=mc.ORANGE, label=r"$|0\rangle$ and replicas")
    axis.legend(loc="upper right")
    mc.panel_label(axis, "b")
    mc.finish(fig, "dot")


def quoted(record):
    """Model parameters of the figure."""
    parameters = record["parameters"]
    return dict(u=(parameters["u"], "g"), epsilon=(parameters["detuning"], "g"),
                fan_amplitude=(record["fan_amplitude"], "g"))


if __name__ == "__main__":
    mc.main(compute, plot)
