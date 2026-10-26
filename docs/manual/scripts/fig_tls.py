"""Proximitized dot as a driven two-level system: Rabi, Autler-Townes, Bloch-Siegert, multiphoton."""

import numpy as np

import manual_common as mc

U = 1.


def model(detuning, pairing):
    """Large-gap (proximitized) dot: U/2(n-1)^2 + eps(n-1) - Gamma_phi(d+_up d+_dn + h.c.)."""
    from qdjj_solver import Hamiltonian, Impurity, create, hermitian_pair
    impurity = Impurity.anderson(U, detuning)
    operator = impurity.operator + hermitian_pair(-pairing*create(0)*create(1))
    return Hamiltonian(operator, 2, (1, -1), observables=impurity.observables)


def periodic(detuning, pairing, amplitude, omega, samples=16):
    from qdjj_solver import PeriodicHamiltonian
    return PeriodicHamiltonian.from_callable(lambda t: model(detuning + amplitude*np.cos(omega*t), pairing),
                                             omega, samples=samples, assume_fixed_coordinates=True)


def pair(detuning, pairing, amplitude, omega, harmonics=12):
    """Quasienergies of the two tracked states (S first, then S*) and the full result."""
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    result = solve_floquet(periodic(detuning, pairing, amplitude, omega), 0, Sector(0, 0),
                           dict(floquet_cutoff=harmonics, states=2))
    return result.quasienergies, result


def splitting(detuning, pairing, amplitude, omega, photons=1):
    from qdjj_solver.common.floquet import quasienergy_distance
    (lower, upper), _ = pair(detuning, pairing, amplitude, omega, harmonics=max(12, 4*photons+8))
    return float(quasienergy_distance(lower, upper, omega))


def resonance(detuning, pairing, amplitude, guess):
    """Frequency of the minimal quasienergy splitting near a one-photon resonance."""
    from scipy.optimize import minimize_scalar
    # The Bloch-Siegert shift is positive and grows as A^2, so search above the bare resonance.
    found = minimize_scalar(lambda omega: splitting(detuning, pairing, amplitude, omega),
                            bounds=(guess-.03, guess+.12), method="bounded", options=dict(xatol=1e-9))
    return float(found.x), float(found.fun)


def compute():
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import rwa_model
    # At epsilon=0 the gate drive is purely transverse in the eigenbasis: the textbook Rabi problem.
    epsilon, pairing = 0., .5
    root = np.hypot(epsilon, pairing)
    omega0 = 2*root
    weak = dict(epsilon=epsilon, pairing=pairing, omega0=omega0, sin_theta=pairing/root)
    amplitude = .1
    sweep = []
    for omega in np.linspace(.85, 1.15, 121):
        (lower, upper), result = pair(epsilon, pairing, amplitude, omega)
        rwa = rwa_model(periodic(epsilon, pairing, amplitude, omega), [0, 1], [0, -1], 0, Sector(0, 0))
        order = np.argsort(rwa.labels)
        sweep.append(dict(omega=float(omega), lower=float(lower), upper=float(upper - omega),
                          rwa=sorted(rwa.quasienergies[order]*[1, 1] - [0, omega]),
                          weight=float(result.static_weights[0, 0])))
    weak.update(amplitude=amplitude, sweep=sweep)
    rows = []
    for amplitude in np.linspace(.02, .5, 25):
        omega, gap = resonance(epsilon, pairing, amplitude, omega0)
        rows.append(dict(amplitude=float(amplitude), resonance=omega, gap=gap))
    weak["resonances"] = rows
    # Strong driving with a weak pairing: k-photon resonances at 2 eps = k omega.
    strong = dict(pairing=.05, omega=1.)
    curves = {}
    for photons in (1, 2, 3):
        values = []
        for amplitude in np.linspace(0., 6., 121):
            values.append(dict(amplitude=float(amplitude),
                               gap=splitting(photons*strong["omega"]/2, strong["pairing"], amplitude, strong["omega"],
                                             photons)))
        curves[str(photons)] = values
    strong["curves"] = curves
    # Quasienergies versus bias at the first zero of J_2(2A/omega).
    from scipy.special import jn_zeros
    amplitude = float(jn_zeros(2, 1)[0]/2)
    fan = []
    for epsilon in np.linspace(0., 2.2, 441):
        (lower, upper), _ = pair(epsilon, strong["pairing"], amplitude, strong["omega"], harmonics=24)
        fan.append(dict(epsilon=float(epsilon), lower=float(lower), upper=float(upper)))
    strong.update(fan_amplitude=amplitude, fan=fan)
    mc.save("tls", dict(u=U, weak=weak, strong=strong), __file__)


def plot():
    from scipy.special import jv
    record = mc.load("tls")
    weak = record["weak"]
    fig, axes = mc.figure(1, 3, height=2.3)
    axis = axes[0]
    sweep = weak["sweep"]
    omega = np.array([r["omega"] for r in sweep])
    base = U/2 - np.hypot(weak["epsilon"], weak["pairing"])
    axis.plot(omega, np.full_like(omega, base), color=mc.MUTED, lw=.8, ls=":")
    axis.plot(omega, base + weak["omega0"] - omega, color=mc.MUTED, lw=.8, ls=":")
    # Labels swap at resonance, where the hybrids are equal mixtures; sort the two branches.
    branches = np.sort([[r["lower"], r["upper"]] for r in sweep], axis=1)
    axis.plot(omega, branches[:, 0], color=mc.BLUE, label="Sambe")
    axis.plot(omega, branches[:, 1], color=mc.BLUE)
    rwa = np.array([r["rwa"] for r in sweep])
    axis.plot(omega, rwa[:, 0], color=mc.ORANGE, ls="--", lw=1.1, label="RWA")
    axis.plot(omega, rwa[:, 1], color=mc.ORANGE, ls="--", lw=1.1)
    axis.set_xlabel(r"drive frequency $\Omega$")
    axis.set_ylabel("rotating-frame quasienergy")
    axis.set_xlim(omega[0], omega[-1])
    axis.legend(loc="lower left")
    mc.panel_label(axis, "a")
    rows = weak["resonances"]
    amplitude = np.array([r["amplitude"] for r in rows])
    axis = axes[1]
    axis.plot(amplitude, [r["gap"] for r in rows], "o", ms=3, mfc="none", mew=.9, color=mc.BLUE, label="Sambe")
    axis.plot(amplitude, amplitude*weak["sin_theta"], color=mc.ORANGE, lw=1.1, label=r"RWA: $A$")
    axis.set_xlabel(r"drive amplitude $A$")
    axis.set_ylabel(r"splitting at resonance")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "b")
    axis = axes[2]
    shift = np.array([r["resonance"] for r in rows]) - weak["omega0"]
    axis.plot(amplitude, shift, "o", ms=3, mfc="none", mew=.9, color=mc.BLUE, label="Sambe")
    axis.plot(amplitude, amplitude**2/(4*weak["omega0"]), color=mc.ORANGE, lw=1.1,
              label=r"Bloch–Siegert $A^2/4\omega_0$")
    axis.set_xlabel(r"drive amplitude $A$")
    axis.set_ylabel(r"resonance shift $\Omega_{\mathrm{res}}-\omega_0$")
    axis.legend(loc="upper left")
    mc.panel_label(axis, "c")
    mc.finish(fig, "tls_weak")

    strong = record["strong"]
    fig, axes = mc.figure(1, 2, height=2.5)
    axis = axes[0]
    fine = np.linspace(0, 6, 400)
    for photons, color in zip(("1", "2", "3"), (mc.BLUE, mc.ORANGE, mc.AQUA), strict=True):
        rows = strong["curves"][photons]
        axis.plot(fine, 2*strong["pairing"]*abs(jv(int(photons), 2*fine/strong["omega"])), color=color, lw=1.1,
                  label=fr"$k={photons}$")
        axis.plot([r["amplitude"] for r in rows][::4], [r["gap"] for r in rows][::4], "o", ms=3, mfc="none",
                  mew=.9, color=color)
    axis.set_xlabel(r"drive amplitude $A/\Omega$")
    axis.set_ylabel(r"splitting at $2ϵ_0=k\Omega$")
    axis.legend(loc="upper right", title=r"lines: $2\Gamma_\phi|J_k(2A/\Omega)|$", title_fontsize=7)
    mc.panel_label(axis, "a")
    axis = axes[1]
    from qdjj_solver.common.floquet import quasienergy_distance
    fan = strong["fan"]
    epsilon = np.array([r["epsilon"] for r in fan])
    gap = quasienergy_distance([r["upper"] for r in fan], [r["lower"] for r in fan], strong["omega"])
    axis.semilogy(epsilon, gap, color=mc.BLUE, lw=1.1, label="Sambe")
    for photons in range(0, 5):
        value = 2*strong["pairing"]*abs(jv(photons, 2*strong["fan_amplitude"]/strong["omega"]))
        axis.plot([photons/2], [max(value, 2e-5)], marker="v" if value < 1e-4 else "o", ms=4, color=mc.ORANGE,
                  ls="none", label=r"$2\Gamma_\phi|J_k(2A/\Omega)|$" if photons == 0 else None)
    axis.set_xlabel(r"detuning $ϵ_0/\Omega$")
    axis.set_ylabel(r"quasienergy splitting")
    axis.set_xlim(epsilon[0], epsilon[-1])
    axis.set_ylim(1e-5, .6)
    axis.text(1., 3e-5, r"$J_2=0$", ha="center", fontsize=7, color=mc.SECONDARY)
    axis.legend(loc="lower right")
    mc.panel_label(axis, "b")
    mc.finish(fig, "tls_strong")


def quoted(record):
    """Parameters of the weak- and strong-driving figures."""
    weak, strong = record["weak"], record["strong"]
    return dict(weak_epsilon=(weak["epsilon"], "g"), weak_pairing=(weak["pairing"], "g"),
                weak_omega0=(weak["omega0"], "g"), weak_amplitude=(weak["amplitude"], "g"),
                strong_pairing=(strong["pairing"], "g"), strong_omega=(strong["omega"], "g"))


if __name__ == "__main__":
    mc.main(compute, plot)
