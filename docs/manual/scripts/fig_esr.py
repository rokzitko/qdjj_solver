"""Electron spin resonance of the Zeeman-split doublet and its Knight-shift renormalization."""

import numpy as np

import manual_common as mc

ESR = dict(u=1., detuning=0., phi=.8, field=.02, pairs=4, bandwidth=10.)


def esr_model(transverse=0., gamma=.3, u=ESR["u"], detuning=ESR["detuning"], field=ESR["field"],
              phi=ESR["phi"], pairs=ESR["pairs"], bandwidth=ESR["bandwidth"]):
    """Reference junction plus a transverse field b S_x on the dot (S_z is then not conserved)."""
    from qdjj_solver import Impurity, annihilate, create, hermitian_pair, make_model
    anderson = Impurity.anderson(u, detuning, field)
    # hermitian_pair gives (b/2)(d+_up d_dn + d+_dn d_up) = b S_x.
    operator = anderson.operator + hermitian_pair(.5*transverse*create(0)*annihilate(1))
    impurity = Impurity(1, operator, metadata=dict(kind="anderson-transverse", u=u, detuning=detuning, field=field))
    bath = mc.bath(pairs=pairs, bandwidth=bandwidth)
    hopping = np.sqrt(gamma/(np.pi*bath.rho*2))*np.eye(2)
    return make_model(impurity, [bath, bath], [hopping, hopping], phases=[-phi/2, phi/2],
                      phase_velocities=[-.5, .5])


def periodic(gamma, amplitude):
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian
    return PeriodicHamiltonian.from_parameters(esr_model, dict(transverse=0., gamma=gamma),
                                               ParameterDrive.cosine("transverse", amplitude), 1., samples=16)


def doublet(gamma):
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import SolverOptions, solve
    return solve(esr_model(gamma=gamma), 2, Sector(1, None), SolverOptions(method="dense", eigenpairs=40))


def rotating(drive, omega, static):
    """Quasienergies of the two doublet states, the upper one referred to the -1 photon replica."""
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import solve_floquet
    # At resonance the replicas form a chain |D_dn, m> - |D_up, m-1>; the Floquet cutoff must hold several links.
    result = solve_floquet(mc.with_omega(drive, omega), 2, Sector(1, None),
                           dict(floquet_cutoff=8, states=2, basis="eigen",
                                basis_states=mc.basis_size(static.energies, 24)), static=static)
    order = np.argsort(result.labels)
    eps = result.quasienergies[order]
    return np.sort([eps[0], eps[1] - omega]), result


def compute():
    from scipy.optimize import minimize_scalar
    from qdjj_solver import Sector
    from qdjj_solver.qp_solver import rwa_model
    sector = Sector(1, None)
    gamma, amplitude = .3, .004
    static = doublet(gamma)
    zeeman = float(static.energies[1] - static.energies[0])
    drive = periodic(gamma, amplitude)
    sweep = []
    for omega in np.linspace(zeeman - .006, zeeman + .006, 121):
        branches, _ = rotating(drive, omega, static)
        rwa = rwa_model(mc.with_omega(drive, omega), [0, 1], [0, -1], 2, sector, static=static)
        order = np.argsort(rwa.labels)
        sweep.append(dict(omega=float(omega), branches=branches.tolist(),
                          rwa=sorted((rwa.quasienergies[order] - [0, omega]).tolist())))
    renormalization = []
    for gamma in np.linspace(.0, 1.2, 13):
        static = doublet(float(gamma))
        spin = float(abs(np.real(static.observables["impurity_spin_z"][0])))
        zeeman = float(static.energies[1] - static.energies[0])
        drive = periodic(float(gamma), amplitude)
        found = minimize_scalar(lambda omega, drive=drive, static=static: np.diff(rotating(drive, omega, static)[0])[0],
                                bounds=(zeeman - .002, zeeman + .002), method="bounded", options=dict(xatol=1e-11))
        renormalization.append(dict(gamma=float(gamma), spin=spin, zeeman=zeeman, rabi=float(found.fun),
                                    resonance=float(found.x)))
    mc.save("esr", dict(parameters=ESR, gamma=.3, amplitude=amplitude, sweep=sweep,
                        renormalization=renormalization), __file__)


def plot():
    record = mc.load("esr")
    fig, axes = mc.figure(1, 2, height=2.5)
    axis = axes[0]
    sweep = record["sweep"]
    omega = np.array([r["omega"] for r in sweep])
    branches = np.array([r["branches"] for r in sweep])
    rwa = np.array([r["rwa"] for r in sweep])
    center = branches[len(branches)//2].mean()
    axis.plot(1e3*omega, 1e3*(branches[:, 0]-center), color=mc.BLUE, label="Sambe")
    axis.plot(1e3*omega, 1e3*(branches[:, 1]-center), color=mc.BLUE)
    axis.plot(1e3*omega, 1e3*(rwa[:, 0]-center), color=mc.ORANGE, ls="--", lw=1.1, label="RWA")
    axis.plot(1e3*omega, 1e3*(rwa[:, 1]-center), color=mc.ORANGE, ls="--", lw=1.1)
    axis.set_xlabel(r"drive frequency $\Omega\times10^3$")
    axis.set_ylabel(r"rotating-frame quasienergy $\times10^3$")
    axis.legend(loc="lower left")
    mc.panel_label(axis, "a")
    axis = axes[1]
    rows = record["renormalization"]
    gamma = np.array([r["gamma"] for r in rows])
    field, amplitude = record["parameters"]["field"], record["amplitude"]
    axis.plot(gamma, [2*r["spin"] for r in rows], color=mc.INK, lw=1.1, label=r"$2|\langle S^z_d\rangle|$ (static)")
    axis.plot(gamma, [r["zeeman"]/field for r in rows], "o", ms=3.5, mfc="none", mew=.9, color=mc.BLUE,
              label=r"Zeeman splitting $/B$")
    axis.plot(gamma, [r["rabi"]/(amplitude/2) for r in rows], "s", ms=3.5, mfc="none", mew=.9, color=mc.ORANGE,
              label=r"Rabi splitting $/(b/2)$")
    axis.set_xlabel(r"hybridization $\Gamma/\Delta$")
    axis.set_ylabel(r"$g_{\mathrm{eff}}/g$")
    axis.set_ylim(.5, 1.05)
    axis.legend(loc="lower left")
    mc.panel_label(axis, "b")
    mc.finish(fig, "esr")


def _scientific(value):
    """One significant digit in LaTeX, for example 2\\times10^{-3}."""
    exponent = int(np.floor(np.log10(value)))
    mantissa = round(value/10**exponent)
    if mantissa == 10:
        mantissa, exponent = 1, exponent + 1
    return f"10^{{{exponent}}}" if mantissa == 1 else f"{mantissa}\\times10^{{{exponent}}}"


def quoted(record):
    parameters, rows = record["parameters"], record["renormalization"]
    field, amplitude = parameters["field"], record["amplitude"]
    # The three determinations of g_eff/g: static polarization, Zeeman splitting / B, Rabi splitting / (b/2).
    estimates = np.array([[2*r["spin"], r["zeeman"]/field, r["rabi"]/(amplitude/2)] for r in rows])
    spread = float(np.max(estimates.max(axis=1) - estimates.min(axis=1)))
    (row,) = [r for r in rows if np.isclose(r["gamma"], record["gamma"])]
    return {
        "gamma": (record["gamma"], ".1f"),
        "amplitude": (amplitude, ".3f"),
        "epsilon": (parameters["detuning"], ".0f"),
        "u": (parameters["u"], ".0f"),
        "phi": (parameters["phi"], ".1f"),
        "field": (field, ".2f"),
        "geff": (2*row["spin"], ".3f"),
        "agreement": (_scientific(spread), None),
        "geff_first": (2*rows[0]["spin"], ".2g"),
        "gamma_first": (rows[0]["gamma"], "g"),
        "geff_last": (2*rows[-1]["spin"], ".2f"),
        "gamma_last": (rows[-1]["gamma"], ".1f"),
    }


if __name__ == "__main__":
    mc.main(compute, plot)
