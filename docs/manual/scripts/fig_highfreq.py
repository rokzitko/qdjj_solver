"""High-frequency limit: van Vleck error scaling and coherent destruction of tunnelling."""

import numpy as np

import manual_common as mc

HIGH = dict(u=1., gamma=.8, phi=.8, detuning=0., pairs=1, bandwidth=2., omega=40.)


def single(phase=0., detuning=.15, u=1.3, gamma=.35, bandwidth=1.5):
    from qdjj_solver import Impurity, cosh_grid, make_model
    bath = cosh_grid(pairs=1, bandwidth=bandwidth)
    t = np.sqrt(gamma/(np.pi*bath.rho))*np.eye(2)
    return make_model(Impurity.anderson(u, detuning), [bath], [t], phases=[phase])


def match(a, b, omega):
    from scipy.optimize import linear_sum_assignment
    from qdjj_solver.common.floquet import quasienergy_distance
    distance = quasienergy_distance(np.asarray(a)[:, None], np.asarray(b)[None, :], omega)
    rows, columns = linear_sum_assignment(distance)
    return float(distance[rows, columns].max())


def compute():
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet
    # (a) Phase and level driven at one harmonic with a relative phase: all three orders are finite.
    drives = (ParameterDrive.cosine("phase", .8), ParameterDrive.cosine("detuning", .5, .7))
    sector = Sector(1, 1)
    scaling = []
    for omega in (10., 14., 20., 28., 40., 56., 80., 113., 160.):
        periodic = PeriodicHamiltonian.from_parameters(single, dict(phase=.5, detuning=.15), drives, omega,
                                                       samples=16)
        exact = solve_floquet(periodic, None, sector, dict(selection="all", floquet_cutoff=4, method="dense"))
        errors = []
        for order in (0, 1, 2):
            energies = solve(periodic.van_vleck(order), None, sector,
                             SolverOptions(method="dense", eigenpairs=exact.basis.dimension)).energies
            errors.append(match(energies, exact.quasienergies, omega))
        scaling.append(dict(omega=omega, errors=errors))
    # (b) Strong phase drive at omega >> many-body width: Gamma -> Gamma J_0(A/4)^2.
    base = mc.junction(**{k: HIGH[k] for k in ("u", "gamma", "phi", "detuning", "pairs", "bandwidth")})
    rows = []
    for amplitude in np.linspace(0., 12., 49):
        periodic = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", amplitude),
                                                       HIGH["omega"], samples=48)
        row = dict(amplitude=float(amplitude))
        for name, sector, count in (("even", Sector(0, 0), 2), ("odd", Sector(1, 1), 1)):
            result = solve_floquet(periodic, 2, sector, dict(selection="all", floquet_cutoff=3, method="dense"))
            row[name] = np.sort(result.quasienergies)[:count].tolist()
        rows.append(row)
    effective = []
    from scipy.special import jv
    for amplitude in np.linspace(0., 12., 241):
        h = reference_model(**{**base, "gamma": HIGH["gamma"]*jv(0, amplitude/4)**2})
        even = solve(h, 2, Sector(0, 0), SolverOptions(eigenpairs=2, method="dense")).energies
        odd = solve(h, 2, Sector(1, 1), SolverOptions(eigenpairs=1, method="dense")).energies
        effective.append(dict(amplitude=float(amplitude), even=even.tolist(), odd=odd.tolist()))
    mc.save("highfreq", dict(scaling=scaling, strong=dict(parameters=HIGH, floquet=rows, effective=effective)),
            __file__)


def plot():
    record = mc.load("highfreq")
    fig, axes = mc.figure(1, 2, height=2.6)
    axis = axes[0]
    omega = np.array([r["omega"] for r in record["scaling"]])
    errors = np.array([r["errors"] for r in record["scaling"]])
    for order, color, marker in zip(range(3), (mc.BLUE, mc.ORANGE, mc.AQUA), "osD", strict=True):
        axis.loglog(omega, errors[:, order], marker=marker, ms=3.5, mfc="none", mew=.9, color=color,
                    label=fr"order {order}")
        reference = errors[-1, order]*(omega/omega[-1])**(-(order+1))
        axis.loglog(omega, reference, color=color, lw=.7, ls=":")
    axis.set_xlabel(r"drive frequency $\Omega$")
    axis.set_ylabel(r"max $|\varepsilon_{\mathrm{vV}}-\varepsilon_{\mathrm{Sambe}}|$")
    axis.legend(loc="lower left", title=r"dotted: $\Omega^{-(n+1)}$", title_fontsize=7)
    mc.panel_label(axis, "a")
    axis = axes[1]
    strong = record["strong"]
    effective = strong["effective"]
    amplitude = np.array([r["amplitude"] for r in effective])
    even = np.array([r["even"] for r in effective])
    odd = np.array([r["odd"][0] for r in effective])
    axis.plot(amplitude, even[:, 0], color=mc.BLUE, lw=1.1, label=r"singlet $S$")
    axis.plot(amplitude, even[:, 1], color=mc.BLUE, lw=1.1, ls="--", label=r"$S^*$")
    axis.plot(amplitude, odd, color=mc.ORANGE, lw=1.1, label=r"doublet $D$")
    floquet = strong["floquet"]
    a = np.array([r["amplitude"] for r in floquet])
    axis.plot(a, [r["even"][0] for r in floquet], "o", ms=3, mfc="none", mew=.9, color=mc.BLUE)
    axis.plot(a, [r["even"][1] for r in floquet], "o", ms=3, mfc="none", mew=.9, color=mc.BLUE)
    axis.plot(a, [r["odd"][0] for r in floquet], "s", ms=3, mfc="none", mew=.9, color=mc.ORANGE)
    zero = 4*2.404825557695773
    axis.axvline(zero, color=mc.GRID, lw=.8, zorder=0)
    axis.text(zero-.15, -.6, r"$J_0(A/4)=0$", ha="right", fontsize=7, color=mc.SECONDARY)
    axis.set_xlabel(r"phase-drive amplitude $A$")
    axis.set_ylabel(r"quasienergy ($\Delta$)")
    axis.set_xlim(0, 12)
    for text, y, color in ((r"$S$", even[0, 0], mc.BLUE), (r"$S^*$", even[0, 1], mc.BLUE),
                           (r"$D$", odd[0], mc.ORANGE)):
        axis.text(.25, y + .035, text, color=color, fontsize=8, va="bottom")
    axis.text(.25, .47, "lines: $\\Gamma\\to\\Gamma J_0(A/4)^2$\nmarkers: Floquet", fontsize=7,
              color=mc.SECONDARY, va="top")
    mc.panel_label(axis, "b")
    mc.finish(fig, "highfreq")


WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def quoted(record):
    """Parameters of panel (b) and the drive amplitude at which the dressed doublet becomes lowest."""
    strong = record["strong"]
    parameters = strong["parameters"]
    amplitude = np.array([r["amplitude"] for r in strong["floquet"]])
    gap = np.array([r["odd"][0] - r["even"][0] for r in strong["floquet"]])
    # Linear interpolation of the first sign change of the doublet-singlet quasienergy difference.
    j = int(np.flatnonzero(gap < 0)[0])
    crossing = amplitude[j-1] - gap[j-1]*(amplitude[j] - amplitude[j-1])/(gap[j] - gap[j-1])
    return dict(u=(parameters["u"], "g"), gamma=(parameters["gamma"], "g"), detuning=(parameters["detuning"], "g"),
                phi=(parameters["phi"], "g"), bandwidth=(parameters["bandwidth"], "g"),
                pairs=(rf"\text{{{WORDS[parameters['pairs']]}}}", None), omega=(parameters["omega"], "g"),
                crossing=(float(crossing), ".1f"))


if __name__ == "__main__":
    mc.main(compute, plot)
