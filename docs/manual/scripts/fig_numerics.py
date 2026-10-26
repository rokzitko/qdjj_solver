"""Convergence tutorial: Floquet cutoff M, static basis K and the discretized continuum."""

import math

import numpy as np

import manual_common as mc

CUTOFF = 2                # QP cutoff of every calculation
OMEGA = .45               # drive frequency of (a) and (b)
STRONG_AMPLITUDE = 1.5    # phase-drive amplitude of (a)
WEAK_AMPLITUDE = .3       # phase-drive amplitude of (b) and (c, d)
NUMBERS = "zero one two three four five six seven eight nine ten".split()


def compute():
    from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet
    sector = Sector(0, 0)
    base = mc.junction()
    full = solve(reference_model(**base), CUTOFF, sector, SolverOptions(method="dense", eigenpairs=546))
    static = full
    # (a) Floquet cutoff at strong drive, in a static basis of at least 60 states (whole multiplets).
    strong = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", STRONG_AMPLITUDE),
                                                 OMEGA, samples=48)
    kept = mc.basis_size(static.energies, 60)
    options = dict(states=2, basis="eigen", basis_states=kept, dense_dimension=6000)
    reference = solve_floquet(strong, CUTOFF, sector, dict(floquet_cutoff=16, **options), static=static)
    cutoff = []
    for harmonics in range(1, 11):
        result = solve_floquet(strong, CUTOFF, sector, dict(floquet_cutoff=harmonics, **options), static=static)
        cutoff.append(dict(harmonics=harmonics, error=abs(result.quasienergies-reference.quasienergies).tolist(),
                           bound=result.residuals.tolist(), boundary=result.boundary_weights.tolist()))
    # (b) Static basis K against the full cutoff-2 Fock space at weaker drive.
    weak = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", WEAK_AMPLITUDE),
                                               OMEGA, samples=32)
    fock = solve_floquet(weak, CUTOFF, sector, dict(floquet_cutoff=4, states=2), static=static)
    basis = []
    for minimum in (4, 8, 16, 32, 64, 128, 256, 400, 546):
        size = mc.basis_size(static.energies, minimum)
        result = solve_floquet(weak, CUTOFF, sector, dict(floquet_cutoff=4, states=2, basis="eigen", basis_states=size,
                                                          dense_dimension=6000), static=static)
        basis.append(dict(basis_states=size, error=abs(result.quasienergies-fock.quasienergies).tolist(),
                          fock_residuals=result.fock_residuals.tolist()))
    # (c, d) The dressed singlet near one-photon resonances with the discretized continuum.
    frequencies = np.linspace(.95, 1.75, 321)
    continuum = {}
    for pairs in (2, 4, 8):
        model = mc.junction(pairs=pairs)
        h = reference_model(**model)
        states = solve(h, CUTOFF, sector, SolverOptions(method="dense", eigenpairs=120))
        odd = solve(h, CUTOFF, Sector(1, 1), SolverOptions(eigenpairs=1))
        size = mc.basis_size(states.energies, int(np.sum(states.energies < states.energies[0] + 2.4)))
        periodic = PeriodicHamiltonian.from_parameters(reference_model, model,
                                                       ParameterDrive.cosine("phi", WEAK_AMPLITUDE), 1., samples=32)
        rows = []
        for omega in frequencies:
            result = solve_floquet(mc.with_omega(periodic, omega), CUTOFF, sector,
                                   dict(floquet_cutoff=3, states=1, basis="eigen", basis_states=size,
                                        dense_dimension=6000, tracking_threshold=0.), static=states)
            rows.append(dict(omega=float(omega), shift=float(result.quasienergies[0]-states.energies[0]),
                             weight=float(result.static_weights[0, 0])))
        continuum[str(pairs)] = dict(rows=rows, basis_states=size, dimension=states.basis.dimension,
                                     threshold=float(odd.energies[0] + 1. - states.energies[0]),
                                     levels=(states.energies[:size]-states.energies[0]).tolist())
    mc.save("numerics", dict(cutoff=cutoff, cutoff_basis_states=kept, basis=basis, continuum=continuum,
                             fock_dimension=static.basis.dimension), __file__)


def plot():
    record = mc.load("numerics")
    fig, axes = mc.figure(2, 2, height=4.6)
    axis = axes[0, 0]
    rows = record["cutoff"]
    m = np.array([r["harmonics"] for r in rows])
    for state, color, name in ((0, mc.BLUE, "S"), (1, mc.ORANGE, "S^*")):
        axis.semilogy(m, [max(r["error"][state], 1e-16) for r in rows], "o-", ms=3.5, mfc="none", mew=.9,
                      color=color, label=fr"error, ${name}$")
        axis.semilogy(m, [r["bound"][state] for r in rows], ls="--", lw=1., color=color, label=fr"bound, ${name}$")
    axis.set_xlabel(r"Floquet cutoff $M$")
    axis.set_ylabel(r"$|\varepsilon_M-\varepsilon_{16}|$")
    axis.legend(loc="lower left", ncols=2)
    mc.panel_label(axis, "a")
    axis = axes[0, 1]
    rows = record["basis"]
    k = np.array([r["basis_states"] for r in rows])
    for state, color, name in ((0, mc.BLUE, "S"), (1, mc.ORANGE, "S^*")):
        axis.loglog(k, [max(r["error"][state], 1e-16) for r in rows], "o-", ms=3.5, mfc="none", mew=.9,
                    color=color, label=fr"error, ${name}$")
        axis.loglog(k, [r["fock_residuals"][state] for r in rows], ls="--", lw=1., color=color,
                    label=fr"Fock residual, ${name}$")
    axis.set_xlabel(r"static basis size $K$")
    axis.set_ylabel(r"$|\varepsilon_K-\varepsilon_{\mathrm{Fock}}|$")
    axis.legend(loc="center left", ncols=1)
    mc.panel_label(axis, "b")
    continuum = record["continuum"]
    for (pairs, data), color in zip(continuum.items(), (mc.BLUE, mc.ORANGE, mc.AQUA), strict=True):
        omega = np.array([r["omega"] for r in data["rows"]])
        axes[1, 0].plot(omega, [1e3*r["shift"] for r in data["rows"]], color=color, lw=1.,
                        label=f"{pairs} ({data['dimension']})")
        axes[1, 1].plot(omega, [r["weight"] for r in data["rows"]], color=color, lw=1.)
    threshold = continuum["4"]["threshold"]
    for axis in axes[1]:
        axis.axvline(threshold, color=mc.MUTED, lw=.8, ls=":")
        axis.set_xlabel(r"drive frequency $\Omega$")
        axis.set_xlim(omega[0], omega[-1])
    axes[1, 0].text(threshold+.01, .97, r"$E_D+\Delta-E_S$", transform=axes[1, 0].get_xaxis_transform(),
                    fontsize=7, color=mc.SECONDARY, va="top")
    axes[1, 0].set_ylabel(r"dressed-singlet shift $(\varepsilon_S-E_S)\times10^3$")
    axes[1, 0].set_ylim(-6, 6)
    axes[1, 0].legend(loc="lower left", title="pairs per lead (dimension)", title_fontsize=7)
    axes[1, 1].set_ylabel(r"static weight on $S$")
    axes[1, 1].set_ylim(.4, 1.02)
    mc.panel_label(axes[1, 0], "c")
    mc.panel_label(axes[1, 1], "d")
    mc.finish(fig, "numerics")


def spelled(*counts):
    """Small counts in words, as the text writes them ("two, four and eight"); \\text works in math mode."""
    words = [NUMBERS[count] for count in counts]
    text = words[-1] if len(words) == 1 else f"{', '.join(words[:-1])} and {words[-1]}"
    return rf"\text{{{text}}}"


def quoted(record):
    """Numbers of the convergence chapter, from the record and the parameters above."""
    basis, continuum = record["basis"], record["continuum"]
    # "About 10^-3 already for K~10": order of magnitude of the larger error at the scanned K nearest to 10.
    small = min(basis, key=lambda row: abs(math.log(row["basis_states"]/10)))
    # "Below 10^-4 when K approaches the full dimension": the largest scanned K short of it.
    large = max((row for row in basis if row["basis_states"] < record["fock_dimension"]),
                key=lambda row: row["basis_states"])
    pairs = sorted(int(key) for key in continuum)
    # Below the lowest one-photon threshold, away from the narrow (two-photon) dips of the static
    # weight, the AC Stark shifts of the different baths agree to within this relative spread.
    limit = min(data["threshold"] for data in continuum.values())
    shifts = {}
    for data in continuum.values():
        for row in data["rows"]:
            if row["omega"] < limit and row["weight"] >= .99:
                shifts.setdefault(round(row["omega"], 6), []).append(row["shift"])
    spread = max((max(s) - min(s))/abs(np.mean(s)) for s in shifts.values() if len(s) == len(continuum))
    dips = sum(row["weight"] < .95 for data in continuum.values() for row in data["rows"] if row["omega"] < limit)
    if not dips:
        raise ValueError("the text mentions avoided crossings below the one-photon threshold")
    return dict(
        cutoff=(CUTOFF, "d"), omega=(OMEGA, "g"), strong_amplitude=(STRONG_AMPLITUDE, "g"),
        weak_amplitude=(WEAK_AMPLITUDE, "g"), cutoff_basis_states=(record["cutoff_basis_states"], "d"),
        fock_dimension=(record["fock_dimension"], "d"), bandwidth=(mc.BATH["bandwidth"], "g"),
        small_basis_error=(f"10^{{{math.floor(math.log10(max(small['error'])))}}}", None),
        large_basis_error=(f"10^{{{math.ceil(math.log10(max(large['error'])))}}}", None),
        threshold=(continuum[str(mc.BATH["pairs"])]["threshold"], ".2f"),
        bath_pairs=(spelled(*pairs), None), bath_count=(spelled(len(pairs)), None),
        below_spread=(100*spread, ".0f"))


if __name__ == "__main__":
    mc.main(compute, plot)
