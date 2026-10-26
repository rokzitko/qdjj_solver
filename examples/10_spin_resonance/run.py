"""Tutorial 10: electron spin resonance of the doublet of a superconducting junction.

    python -B examples/10_spin_resonance/run.py                  # quick profile, a few seconds
    python -B examples/10_spin_resonance/run.py --profile full   # the complete calculation, under a minute

A static field B S^z_d splits the doublet into D_dn and D_up; a transverse drive b cos(Omega t) S^x_d
flips the spin. The coupling to the leads partially screens the impurity spin, which renormalizes
the g-factor. The script computes (1) g_eff/g from the static spin polarization, the Zeeman
splitting and the Floquet Rabi splitting at resonance, (2) the probes of the driven doublet: the
Mollow pattern of a transverse probe S^x_d and the line of a longitudinal probe S^z_d at the Rabi
frequency, and (3) g_eff/g versus the hybridization. Results go to
runs/examples/10_spin_resonance/<profile>/ (results.json, and spin_resonance.png if matplotlib is installed).
"""

import argparse
import dataclasses
import json
import os
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np  # noqa: E402

from qdjj_solver import (Impurity, ParameterDrive, PeriodicHamiltonian, Sector, annihilate, cosh_grid,  # noqa: E402
                         create, hermitian_pair, make_model)
from qdjj_solver.qp_solver import SolverOptions, solve, solve_floquet  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The doublet of the physics manual's ESR chapter: Delta = 1, four bath level pairs per lead, half filling.
GAMMA = .3          # total hybridization
FIELD = .02         # static Zeeman field B (energy units, g mu_B = 1)
DRIVE = .004        # transverse drive amplitude b
HARMONICS = 8       # Floquet cutoff M; the resonant replica chain needs M >~ 5
BASIS = 24          # minimal size K of the static eigenbasis
WIDTH = .006        # half width of the plotted lines, in units of the Zeeman splitting
SECTOR = Sector(1, None)    # odd parity, all S_z: the transverse drive does not conserve S_z
PROFILES = {
    # drive frequencies Omega/Z of the probe sweep, hybridizations of the g-factor scan
    "quick": dict(sweep=np.linspace(.8, 1.2, 5), gammas=np.array([0., .3, 1.2])),
    "full": dict(sweep=np.linspace(.8, 1.2, 81), gammas=np.round(np.linspace(0., 1.2, 13), 10)),
}


def esr_model(transverse=0., gamma=GAMMA, u=1., detuning=0., field=FIELD, phi=.8, pairs=4, bandwidth=10.):
    """Reference junction with a Zeeman field B S^z_d and a transverse field b S^x_d on the dot."""
    bath = cosh_grid(pairs=pairs, delta=1., bandwidth=bandwidth)
    anderson = Impurity.anderson(u, detuning, field)       # U/2 (n-1)^2 + eps (n-1) + B S^z_d
    # hermitian_pair(X) = X + X^dagger, so this adds (b/2)(d+_up d_dn + d+_dn d_up) = b S^x_d.
    operator = anderson.operator + hermitian_pair(.5*transverse*create(0)*annihilate(1))
    impurity = Impurity(1, operator, metadata=dict(kind="anderson-transverse"))
    hopping = np.sqrt(gamma/(np.pi*bath.rho*2))*np.eye(2)  # gamma is the total hybridization of both leads
    return make_model(impurity, [bath, bath], [hopping, hopping], phases=[-phi/2, phi/2],
                      phase_velocities=[-.5, .5])


def basis_size(energies, minimum):
    """Smallest K >= minimum that does not split a degenerate static multiplet."""
    gaps = np.flatnonzero(np.diff(energies) > 1e-7) + 1
    gaps = gaps[gaps >= minimum]
    return int(gaps[0]) if len(gaps) else len(energies)


def transverse_drive(gamma):
    """The drive b cos(Omega t) S^x_d at hybridization gamma, sampled at Omega = 1."""
    return PeriodicHamiltonian.from_parameters(esr_model, dict(transverse=0., gamma=gamma),
                                               ParameterDrive.cosine("transverse", DRIVE), 1., samples=8)


def driven(periodic, omega, static):
    """All Floquet states of the driven doublet in the static basis, at drive frequency omega."""
    # The harmonics of p0 + A cos(Omega t) do not depend on Omega, so only the frequency is replaced.
    periodic = dataclasses.replace(periodic, omega=omega)
    return solve_floquet(periodic, 2, SECTOR,
                         dict(floquet_cutoff=HARMONICS, selection="all", method="dense", basis="eigen",
                              basis_states=basis_size(static.energies, BASIS)), static=static)


def doublet_states(states):
    """Indices of the two Floquet states with the largest weight on the static doublet D_dn, D_up."""
    return np.argsort(states.static_weights[:, 0] + states.static_weights[:, 1])[-2:]


def splitting(states):
    """Rotating-frame splitting w of the two dressed doublet states: their quasienergy difference modulo Omega."""
    first, second = doublet_states(states)
    difference = states.quasienergies[second] - states.quasienergies[first]
    return float(abs((difference + states.omega/2) % states.omega - states.omega/2))


def resonance(periodic, static, zeeman, rabi):
    """Rabi splitting and resonance frequency from a parabola w^2 = a (Omega - Omega_r)^2 + delta^2 (a = 1 in the RWA)."""
    omegas = zeeman + rabi*np.linspace(-1., 1., 5)
    squares = [splitting(driven(periodic, float(omega), static))**2 for omega in omegas]
    a, b, c = np.polyfit(omegas - zeeman, squares, 2)
    return float(np.sqrt(c - b*b/(4*a))), float(zeeman - b/(2*a)), float(a)


def probes(states, zeeman):
    """Mollow pattern of S^x_d and the Rabi line of S^z_d out of the dressed lower Zeeman level.

    Exactly at the resonance both dressed states carry equal weights of D_dn, and the continuation of D_dn
    is not defined; at Omega = Z, slightly below the shifted resonance, it is.
    """
    initial = states.dressed_state(0)
    sx = states.matrix_elements("spin_x").lines(initial=initial)
    sz = states.matrix_elements("spin_z").lines(initial=initial)
    sx = [(f, s) for f, s in zip(sx.frequencies, sx.strengths, strict=True) if abs(f) < 2*zeeman]
    sz = [(f, s) for f, s in zip(sz.frequencies, sz.strengths, strict=True) if abs(f) < zeeman/2]
    return dict(spin_x=dict(frequencies=[float(f) for f, _ in sx], strengths=[float(s) for _, s in sx]),
                spin_z=dict(frequencies=[float(f) for f, _ in sz], strengths=[float(s) for _, s in sz]))


def mollow(record, omega, w):
    """Strengths of the emission sideband at -(Omega - w), the elastic lines at +-Omega and the absorption
    sideband at Omega + w, for a record of ``probes``; w is the rotating-frame splitting."""
    f, s = np.array(record["frequencies"]), np.array(record["strengths"])
    def near(center):
        return float(s[abs(f - center) < w/2].max())
    return near(-(omega - w)), (near(omega) + near(-omega))/2, near(omega + w)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    parser.add_argument("--output", type=Path, help="output directory (default: runs/examples/10_spin_resonance/<profile>)")
    args = parser.parse_args(argv)
    settings = PROFILES[args.profile]
    output = args.output or ROOT/"runs"/"examples"/HERE.name/args.profile
    output.mkdir(parents=True, exist_ok=True)

    # 1. Static doublet and Floquet resonance for each hybridization (computed once per Gamma).
    cache = {}
    def doublet(gamma):
        if gamma not in cache:
            static = solve(esr_model(gamma=gamma), 2, SECTOR, SolverOptions(method="dense", eigenpairs=60))
            zeeman = float(static.energies[1] - static.energies[0])        # E_{D_up} - E_{D_dn}
            spin = float(abs(np.real(static.observables["impurity_spin_z"][0])))   # lower level D_dn
            transverse = float(abs(static.matrix_elements("spin_x")[0, 1]))       # |<D_dn|S^x_d|D_up>|
            periodic = transverse_drive(gamma)
            rabi, omega_r, curvature = resonance(periodic, static, zeeman, DRIVE*transverse)
            cache[gamma] = dict(static=static, periodic=periodic, zeeman=zeeman, spin=spin, transverse=transverse,
                                rabi=rabi, omega_r=omega_r, curvature=curvature)
        return cache[gamma]
    main_point = doublet(GAMMA)
    static, periodic, zeeman, spin, transverse = (main_point[key] for key in
                                                  ("static", "periodic", "zeeman", "spin", "transverse"))
    rabi, omega_r = main_point["rabi"], main_point["omega_r"]
    rabi_rwa = DRIVE*transverse                 # RWA Rabi splitting b |<D_up|S^x_d|D_dn>| = (g_eff/g) b/2
    print(f"Gamma = {GAMMA}: Z = {zeeman:.6f}; g_eff/g from 2|<S^z_d>| = {2*spin:.6f}, Z/B = {zeeman/FIELD:.6f}, "
          f"2|<D_dn|S^x_d|D_up>| = {2*transverse:.6f}")
    print(f"Rabi splitting: RWA {rabi_rwa:.6g}, Floquet {rabi:.6g}; resonance shift Omega_r - Z = "
          f"{omega_r - zeeman:.4g} (delta^2/4Z = {rabi**2/(4*zeeman):.4g})")

    # 2. Probes of the driven doublet: versus drive frequency, and at Omega = Z.
    sweep = []
    for ratio in settings["sweep"]:
        omega = float(ratio*zeeman)
        states = driven(periodic, omega, static)
        lines = states.matrix_elements("spin_x").lines(initial=states.dressed_state(0))
        keep = abs(abs(lines.frequencies) - omega) < zeeman/2
        sweep.append(dict(omega=omega, frequencies=lines.frequencies[keep].tolist(),
                          strengths=lines.strengths[keep].tolist()))
    at_zeeman = driven(periodic, zeeman, static)
    spectra = probes(at_zeeman, zeeman)
    lower, elastic, upper = mollow(spectra["spin_x"], zeeman, splitting(at_zeeman))
    (sz_frequency,), (sz_line,) = spectra["spin_z"]["frequencies"], spectra["spin_z"]["strengths"]
    fock_residual = float(at_zeeman.fock_residuals[doublet_states(at_zeeman)].max())

    # 3. The effective g-factor versus hybridization: three independent determinations.
    scan = []
    for gamma in settings["gammas"]:
        point = doublet(float(gamma))
        scan.append(dict(gamma=float(gamma), spin=2*point["spin"], zeeman=point["zeeman"]/FIELD,
                         rabi=point["rabi"]/(DRIVE/2)))
        print(f"Gamma = {gamma:.2f}: 2|<S^z_d>| = {scan[-1]['spin']:.5f}, Z/B = {scan[-1]['zeeman']:.5f}, "
              f"delta/(b/2) = {scan[-1]['rabi']:.5f}")
    spread = max(max(r["spin"], r["zeeman"], r["rabi"]) - min(r["spin"], r["zeeman"], r["rabi"]) for r in scan)

    values = dict(
        zeeman=zeeman, geff=2*spin, geff_zeeman=zeeman/FIELD, geff_transverse=2*transverse,
        rabi_rwa=rabi_rwa, rabi=rabi, geff_rabi=rabi/(DRIVE/2), rabi_reduction=1 - rabi/rabi_rwa,
        curvature=main_point["curvature"], bloch_siegert=omega_r - zeeman,
        bloch_siegert_estimate=rabi**2/(4*zeeman),
        sx_lower=lower, sx_elastic=elastic, sx_upper=upper, sx_rwa=transverse**2/4,
        sx_asymmetry=(upper - lower)/(upper + lower), rabi_over_zeeman=rabi/zeeman,
        sz_frequency=sz_frequency, sz_line=sz_line, sz_rwa=transverse**2,
        geff_last=scan[-1]["spin"], geff_spread=spread, fock_residual=fock_residual,
    )
    for key, value in values.items():
        print(f"{key:24s} {value:.6g}")
    record = dict(tutorial=HERE.name, profile=args.profile, values=values, field=FIELD, drive=DRIVE,
                  omega_resonance=omega_r, sweep=sweep, spectra=spectra, scan=scan)
    (output/"results.json").write_text(json.dumps(record, indent=1) + "\n")
    plot(record, output/"spin_resonance.png")
    print(f"Saved {output}")


def plot(record, path):
    """Four panels: S^x lines versus drive frequency, S^x and S^z spectra at Omega = Z, g_eff/g versus Gamma."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping the figure")
        return
    from qdjj_solver import line_shape
    values = record["values"]
    zeeman, rabi = values["zeeman"], values["rabi"]
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 6), layout="constrained")
    axes = axes.ravel()

    x, y, s = [], [], []
    for row in record["sweep"]:
        x += [row["omega"]/zeeman]*len(row["frequencies"])
        y += [abs(f)/zeeman for f in row["frequencies"]]
        s += row["strengths"]
    s = np.array(s)
    order = np.argsort(s)
    axes[0].scatter(np.array(x)[order], np.array(y)[order], c=s[order]/s.max(), cmap="Blues", vmin=0, vmax=1,
                    s=8, lw=0)
    detuning = np.linspace(-.2, .2, 200)*zeeman
    for sign in (-1, 1):
        axes[0].plot(1 + detuning/zeeman, 1 + (detuning + sign*np.hypot(detuning, rabi))/zeeman, "--",
                     color="#eb6834", lw=.8)
    axes[0].set(xlabel=r"drive frequency $\Omega/Z$", ylabel=r"$S^x_d$ probe frequency $|\omega|/Z$", ylim=(.7, 1.3))

    def rate(grid, lines):
        f, st = np.array(lines["frequencies"]), np.array(lines["strengths"])
        return sum(a*(line_shape(grid - b, WIDTH*zeeman) + line_shape(-grid - b, WIDTH*zeeman))
                   for b, a in zip(f, st, strict=True))
    spectra = record["spectra"]
    grid = np.linspace(.8, 1.2, 801)*zeeman
    axes[1].plot(grid/zeeman, rate(grid, spectra["spin_x"]), color="#2a78d6")
    for value in (zeeman - rabi, zeeman, zeeman + rabi):
        axes[1].axvline(value/zeeman, color="gray", lw=.6, ls=":", zorder=0)
    axes[1].set(xlabel=r"probe frequency $\omega/Z$", ylabel=r"$S^x_d$ probe (rate), $\Omega=Z$")

    grid = np.linspace(0., .3, 601)*zeeman
    axes[2].plot(grid/zeeman, rate(grid, spectra["spin_z"]), color="#2a78d6")
    axes[2].axvline(rabi/zeeman, color="gray", lw=.6, ls=":", zorder=0)
    axes[2].set(xlabel=r"probe frequency $\omega/Z$", ylabel=r"$S^z_d$ probe (rate), $\Omega=Z$")

    scan = record["scan"]
    gamma = [r["gamma"] for r in scan]
    axes[3].plot(gamma, [r["spin"] for r in scan], color="gray", label=r"$2\vert\langle S^z_d\rangle\vert$ (static)")
    axes[3].plot(gamma, [r["zeeman"] for r in scan], "o", ms=5, mfc="none", color="#2a78d6",
                 label=r"Zeeman splitting $/B$")
    axes[3].plot(gamma, [r["rabi"] for r in scan], "s", ms=3.5, mfc="none", color="#eb6834",
                 label=r"Floquet Rabi splitting $/(b/2)$")
    axes[3].set(xlabel=r"hybridization $\Gamma/\Delta$", ylabel=r"$g_{\mathrm{eff}}/g$", ylim=(.5, 1.05))
    axes[3].legend(frameon=False, fontsize=8, loc="lower left")
    for axis, label in zip(axes, "abcd", strict=True):
        axis.set_title(f"({label})", loc="left", fontsize=10)
    fig.savefig(path, dpi=100)
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
