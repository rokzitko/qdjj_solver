"""Static subgap spectrum of the manual's reference junction (chapter on the model)."""

import numpy as np

import manual_common as mc


def compute():
    from qdjj_solver import Sector, reference_model
    from qdjj_solver.qp_solver import SolverOptions, solve
    results = {}
    scans = dict(phi=np.linspace(0., 2*np.pi, 61), detuning=np.linspace(-1.2, 1.2, 49))
    for name, values in scans.items():
        rows = []
        for value in values:
            h = reference_model(**mc.junction(**{name: float(value)}))
            even = solve(h, 2, Sector(0, 0), SolverOptions(eigenpairs=2))
            odd = solve(h, 2, Sector(1, 1), SolverOptions(eigenpairs=1))
            rows.append(dict(value=float(value), singlet=float(even.energies[0]), excited=float(even.energies[1]),
                             excited_qp0=float(even.qp_weights[1, 0]), doublet=float(odd.energies[0]),
                             current=float(np.real(even.observables["phase_derivative"][0])),
                             doublet_current=float(np.real(odd.observables["phase_derivative"][0]))))
        results[name] = rows
    mc.save("static", dict(parameters={k: v for k, v in mc.junction().items() if k != "bath"},
                           bath=mc.BATH, cutoff=2, scans=results), __file__)


def plot():
    record = mc.load("static")
    fig, axes = mc.figure(1, 2, height=2.4)
    for axis, (name, label) in zip(axes, (("phi", r"phase $\phi/\pi$"), ("detuning", r"detuning $ϵ/\Delta$")),
                                   strict=True):
        rows = record["scans"][name]
        x = np.array([r["value"] for r in rows])/(np.pi if name == "phi" else 1.)
        singlet, excited, doublet = (np.array([r[k] for r in rows]) for k in ("singlet", "excited", "doublet"))
        ground = np.minimum(singlet, doublet)
        # Lowest threshold across parities; a fixed-sector probe can have a higher threshold.
        axis.axhspan(1., 2.2, color="#f0efec", lw=0, zorder=0)
        axis.plot(x, singlet-ground, color=mc.BLUE, label=r"singlet $S$")
        axis.plot(x, doublet-ground, color=mc.ORANGE, label=r"doublet $D$")
        axis.plot(x, excited-ground, color=mc.BLUE, ls="--", label=r"excited singlet $S^*$")
        axis.set_xlabel(label)
        axis.set_ylim(-.05, 2.0)
        axis.set_xlim(x[0], x[-1])
        axis.text(.03, .93, "lowest continuum threshold", transform=axis.transAxes, color=mc.SECONDARY, fontsize=7.5, va="top")
    axes[0].set_ylabel(r"$E-E_{\mathrm{ground}}$  ($\Delta$)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncols=3)
    mc.panel_label(axes[0], "a")
    mc.panel_label(axes[1], "b")
    mc.finish(fig, "static")


WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")


def quoted(record):
    """Reference parameters and the excitation energies above the singlet at the reference point."""
    parameters, bath = record["parameters"], record["bath"]
    # The detuning scan runs at the reference phase and passes through the reference detuning.
    row = min(record["scans"]["detuning"], key=lambda r: abs(r["value"] - parameters["detuning"]))
    assert abs(row["value"] - parameters["detuning"]) < 1e-9, "reference detuning missing from the scan"
    return dict(u=(parameters["u"], "g"), gamma=(parameters["gamma"], "g"), detuning=(parameters["detuning"], "g"),
                phi=(parameters["phi"], "g"), bandwidth=(bath["bandwidth"], "g"), delta=(bath["delta"], "g"),
                pairs=(rf"\text{{{WORDS[bath['pairs']]}}}", None), cutoff=(record["cutoff"], "d"),
                doublet_energy=(row["doublet"] - row["singlet"], ".3f"),
                excited_energy=(row["excited"] - row["singlet"], ".3f"))


if __name__ == "__main__":
    mc.main(compute, plot)
