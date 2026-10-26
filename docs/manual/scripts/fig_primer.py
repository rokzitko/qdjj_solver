"""Introductory analytic plots and original schematics; compute uses no impurity solver.

The tracked primer record stores analytic curves and their assumptions. The current
panel reuses static.json, preserving its finite-bath provenance. Schematic geometry
is defined here and has no quantitative interpretation. plot never recomputes curves.
"""

import numpy as np

import manual_common as mc


def contact(phi, transmission):
    """Positive ABS energy / Delta and even-ground current / (e Delta / hbar).

    Equal-gap, spin-degenerate short channel; exclude the ambiguous derivative
    at the perfectly transmitting zero-energy crossing from current plots.
    """
    phi = np.asarray(phi)
    energy = np.sqrt(1. - transmission*np.sin(phi/2.)**2)
    current = np.divide(transmission*np.sin(phi), 2.*energy,
                        out=np.full_like(energy, np.nan), where=energy > 1e-12)
    return energy, current


def atomic(detuning, u, pairing):
    """Even energies and lower-even charge in the manual's shifted energy zero."""
    detuning = np.asarray(detuning)
    radius = np.sqrt(detuning**2 + abs(pairing)**2)
    charge = 1. - np.divide(detuning, radius, out=np.zeros_like(radius), where=radius > 0.)
    # Charge at the fully degenerate even block is basis dependent, not uniquely one.
    charge = np.where(radius > 0., charge, np.nan)
    return u/2. - radius, u/2. + radius, charge


def compute():
    phi = np.linspace(0., 2.*np.pi, 161)
    transmissions = [.2, .7, .95]
    channels = []
    for transmission in transmissions:
        energy, current = contact(phi, transmission)
        channels.append(dict(transmission=transmission, energy=energy, current=current))
    detuning = np.linspace(-1.5, 1.5, 161)
    u, pairing = 2., .6
    singlet, excited, charge = atomic(detuning, u, pairing)
    alpha = np.unique(np.r_[np.linspace(0., 3., 161), 1.])  # resolve the exact parity crossing
    thermal_x = np.linspace(-6., 6., 161)
    energy = np.linspace(-3., 3., 601)
    dos = np.zeros_like(energy)
    outside = abs(energy) > 1.001
    dos[outside] = abs(energy[outside])/np.sqrt(energy[outside]**2 - 1.)
    distance = np.linspace(0., 5., 101)
    envelopes = [dict(energy=e, amplitude=np.exp(-distance*np.sqrt(1.-e**2))) for e in (.2, .9)]
    mc.save("primer", dict(
        kind="analytic teaching examples; schematic drawings are not experimental data",
        assumptions=dict(contact="short, spin degenerate, equal gaps, noninteracting, even ground state",
                         atomic="wide band before infinite gap; zero field; shifted dot energy zero",
                         ysr="classical spin; flat wide band; no potential scattering; alpha=pi*rho0*J*S/2",
                         thermal="one singlet and one doublet; equilibrium between both parities"),
        contact=dict(phi=phi, channels=channels),
        atomic=dict(u=u, pairing=pairing, detuning=detuning, singlet=singlet, excited=excited,
                    singlet_charge=charge, ground_charge=np.where(singlet < 0., charge, 1.)),
        ysr=dict(alpha=alpha, signed_energy=(1.-alpha**2)/(1.+alpha**2)),
        thermal=dict(x=thermal_x, doublet_probability=2./(2.+np.exp(thermal_x))),
        gap=dict(energy=energy, dos=dos, distance=distance, envelopes=envelopes),
    ), __file__, compact=True)


def schematic(columns=1, height=2.5):
    fig, axes = mc.figure(1, columns, height=height, squeeze=False)
    for axis in axes.flat:
        axis.set(xlim=(0., 1.), ylim=(0., 1.))
        axis.axis("off")
    return fig, axes[0]


def node(axis, xy, text, width=.28, height=.18, color=mc.BLUE, fontsize=8.5):
    from matplotlib.patches import FancyBboxPatch
    x, y = xy
    axis.add_patch(FancyBboxPatch((x-width/2., y-height/2.), width, height,
                                boxstyle="round,pad=0.015,rounding_size=0.02",
                                edgecolor=color, facecolor=color+"15", linewidth=1.))
    axis.text(x, y, text, ha="center", va="center", color=mc.INK, fontsize=fontsize)


def arrow(axis, start, end, color=mc.SECONDARY, both=False, style="-"):
    axis.annotate("", xy=end, xytext=start,
                  arrowprops=dict(arrowstyle="<->" if both else "->", color=color,
                                  lw=1.3, linestyle=style, shrinkA=0., shrinkB=0.))


def level(axis, y, label, x=.35, width=.35, color=mc.BLUE):
    axis.plot([x-width/2., x+width/2.], [y, y], color=color, lw=2.)
    axis.text(x+width/2.+.035, y, label, va="center", fontsize=8.5, color=color)


def save(fig, name):
    import matplotlib.pyplot as plt
    mc.finish(fig, "primer_"+name)
    plt.close(fig)


def plot_devices():
    fig, axes = schematic(3, height=2.15)
    for axis, title in zip(axes, ("Short contact / QPC", "Quantum-dot junction", "Magnetic adatom"), strict=True):
        axis.set_title(title)
    for axis in axes[:2]:
        node(axis, (.14, .62), "S", width=.20, height=.30)
        node(axis, (.86, .62), "S", width=.20, height=.30)
    axes[0].plot([.25, .4, .6, .75], [.74, .66, .66, .74], color=mc.BLUE, lw=2.)
    axes[0].plot([.25, .4, .6, .75], [.50, .58, .58, .50], color=mc.BLUE, lw=2.)
    axes[0].text(.5, .35, r"channels $\tau_n$"+"\nphase difference "+r"$\phi$", ha="center")
    node(axes[1], (.5, .62), "dot", width=.23, color=mc.ORANGE)
    for x in (.31, .69):
        axes[1].plot([x, x], [.48, .76], color=mc.SECONDARY, lw=2.)
    arrow(axes[1], (.5, .20), (.5, .48), color=mc.ORANGE)
    axes[1].text(.5, .10, r"gate: $ϵ$, barriers: $\Gamma_l$", ha="center", fontsize=8.)
    node(axes[2], (.5, .45), "superconducting surface", width=.85, height=.20, fontsize=7.5)
    axes[2].scatter([.5], [.66], s=170, color=mc.ORANGE)
    arrow(axes[2], (.5, .64), (.5, .82), color=mc.INK)
    arrow(axes[2], (.78, .94), (.57, .73))
    axes[2].text(.5, .16, "local spin; STM probe", ha="center", fontsize=8.)
    save(fig, "devices")


def plot_gap(record):
    data = record["gap"]
    fig, axes = mc.figure(1, 3, height=2.15)
    axes[0].plot(data["energy"], data["dos"], color=mc.BLUE)
    axes[0].axvspan(-1., 1., color=mc.GRID, alpha=.5)
    axes[0].set(xlabel=r"$E/\Delta$", ylabel=r"$\rho/\rho_0$", ylim=(0., 4.), title="Bulk DOS")
    axes[1].axvspan(-1., 1., color=mc.GRID, alpha=.5)
    axes[1].vlines([-.45, .45], [0., 0.], [.5, 1.], colors=[mc.AQUA, mc.ORANGE], lw=2.)
    axes[1].set(xlabel=r"$E/\Delta$", xlim=(-1.4, 1.4), ylim=(0., 1.2), yticks=[], title="Local poles")
    axes[1].text(-.45, .55, "remove", ha="center", fontsize=7.)
    axes[1].text(.45, 1.05, "add", ha="center", fontsize=7.)
    for envelope, color in zip(data["envelopes"], (mc.BLUE, mc.ORANGE), strict=True):
        axes[2].plot(data["distance"], envelope["amplitude"], color=color,
                     label=rf"$E/\Delta={envelope['energy']}$")
    axes[2].set(xlabel=r"$x\Delta/(\hbar v_F)$", ylabel="envelope", title="Decay into lead")
    axes[2].legend(loc="upper right")
    for axis, label in zip(axes, "abc", strict=True):
        mc.panel_label(axis, label)
    save(fig, "gap")


def plot_reflection():
    fig, axes = schematic(2, height=2.3)
    axes[0].set_title("One interface")
    node(axes[0], (.80, .54), "S\ncondensate", width=.32, height=.60)
    axes[0].text(.24, .91, "normal region", ha="center")
    arrow(axes[0], (.06, .68), (.61, .68), color=mc.BLUE)
    arrow(axes[0], (.61, .39), (.06, .39), color=mc.AQUA)
    axes[0].text(.28, .75, "electron", color=mc.BLUE, ha="center")
    axes[0].text(.28, .27, "hole", color=mc.AQUA, ha="center")
    axes[0].text(.5, .08, "pair transferred into S", ha="center")
    axes[1].set_title("Two interfaces: a bound orbit")
    node(axes[1], (.10, .54), r"$S_L$", width=.15, height=.60)
    node(axes[1], (.90, .54), r"$S_R$", width=.15, height=.60)
    arrow(axes[1], (.23, .70), (.77, .70), color=mc.BLUE)
    arrow(axes[1], (.77, .37), (.23, .37), color=mc.AQUA)
    arrow(axes[1], (.78, .67), (.78, .40), color=mc.ORANGE)
    arrow(axes[1], (.22, .40), (.22, .67), color=mc.ORANGE)
    axes[1].text(.5, .54, "coherent phase\naccumulation", ha="center", va="center")
    axes[1].text(.5, .08, r"$\phi=\phi_R-\phi_L$", ha="center")
    save(fig, "reflection")


def plot_occupations():
    fig, axes = schematic(2, height=2.5)
    for axis in axes:
        level(axis, .16, r"$|g\rangle$  even", x=.24, width=.28)
        level(axis, .50, r"$|o_\uparrow\rangle, |o_\downarrow\rangle$", x=.24, width=.28, color=mc.ORANGE)
        level(axis, .84, r"$|e\rangle$  even", x=.24, width=.28)
    axes[0].set_title("Single-electron probe")
    arrow(axes[0], (.23, .19), (.23, .47), color=mc.ORANGE, both=True)
    axes[0].text(.70, .30, r"cost $E_A$", ha="center", color=mc.ORANGE)
    axes[1].set_title("Pair excitation")
    arrow(axes[1], (.23, .19), (.23, .81), color=mc.BLUE, both=True)
    axes[1].text(.75, .66, r"cost $2E_A$", ha="center", color=mc.BLUE)
    for axis in axes:
        axis.text(.5, .01, "quasiparticle occupations: 0, 1, 2", ha="center", fontsize=7.5)
    save(fig, "occupations")


def plot_contact(record):
    data = record["contact"]
    fig, axes = mc.figure(1, 2, height=2.4)
    for channel, color in zip(data["channels"], (mc.BLUE, mc.AQUA, mc.ORANGE), strict=True):
        x = np.asarray(data["phi"])/np.pi
        label = rf"$\tau={channel['transmission']}$"
        axes[0].plot(x, channel["energy"], color=color, label=label)
        axes[1].plot(x, channel["current"], color=color, label=label)
    axes[0].set(ylabel=r"$E_A/\Delta$", ylim=(0., 1.05))
    axes[1].set(ylabel=r"$I_g/(e\Delta/\hbar)$")
    for axis in axes:
        axis.set(xlabel=r"$\phi/\pi$", xlim=(0., 2.))
    axes[0].legend(loc="lower left", ncols=1)
    save(fig, "contact")


def plot_dot():
    fig, (axis,) = schematic(height=2.35)
    node(axis, (.15, .75), r"$|0\rangle$"+"\n"+r"$U/2-ϵ$", width=.22, height=.23)
    node(axis, (.85, .75), r"$|2\rangle$"+"\n"+r"$U/2+ϵ$", width=.22, height=.23)
    node(axis, (.5, .22), r"$|\uparrow\rangle,|\downarrow\rangle$"+"\n"+r"$\pm B/2$",
         width=.28, height=.23, color=mc.ORANGE)
    arrow(axis, (.29, .75), (.71, .75), color=mc.BLUE, both=True)
    axis.text(.5, .86, r"pair mixing $\Gamma_\phi$", ha="center", color=mc.BLUE)
    arrow(axis, (.23, .60), (.37, .35), color=mc.ORANGE, both=True)
    arrow(axis, (.77, .60), (.63, .35), color=mc.ORANGE, both=True)
    axis.text(.16, .35, "one electron", ha="center", color=mc.ORANGE)
    axis.text(.84, .35, "one electron", ha="center", color=mc.ORANGE)
    save(fig, "dot")


def plot_atomic(record):
    data = record["atomic"]
    x = data["detuning"]
    fig, axes = mc.figure(1, 2, height=2.45)
    axes[0].plot(x, data["singlet"], color=mc.BLUE, label=r"$S$")
    axes[0].plot(x, data["excited"], color=mc.BLUE, ls="--", label=r"$S^*$")
    axes[0].axhline(0., color=mc.ORANGE, label=r"$D$")
    axes[0].set(ylabel="energy (illustrative unit)")
    axes[0].legend(ncols=3, loc="upper center")
    axes[1].plot(x, data["singlet_charge"], color=mc.BLUE, ls="--", label="lower even state")
    ground = np.asarray(data["ground_charge"])
    # Do not draw a false interpolating state across a first-order charge jump.
    jump = np.r_[False, np.abs(np.diff(ground)) > .1]
    axes[1].plot(x, np.where(jump, np.nan, ground), color=mc.INK, label="global ground state")
    axes[1].set(ylabel=r"$\langle n_d\rangle$", ylim=(-.1, 2.1))
    axes[1].legend(loc="upper right")
    for axis, label in zip(axes, "ab", strict=True):
        axis.set(xlabel="detuning ϵ (illustrative unit)", xlim=(x[0], x[-1]))
        mc.panel_label(axis, label)
    save(fig, "atomic")


def plot_current():
    data = mc.load("static")
    rows = data["scans"]["phi"]
    x = np.array([r["value"] for r in rows])/np.pi
    singlet_ground = np.array([r["singlet"] < r["doublet"] for r in rows])
    fig, axis = mc.figure(height=2.1)
    for key, mask, color, label in (("current", singlet_ground, mc.BLUE, "singlet"),
                                    ("doublet_current", ~singlet_ground, mc.ORANGE, "doublet")):
        current = np.array([r[key] for r in rows])
        axis.plot(x, current, color=color, ls="--", lw=1., label=label+" branch")
        axis.plot(x, np.where(mask, current, np.nan), color=color, lw=2.5)
    axis.set(xlabel=r"$\phi/\pi$", ylabel=r"$I/(2e\Delta/\hbar)$", xlim=(0., 2.))
    axis.legend(loc="upper right", title="thick: global ground state", fontsize=8.)
    save(fig, "current")


def plot_screening():
    fig, axes = schematic(3, height=2.5)
    titles = ("Charge screening", "Proximity pairing", "Kondo screening")
    for axis, title in zip(axes, titles, strict=True):
        axis.set_title(title)
    node(axes[0], (.5, .58), "+", width=.18, height=.18, color=mc.ORANGE, fontsize=14.)
    for theta in np.linspace(0., 2.*np.pi, 6, endpoint=False):
        axes[0].text(.5+.32*np.cos(theta), .58+.25*np.sin(theta), "−", color=mc.BLUE,
                     ha="center", va="center", fontsize=15.)
    axes[0].text(.5, .10, "charge redistribution", ha="center", fontsize=8.)
    node(axes[1], (.23, .60), r"$|0\rangle$", width=.22)
    node(axes[1], (.77, .60), r"$|2\rangle$", width=.22)
    arrow(axes[1], (.38, .60), (.62, .60), color=mc.BLUE, both=True)
    axes[1].text(.5, .30, r"$u|0\rangle+v|2\rangle$", ha="center")
    axes[1].text(.5, .10, "local even singlet", ha="center", fontsize=8.)
    node(axes[2], (.5, .76), "dot spin", width=.48, color=mc.ORANGE)
    node(axes[2], (.5, .36), "many-electron\nbath state", width=.75, height=.22)
    arrow(axes[2], (.5, .65), (.5, .49), color=mc.AQUA, both=True)
    axes[2].text(.5, .10, "entangled total singlet", ha="center", fontsize=8.)
    save(fig, "screening")


def plot_scales():
    fig, axes = schematic(2, height=2.5)
    for axis, title, tk, gap in zip(axes, (r"$k_BT_K\gg\Delta$", r"$k_BT_K\ll\Delta$"),
                                     (.62, .25), (.25, .62), strict=True):
        axis.set_title(title)
        arrow(axis, (.16, .90), (.16, .10))
        axis.text(.09, .5, "lower energy", rotation=90, va="center", ha="center")
        axis.axhspan(.08, gap, xmin=.25, xmax=.95, color=mc.GRID, alpha=.6)
        level(axis, gap, r"$\Delta$", x=.47, width=.38, color=mc.BLUE)
        level(axis, tk, r"$k_BT_K$", x=.47, width=.38, color=mc.ORANGE)
        axis.text(.5, .96, "normal-state high-energy scale", ha="center", fontsize=7.5)
    save(fig, "scales")


def plot_ysr(record):
    data = record["ysr"]
    x, signed = np.asarray(data["alpha"]), np.asarray(data["signed_energy"])
    fig, axes = mc.figure(1, 2, height=2.45)
    axes[0].plot(x, signed, color=mc.BLUE, label="signed branch")
    axes[0].plot(x, -signed, color=mc.ORANGE, ls="--", label="partner")
    axes[0].set(ylabel=r"$E/\Delta$", ylim=(-1.05, 1.05))
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside upper center", ncols=2)
    axes[1].plot(x, abs(signed), color=mc.BLUE)
    axes[1].set(ylabel=r"excitation cost $|E|/\Delta$", ylim=(-.05, 1.05))
    for axis, label in zip(axes, "ab", strict=True):
        axis.axvline(1., color=mc.MUTED, ls=":")
        axis.set(xlabel=r"exchange strength $\alpha$", xlim=(0., 3.))
        mc.panel_label(axis, label)
    save(fig, "ysr")


def plot_thermal(record):
    data = record["thermal"]
    x = np.asarray(data["x"])
    fig, axes = mc.figure(1, 2, height=2.2)
    axes[0].plot(x, x, color=mc.ORANGE, label=r"$E_D-E_S$")
    axes[0].axhline(0., color=mc.BLUE, label="singlet reference")
    axes[0].set(ylabel=r"energy / $k_BT$")
    axes[0].legend(loc="upper left", fontsize=7.)
    axes[1].plot(x, data["doublet_probability"], color=mc.ORANGE)
    axes[1].plot([0.], [2./3.], "o", color=mc.INK, ms=4.)
    axes[1].annotate(r"$2/3$", (0., 2./3.), xytext=(1.2, .80),
                     arrowprops=dict(arrowstyle="->", color=mc.SECONDARY))
    axes[1].set(ylabel=r"doublet population $P_D$", ylim=(0., 1.05))
    for axis in axes:
        axis.set(xlabel=r"$\delta_{SD}/(k_BT)$", xlim=(-6., 6.))
    save(fig, "thermal")


def plot_probes():
    fig, axes = schematic(2, height=2.5)
    axes[0].set_title("Electron tunnelling")
    axes[1].set_title("Microwave coupling")
    node(axes[0], (.18, .65), "normal\nprobe", width=.28, color=mc.ORANGE)
    node(axes[0], (.77, .65), "dot + S", width=.35)
    arrow(axes[0], (.37, .65), (.55, .65), color=mc.ORANGE, both=True)
    axes[0].text(.5, .38, r"$d_\sigma,\ d_\sigma^\dagger$"+"\nchange subsystem parity", ha="center", va="center")
    axes[0].text(.5, .08, r"detect $dI/dV$", ha="center")
    node(axes[1], (.18, .65), "ac gate\nor flux", width=.28, color=mc.AQUA)
    node(axes[1], (.77, .65), "junction", width=.35)
    arrow(axes[1], (.37, .65), (.55, .65), color=mc.AQUA)
    axes[1].text(.5, .38, r"$\partial H/\partial\lambda$"+"\npreserves total parity", ha="center", va="center")
    axes[1].text(.5, .08, "detect absorption / resonator shift", ha="center", fontsize=8.)
    save(fig, "probes")


def plot_spectroscopy(record):
    data = record["atomic"]
    x, signed = np.asarray(data["detuning"]), -np.asarray(data["singlet"])
    fig, axes = mc.figure(1, 2, height=2.25)
    axes[0].plot(x, signed, color=mc.ORANGE)
    axes[0].axhline(0., color=mc.MUTED, lw=.8)
    axes[0].set(ylabel=r"$E_D-E_S$", title="Signed energy difference")
    for sign, color in ((1., mc.ORANGE), (-1., mc.AQUA)):
        axes[1].plot(x, sign*abs(signed), color=color)
    axes[1].set(ylabel="probe energy eV", title="Lowest transition positions")
    for axis, label in zip(axes, "ab", strict=True):
        axis.set(xlabel="detuning ϵ (illustrative unit)", xlim=(x[0], x[-1]))
        mc.panel_label(axis, label)
    save(fig, "spectroscopy")


def plot_qubits():
    fig, axes = schematic(2, height=2.5)
    axes[0].set_title("Andreev-level / pair qubit")
    level(axes[0], .20, r"$|g\rangle$", x=.24, width=.32)
    level(axes[0], .80, r"$|e\rangle$", x=.24, width=.32)
    arrow(axes[0], (.24, .23), (.24, .77), color=mc.BLUE, both=True)
    axes[0].text(.70, .51, "even parity\n"+r"$hf=2E_A$", ha="center", va="center")
    axes[1].set_title("Andreev spin qubit")
    level(axes[1], .34, r"$|D_\downarrow\rangle$", x=.24, width=.32, color=mc.ORANGE)
    level(axes[1], .66, r"$|D_\uparrow\rangle$", x=.24, width=.32, color=mc.ORANGE)
    arrow(axes[1], (.24, .37), (.24, .63), color=mc.ORANGE, both=True)
    axes[1].text(.70, .12, "odd parity\nspin / pseudospin splitting", ha="center", fontsize=8.)
    axes[0].text(.5, .04, "non-interacting channel limit", ha="center", fontsize=7.5)
    save(fig, "qubits")


def plot_circuit():
    fig, (axis,) = schematic(height=2.9)
    # Upper and lower electrical nodes are joined by three parallel branches:
    # the two junctions (enclosing flux) and a shunt capacitor. Interrupt wires
    # at each component; a line through a component would incorrectly short it.
    for y in (.38, .85):
        axis.plot([.07, .58], [y, y], color=mc.SECONDARY)
    for x in (.26, .58):
        axis.plot([x, x], [.38, .48], color=mc.SECONDARY)
        axis.plot([x, x], [.72, .85], color=mc.SECONDARY)
    node(axis, (.26, .60), "weak\nlink", width=.15, height=.21, color=mc.ORANGE)
    node(axis, (.58, .60), "reference\njunction", width=.16, height=.21)
    axis.text(.42, .60, r"flux $\Phi$", ha="center", va="center")
    axis.plot([.07, .07], [.38, .56], color=mc.SECONDARY)
    axis.plot([.07, .07], [.64, .85], color=mc.SECONDARY)
    for y in (.56, .64):
        axis.plot([.015, .125], [y, y], color=mc.BLUE, lw=1.5)
    axis.text(.07, .28, "C", ha="center", color=mc.BLUE)
    node(axis, (.87, .60), "readout\nresonator", width=.20, height=.23, color=mc.AQUA)
    arrow(axis, (.685, .60), (.745, .60), color=mc.AQUA, both=True)
    axis.text(.83, .95, "microwave in / out", ha="center", fontsize=8.)
    arrow(axis, (.87, .87), (.87, .75), color=mc.AQUA, both=True)
    for x, text, color in ((.16, "prepare", mc.ORANGE), (.5, "drive", mc.BLUE), (.84, "read out", mc.AQUA)):
        node(axis, (x, .12), text, width=.22, height=.12, color=color)
    arrow(axis, (.29, .12), (.37, .12))
    arrow(axis, (.63, .12), (.71, .12))
    save(fig, "circuit")


def plot_methods():
    fig, axes = schematic(3, height=2.65)
    axes[0].set_title("Continuum reservoir")
    for x in np.linspace(.16, .84, 15):
        axes[0].plot([x, x], [.48, .88], color=mc.BLUE, alpha=.6)
    node(axes[0], (.5, .24), "dot", width=.24, color=mc.ORANGE)
    arrow(axes[0], (.5, .37), (.5, .45), both=True)
    axes[0].text(.5, .03, "discretization is an approximation", ha="center", fontsize=7.)
    axes[1].set_title("Finite star")
    dot = (.5, .55)
    for theta in np.linspace(0., 2.*np.pi, 5, endpoint=False):
        point = (.5+.36*np.cos(theta), .55+.30*np.sin(theta))
        axes[1].plot([dot[0], point[0]], [dot[1], point[1]], color=mc.SECONDARY, lw=1.)
        axes[1].scatter(*point, color=mc.BLUE, s=70, zorder=3)
    axes[1].scatter(*dot, color=mc.ORANGE, s=110, zorder=4)
    axes[1].text(.5, .03, "ED / QP / MPS representations", ha="center", fontsize=7.)
    axes[2].set_title("Chain coordinates")
    x = np.linspace(.10, .90, 6)
    axes[2].plot(x, np.full(6, .55), color=mc.SECONDARY)
    axes[2].scatter(x, np.full(6, .55), color=[mc.ORANGE]+[mc.BLUE]*5, s=75, zorder=3)
    axes[2].text(.5, .28, "unitary bath transformation\n≠ physical lead geometry", ha="center", fontsize=8.)
    axes[2].text(.5, .03, "NRG: logarithmic scales + truncation", ha="center", fontsize=7.)
    save(fig, "methods")


def plot():
    record = mc.load("primer")
    for draw in (plot_devices, plot_reflection, plot_occupations, plot_dot, plot_current,
                 plot_screening, plot_scales, plot_probes, plot_qubits, plot_circuit, plot_methods):
        draw()
    for draw in (plot_gap, plot_contact, plot_atomic, plot_ysr, plot_thermal, plot_spectroscopy):
        draw(record)


def quoted(record):
    return dict(atomic_u=(record["atomic"]["u"], "g"),
                atomic_pairing=(record["atomic"]["pairing"], "g"))


if __name__ == "__main__":
    mc.main(compute, plot)
