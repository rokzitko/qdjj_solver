"""Named probe operators for matrix elements and spectra, in canonical coordinates.

``current`` is the phase derivative dH/dphi of the model (the supercurrent is
(2e/hbar) times its expectation value), and ``charge`` and ``spin_z`` are the
impurity observables. The spin-flip and single-electron probes are physical
electron operators of one impurity orbital, transformed with the model's
electron-to-canonical map.
"""

from __future__ import annotations

from .algebra import annihilate, create, is_integer

OBSERVABLE_PROBES = {"current": "phase_derivative", "charge": "impurity_charge", "spin_z": "impurity_spin_z"}
ELECTRON_PROBES = ("spin_x", "spin_y", "spin_plus", "spin_minus", "d_up", "d_dn", "d_up_dagger", "d_dn_dagger",
                   "charge", "spin_z")
PROBES = ("current", "charge", "spin_z", "spin_x", "spin_y", "spin_plus", "spin_minus", "d_up", "d_dn",
          "d_up_dagger", "d_dn_dagger")
# Probe name -> (operator, typical use). GitHub Markdown with $...$ mathematics; rendered by
# tools/gen_docs.py in docs/reference/observables.md.
PROBE_DOCS = {
    "current": ("$\\partial H/\\partial\\phi$ (observable `phase_derivative`), with its drive harmonics",
                "Flux or phase probe; supercurrent"),
    "charge": ("$n_d$ (observable `impurity_charge`)", "Gate probe"),
    "spin_z": ("$S^z_d$ (observable `impurity_spin_z`)", "Longitudinal spin probe"),
    "spin_x": ("$S^x_d=(S^+_d+S^-_d)/2$", "Transverse spin resonance"),
    "spin_y": ("$S^y_d=(S^+_d-S^-_d)/(2\\mathrm{i})$", "Transverse spin resonance"),
    "spin_plus": ("$S^+_d=d^\\dagger_\\uparrow d_\\downarrow$", "Spin raising"),
    "spin_minus": ("$S^-_d=d^\\dagger_\\downarrow d_\\uparrow$", "Spin lowering"),
    "d_up": ("$d_\\uparrow$", "Electron removal (tunnelling spectroscopy)"),
    "d_dn": ("$d_\\downarrow$", "Electron removal (tunnelling spectroscopy)"),
    "d_up_dagger": ("$d^\\dagger_\\uparrow$", "Electron addition (tunnelling spectroscopy)"),
    "d_dn_dagger": ("$d^\\dagger_\\downarrow$", "Electron addition (tunnelling spectroscopy)"),
}


def electron_probe(name, orbital=0):
    """The physical electron operator of a built-in probe for one impurity orbital."""
    up, down = 2*orbital, 2*orbital + 1
    plus = create(up)*annihilate(down)
    minus = create(down)*annihilate(up)
    operators = {"d_up": lambda: annihilate(up), "d_dn": lambda: annihilate(down),
                 "d_up_dagger": lambda: create(up), "d_dn_dagger": lambda: create(down),
                 "spin_plus": lambda: plus, "spin_minus": lambda: minus,
                 "spin_x": lambda: (plus + minus)*.5, "spin_y": lambda: (plus - minus)*(-.5j),
                 "spin_z": lambda: (create(up)*annihilate(up) - create(down)*annihilate(down))*.5,
                 "charge": lambda: create(up)*annihilate(up) + create(down)*annihilate(down)}
    if name not in operators:
        raise ValueError(f"{name!r} is not an electron probe")
    return operators[name]()


def _check(name, orbital, nimp):
    if not isinstance(name, str):
        raise ValueError("a probe name must be a string")
    if orbital is not None and (not is_integer(orbital) or not 0 <= orbital < nimp//2):
        raise ValueError(f"orbital must be an impurity orbital index between 0 and {nimp//2 - 1}")


def probe_operator(hamiltonian, name, orbital=None):
    """A probe as a FermionOperator in the canonical coordinates of ``hamiltonian``.

    Built-in names are listed in ``PROBES``; any other name selects a model
    observable. ``current``, and ``charge`` and ``spin_z`` without ``orbital``,
    are model observables. The spin-flip and single-electron probes, and
    ``charge`` and ``spin_z`` with ``orbital``, act on one impurity orbital
    (default 0) through the model's electron map, which requires an uncompressed
    model.
    """
    _check(name, orbital, hamiltonian.nimp)
    if name == "current" and orbital is not None:
        raise ValueError("the current probe is a property of the junction, not of an orbital")
    if name in OBSERVABLE_PROBES and orbital is None:
        key = OBSERVABLE_PROBES[name]
        if key not in hamiltonian.observables:
            raise ValueError(f"probe {name!r} needs the model observable {key!r}"
                             + (" (set phase_velocities)" if name == "current" else ""))
        return hamiltonian.observables[key]
    if name in ELECTRON_PROBES:
        return hamiltonian.physical_operator(electron_probe(name, orbital or 0))
    if name in hamiltonian.observables and orbital is None:
        return hamiltonian.observables[name]
    raise ValueError(f"unknown probe {name!r}; built-in probes: {', '.join(PROBES)}")


def periodic_probe(periodic, name, orbital=None):
    """Harmonics {k: operator} of a probe under the drive of a PeriodicHamiltonian.

    Observables depend on the driven parameters (the current under a phase drive,
    for instance) and use their sampled harmonics; electron operators are time
    independent in the fixed canonical coordinates.
    """
    _check(name, orbital, periodic.nimp)
    key = OBSERVABLE_PROBES.get(name, name) if orbital is None else None
    if key is not None and key in periodic.reference.observables:
        if key not in periodic.observables:
            raise ValueError(f"observable {key!r} was not sampled; build the PeriodicHamiltonian with observables=True")
        return {k: v for k, v in periodic.observables[key].items() if v.terms}
    return {0: probe_operator(periodic.reference, name, orbital)}


__all__ = ["PROBES", "electron_probe", "periodic_probe", "probe_operator"]
