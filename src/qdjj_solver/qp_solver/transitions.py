"""Probe matrix elements between Floquet states of one drive.

With u_{a,m} the harmonics of the Floquet modes and X(t)=sum_k X_k exp(ik omega t),
X^(q)_ab = sum_m <u_{a,m-q}|w_{b,m}>,  w_{b,m} = sum_k X_k u_{b,m-k}.
The probe acts on the ket modes in their Fock basis and maps them into the Fock
basis of the bra states (the sectors may differ); the result is contracted with
the bra modes in their own representation.
"""

from __future__ import annotations

import numpy as np
from threadpoolctl import threadpool_limits

from ..common.algebra import FermionOperator, is_integer
from ..common.problem import check_transition
from ..common.spectra import FloquetTransitions
from .solver import ProjectedOperator, same_basis


def _probe_harmonics(periodic, probe, orbital):
    if isinstance(probe, str):
        harmonics, name = periodic.probe(probe, orbital), probe
    elif isinstance(probe, FermionOperator):
        if orbital is not None:
            raise ValueError("orbital applies to named probes only")
        harmonics, name = {0: probe}, "operator"
    elif isinstance(probe, dict):
        if orbital is not None:
            raise ValueError("orbital applies to named probes only")
        if not all(is_integer(k) and isinstance(v, FermionOperator) for k, v in probe.items()):
            raise ValueError("probe harmonics must map integers to FermionOperators")
        harmonics, name = probe, "harmonics"
    else:
        raise ValueError("the probe must be a name, a FermionOperator or a dict of harmonics")
    harmonics = {int(k): v for k, v in harmonics.items() if v.terms}
    if not harmonics:
        raise ValueError("the probe is zero")
    return harmonics, name


def _check_pair(bra, ket):
    if bra.periodic is not ket.periodic and bra.periodic.fingerprint() != ket.periodic.fingerprint():
        raise ValueError("the Floquet states belong to different drives (PeriodicHamiltonians)")
    def gauge(result):
        return result.metadata.get("options", {}).get("drop_scalar_harmonics")
    if gauge(bra) != gauge(ket):
        raise ValueError("the Floquet states use different scalar gauges (drop_scalar_harmonics)")
    if bra.basis.spins != ket.basis.spins or bra.basis.nimp != ket.basis.nimp:
        raise ValueError("the Floquet states use different mode labels")


def _indices(selection, count, name):
    if selection is None:
        return np.arange(count)
    indices = np.atleast_1d(np.asarray(selection, dtype=int))
    if indices.ndim != 1 or len(indices) == 0 or np.any(indices < 0) or np.any(indices >= count) \
            or len(set(indices.tolist())) != len(indices):
        raise ValueError(f"{name} must list distinct state indices")
    return indices


def floquet_matrix_elements(bra, probe, ket=None, *, orbital=None, ket_states=None, bra_states=None,
                            max_memory_gib=None) -> FloquetTransitions:
    """Probe matrix elements X^(q) from the ket states (initial) to the bra states (final) of one drive.

    See ``FloquetStates.matrix_elements`` for the arguments; ``max_memory_gib``
    defaults to the setting of the bra calculation.
    """
    from .floquet import FloquetOptions, FloquetStates
    same = ket is None or ket is bra
    ket = bra if ket is None else ket
    if not isinstance(bra, FloquetStates) or not isinstance(ket, FloquetStates):
        raise ValueError("matrix elements need FloquetStates on both sides")
    _check_pair(bra, ket)
    periodic = bra.periodic
    harmonics, name = _probe_harmonics(periodic, probe, orbital)
    connecting = 0  # number of probe harmonics with a term from the ket to the bra sector
    for operator in harmonics.values():
        try:
            check_transition(operator, periodic.spins, periodic.eta_labels, ket.basis.sector, bra.basis.sector)
            connecting += 1
        except ValueError:
            continue
    if not connecting:
        raise ValueError("no harmonic of the probe connects the ket sector to the bra sector")
    bra_idx = _indices(bra_states, len(bra.quasienergies), "bra_states")
    ket_idx = _indices(ket_states, len(ket.quasienergies), "ket_states")
    options = bra.metadata.get("options", {})
    gib = options.get("max_memory_gib", FloquetOptions.max_memory_gib) if max_memory_gib is None else max_memory_gib
    ket_modes = ket.fock_modes()[:, :, ket_idx]
    ket_harmonics = np.asarray(ket.harmonic_indices)
    probe_harmonics = sorted(harmonics)
    window = np.arange(ket_harmonics[0] + probe_harmonics[0], ket_harmonics[-1] + probe_harmonics[-1] + 1)
    bra_dimension = bra.basis.dimension
    estimate = 16*(len(window)*bra_dimension*len(ket_idx) + ket_modes.size)
    if estimate > gib*1024**3:
        raise MemoryError(f"the probe images need {estimate/1024**3:.3g} GiB, more than "
                          f"max_memory_gib={gib:g}; restrict ket_states")
    with threadpool_limits(limits=options.get("threads", FloquetOptions.threads)):
        images = np.zeros((len(window), bra_dimension, len(ket_idx)), complex)
        flat = ket_modes.transpose(1, 0, 2).reshape(ket_modes.shape[1], -1)
        shared = same_basis(bra.basis, ket.basis)
        for k in probe_harmonics:
            projected = ProjectedOperator(harmonics[k], ket.basis)
            if shared:
                result = projected.sparse_matrix() @ flat
            else:
                result = np.column_stack([projected.action_between(np.ascontiguousarray(v), bra.basis)
                                          for v in flat.T])
            result = result.reshape(bra_dimension, len(ket_harmonics), len(ket_idx))
            for position, m in enumerate(ket_harmonics):
                images[m + k - window[0]] += result[:, position, :]
        total = np.sum(np.abs(images)**2, axis=(0, 1))
        if bra.representation == "static":
            vectors = np.asarray(bra.static.vectors[:, :bra.modes.shape[1]], complex)
            images = np.einsum("Dd,mDn->mdn", vectors.conj(), images)
        bra_modes = bra.modes[:, :, bra_idx]
        bra_harmonics = np.asarray(bra.harmonic_indices)
        # F(theta) = <u_bra(theta)|w(theta)> has len(photons) Fourier components; as many equally
        # spaced samples give every X^(q) exactly by a discrete Fourier transform, without aliasing.
        photons = np.arange(window[0] - bra_harmonics[-1], window[-1] - bra_harmonics[0] + 1)
        samples = len(photons)
        theta = 2*np.pi*np.arange(samples)/samples
        values = np.empty((samples, len(bra_idx), len(ket_idx)), complex)
        for j, angle in enumerate(theta):
            bra_t = np.tensordot(np.exp(1j*bra_harmonics*angle), bra_modes, axes=1)
            image_t = np.tensordot(np.exp(1j*window*angle), images, axes=1)
            values[j] = bra_t.conj().T @ image_t
        transform = np.fft.fft(values, axis=0)/samples
        elements = transform[np.mod(photons, samples)]
    default = ket.dressed_state(0)
    if default not in ket_idx:
        default = int(ket_idx[0])
    split = bool(ket.metadata.get("multiplet_split")) or bool(bra.metadata.get("multiplet_split"))
    metadata = dict(probe_harmonics=probe_harmonics, connecting_harmonics=connecting, bra_sector=bra.metadata.get("sector"),
                    ket_sector=ket.metadata.get("sector"), bra_method=bra.method, ket_method=ket.method,
                    multiplet_split=split, convention="X^(q)_ab = sum_m <u_{a,m-q}|(X u_b)_m>; "
                    "line frequency eps_a - eps_b + q omega")
    return FloquetTransitions(bra.omega, photons, elements, bra.quasienergies[bra_idx], ket.quasienergies[ket_idx],
                              bra.labels[bra_idx], ket.labels[ket_idx], bra_idx, ket_idx, total, same,
                              bool(bra.metadata.get("complete")) and bra_states is None, default, name, metadata)


__all__ = ["floquet_matrix_elements"]
