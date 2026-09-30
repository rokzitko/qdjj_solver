"""Symmetry sectors and truncated Fock spaces."""

from __future__ import annotations

import numpy as np

from . import _core
from ..common.problem import Sector as Sector, DEFAULT_SECTOR


def dimension(nimp, spins, cutoff, sector, eta_labels=None):
    """Exact integer dimension by counting occupations, without forming kets."""
    if (not isinstance(nimp, (int, np.integer)) or not 0 <= nimp <= len(spins)
            or any(s not in (-1, 1) for s in spins)):
        raise ValueError("invalid impurity size or spin labels")
    if not isinstance(cutoff, (int, np.integer)) or not 0 <= cutoff <= len(spins) - nimp:
        raise ValueError("QP cutoff must lie between zero and the bath mode count")
    if eta_labels is not None and (len(eta_labels) != len(spins) or any(e not in (-1, 1) for e in eta_labels)):
        raise ValueError("invalid eta labels")
    if sector.eta is not None and eta_labels is None:
        raise ValueError("this representation has no diagonal eta symmetry")
    counts = {(0, 0, 0, 1, 0): 1}
    for i, spin in enumerate(spins):
        updated = counts.copy()
        dq = int(i >= nimp)
        e = eta_labels[i] if sector.eta is not None else 1
        for (q, p, sz, eta, particles), count in counts.items():
            if q + dq > cutoff:
                continue
            if sector.particle_number is not None and particles+1 > sector.particle_number:
                continue
            key = (q + dq, 1 - p, sz + spin if sector.twice_sz is not None else 0, eta * e,
                   particles+1 if sector.particle_number is not None else 0)
            updated[key] = updated.get(key, 0) + count
        counts = updated
    return sum(count for (q, p, sz, eta, particles), count in counts.items()
               if p == sector.parity and (sector.twice_sz is None or sz == sector.twice_sz)
               and (sector.eta is None or eta == sector.eta)
               and (sector.particle_number is None or particles == sector.particle_number))


class FockBasis:
    def __init__(self, nimp, spins, cutoff, sector=DEFAULT_SECTOR, eta_labels=None,
                 max_dimension=5_000_000):
        self.nimp, self.spins, self.cutoff = nimp, tuple(spins), cutoff
        self.sector = sector
        eta_labels = None if eta_labels is None else tuple(eta_labels)
        self.dimension = dimension(nimp, self.spins, cutoff, sector, eta_labels)
        if self.dimension == 0:
            raise ValueError("empty symmetry sector")
        if self.dimension > max_dimension:
            raise MemoryError(f"sector dimension {self.dimension:,} exceeds max_dimension={max_dimension:,}")
        self.core = _core.Basis(nimp, list(self.spins), cutoff, sector.parity,
                                _core_no_sz if sector.twice_sz is None else sector.twice_sz,
                                [] if eta_labels is None else list(eta_labels),
                                0 if sector.eta is None else sector.eta,
                                int(max_dimension), -1 if sector.particle_number is None else sector.particle_number)
        if self.core.dimension != self.dimension:
            raise RuntimeError("Fock-space counting and enumeration disagree")

    @property
    def bytes(self):
        return self.core.bytes

    def occupations(self):
        """Copy occupations as rows of 64-bit words (least-significant word first)."""
        return self.core.occupations()

    def qp_counts(self):
        return self.core.qp_counts()


_core_no_sz = 1_000_000_000
