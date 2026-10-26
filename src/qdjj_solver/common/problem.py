"""Finite fermionic problems and shared symmetry conventions."""
from dataclasses import dataclass, field
import numpy as np
from .algebra import FermionOperator


@dataclass(frozen=True)
class Sector:
    """Symmetry sector: fermion parity (0 even, 1 odd), twice S_z, and optionally the eta parity
    and the particle number of the canonical modes; ``None`` leaves a quantum number unrestricted."""
    parity: int = 1
    """Total fermion parity: zero is even and one is odd; it includes impurity and bath."""
    twice_sz: int | None = 1
    """Twice total spin projection, an integer compatible with parity; ``None`` leaves it unrestricted."""
    eta: int | None = None
    """Optional eta parity, +1 or -1; requires declared eta labels and a conserving Hamiltonian."""
    particle_number: int | None = None
    """Optional nonnegative total occupation of the chosen canonical modes, compatible with parity and spin."""

    def __post_init__(self):
        if self.parity not in (0, 1):
            raise ValueError("parity must be 0 (even) or 1 (odd)")
        if self.twice_sz is not None:
            if not isinstance(self.twice_sz, (int, np.integer)):
                raise ValueError("twice_sz must be an integer or None")
            if self.twice_sz % 2 != self.parity:
                raise ValueError("fermion parity and twice_sz are incompatible")
        if self.eta not in (None, -1, 1):
            raise ValueError("eta must be +1, -1, or None")
        if self.particle_number is not None:
            if (not isinstance(self.particle_number, (int, np.integer)) or self.particle_number < 0
                    or self.particle_number % 2 != self.parity):
                raise ValueError("particle_number must be nonnegative and compatible with parity")
            if self.twice_sz is not None and abs(self.twice_sz) > self.particle_number:
                raise ValueError("particle_number and twice_sz are incompatible")


DEFAULT_SECTOR = Sector()


def check_operator_sector(operator, spins, eta_labels, sector):
    """Require every term of a (not necessarily Hermitian) operator to conserve a sector."""
    if sector.eta is not None and eta_labels is None:
        raise ValueError("eta reduction is not defined for this Hamiltonian")
    if sector.particle_number is not None and sector.particle_number > len(spins):
        raise ValueError("particle_number exceeds the mode count")
    for key in operator.terms:
        if sector.particle_number is not None and sum(1 if op > 0 else -1 for op in key):
            raise ValueError("Hamiltonian does not conserve particle number")
        if sector.twice_sz is not None:
            change = sum((1 if op > 0 else -1) * spins[abs(op)-1] for op in key)
            if change:
                raise ValueError("Hamiltonian does not conserve S_z; choose twice_sz=None")
        if sector.eta is not None:
            if np.prod([eta_labels[abs(op)-1] for op in key], dtype=int) != 1:
                raise ValueError("Hamiltonian does not conserve the requested eta parity")


def check_transition(operator, spins, eta_labels, source, destination, tolerance=1e-12):
    """Require some term of an operator to map the source sector into the destination.

    Terms that leave the destination sector are dropped by projected actions, as they
    should be; an operator without any connecting term indicates a wrong pairing of
    sectors. An unrestricted quantum number (None) on either side is not checked. Terms
    below ``tolerance`` relative to the largest coefficient (rounding noise of coordinate
    transformations) do not count as connections.
    """
    if not operator.terms:
        raise ValueError("the operator is zero")
    if (source.eta is not None or destination.eta is not None) and eta_labels is None:
        raise ValueError("eta reduction is not defined for this Hamiltonian")
    scale = max(abs(v) for v in operator.terms.values())
    for key, value in operator.terms.items():
        if abs(value) <= tolerance*scale:
            continue
        if (source.parity + len(key)) % 2 != destination.parity:
            continue
        if source.twice_sz is not None and destination.twice_sz is not None:
            change = sum((1 if op > 0 else -1) * spins[abs(op)-1] for op in key)
            if source.twice_sz + change != destination.twice_sz:
                continue
        if source.eta is not None and destination.eta is not None:
            if source.eta * np.prod([eta_labels[abs(op)-1] for op in key], dtype=int) != destination.eta:
                continue
        if source.particle_number is not None and destination.particle_number is not None:
            if source.particle_number + sum(1 if op > 0 else -1 for op in key) != destination.particle_number:
                continue
        return
    raise ValueError("no term of the operator connects the source sector to the destination sector")


@dataclass
class Hamiltonian:
    """A finite, Hermitian, parity-conserving Hamiltonian in canonical fermionic coordinates.

    The first ``nimp`` modes are impurity spin orbitals; ``spins`` holds +1 or -1 for every mode and
    ``eta_labels`` optional eta-parity labels. ``observables`` maps names to operators in the same
    coordinates; ``metadata`` records the construction, including the electron-to-canonical map that
    ``physical_operator`` uses.
    """

    operator: FermionOperator
    """Finite Hermitian, parity-conserving polynomial in the declared canonical modes."""
    nimp: int
    """Number of impurity spin modes, occupying the first mode indices; zero is allowed."""
    spins: tuple[int, ...]
    """Twice-spin labels +1 or -1 for every canonical mode, in basis order."""
    eta_labels: tuple[int, ...] | None = None
    """Optional +1/-1 eta label per mode; required when imposing eta parity."""
    observables: dict[str, FermionOperator] = field(default_factory=dict)
    """Named operators in the same canonical coordinates as the Hamiltonian."""
    metadata: dict = field(default_factory=dict)
    """Construction, energy-reference and coordinate information; may include user entries."""

    def __post_init__(self):
        self.spins = tuple(self.spins)
        if not 0 <= self.nimp <= len(self.spins) or any(s not in (-1, 1) for s in self.spins):
            raise ValueError("invalid impurity dimension or mode spin labels")
        if self.eta_labels is not None:
            self.eta_labels = tuple(self.eta_labels)
            if len(self.eta_labels) != len(self.spins) or any(e not in (-1, 1) for e in self.eta_labels):
                raise ValueError("invalid eta labels")
        for key in self.operator.terms:
            if any(abs(op) > len(self.spins) for op in key):
                raise ValueError("Hamiltonian acts on an undeclared fermionic mode")
            if len(key) % 2:
                raise ValueError("a physical Hamiltonian must preserve fermion parity")
        scale = max((abs(v) for v in self.operator.terms.values()), default=1.)
        if self.operator.hermiticity_error() > 1e-12 * scale:
            raise ValueError("Hamiltonian is not Hermitian")

    @property
    def bath_modes(self):
        """Number of declared non-impurity spin modes, including retained dark modes."""
        return len(self.spins) - self.nimp

    def record(self):
        """Full operator, mode labels, observables and metadata as a JSON-ready record."""
        return dict(operator=self.operator.to_records(), nimp=self.nimp, spins=list(self.spins),
                    eta_labels=None if self.eta_labels is None else list(self.eta_labels),
                    observables={k: v.to_records() for k, v in self.observables.items()}, metadata=self.metadata)

    @classmethod
    def from_record(cls, record):
        """Reconstruct and validate a Hamiltonian record, preserving its additive constant."""
        return cls(FermionOperator.from_records(record["operator"]), record["nimp"], tuple(record["spins"]),
                   record.get("eta_labels"),
                   {k: FermionOperator.from_records(v) for k, v in record.get("observables", {}).items()},
                   record.get("metadata", {}))

    def check_sector(self, sector):
        """Require every term of the Hamiltonian to conserve the sector's quantum numbers."""
        check_operator_sector(self.operator, self.spins, self.eta_labels, sector)

    def fingerprint(self):
        """SHA256 of the full canonical record, including observables and metadata."""
        from .io import fingerprint
        return fingerprint(self.record())

    def coordinate_fingerprint(self):
        """A hash that identifies the canonical coordinates, independent of the numerical backend.

        User-built Hamiltonians can declare their basis in ``metadata['coordinates']``.
        Without this declaration the hash covers the entire Hamiltonian, so that states
        computed in different bases are never combined by accident.
        """
        from .io import fingerprint
        coordinates = self.metadata.get("coordinates")
        if coordinates is None:
            return self.fingerprint()
        return fingerprint(dict(coordinates=coordinates, nimp=self.nimp,
                                spins=self.spins, eta_labels=self.eta_labels))

    def probe(self, name, orbital=None):
        """A named probe operator in canonical coordinates (see ``common.probes``)."""
        from .probes import probe_operator
        return probe_operator(self, name, orbital)

    def physical_operator(self, operator):
        """Transform a physical electron operator into the model's canonical coordinates.

        The electron modes are ordered as the impurity modes, then the reservoir
        modes by reservoir, level and spin. Arbitrary products require the complete
        bath map, so this method rejects models with paired-mode compression, whose
        discarded dark modes are frozen in their vacuum; the model's built-in
        observables already include the contractions required by that restriction.
        """
        coordinates = self.metadata.get("coordinates", {})
        if coordinates.get("restriction") != "none":
            raise ValueError("physical operators require an uncompressed declared electron map")
        maps = coordinates.get("electron_annihilators")
        if maps is None:
            raise ValueError("this Hamiltonian has no declared electron map")
        if any(abs(op) > len(maps) for key in operator.terms for op in key):
            raise ValueError("physical operator acts on an undeclared electron mode")
        return operator.substitute([FermionOperator.from_records(m) for m in maps])
