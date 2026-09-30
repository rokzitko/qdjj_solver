"""Finite fermionic problems and shared symmetry conventions."""
from dataclasses import dataclass, field
import numpy as np
from .algebra import FermionOperator


@dataclass(frozen=True)
class Sector:
    """Parity (0 even, 1 odd), twice S_z, optional eta and mode particle number."""
    parity: int = 1
    twice_sz: int | None = 1
    eta: int | None = None
    particle_number: int | None = None

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


@dataclass
class Hamiltonian:
    operator: FermionOperator
    nimp: int
    spins: tuple[int, ...]
    eta_labels: tuple[int, ...] | None = None
    observables: dict[str, FermionOperator] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

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
        return len(self.spins) - self.nimp

    def record(self):
        return dict(operator=self.operator.to_records(), nimp=self.nimp, spins=list(self.spins),
                    eta_labels=None if self.eta_labels is None else list(self.eta_labels),
                    observables={k: v.to_records() for k, v in self.observables.items()}, metadata=self.metadata)

    @classmethod
    def from_record(cls, record):
        return cls(FermionOperator.from_records(record["operator"]), record["nimp"], tuple(record["spins"]),
                   record.get("eta_labels"),
                   {k: FermionOperator.from_records(v) for k, v in record.get("observables", {}).items()},
                   record.get("metadata", {}))

    def check_sector(self, sector):
        if sector.eta is not None and self.eta_labels is None:
            raise ValueError("eta reduction is not defined for this Hamiltonian")
        if sector.particle_number is not None and sector.particle_number > len(self.spins):
            raise ValueError("particle_number exceeds the mode count")
        for key in self.operator.terms:
            if sector.particle_number is not None and sum(1 if op > 0 else -1 for op in key):
                raise ValueError("Hamiltonian does not conserve particle number")
            if sector.twice_sz is not None:
                change = sum((1 if op > 0 else -1) * self.spins[abs(op)-1] for op in key)
                if change:
                    raise ValueError("Hamiltonian does not conserve S_z; choose twice_sz=None")
            if sector.eta is not None:
                if np.prod([self.eta_labels[abs(op)-1] for op in key], dtype=int) != 1:
                    raise ValueError("Hamiltonian does not conserve the requested eta parity")

    def fingerprint(self):
        from .io import fingerprint
        return fingerprint(self.record())

    def coordinate_fingerprint(self):
        """Identify the many-body coordinates, independently of numerical backend.

        User-built Hamiltonians can supply ``metadata['coordinates']`` as an
        explicit basis declaration. Without it we conservatively identify the
        entire Hamiltonian, preventing accidental transitions between bases.
        """
        from .io import fingerprint
        coordinates = self.metadata.get("coordinates")
        if coordinates is None:
            return self.fingerprint()
        return fingerprint(dict(coordinates=coordinates, nimp=self.nimp,
                                spins=self.spins, eta_labels=self.eta_labels))

    def physical_operator(self, operator):
        """Transform an electron operator using the model's fixed canonical map.

        Electron order is impurity modes, then reservoir/level/spin modes.
        Arbitrary products require the complete bath map, so this method rejects
        dark-vacuum compression; the model's built-in observables already include
        the contractions required by that restriction.
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
