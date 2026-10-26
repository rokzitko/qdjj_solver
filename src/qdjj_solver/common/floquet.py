"""Periodically driven Hamiltonians in fixed canonical coordinates.

Fourier convention: H(t)=sum_k H_k exp(i k omega t), with H_{-k}=H_k^dagger and
time origin t=0; omega is the drive frequency. The harmonics are obtained from
equally spaced samples of a model over one period. The canonical coordinates must
not depend on time: a moving quasiparticle basis would add the gauge term
i (dW/dt) W^dagger of the time-dependent transformation W(t).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from .algebra import FermionOperator, is_integer
from .problem import Hamiltonian, check_operator_sector

FOURIER_CONVENTION = "H(t)=sum_k H_k exp(i k omega t); H_{-k}=H_k^dagger; t0=0"

_MOVING_COORDINATES = (
    "the canonical coordinates change during the drive period. Not supported are drives of the gap "
    "delta or of the bath grid, a phase drive with paired-mode compression, a detuning or rho_wn drive "
    "with symmetry=True (eta coordinates), and driven phases or direct hopping with "
    "bath_reference='coupled'; build the model with symmetry=False, compress=False and "
    "bath_reference='isolated'")


@dataclass(frozen=True)
class Harmonic:
    """One cosine component A cos(n omega t + theta) of a parameter modulation."""

    order: int
    """Positive integer multiple of the fundamental drive frequency."""
    amplitude: float
    """Finite real cosine amplitude in the driven parameter's units."""
    phase: float = 0.
    """Finite cosine phase offset in radians."""

    def __post_init__(self):
        if not is_integer(self.order) or self.order < 1:
            raise ValueError("harmonic order must be a positive integer")
        if not (np.isreal(self.amplitude) and np.isreal(self.phase)
                and np.isfinite(self.amplitude) and np.isfinite(self.phase)):
            raise ValueError("harmonic amplitude and phase must be real and finite")

    def record(self):
        """JSON-ready order, amplitude and phase for this cosine component."""
        return dict(order=int(self.order), amplitude=float(self.amplitude), phase=float(self.phase))


@dataclass(frozen=True)
class ParameterDrive:
    """p(t)=p0+sum_n A_n cos(n omega t+theta_n) for one builder keyword.

    ``index`` selects one element of an array keyword, such as ``phases`` of
    :func:`make_model`. Amplitudes carry the parameter's units (radians for phases).
    """

    parameter: str
    """Nonempty keyword name in the model builder's undriven parameter mapping."""
    harmonics: tuple[Harmonic, ...]
    """Nonempty sequence of cosine components with distinct positive orders."""
    index: int | None = None
    """Optional nonnegative index of an array parameter; ``None`` selects a scalar parameter."""

    def __post_init__(self):
        if not isinstance(self.parameter, str) or not self.parameter:
            raise ValueError("the driven parameter must be a keyword name")
        harmonics = tuple(self.harmonics)
        if not harmonics or not all(isinstance(h, Harmonic) for h in harmonics):
            raise ValueError("a parameter drive needs at least one Harmonic")
        if len({h.order for h in harmonics}) != len(harmonics):
            raise ValueError("combine repeated harmonic orders into one Harmonic")
        if self.index is not None and (not is_integer(self.index) or self.index < 0):
            raise ValueError("the parameter index must be a nonnegative integer")
        object.__setattr__(self, "harmonics", harmonics)

    @classmethod
    def cosine(cls, parameter, amplitude, phase=0., order=1, index=None):
        """A single harmonic: p(t)=p0+amplitude*cos(order*omega*t+phase)."""
        return cls(parameter, (Harmonic(order, amplitude, phase),), index)

    @property
    def max_order(self):
        """Highest cosine order in this parameter modulation."""
        return max(h.order for h in self.harmonics)

    def value(self, t, omega):
        """Displacement p(t)-p0."""
        return sum(h.amplitude*np.cos(h.order*omega*np.asarray(t, float) + h.phase) for h in self.harmonics)

    def record(self):
        """JSON-ready parameter name, optional array index and harmonic records."""
        return dict(parameter=self.parameter, index=self.index, harmonics=[h.record() for h in self.harmonics])


def fold_quasienergies(energies, omega, center=0.):
    """Map quasienergies into the zone [center-omega/2, center+omega/2)."""
    energies = np.asarray(energies, float)
    return center + np.mod(energies - center + omega/2, omega) - omega/2


def quasienergy_distance(a, b, omega):
    """Distance between quasienergies modulo omega."""
    difference = np.mod(np.asarray(a, float) - np.asarray(b, float), omega)
    return np.minimum(difference, omega - difference)


def _commutator(a, b):
    return a*b - b*a


def without_scalar(operator):
    """The operator without its c-number term."""
    return FermionOperator({key: value for key, value in operator.terms.items() if key})


def _time_average(omega, operator, nimp, spins, eta_labels, observables, coordinates):
    """H_0 as a Hamiltonian, with the time averages (k=0 harmonics) of the observables."""
    metadata = dict(floquet=dict(kind="time-average", omega=omega))
    if coordinates is not None:
        metadata["coordinates"] = coordinates
    averages = {name: values[0] for name, values in observables.items() if 0 in values}
    return Hamiltonian(operator, nimp, spins, eta_labels, averages, metadata)


@dataclass
class PeriodicHamiltonian:
    """Hermitian H(t+T)=H(t), stored as harmonics H_k (k>=0) in fixed canonical coordinates.

    ``reference`` is the static Hamiltonian whose eigenstates label the dressed
    states; it defaults to the time average H_0. ``observables[name][k]`` holds the
    harmonics (both signs of k) of the model observables, used for probe matrix
    elements between Floquet states.
    """

    omega: float
    """Positive finite fundamental angular frequency in energy units with hbar=1."""
    harmonics: dict[int, FermionOperator]
    """Nonnegative Fourier indices mapped to operators; negative harmonics are Hermitian conjugates."""
    nimp: int
    """Number of impurity spin modes, identical for every harmonic and the reference."""
    spins: tuple[int, ...]
    """Twice-spin label +1 or -1 per fixed canonical mode."""
    eta_labels: tuple[int, ...] | None = None
    """Optional eta label per mode; fixed across all harmonics."""
    reference: Hamiltonian | None = None
    """Static Hamiltonian defining tracked state labels; ``None`` uses the time average."""
    observables: dict[str, dict[int, FermionOperator]] = field(default_factory=dict)
    """Named operator harmonics, allowing both positive and negative Fourier indices."""
    metadata: dict = field(default_factory=dict)
    """Coordinate declarations, drive settings, sampling diagnostics and optional user entries."""

    def __post_init__(self):
        if not np.isreal(self.omega) or not np.isfinite(self.omega) or self.omega <= 0:
            raise ValueError("the drive frequency omega must be positive and finite")
        self.omega = float(self.omega)
        self.spins = tuple(self.spins)
        self.eta_labels = None if self.eta_labels is None else tuple(self.eta_labels)
        harmonics = {}
        for k, operator in dict(self.harmonics).items():
            if not is_integer(k) or k < 0:
                raise ValueError("store harmonics for k >= 0 only; H_{-k} is H_k^dagger")
            if not isinstance(operator, FermionOperator):
                raise ValueError("harmonics must be FermionOperators")
            if k == 0 or operator.terms:
                harmonics[int(k)] = operator
        harmonics.setdefault(0, FermionOperator())
        self.harmonics = dict(sorted(harmonics.items()))
        for k, operator in self.harmonics.items():
            for key in operator.terms:
                if any(abs(op) > len(self.spins) for op in key):
                    raise ValueError(f"harmonic {k} acts on an undeclared fermionic mode")
                if len(key) % 2:
                    raise ValueError(f"harmonic {k} does not preserve fermion parity")
        observables = {}
        for name, values in dict(self.observables).items():
            if not all(is_integer(k) and isinstance(v, FermionOperator) for k, v in values.items()):
                raise ValueError(f"observable {name!r} needs integer harmonics of FermionOperators")
            observables[name] = {int(k): v for k, v in sorted(values.items())}
        self.observables = observables
        if self.reference is None:
            self.reference = self._average(self.metadata.get("coordinates"))
        elif not isinstance(self.reference, Hamiltonian):
            raise ValueError("the reference must be a Hamiltonian")
        elif ((self.reference.nimp, self.reference.spins, self.reference.eta_labels)
              != (self.nimp, self.spins, self.eta_labels)):
            raise ValueError("the reference Hamiltonian has different mode labels")
        else:
            # Validates the modes, parity and Hermiticity of the time average.
            self._average(None)

    def _average(self, coordinates):
        return _time_average(self.omega, self.harmonics[0], self.nimp, self.spins, self.eta_labels,
                             self.observables, coordinates)

    @property
    def period(self):
        """Drive period 2*pi/omega in inverse Hamiltonian-energy units."""
        return 2*np.pi/self.omega

    @property
    def max_harmonic(self):
        """Largest stored nonnegative Hamiltonian harmonic index."""
        return max(self.harmonics)

    def harmonic(self, k):
        """H_k for any integer k, with H_{-k}=H_k^dagger."""
        if not is_integer(k):
            raise ValueError("harmonic index must be an integer")
        if k >= 0:
            return self.harmonics.get(int(k), FermionOperator())
        return self.harmonics.get(-int(k), FermionOperator()).dagger()

    def operator_at(self, t):
        """H(t) as a FermionOperator."""
        result = self.harmonics[0]
        for k, operator in self.harmonics.items():
            if k:
                term = operator*complex(np.exp(1j*k*self.omega*t))
                result = result + term + term.dagger()
        return result

    def time_average(self) -> Hamiltonian:
        """The time average H_0 as a Hamiltonian, with the time-averaged observables."""
        return self._average(self.reference.metadata.get("coordinates"))

    def check_sector(self, sector):
        """Require every harmonic to conserve the sector; a projection would silently drop terms."""
        check_operator_sector(self.reference.operator, self.spins, self.eta_labels, sector)
        for operator in self.harmonics.values():
            check_operator_sector(operator, self.spins, self.eta_labels, sector)

    def coordinate_fingerprint(self):
        """Coordinate identity of the static reference; independent of Fourier truncation."""
        return self.reference.coordinate_fingerprint()

    def probe(self, name, orbital=None):
        """Harmonics {k: operator} of a named probe under this drive (see ``common.probes``)."""
        from .probes import periodic_probe
        return periodic_probe(self, name, orbital)

    def record(self):
        """JSON-ready periodic Hamiltonian; harmonic keys are strings and operators are term lists."""
        return dict(omega=self.omega, nimp=self.nimp, spins=list(self.spins),
                    eta_labels=None if self.eta_labels is None else list(self.eta_labels),
                    harmonics={str(k): op.to_records() for k, op in self.harmonics.items()},
                    reference=self.reference.record(),
                    observables={name: {str(k): op.to_records() for k, op in values.items()}
                                 for name, values in self.observables.items()},
                    metadata=self.metadata)

    @classmethod
    def from_record(cls, record):
        """Reconstruct a periodic-Hamiltonian record, including reference and observable harmonics."""
        return cls(record["omega"], {int(k): FermionOperator.from_records(v) for k, v in record["harmonics"].items()},
                   record["nimp"], tuple(record["spins"]), record.get("eta_labels"),
                   Hamiltonian.from_record(record["reference"]),
                   {name: {int(k): FermionOperator.from_records(v) for k, v in values.items()}
                    for name, values in record.get("observables", {}).items()},
                   record.get("metadata", {}))

    def fingerprint(self):
        """SHA256 of the full periodic record, used to check transition compatibility."""
        from .io import fingerprint
        return fingerprint(self.record())

    def van_vleck(self, order=2, *, tolerance=1e-14, max_harmonic=None) -> Hamiltonian:
        """High-frequency (van Vleck) effective Hamiltonian through O(omega^-order).

        Commutators are formed in the full Fock space, before any QP truncation.
        The expansion requires omega to exceed the many-body spectral width; its
        cost grows with products of term counts. Observables are copied from the
        time average without the O(1/omega) kick-operator correction.
        """
        if order not in (0, 1, 2):
            raise ValueError("van Vleck order must be 0, 1 or 2")
        top = self.max_harmonic if max_harmonic is None else max_harmonic
        if not is_integer(top) or top < 0:
            raise ValueError("max_harmonic must be a nonnegative integer")
        # Constant parts of H_k (k != 0) commute with everything.
        h = {k: (self.harmonic(k) if k == 0 else without_scalar(self.harmonic(k)))
             for k in range(-top, top+1)}
        h = {k: v for k, v in h.items() if k == 0 or v.terms}
        result = self.harmonics[0]
        if order >= 1:
            first = FermionOperator()
            for m in range(1, top+1):
                if m in h:
                    first += _commutator(h[m], h[-m])/m
            result = result + first/self.omega
        if order >= 2:
            second = FermionOperator()
            nonzero = [k for k in h if k]
            for m in nonzero:
                second += _commutator(_commutator(h[-m], h[0]), h[m])/(2*m*m)
                for mp in nonzero:
                    if mp != m and mp-m in h:
                        second += _commutator(_commutator(h[-mp], h[mp-m]), h[m])/(3*m*mp)
            result = result + second/self.omega**2
        exact = (result + result.dagger())*.5
        cleaned = exact.cleaned(tolerance)
        average = self.time_average()
        metadata = dict(average.metadata)
        metadata["floquet_effective"] = dict(
            method="van Vleck", order=order, omega=self.omega, max_harmonic=top,
            cleanup_bound=float(sum(abs(v) for v in (exact-cleaned).terms.values())),
            periodic_sha256=self.fingerprint(),
            observables="time-averaged lab-frame operators without kick-operator corrections")
        return Hamiltonian(cleaned, self.nimp, self.spins, self.eta_labels, average.observables, metadata)

    @classmethod
    def from_parameters(cls, builder, base: Mapping, drives, omega, *, samples=64,
                        coefficient_tolerance=1e-14, sampling_tolerance=1e-9, verify=True,
                        observables=True) -> PeriodicHamiltonian:
        """Harmonics of a model constructor whose keyword arguments are driven periodically.

        ``base`` holds the undriven keywords; each :class:`ParameterDrive`
        displaces one of them. ``builder`` is evaluated at ``samples`` equally
        spaced times of one period, and the operator coefficients are Fourier
        transformed. With ``verify=True`` the Fourier series is compared with the
        model at the midpoints between samples and must agree within
        ``sampling_tolerance``. The undriven model ``builder(**base)`` labels the
        dressed states. For the reference junction use, e.g.,
        ``builder=reference_model, base=dict(bath=bath, phi=.8, symmetry=False, compress=False)``.
        """
        base = dict(base)
        drives = (drives,) if isinstance(drives, ParameterDrive) else tuple(drives)
        if not drives or not all(isinstance(d, ParameterDrive) for d in drives):
            raise ValueError("provide one or more ParameterDrive objects")
        if len({(d.parameter, d.index) for d in drives}) != len(drives):
            raise ValueError("each parameter (and index) may be driven only once")
        for drive in drives:
            if drive.parameter not in base:
                raise ValueError(f"driven parameter {drive.parameter!r} needs an undriven value in base")
            value = np.asarray(base[drive.parameter])
            if drive.index is None:
                if value.ndim or not np.isrealobj(value) or not np.isfinite(value):
                    raise ValueError(f"driven parameter {drive.parameter!r} must be a real scalar; "
                                     "use index= for array parameters")
            elif value.ndim != 1 or not np.isrealobj(value) or drive.index >= len(value):
                raise ValueError(f"index {drive.index} does not select an element of {drive.parameter!r}")

        def parameters(t):
            result = dict(base)
            for drive in drives:
                if drive.index is None:
                    result[drive.parameter] = float(base[drive.parameter]) + float(drive.value(t, omega))
                else:
                    array = np.array(result[drive.parameter], float)
                    array[drive.index] += float(drive.value(t, omega))
                    result[drive.parameter] = array
            return result

        reference = _sample(lambda t: builder(**base), 0.)
        scalars = {key: value for key, value in base.items()
                   if isinstance(value, (bool, int, float, str, np.integer, np.floating))}
        extra = dict(drives=[d.record() for d in drives], base=scalars,
                     builder=getattr(builder, "__qualname__", repr(builder)))
        return cls._from_samples(lambda t: builder(**parameters(t)), omega, samples, reference, False,
                                 coefficient_tolerance, sampling_tolerance, verify, observables, extra)

    @classmethod
    def from_callable(cls, function, omega, *, samples=64, reference=None, assume_fixed_coordinates=False,
                      coefficient_tolerance=1e-14, sampling_tolerance=1e-9, verify=True,
                      observables=True) -> PeriodicHamiltonian:
        """Harmonics of an arbitrary periodic ``function(t) -> Hamiltonian``.

        Without ``reference`` the time average labels the dressed states. The
        sampled Hamiltonians must declare their canonical coordinates in
        ``metadata['coordinates']``, and these must not change during the period.
        ``assume_fixed_coordinates=True`` accepts Hamiltonians without this
        declaration; only their mode labels are then compared.
        """
        extra = dict(builder=getattr(function, "__qualname__", repr(function)))
        return cls._from_samples(function, omega, samples, reference, assume_fixed_coordinates,
                                 coefficient_tolerance, sampling_tolerance, verify, observables, extra)

    @classmethod
    def _from_samples(cls, function, omega, samples, reference, assume_fixed, coefficient_tolerance,
                      sampling_tolerance, verify, observables, extra):
        if not np.isreal(omega) or not np.isfinite(omega) or omega <= 0:
            raise ValueError("the drive frequency omega must be positive and finite")
        if not is_integer(samples) or samples < 3:
            raise ValueError("at least three samples per period are required")
        for name, value in (("coefficient_tolerance", coefficient_tolerance),
                            ("sampling_tolerance", sampling_tolerance)):
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be nonnegative and finite")
        omega, samples = float(omega), int(samples)
        period = 2*np.pi/omega
        sampled = [_sample(function, j*period/samples) for j in range(samples)]
        template = sampled[0] if reference is None else reference
        check = _coordinate_check(template, assume_fixed)
        for h in sampled:
            check(h)
        keys = sorted(set().union(*(h.operator.terms for h in sampled)), key=lambda k: (len(k), k))
        coefficients = np.array([[h.operator.terms.get(key, 0) for key in keys] for h in sampled], complex)
        transform = np.fft.fft(coefficients, axis=0)/samples
        top = (samples-1)//2
        harmonics, l1, cleanup = {}, [], 0.
        for k in range(top+1):
            exact = FermionOperator(dict(zip(keys, transform[k], strict=True)))
            if k == 0:
                exact = (exact + exact.dagger())*.5
            cleaned = exact.cleaned(coefficient_tolerance)
            cleanup += sum(abs(v) for v in (exact-cleaned).terms.values())
            l1.append(float(sum(abs(v) for v in cleaned.terms.values())))
            if k == 0 or cleaned.terms:
                harmonics[k] = cleaned
        hermiticity = 0.
        for k in range(1, top+1):
            negative = FermionOperator(dict(zip(keys, transform[samples-k], strict=True)))
            difference = negative - FermionOperator(dict(zip(keys, transform[k], strict=True))).dagger()
            hermiticity = max(hermiticity, max(map(abs, difference.terms.values()), default=0.))
        observable_harmonics = {}
        if observables:
            names = set(sampled[0].observables)
            if any(set(h.observables) != names for h in sampled):
                raise ValueError("the sampled models declare different observables")
            for name in sorted(names):
                okeys = sorted(set().union(*(h.observables[name].terms for h in sampled)), key=lambda k: (len(k), k))
                values = np.array([[h.observables[name].terms.get(key, 0) for key in okeys] for h in sampled], complex)
                values = np.fft.fft(values, axis=0)/samples
                result = {}
                for k in range(-top, top+1):
                    operator = FermionOperator(dict(zip(okeys, values[k % samples], strict=True)))
                    operator = operator.cleaned(coefficient_tolerance)
                    if operator.terms:
                        result[k] = operator
                observable_harmonics[name] = result
        coordinates = template.metadata.get("coordinates")
        if reference is None:
            reference = _time_average(omega, harmonics[0], template.nimp, template.spins, template.eta_labels,
                                      observable_harmonics, coordinates)
        scalars = {str(k): [complex(op.terms.get((), 0)).real, complex(op.terms.get((), 0)).imag]
                   for k, op in harmonics.items() if k and () in op.terms}
        diagnostics = dict(fourier_convention=FOURIER_CONVENTION, time_origin=0., samples=samples,
                           max_harmonic=top, retained_harmonics=sorted(harmonics), harmonic_l1=l1,
                           hermiticity_error=float(hermiticity), cleanup_bound=float(cleanup),
                           coefficient_tolerance=coefficient_tolerance, scalar_harmonics=scalars,
                           sample_cleanup_bound=float(max(h.metadata.get("coefficient_cleanup_bound", 0.)
                                                          for h in sampled)))
        metadata = dict(extra)
        metadata["floquet"] = diagnostics
        result = cls(omega, harmonics, template.nimp, template.spins, template.eta_labels, reference,
                     observable_harmonics, metadata)
        if verify:
            # The midpoints between samples test aliasing, harmonic truncation and cleanup at once.
            error = 0.
            for j in range(samples):
                t = (j+.5)*period/samples
                exact = _sample(function, t)
                check(exact)
                error = max(error, sum(abs(v) for v in (exact.operator - result.operator_at(t)).terms.values()))
            diagnostics.update(interpolation_error=float(error), sampling_tolerance=sampling_tolerance)
            if error > sampling_tolerance:
                raise ValueError(f"Fourier interpolation error {error:.3e} exceeds "
                                 f"sampling_tolerance={sampling_tolerance:.3e}; "
                                 "increase samples")
        return result


def _sample(function, t):
    try:
        h = function(t)
    except (ValueError, TypeError, ArithmeticError) as error:
        raise ValueError(f"model construction failed at t={t:.6g}: {error}") from error
    if not isinstance(h, Hamiltonian):
        raise TypeError("a periodic model must return a Hamiltonian")
    return h


def _coordinate_check(template, assume_fixed):
    labels = (template.nimp, template.spins, template.eta_labels)
    if not assume_fixed and "coordinates" not in template.metadata:
        raise ValueError("the model does not declare metadata['coordinates']; declare its canonical "
                         "coordinates or pass assume_fixed_coordinates=True")
    expected = None if assume_fixed else template.coordinate_fingerprint()

    def check(h):
        if (h.nimp, h.spins, h.eta_labels) != labels:
            raise ValueError("mode labels differ: " + _MOVING_COORDINATES)
        if expected is not None and ("coordinates" not in h.metadata or h.coordinate_fingerprint() != expected):
            raise ValueError(_MOVING_COORDINATES)
    return check


__all__ = ["FOURIER_CONVENTION", "Harmonic", "ParameterDrive", "PeriodicHamiltonian",
           "fold_quasienergies", "quasienergy_distance"]
