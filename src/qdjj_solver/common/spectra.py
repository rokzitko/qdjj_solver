"""Transition tables, spectral lines and line shapes for weak probes of driven states.

For Floquet states psi(t)=exp(-i eps t) sum_m exp(i m Omega t) u_m and a probe
X(t)=sum_k X_k exp(i k Omega t), with Omega the drive frequency (attribute ``omega``),
<psi_a(t)|X(t)|psi_b(t)> = sum_q X^(q)_ab exp[i(eps_a - eps_b + q Omega) t],
X^(q)_ab = sum_{m,k} <u_{a,m+k-q}|X_k|u_{b,m}>.
A line at omega_ab^(q) = eps_a - eps_b + q Omega carries the strength |X^(q)_ab|^2;
a positive frequency absorbs a probe quantum, a negative one corresponds to
stimulated emission into the probe. Spectra are functions of the probe frequency omega.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

SHAPES = ("lorentzian", "gaussian")
KINDS = ("signed", "rate", "net")


def line_shape(x, width, shape="lorentzian"):
    """Unit-area line shape; ``width`` is the half width at half maximum."""
    if not np.isfinite(width) or width <= 0:
        raise ValueError("the line width must be positive and finite")
    x = np.asarray(x, float)
    if shape == "lorentzian":
        return width/np.pi/(x*x + width*width)
    if shape == "gaussian":
        sigma = width/np.sqrt(2*np.log(2))
        return np.exp(-x*x/(2*sigma*sigma))/(sigma*np.sqrt(2*np.pi))
    raise ValueError(f"shape must be one of: {', '.join(SHAPES)}")


def degenerate_groups(values, omega=None, tolerance=1e-8):
    """Assign common labels to (quasi)energies that agree within ``tolerance`` (modulo ``omega``)."""
    values = np.asarray(values, float)
    if values.ndim != 1:
        raise ValueError("degenerate groups need a one-dimensional array of energies")
    if len(values) == 0:
        return np.zeros(0, dtype=int)
    folded = values if omega is None else np.mod(values, omega)
    order = np.argsort(folded, kind="stable")
    labels = np.empty(len(values), dtype=int)
    group = 0
    labels[order[0]] = 0
    for previous, current in zip(order[:-1], order[1:], strict=True):
        if folded[current] - folded[previous] > tolerance:
            group += 1
        labels[current] = group
    if omega is not None and group and folded[order[0]] + omega - folded[order[-1]] <= tolerance:
        labels[labels == group] = labels[order[0]]
    return labels


@dataclass
class SpectralLines:
    """Discrete lines omega_l with strengths s_l, weighted by the initial populations.

    ``final`` and ``initial`` list the state indices merged into each line (degenerate
    groups), ``photons`` the photon number q, and ``elastic`` is true for lines between
    states of one group at q != 0 (the central line of a Mollow triplet).
    """

    frequencies: np.ndarray
    """Sorted line frequencies; positive values correspond to absorption of a probe quantum."""
    strengths: np.ndarray
    """Line strengths, weighted by the initial populations."""
    final: list
    """Final-state indices merged into each line (a degenerate group)."""
    initial: list
    """Initial-state indices that contribute to each line."""
    photons: np.ndarray
    """Photon number q of each line."""
    elastic: np.ndarray
    """True for elastic lines, between states of one degenerate group at q != 0."""
    metadata: dict = field(default_factory=dict)
    """Probe, initial states, weights and merging options."""

    def __len__(self):
        return len(self.frequencies)

    def spectrum(self, frequencies, width, *, shape="lorentzian", kind="signed"):
        """Broadened spectrum S(omega) = sum_l s_l L(omega - omega_l) with a unit-area line shape L.

        ``kind="signed"`` returns S(omega); ``kind="rate"`` returns S(omega)+S(-omega),
        proportional to the golden-rule transition rate out of the initial states at
        probe frequency omega>0; ``kind="net"`` returns S(omega)-S(-omega), absorption
        minus stimulated emission. ``width`` is the half width at half maximum.
        """
        if kind not in KINDS:
            raise ValueError(f"kind must be one of: {', '.join(KINDS)}")
        omega = np.asarray(frequencies, float)
        def signed(x):
            return np.sum(self.strengths[None, :]*line_shape(x[:, None] - self.frequencies[None, :], width, shape),
                          axis=1) if len(self) else np.zeros(len(x))
        if kind == "signed":
            return signed(omega)
        if kind == "rate":
            return signed(omega) + signed(-omega)
        return signed(omega) - signed(-omega)

    def record(self):
        """JSON-ready line positions, strengths, contributing state indices and metadata."""
        return dict(frequencies=self.frequencies.tolist(), strengths=self.strengths.tolist(),
                    final=[list(map(int, f)) for f in self.final], initial=[list(map(int, i)) for i in self.initial],
                    photons=self.photons.tolist(), elastic=self.elastic.tolist(), metadata=self.metadata)


@dataclass
class FloquetTransitions:
    """Probe matrix elements X^(q) between two sets of Floquet states of one drive.

    ``elements[i, a, b]`` is X^(q) for ``photons[i]``, bra state ``bra_indices[a]``
    (final) and ket state ``ket_indices[b]`` (initial) of the two FloquetStates.
    ``total[b]`` is ``sum_m ||(X u_b)_m||^2``, the right-hand side of the sum rule
    ``sum_a sum_q |X^(q)_ab|^2 = total[b]``, which holds for a complete set of bra states.
    """

    omega: float
    """Drive frequency."""
    photons: np.ndarray
    """Photon numbers q of the stored harmonics."""
    elements: np.ndarray
    """``elements[i, a, b]`` is X^(q) for ``photons[i]``, final state a and initial state b."""
    bra_quasienergies: np.ndarray
    """Quasienergies of the final (bra) states."""
    ket_quasienergies: np.ndarray
    """Quasienergies of the initial (ket) states."""
    bra_labels: np.ndarray
    """Static labels of the final states."""
    ket_labels: np.ndarray
    """Static labels of the initial states."""
    bra_indices: np.ndarray
    """Indices of the final states in their ``FloquetStates``."""
    ket_indices: np.ndarray
    """Indices of the initial states in their ``FloquetStates``."""
    total: np.ndarray
    """Time-averaged squared norm of X(t) applied to each initial state: the right-hand side of the sum rule."""
    same_states: bool
    """True when bra and ket come from the same ``FloquetStates`` object."""
    complete: bool
    """True when the bra states form the complete set of their representation."""
    default_initial: int
    """Initial state used by ``lines()`` without arguments: the dressed ground state of the ket sector,
    ``dressed_state(0)``, or the first ket state if that state was not selected."""
    probe: str
    """Name of the probe."""
    metadata: dict = field(default_factory=dict)
    """Probe harmonics, sectors, methods and the matrix-element convention."""

    @property
    def frequencies(self):
        """Line frequencies eps_a - eps_b + q Omega, with the shape of ``elements``."""
        return (self.bra_quasienergies[None, :, None] - self.ket_quasienergies[None, None, :]
                + self.photons[:, None, None]*self.omega)

    @property
    def strengths(self):
        """Line strengths ``|X^(q)_ab|^2``, with the shape of ``elements``."""
        return np.abs(self.elements)**2

    @property
    def captured(self):
        """Total transition strength from each initial state to the computed final (bra) states."""
        return self.strengths.sum(axis=(0, 1))

    @property
    def completeness(self):
        """``captured/total``: the fraction of the sum rule exhausted by the computed bra states."""
        return self.captured/np.maximum(self.total, np.finfo(float).tiny)

    def element(self, q):
        """The matrix X^(q), shape (bra states, ket states); zero if q is not stored."""
        positions = np.flatnonzero(self.photons == q)
        if not len(positions):
            return np.zeros(self.elements.shape[1:], complex)
        return self.elements[positions[0]]

    def lines(self, initial=None, weights=None, *, merge_degenerate=True, elastic=True, threshold=None,
              tolerance=1e-8):
        """Spectral lines originating from the specified initial states (indices of the ket result).

        Default: ``default_initial`` with unit weight. Merging distributes each initial
        weight uniformly over its degenerate group and sums the strengths of degenerate
        final states, which makes the lines independent of the basis chosen within
        degenerate multiplets. Terms between states of one group at q=0 are expectation
        values, not transitions, and are dropped. Lines weaker than ``threshold`` (default:
        1e-12 of the strongest line, i.e. numerical zeros) are omitted.
        """
        if merge_degenerate and self.metadata.get("multiplet_split"):
            raise ValueError("the computed states split a degenerate static multiplet; increase states or basis_states")
        initial = [self.default_initial] if initial is None else list(np.atleast_1d(initial))
        positions = []
        for state in initial:
            found = np.flatnonzero(self.ket_indices == state)
            if not len(found):
                raise ValueError(f"initial state {state} is not among the ket states")
            positions.append(int(found[0]))
        weights = np.full(len(positions), 1./len(positions)) if weights is None else np.asarray(weights, float)
        if weights.shape != (len(positions),) or np.any(weights < 0) or not np.all(np.isfinite(weights)):
            raise ValueError("give one nonnegative finite weight per initial state")
        population = np.zeros(len(self.ket_indices))
        ket_groups = (degenerate_groups(self.ket_quasienergies, self.omega, tolerance) if merge_degenerate
                      else np.arange(len(self.ket_indices)))
        bra_groups = (degenerate_groups(self.bra_quasienergies, self.omega, tolerance) if merge_degenerate
                      else np.arange(len(self.bra_indices)))
        for position, weight in zip(positions, weights, strict=True):
            members = np.flatnonzero(ket_groups == ket_groups[position])
            population[members] += weight/len(members)
        # Bra and ket from one set of states: identify their groups through the state indices.
        identical = None
        if self.same_states:
            lookup = {int(state): int(g) for state, g in zip(self.bra_indices, bra_groups, strict=True)}
            identical = np.array([lookup.get(int(state), -1) for state in self.ket_indices])
        merged = {}
        frequencies, strengths = self.frequencies, self.strengths
        for b in np.flatnonzero(population > 0):
            for i, q in enumerate(self.photons):
                for a in range(len(self.bra_indices)):
                    value = population[b]*strengths[i, a, b]
                    if value == 0:
                        continue
                    same_group = identical is not None and identical[b] == bra_groups[a]
                    if same_group and q == 0:
                        continue
                    if same_group and not elastic:
                        continue
                    key = (int(bra_groups[a]), int(ket_groups[b]), int(q))
                    entry = merged.setdefault(key, dict(strength=0., frequency=[], final=set(), initial=set(),
                                                        elastic=bool(same_group), photons=int(q)))
                    entry["strength"] += value
                    entry["frequency"].append(frequencies[i, a, b])
                    entry["final"].add(int(self.bra_indices[a]))
                    entry["initial"].add(int(self.ket_indices[b]))
        if threshold is None:
            threshold = 1e-12*max((e["strength"] for e in merged.values()), default=0.)
        rows = sorted(((float(np.mean(e["frequency"])), e) for e in merged.values() if e["strength"] > threshold),
                      key=lambda row: row[0])
        return SpectralLines(np.array([f for f, _ in rows]), np.array([e["strength"] for _, e in rows]),
                             [tuple(sorted(e["final"])) for _, e in rows], [tuple(sorted(e["initial"])) for _, e in rows],
                             np.array([e["photons"] for _, e in rows], dtype=int),
                             np.array([e["elastic"] for _, e in rows], dtype=bool),
                             dict(probe=self.probe, initial=[int(s) for s in initial], weights=weights.tolist(),
                                  merge_degenerate=merge_degenerate))

    def spectrum(self, frequencies, width, *, shape="lorentzian", kind="signed", **line_options):
        """Broadened spectrum of ``lines(**line_options)``; see ``SpectralLines.spectrum``."""
        return self.lines(**line_options).spectrum(frequencies, width, shape=shape, kind=kind)


__all__ = ["FloquetTransitions", "SpectralLines", "degenerate_groups", "line_shape"]
