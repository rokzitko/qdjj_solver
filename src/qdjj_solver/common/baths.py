"""Positive reservoir quadratures, surrogate fits and Padé chain expansions."""

from __future__ import annotations

from dataclasses import dataclass, field
from time import perf_counter

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import least_squares
from scipy.special import roots_legendre


@dataclass(frozen=True)
class DiscreteBath:
    """Spin-degenerate normal levels and contact weights for one BCS reservoir.

    ``weights`` represents rho d xi, not the quasiparticle density of states.
    Each normal level is a spinful orbital. All QP excitation energies are
    positive. Surrogate contact weights need not sum to one at finite size.
    """

    xi: np.ndarray
    weights: np.ndarray
    delta: float = 1.
    bandwidth: float = 100.
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        xi, weights = np.array(self.xi, dtype=float, copy=True), np.array(self.weights, dtype=float, copy=True)
        if xi.ndim != 1 or xi.size == 0 or xi.shape != weights.shape:
            raise ValueError("xi and weights must be nonempty one-dimensional arrays of equal size")
        if not np.all(np.isfinite(xi)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
            raise ValueError("bath energies must be finite and contact weights strictly positive")
        if not np.isfinite(self.delta) or self.delta <= 0 or not np.isfinite(self.bandwidth) or self.bandwidth <= 0:
            raise ValueError("gap and half bandwidth must be positive and finite")
        xi.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "xi", xi)
        object.__setattr__(self, "weights", weights)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def levels(self):
        return self.xi.size

    @property
    def rho(self):
        return 1 / (2*self.bandwidth)

    @property
    def energies(self):
        return np.hypot(self.xi, self.delta)

    @property
    def paired(self):
        """Exact mirrored coordinates and weights, as required for exact mode reduction."""
        return (np.array_equal(self.xi, -self.xi[::-1])
                and np.array_equal(self.weights, self.weights[::-1]))

    def record(self):
        return dict(xi=self.xi.tolist(), weights=self.weights.tolist(), delta=self.delta,
                    bandwidth=self.bandwidth, metadata=self.metadata)

    @classmethod
    def from_record(cls, record):
        return cls(**record)


def cosh_grid(pairs: int, delta=1., bandwidth=100.) -> DiscreteBath:
    """Gauss--Legendre quadrature in xi=±Delta sinh(t), E=Delta cosh(t).

    There are ``2*pairs`` signed normal levels *per reservoir*. Coarse rules
    retain their actual quadrature weights; silently normalizing their sum
    would change the discretization. The zeroth-moment error is recorded.
    """
    if not isinstance(pairs, (int, np.integer)) or pairs < 1:
        raise ValueError("pairs must be a positive integer")
    if not np.isfinite(delta) or not np.isfinite(bandwidth) or delta <= 0 or bandwidth <= 0:
        raise ValueError("gap and half bandwidth must be positive and finite")
    points, weights = np.polynomial.legendre.leggauss(pairs)
    upper = np.arcsinh(bandwidth/delta)
    t = upper*(points + 1)/2
    positive = delta*np.sinh(t)
    weights = weights*upper/2*delta*np.cosh(t)/(2*bandwidth)
    return DiscreteBath(np.r_[-positive[::-1], positive], np.r_[weights[::-1], weights],
                        delta, bandwidth, dict(kind="cosh-grid", positive_pairs=int(pairs),
                                               measure_error=float(2*weights.sum() - 1)))


def hybridization_g(frequency, delta=1., bandwidth=100.):
    frequency = np.asarray(frequency, dtype=float)
    a = np.hypot(frequency, delta)
    return 2/np.pi*np.arctan(bandwidth/a)/a


def discrete_g(bath, frequency):
    frequency = np.asarray(frequency, dtype=float)
    return np.sum((bath.weights/(np.pi*bath.rho)) /
                  (frequency[..., None]**2 + bath.energies**2), axis=-1)


def _measure_jacobi(nodes, weights, size):
    """Jacobi recurrence of a positive discrete measure, with reorthogonalization."""
    vector = np.sqrt(weights/weights.sum())
    history = np.empty((size, len(nodes)))
    diagonal, off = [], []
    previous, beta = np.zeros_like(vector), 0.
    for j in range(size):
        history[j] = vector
        residual = nodes*vector-beta*previous
        alpha = float(vector@residual)
        diagonal.append(alpha)
        if j == size-1:
            break
        residual -= alpha*vector
        for _ in range(2):
            residual -= (history[:j+1]@residual)@history[:j+1]
        beta = float(np.linalg.norm(residual))
        if beta == 0:
            raise ValueError("chain expansion lost numerical rank")
        off.append(beta)
        previous, vector = vector, residual/beta
    return np.array(diagonal), np.array(off)


def chain_expansion(levels: int, delta=1., bandwidth=100., *, wide_band=False) -> DiscreteBath:
    """Bobok--Frk--Pokorny--Zonda ChE(L), returned in exact mirrored star coordinates.

    ``levels`` is the number of spinful sites in one reservoir, not a QP cutoff.
    The default targets a finite flat band. ``wide_band=True`` uses Eq. (17)'s
    analytic infinite-band coefficients; then ``bandwidth`` only specifies the
    density convention rho=1/(2D), and cancels from physical hybridizations.

    Finite-band Padé coefficients are obtained from equivalent positive Gauss
    (even L) or endpoint-Radau (odd L) quadrature, avoiding a Hankel moment solve.
    Coarse contact weights are NOT normalized to unit normal-state measure.
    """
    if (not isinstance(levels, (int, np.integer)) or isinstance(levels, (bool, np.bool_))
            or levels < 1):
        raise ValueError("levels must be a positive integer")
    if not isinstance(wide_band, (bool, np.bool_)):
        raise ValueError("wide_band must be boolean")
    if not np.all(np.isfinite([delta, bandwidth])) or delta <= 0 or bandwidth <= 0:
        raise ValueError("gap and half bandwidth must be positive and finite")
    levels = int(levels)
    if wide_band:
        ell = np.arange(1, levels, dtype=float)
        h = np.r_[float(levels), (levels**2-ell**2)/(4*ell**2-1)]
        nodes, vectors = eigh_tridiagonal(np.zeros(levels), np.sqrt(h[1:]))
        residues = h[0]*vectors[0]**2
        # The bipartite chain has exact +/- eigenpairs. Preserve that identity
        # in the stored coordinates despite roundoff in the eigendecomposition.
        nodes = (nodes-nodes[::-1])/2
        residues = (residues+residues[::-1])/2
        order = 0
    else:
        angle = float(np.arctan2(bandwidth, delta))
        order = max(64, 4*levels)
        points, measure = roots_legendre(order)
        theta = angle*(points+1)/2
        measure = measure*angle/np.pi
        # z in [0,1] stays well scaled even for D << Delta. The Stieltjes
        # measure is (2/pi)dtheta; its kernel is 1/[1+x*cos(theta)**2].
        scale = np.sin(angle)**2
        z = (np.sin(theta)/np.sin(angle))**2
        pairs, odd = divmod(levels, 2)
        if pairs:
            quadrature_measure = z*measure if odd else measure
            diagonal, off = _measure_jacobi(z, quadrature_measure, pairs)
            nodes, vectors = eigh_tridiagonal(diagonal, off)
            masses = quadrature_measure.sum()*vectors[0]**2
            if odd:
                masses /= nodes
            r = scale*nodes
            positive = np.sqrt(r/(1-r))
            residues = masses/(2*(1-r))
        else:
            positive, residues, masses = np.array([]), np.array([]), np.array([])
        center = np.array([measure.sum()-masses.sum()]) if odd else np.array([])
        nodes = np.r_[-positive[::-1], np.zeros(odd), positive]
        residues = np.r_[residues[::-1], center, residues]
        if np.any(residues <= 0) or not np.all(np.isfinite(nodes)):
            raise ValueError("chain expansion is not resolved at these energy scales")
        # Recover the chain for diagnostics and independent chain-space solves.
        _, hopping = _measure_jacobi(nodes, residues, levels)
        h = np.r_[residues.sum(), hopping**2]
    weights = np.pi/2*(delta/bandwidth)*residues
    return DiscreteBath(delta*nodes, weights, delta, bandwidth, dict(
        kind="chain-expansion", levels=levels, wide_band=bool(wide_band),
        target="wide-band" if wide_band else "finite-band", chain_h=h.tolist(),
        coefficient_quadrature_order=order, measure_error=float(weights.sum()-1),
        citation="Bobok, Frk, Pokorny and Zonda, Phys. Rev. B 112, 205418 (2025), 10.1103/mxsl-fc96"))


def fit_surrogate(levels: int, delta=1., bandwidth=100., frequency_cutoff=100.,
                  frequency_min=1e-3, frequency_points=1000, starts=4, seed=1729,
                  relative_weight=0., *, max_nfev=5000) -> DiscreteBath:
    """Positive, particle-hole symmetric least-squares fit of the BCS g function.

    The default ``relative_weight=0`` is the absolute least-squares criterion
    on a logarithmic frequency mesh used by Baran, Frost and Paaske (2023).
    A nonzero exponent changes the criterion and is therefore recorded.
    ``frequency_min`` and ``frequency_cutoff`` are in the Hamiltonian's units.
    ``max_nfev`` bounds function evaluations per start, not the fitting accuracy.
    """
    if not isinstance(levels, (int, np.integer)) or levels < 1:
        raise ValueError("levels must be a positive integer")
    if starts < 1 or frequency_points < 10 or not 0 < frequency_min < frequency_cutoff:
        raise ValueError("invalid surrogate fitting mesh or number of starts")
    if (not isinstance(max_nfev, (int, np.integer)) or isinstance(max_nfev, (bool, np.bool_))
            or max_nfev < 1):
        raise ValueError("max_nfev must be a positive integer")
    if delta <= 0 or bandwidth <= 0 or not np.all(np.isfinite([delta, bandwidth, frequency_min, frequency_cutoff, relative_weight])):
        raise ValueError("invalid energy scales")
    start_time = perf_counter()
    omega = np.geomspace(frequency_min, frequency_cutoff, frequency_points)
    target = hybridization_g(omega, delta, bandwidth)
    scaling = target**relative_weight * target.max()**(1-relative_weight)
    pairs, odd = levels//2, levels % 2
    multiplicity = np.r_[1. if odd else [], np.full(pairs, 2.)]
    count = pairs + odd

    def evaluate(parameters, jacobian=False):
        residues = np.exp(parameters[:count])
        positive = np.exp(parameters[count:])
        xi = np.r_[0. if odd else [], positive]
        den = omega[:, None]**2 + delta**2 + xi[None, :]**2
        terms = multiplicity*residues/den
        if jacobian:
            return np.c_[terms, -2*terms[:, odd:]*positive**2/den[:, odd:]]/scaling[:, None]
        return (terms.sum(axis=1) - target)/scaling

    rng = np.random.default_rng(seed)
    fits = []
    for trial in range(starts):
        positive = np.geomspace(0.5*delta, min(bandwidth, frequency_cutoff)*0.75, max(pairs, 2))[:pairs]
        if pairs == 1:
            positive = np.array([1.5*delta])
        positive *= np.exp(rng.normal(0, .35, pairs)) if trial else 1
        residues = np.r_[delta/2 if odd else [], np.maximum(delta/2, positive)]
        initial = np.log(np.r_[residues, positive])
        low = np.log(np.full(initial.size, 1e-10*min(delta, bandwidth)))
        high = np.log(np.full(initial.size, 1e3*max(delta, bandwidth, frequency_cutoff)))
        fit = least_squares(evaluate, initial, jac=lambda p: evaluate(p, True), bounds=(low, high),
                            ftol=2e-13, xtol=2e-13, gtol=2e-13, max_nfev=int(max_nfev), x_scale="jac")
        fits.append(fit)
    trials = []
    for f in fits:
        with np.errstate(over="ignore", invalid="ignore"):
            cost = float(np.dot(f.fun, f.fun))
        trials.append(dict(success=bool(f.success), status=int(f.status), message=str(f.message),
                           evaluations=int(f.nfev), optimality=float(f.optimality) if np.isfinite(f.optimality) else None,
                           cost=cost if np.isfinite(cost) else None))
    successful = [(i, f) for i, f in enumerate(fits)
                  if f.success and np.all(np.isfinite(f.x)) and trials[i]["cost"] is not None]
    if not successful:
        details = "; ".join(f"start {i}: status={t['status']}, nfev={t['evaluations']}, {t['message']}"
                            for i, t in enumerate(trials))
        raise RuntimeError(f"surrogate fit did not converge (max_nfev={max_nfev}; {details})")
    selected, fit = min(successful, key=lambda item: trials[item[0]]["cost"])
    residues = np.exp(fit.x[:count])
    positive = np.exp(fit.x[count:])
    xi = np.r_[-positive, 0. if odd else [], positive]
    gamma = np.r_[residues[odd:], residues[:1] if odd else [], residues[odd:]]
    order = np.argsort(xi)
    weights = np.pi/(2*bandwidth)*gamma[order]
    metadata = dict(kind="surrogate", levels=int(levels), frequency_min=frequency_min,
                    frequency_cutoff=frequency_cutoff, frequency_points=frequency_points,
                    relative_weight=relative_weight, starts=starts, seed=seed,
                    max_nfev=int(max_nfev), selected_start=int(selected), optimizer_trials=trials,
                    fit_seconds=perf_counter()-start_time, cost=trials[selected]["cost"],
                    fit_optimality=float(fit.optimality), evaluations=int(fit.nfev),
                    measure_error=float(weights.sum()-1),
                    citation="Baran, Frost and Paaske, Phys. Rev. B 108, L220506 (2023)")
    bath = DiscreteBath(xi[order], weights, delta, bandwidth, metadata)
    metadata["max_relative_fit_error"] = float(np.max(abs(discrete_g(bath, omega)/target - 1)))
    return DiscreteBath(bath.xi, bath.weights, delta, bandwidth, metadata)
