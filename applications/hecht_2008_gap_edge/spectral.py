"""Numerical spectral constructions for the Hecht gap-edge study.

The quadratic resolvent eliminates the declared DiscreteBath exactly. The
interacting calculation uses native QP operators and a reorthogonalized Lanczos
measure, including electron addition AND removal in both spin channels.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal

from qdjj_solver import DiscreteBath, Sector, annihilate, create, reference_model
from qdjj_solver.qp_solver.basis import FockBasis
from qdjj_solver.qp_solver.solver import ProjectedOperator
from applications._common import parity_states


def logarithmic_bath(cells, bandwidth, order=4):
    """Positive composite Gauss quadrature of rho dxi, including the central cell."""
    if cells < 2 or order < 1 or bandwidth <= 1e-9:
        raise ValueError("positive order, at least two cells and bandwidth > 1e-9 required")
    edges = np.r_[0., np.geomspace(1e-9, bandwidth, cells)]
    nodes, weights = np.polynomial.legendre.leggauss(order)
    mid, half = (edges[1:]+edges[:-1])/2, np.diff(edges)/2
    xi = (mid[:, None]+half[:, None]*nodes).ravel()
    measure = (half[:, None]*weights/(2*bandwidth)).ravel()
    return DiscreteBath(np.r_[-xi[::-1], xi], np.r_[measure[::-1], measure],
                        delta=1., bandwidth=bandwidth,
                        metadata=dict(kind="composite-log-gauss", cells=cells, order=order,
                                      central_edge=1e-9, measure_error=float(2*measure.sum()-1)))


def quadratic_green(z, gamma, epsilon=0., *, bath=None, bandwidth=np.inf):
    """Retarded dot G_up for one real-gap reservoir in Delta=1 units."""
    z = np.asarray(z, complex)
    if np.any(z.imag <= 0):
        raise ValueError("retarded resolvent requires Im(z)>0")
    if bath is None:
        root = np.sqrt(1-z*z)
        g = 1/root if np.isinf(bandwidth) else 2/np.pi*np.arctan(bandwidth/root)/root
    else:
        if not bath.paired or bath.delta != 1.:
            raise ValueError("this reduction requires a mirrored bath with Delta=1")
        g = np.empty(z.size, complex)
        flat = z.ravel()
        for start in range(0, z.size, 32):
            stop = start+32
            # (1-z)*(1+z) preserves more gap-edge precision than 1-z**2.
            den = bath.xi**2+(1-flat[start:stop, None])*(1+flat[start:stop, None])
            g[start:stop] = np.sum(bath.weights/(np.pi*bath.rho)/den, axis=1)
        g = g.reshape(z.shape)
    s = gamma*g
    diagonal = z*(1+s)
    # Factored determinant avoids subtracting two large, almost equal squares.
    determinant = (z+(z-1)*s)*(z+(z+1)*s)-epsilon**2
    return (diagonal+epsilon)/determinant


def analytic_continuum(offset, gamma, epsilon=0.):
    """Eq. (28) implies Eq. (30) with (2 Gamma omega rho)^2 in its denominator."""
    offset = np.asarray(offset, float)
    if np.any(offset <= 0):
        raise ValueError("positive distance from the gap edge required")
    omega = 1+offset
    rho = omega/np.sqrt(offset*(2+offset))
    return gamma*rho/np.pi*((omega+epsilon)**2+gamma**2)/(
        (omega**2-epsilon**2-gamma**2)**2+(2*gamma*omega*rho)**2)


def lanczos_measure(matrix, source, steps, max_memory_gib=1.):
    """Positive quadrature of <source|delta(E-H)|source>, with two-pass reorthogonalization."""
    source = np.real_if_close(np.asarray(source))
    weight = float(np.vdot(source, source).real)
    if weight < 1e-28:
        return np.array([]), np.array([]), dict(steps=0, weight=weight, terminal_beta=0.)
    steps = min(steps, len(source))
    if steps < 1:
        raise ValueError("positive Lanczos step count required")
    dtype = np.result_type(matrix.dtype, source.dtype)
    if steps*len(source)*np.dtype(dtype).itemsize > max_memory_gib*1024**3:
        raise MemoryError("Lanczos history exceeds the declared memory limit")
    history = np.empty((steps, len(source)), dtype=dtype)
    q, previous, beta = source/np.sqrt(weight), np.zeros_like(source), 0.
    diagonal, off = [], []
    threshold = 64*np.finfo(float).eps*max(1., float(np.max(abs(matrix.data))))
    for j in range(steps):
        history[j] = q
        w = matrix@q-beta*previous
        alpha = float(np.vdot(q, w).real)
        w -= alpha*q
        for _ in range(2):
            w -= (history[:j+1].conj()@w)@history[:j+1]
        beta_new = float(np.linalg.norm(w))
        diagonal.append(alpha)
        if beta_new < threshold or j == steps-1:
            beta = beta_new
            break
        off.append(beta_new)
        previous, q, beta = q, w/beta_new, beta_new
    energies, vectors = eigh_tridiagonal(np.array(diagonal), np.array(off))
    return energies, weight*abs(vectors[0])**2, dict(steps=len(diagonal), weight=weight,
                                                  terminal_beta=beta)


def fock_measure(bath, u, gamma, detuning, settings, steps):
    """Spin-averaged zero-temperature spectral measure, retaining both doublet projections by SU(2)."""
    h = reference_model(bath, u=u, gamma=gamma, detuning=detuning, geometry="single",
                        symmetry=False, compress=False, coefficient_tolerance=0.)
    states = parity_states(h, settings)
    parity = int(states[1].energies[0] < states[0].energies[0])
    ground = states[parity]
    energy = float(ground.energies[0])
    cutoff = ground.basis.cutoff
    destinations, rows, diagnostics = {}, [], []
    for spin, sz in ((0, 1), (1, -1)):
        for addition in (True, False):
            target_sz = parity+(sz if addition else -sz)
            if target_sz not in destinations:
                basis = FockBasis(h.nimp, h.spins, cutoff, Sector(1-parity, target_sz))
                matrix = ProjectedOperator(h.operator, basis).sparse_matrix()
                destinations[target_sz] = basis, matrix
            basis, matrix = destinations[target_sz]
            op = h.physical_operator(create(spin) if addition else annihilate(spin))
            source = ProjectedOperator(op, ground.basis).action_between(ground.vectors[:, 0], basis)
            energies, weights, info = lanczos_measure(matrix, source, steps)
            diagnostics.append(dict(spin=spin, addition=addition, dimension=basis.dimension, **info))
            for value, weight in zip(energies, weights, strict=True):
                rows.append(dict(omega=float((value-energy)*(1 if addition else -1)),
                                 weight=float(weight/2), spin=spin, addition=int(addition)))
    return rows, states, diagnostics


def broaden_continuum(offset, poles, width):
    """Normalized log-Gaussian in omega-Delta; keep subgap poles separate."""
    offset = np.asarray(offset, float)
    if np.any(offset <= 0) or width <= 0:
        raise ValueError("positive offsets and broadening width required")
    result = np.zeros_like(offset)
    for pole in poles:
        distance = pole["omega"]-1
        if distance > 0:
            result += pole["weight"]*np.exp(-(np.log(offset/distance)/width)**2)/(
                np.sqrt(np.pi)*width*offset)
    return result
