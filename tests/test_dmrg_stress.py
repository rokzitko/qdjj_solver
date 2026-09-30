"""One bounded nightly check of a long, reordered, ragged DMRG layout."""

from math import comb

import numpy as np
from numpy.testing import assert_allclose
import pytest
from threadpoolctl import threadpool_limits

pytest.importorskip("tenpy")

from qdjj_solver.common import Hamiltonian, Sector, annihilate, create, hermitian_pair, number
from qdjj_solver.dmrg_solver import SolverOptions, solve


pytestmark = [pytest.mark.dmrg, pytest.mark.slow, pytest.mark.convergence, pytest.mark.timeout(90)]


def test_odd_reordered_chain_against_analytic_spectrum():
    # Hopping sqrt((j+1)*(N-1-j)) realizes -2*t*J_x. The lowest two
    # eigenvectors are a binomial amplitude and its linear Krawtchouk polynomial.
    nmodes, seed = 41, 2207
    rng = np.random.default_rng(seed)
    gauge = np.exp(1j*rng.uniform(-np.pi, np.pi, nmodes))
    onsite, constant = 4.+rng.uniform(0., .2), 7.+rng.uniform(0., 1.)
    hopping = rng.uniform(.035, .055)
    positions = np.arange(nmodes)
    strengths = hopping*np.sqrt(np.arange(1, nmodes)*np.arange(nmodes-1, 0, -1))
    energies = constant+onsite+hopping*(2*np.arange(2)-(nmodes-1))
    ground = np.sqrt(np.array([comb(nmodes-1, j) for j in positions], dtype=float)/2.**(nmodes-1))
    vectors = gauge[:, None]*ground[:, None]*np.column_stack((
        np.ones(nmodes), (nmodes-1-2*positions)/np.sqrt(nmodes-1)))
    edges = -strengths*gauge[:-1]*gauge[1:].conj()
    matrix = (constant+onsite)*np.eye(nmodes, dtype=complex)
    matrix += np.diag(edges, 1)+np.diag(edges.conj(), -1)
    assert_allclose(matrix @ vectors, vectors*energies, atol=1e-12, rtol=0.)
    assert_allclose(vectors.conj().T @ vectors, np.eye(2), atol=1e-12, rtol=0.)

    operator = constant+sum(onsite*number(j) for j in range(nmodes))
    operator += sum(hermitian_pair(edge*create(j)*annihilate(j+1))
                    for j, edge in enumerate(edges))
    left, right = nmodes//2-2, nmodes//2+1
    hamiltonian = Hamiltonian(
        operator, 1, tuple(1 if j % 2 == 0 else -1 for j in range(nmodes)),
        observables={"center_density": number(nmodes//2),
                     "coherence": create(left)*annihilate(right),
                     "position": sum(j/(nmodes-1)*number(j) for j in range(nmodes))})
    options = SolverOptions(eigenpairs=2, chi_max=4, group_size=2,
                            mode_order=tuple(reversed(range(nmodes))),
                            seed=seed, seed_trials=1, excitation_operator="position",
                            min_sweeps=8, max_sweeps=28, mixer_sweeps=4,
                            energy_tolerance=1e-12, residual_tolerance=2e-8,
                            require_convergence=True, threads=1)
    with threadpool_limits(limits=1):
        result = solve(hamiltonian, sector=Sector(1, None, particle_number=1), options=options)
    assert len(result.prepared.groups[-1]) == 1
    assert_allclose(result.energies, energies, atol=3e-10, rtol=0.)
    assert_allclose(result.overlaps(), np.eye(2), atol=2e-9, rtol=0.)
    assert_allclose(result.residuals, 0., atol=2e-8, rtol=0.)
    assert_allclose(result.observables["center_density"], abs(vectors[nmodes//2])**2,
                    atol=2e-7, rtol=0.)
    assert_allclose(result.observables["coherence"], vectors[left].conj()*vectors[right],
                    atol=2e-7, rtol=0.)
    assert_allclose(result.observables["position"], positions/(nmodes-1) @ abs(vectors)**2,
                    atol=2e-7, rtol=0.)
