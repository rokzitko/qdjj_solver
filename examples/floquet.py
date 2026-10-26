"""Floquet quasienergies of a phase-driven junction at a subgap resonance.

Requires only the base installation. The phase phi(t)=phi0+A cos(omega t) acts
on the hoppings, i.e. a voltage drive with V=hbar (dphi/dt)/2e, where V=V_L-V_R is the
electric potential of the left lead relative to the right one. The frequency is
tuned to the lowest even-parity excitation, so the drive splits the two dressed
states. Energies and omega are in units of Delta=1.
"""

import numpy as np

from qdjj_solver import ParameterDrive, PeriodicHamiltonian, Sector, cosh_grid, reference_model, solve
from qdjj_solver.common.floquet import quasienergy_distance
from qdjj_solver.qp_solver import ProjectedOperator, floquet_propagator, rwa_model, solve_floquet


def main():
    # symmetry=False and compress=False keep the canonical coordinates fixed during the period.
    # At detuning=0 the 0-1 phase-drive element vanishes: both states have opposite eta parity.
    base = dict(bath=cosh_grid(pairs=2, delta=1., bandwidth=5.), u=1.6, gamma=.3, phi=.8,
                detuning=.2, symmetry=False, compress=False)
    cutoff, sector = 2, Sector(parity=0, twice_sz=0)
    undriven = reference_model(**base)
    static = solve(undriven, cutoff, sector, {"method": "dense", "eigenpairs": 40})
    current = ProjectedOperator(undriven.observables["phase_derivative"], static.basis)
    element = abs(current.expectation(static.vectors[:, 0], static.vectors[:, 1]))
    omega, amplitude = static.energies[1] - static.energies[0], .05
    print(f"Static even states E0={static.energies[0]:.6f}, E1={static.energies[1]:.6f}; "
          f"|<0|dH/dphi|1>|={element:.6f}")
    print(f"Drive phi(t)=phi0+{amplitude} cos(omega t), omega=E1-E0={omega:.6f}")
    periodic = PeriodicHamiltonian.from_parameters(reference_model, base, ParameterDrive.cosine("phi", amplitude), omega)

    # A static-eigenbasis cut must not split a degenerate multiplet.
    gaps = np.flatnonzero(np.diff(static.energies) > 1e-6) + 1
    kept = int(gaps[gaps >= 16][0])
    fock = solve_floquet(periodic, cutoff, sector, dict(floquet_cutoff=4, states=2), static=static)
    eigen = solve_floquet(periodic, cutoff, sector, dict(floquet_cutoff=4, states=2, basis="eigen", basis_states=kept),
                          static=static)
    stepped = floquet_propagator(periodic, cutoff, sector, dict(floquet_cutoff=4, states=2, basis="eigen",
                                                                basis_states=kept, steps=128), static=static)
    rwa = rwa_model(periodic, [0, 1], [0, -1], cutoff, sector, static=static)

    print(f"Fock dimension {fock.basis.dimension}; static basis K={kept}; Floquet cutoff M=4")
    print("method                      eps_0          eps_1         splitting   weights on (0, 1)")
    for name, result in (("Sambe, Fock basis", fock), ("Sambe, static basis", eigen),
                         ("propagator, static basis", stepped), ("two-level RWA", rwa)):
        order = np.argsort(result.labels)
        eps = result.quasienergies[order]
        splitting = quasienergy_distance(eps[0], eps[1], omega)
        weights = np.diag(result.static_weights[order][:, :2])
        print(f"{name:<25} {eps[0]: .9f}  {eps[1]: .9f}  {splitting:.3e}   {weights[0]:.3f}, {weights[1]:.3f}")
    splitting = quasienergy_distance(*fock.quasienergies, omega)
    print(f"Weak-drive estimate A|<0|dH/dphi|1>| = {amplitude*element:.3e}")

    # The RWA splitting is the dipole estimate to O(A^2); counter-rotating terms shift Sambe slightly.
    np.testing.assert_allclose(quasienergy_distance(*rwa.quasienergies, omega), amplitude*element, rtol=1e-2)
    np.testing.assert_allclose(splitting, amplitude*element, rtol=5e-2)
    # The truncated static basis is bounded by its residual in the full cutoff Fock space.
    assert np.all(abs(eigen.quasienergies - fock.quasienergies) <= eigen.fock_residuals)
    # Sambe and the one-period propagator agree in the same basis.
    np.testing.assert_allclose(stepped.quasienergies, eigen.quasienergies, rtol=0, atol=1e-8)
    print(f"Largest Fock-space residual of the K={kept} result: {np.max(eigen.fock_residuals):.2e}")
    print("Convergence in M, K, QP cutoff and bath discretization must still be checked separately.")


if __name__ == "__main__":
    main()
