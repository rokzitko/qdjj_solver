"""Parity-changing excitation and current, with an explicit finite-bath QP ladder.

Requires only the base installation. All spin projections are retained in each
parity sector; the final small full-space calculation uses dense diagonalization.
"""

import numpy as np

from qdjj_solver import Sector, cosh_grid, reference_model, solve


def main():
    bath = cosh_grid(pairs=1, delta=1., bandwidth=5.)
    parameters = dict(u=1.6, gamma=.3, phi=.8, compress=False)
    h = reference_model(bath, **parameters)
    sectors = [Sector(parity=p, twice_sz=None) for p in (0, 1)]

    previous = np.full(2, np.inf)
    print("All spin projections retained; energies in units of Delta=1.")
    print("cutoff       E_even           E_odd        E_odd-E_even     max residual")
    for cutoff in (0, 1, 2, 4, None):
        results = [solve(h, cutoff=cutoff, sector=sector,
                         options={"method": "dense" if cutoff is None else "auto"})
                   for sector in sectors]
        energies = np.array([result.energies[0] for result in results])
        residual = max(result.residuals[0] for result in results)
        label = "full" if cutoff is None else str(cutoff)
        print(f"{label:>6}  {energies[0]: .12f}  {energies[1]: .12f}"
              f"  {energies[1]-energies[0]: .12f}    {residual:.2e}")
        # Each sector's lowest energy is variational under nested cutoffs.
        assert np.all(energies <= previous + 1e-10)
        previous = energies

    ground_parity = int(np.argmin(energies))
    parity_gap = energies[1-ground_parity]-energies[ground_parity]
    ground = results[ground_parity]
    current = float(ground.observables["phase_derivative"][0].real)
    print(f"Full finite-model ground parity: {('even', 'odd')[ground_parity]}")
    print(f"Parity-changing excitation energy / Delta: {parity_gap:.12f}")
    print(f"Ground-state current / (2 e Delta / hbar): {current:.12f}")

    # Independent finite-difference check of the current in this smooth sector.
    # Rebuild the physical model at each phase, retaining the same bath/reference.
    step = 1e-5
    shifted = []
    for sign in (-1, 1):
        shifted_h = reference_model(bath, **(parameters | {"phi": parameters["phi"]+sign*step}))
        state = solve(shifted_h, cutoff=None, sector=sectors[ground_parity], options={"method": "dense"})
        shifted.append(state.energies[0])
    finite_difference = (shifted[1]-shifted[0])/(2*step)
    np.testing.assert_allclose(current, finite_difference, rtol=0, atol=2e-8)
    print("Current agrees with the finite-difference energy derivative.")
    print("Bath convergence remains to be checked by refining nodes at fixed physical parameters.")


if __name__ == "__main__":
    main()
