"""Charge and spin convergence of a doublet in a small, fixed finite bath.

Run with the base installation: python -B examples/qp_convergence.py.
See docs/qp_convergence.md for the larger reference study and conventions.
"""

import numpy as np

from qdjj_solver import Sector, annihilate, cosh_grid, create, number, reference_model, solve


def doublet_values(energy, charge, double_occupancy, spin_z):
    """Physical values for a single-orbital, SU(2) doublet with S_z=+1/2.

    The correlation is inferred from angular-momentum addition, not measured
    independently. P_n counts dot electrons; it is not a bath-QP probability.
    Use the constructor's physical observables, including in the eta basis.
    """
    p1 = charge-2*double_occupancy
    qd = 2*spin_z
    return dict(energy=energy, q_d=qd, P_0=1-charge+double_occupancy,
                P_1=p1, P_2=double_occupancy, C_d_bath=.75*(qd-p1))


def values_from_state(state):
    obs = state.observables
    names = ("impurity_charge", "double_occupancy_0", "impurity_spin_z")
    values = [complex(obs[name][0]) for name in names]
    if max(abs(v.imag) for v in values) > 1e-10:
        raise ValueError("physical charge and spin expectations must be real")
    return doublet_values(float(state.energies[0]), *(v.real for v in values))


def add_spin_observables(h):
    """Direct physical-electron spin operators, for an uncompressed single dot.

    The bath spin is the sum over ALL lead levels, with no contact weights.
    Multiply complete operators before transformation and QP projection.
    """
    if h.nimp != 2:
        raise ValueError("this example requires a single impurity orbital")
    modes = len(h.metadata["coordinates"]["electron_annihilators"])
    zd = (number(0)-number(1))/2
    pd = create(0)*annihilate(1)
    zb = sum((number(i)-number(i+1))/2 for i in range(2, modes, 2))
    pb = sum(create(i)*annihilate(i+1) for i in range(2, modes, 2))
    correlation = zd*zb+(pd*pb.dagger()+pd.dagger()*pb)/2
    z, p = zd+zb, pd+pb
    total_squared = z*z+(p*p.dagger()+p.dagger()*p)/2
    h.observables["dot_bath_spin"] = h.physical_operator(correlation)
    h.observables["total_spin_squared"] = h.physical_operator(total_squared)


def main():
    # D=5, rather than the large study's D=100: this is a finite-model control.
    bath = cosh_grid(pairs=1, delta=1., bandwidth=5.)
    h = reference_model(bath, u=2., gamma=.4, phi=3*np.pi/5, compress=False)
    add_spin_observables(h)
    sector = Sector(1, 1)
    print("Lowest odd state, S_z=+1/2; Delta=1, U=2, Gamma=0.4, D=5, phi=3*pi/5.")
    print("P_n counts dot electrons. C=<S_dot . (S_L+S_R)> is evaluated directly.")
    print("   Q           E          q_d          P_0          P_1          P_2            C")
    previous = np.inf
    rows = []
    for cutoff in (0, 1, 2, 3, 4, 5, 6, None):
        state = solve(h, cutoff=cutoff, sector=sector, options={"method": "dense"})
        row = values_from_state(state)
        direct = float(state.observables["dot_bath_spin"][0].real)
        np.testing.assert_allclose(state.observables["total_spin_squared"][0], .75, atol=2e-11, rtol=0)
        np.testing.assert_allclose(direct, row["C_d_bath"], atol=2e-11, rtol=0)
        np.testing.assert_allclose(row["P_0"], row["P_2"], atol=2e-11, rtol=0)
        assert row["energy"] <= previous+2e-11
        assert max(state.residuals) < 2e-10
        previous = row["energy"]
        label = "full" if cutoff is None else str(cutoff)
        displayed = dict(row, C_d_bath=direct)
        print(f"{label:>4} " + " ".join(f"{value:12.8f}" for value in displayed.values()))
        rows.append(row)
    error = {key: abs(rows[2][key]-rows[-1][key]) for key in rows[-1]}
    print("Absolute errors at Q=2 against full finite-bath ED:")
    print("  " + ", ".join(f"{key}={value:.3e}" for key, value in error.items()))
    even = solve(h, cutoff=None, sector=Sector(0, 0), options={"method": "dense"})
    print(f"Full-bath E_even-E_doublet = {even.energies[0]-rows[-1]['energy']:.10f}")
    print("Direct spin correlation agrees with 3*(q_d-P_1)/4 at every cutoff.")
    print("Bath convergence is a separate check; see docs/qp_convergence.md.")


if __name__ == "__main__":
    main()
