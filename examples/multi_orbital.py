"""Compare both backends for a small complex two-orbital impurity."""
import numpy as np

from qdjj_solver import (DiscreteBath, Impurity, Sector, annihilate, create,
                         hermitian_pair, make_model, number, solve)


def main():
    op = sum((-.3+.1*i)*number(i) for i in range(4))
    op += 1.2*number(0)*number(1)+.9*number(2)*number(3)
    op += hermitian_pair(.08j*create(0)*annihilate(3))
    op += hermitian_pair(.06*create(0)*create(1)*annihilate(3)*annihilate(2))
    bath = DiscreteBath([.4], [1.], delta=.8, bandwidth=3.)
    hopping = np.array([[.3, .07j, .1, 0.], [0., .27, .04j, .12]])
    h = make_model(Impurity(2, op), [bath], [hopping], phases=[.3])
    sector = Sector(1, None)
    qp = solve(h, sector=sector, options={"eigenpairs": 3})
    mps = solve(h, sector=sector, backend="dmrg",
                options={"eigenpairs": 3, "chi_max": 32, "require_convergence": True})
    print("Full finite-Fock-space energies:", qp.energies)
    print("MPS energies:", mps.energies)
    print("Maximum difference:", np.max(abs(qp.energies-mps.energies)))
    print("Physical spin matrix:", mps.matrix_elements(h.observables["impurity_spin_z"]))


if __name__ == "__main__":
    main()
