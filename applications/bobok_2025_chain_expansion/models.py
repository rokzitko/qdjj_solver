"""Published one-dot junction and coherent common-lead double-dot geometries."""

import numpy as np

from qdjj_solver import (Hamiltonian, Impurity, annihilate, create, hermitian_pair,
                         make_model, number, reference_model)


def junction(bath, u, gamma_per_lead, phi, *, compress=True):
    return reference_model(bath, u=u, gamma=2*gamma_per_lead, phi=phi,
                            detuning=0., symmetry=compress, compress=compress)


def gal_junction(u, gamma_per_lead, phi):
    """Dot-only GAL at half filling, without the additional GAL+C band term."""
    gamma = 2*gamma_per_lead
    nu = 1/(1+gamma)
    pair = create(0)*create(1)
    op = Impurity.anderson(u*nu**2).operator-hermitian_pair(nu*gamma*np.cos(phi/2)*pair)
    derivative = hermitian_pair(nu*gamma/2*np.sin(phi/2)*pair)
    return Hamiltonian(op, 2, (1, -1), observables={"phase_derivative": derivative},
                        metadata=dict(kind="GAL", u=u, gamma=gamma, phi=phi, delta=1.))


def spin_operators(orbital):
    up, down = 2*orbital, 2*orbital+1
    return (number(up)-number(down))/2, create(up)*annihilate(down)


def spin_square(orbitals):
    z = sum(spin_operators(j)[0] for j in range(orbitals))
    raising = sum(spin_operators(j)[1] for j in range(orbitals))
    return z*z+(raising*raising.dagger()+raising.dagger()*raising)/2


def double_dot(bath, u, gamma_per_dot, *, exchange_basis=False):
    """zeta=1: one coherent bath, Gamma_ij=gamma_per_dot for every i,j."""
    charges = [number(2*j)+number(2*j+1) for j in (0, 1)]
    op = sum(u/2*(n-1)*(n-1) for n in charges)
    z0, plus0 = spin_operators(0)
    z1, plus1 = spin_operators(1)
    correlation = z0*z1+hermitian_pair(plus0*plus1.dagger())/2
    pairing = sum(hermitian_pair(create(2*j)*create(2*j+1))/2 for j in (0, 1))/2
    observables = {"spin_correlation": correlation, "pairing": pairing}
    contact = np.sqrt(gamma_per_dot/(np.pi*bath.rho))*np.tile(np.eye(2), (1, 2))
    if exchange_basis:
        # Exact dot-only rotation. Dark impurity occupations remain dynamical:
        # local interactions generate exchange and pair hopping between them.
        maps = [(annihilate(s)+annihilate(s+2))/np.sqrt(2) for s in (0, 1)]
        maps += [(annihilate(s)-annihilate(s+2))/np.sqrt(2) for s in (0, 1)]
        observables.update({f"double_occupancy_{j}": number(2*j)*number(2*j+1) for j in (0, 1)})
        op = op.substitute(maps)
        observables = {name: value.substitute(maps) for name, value in observables.items()}
        contact = np.c_[np.sqrt(2)*contact[:, :2], np.zeros((2, 2))]
    impurity = Impurity(2, op, observables=observables)
    h = make_model(impurity, [bath], [contact], phases=[0.])
    if exchange_basis:
        h.eta_labels = (1, 1, -1, -1)+(1,)*h.bath_modes
        # Physical electron operators still refer to the ORIGINAL two dots.
        h.metadata["coordinates"]["electron_annihilators"][:4] = [m.to_records() for m in maps]
        h.metadata["coordinates"]["impurity_basis"] = "symmetric-antisymmetric"
    h.observables["total_spin_squared"] = h.physical_operator(spin_square(2+bath.levels)).cleaned(1e-13)
    h.metadata.update(kind="common-lead-DQD", u=u, gamma_per_dot=gamma_per_dot, zeta=1.,
                      hopping=0., capacitance=0., detuning=0., exchange_basis=exchange_basis)
    return h
