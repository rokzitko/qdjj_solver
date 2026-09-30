"""Interacting impurities and BCS reservoirs, with explicit physical conventions."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .algebra import FermionOperator, annihilate, create, hermitian_pair, number
from .baths import DiscreteBath
from .problem import Hamiltonian


@dataclass
class Impurity:
    """Untruncated spinful impurity; modes are (orbital 0 up, down, 1 up, ...)."""

    orbitals: int
    operator: FermionOperator
    observables: dict[str, FermionOperator] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.orbitals, int) or self.orbitals < 1:
            raise ValueError("impurity orbitals must be a positive integer")
        self.observables = dict(self.observables)
        for name, observable in self.observables.items():
            if any(abs(op) > 2*self.orbitals for key in observable.terms for op in key):
                raise ValueError(f"impurity observable {name!r} acts on an undeclared fermionic mode")
        self.observables.setdefault("impurity_charge", sum(number(i) for i in range(2*self.orbitals)))
        self.observables.setdefault("impurity_spin_z", sum((1-2*(i%2))*number(i)/2 for i in range(2*self.orbitals)))
        for a in range(self.orbitals):
            self.observables.setdefault(f"double_occupancy_{a}", number(2*a)*number(2*a+1))
        Hamiltonian(self.operator, 2*self.orbitals, (1, -1)*self.orbitals)

    @classmethod
    def anderson(cls, u=2., detuning=0., field=0.):
        """H_imp=U(n-1)^2/2+detuning(n-1)+field S_z (single orbital)."""
        if not np.all(np.isfinite([u, detuning, field])):
            raise ValueError("impurity parameters must be finite")
        n = number(0) + number(1)
        op = u/2*(n-1)*(n-1) + detuning*(n-1) + field/2*(number(0)-number(1))
        return cls(1, op, metadata=dict(kind="anderson", u=u, detuning=detuning, field=field))

    @classmethod
    def from_integrals(cls, one_body, interaction=None, constant=0.):
        """H=sum h_ab d†_a d_b + (1/4)sum U_abcd d†_a d†_b d_d d_c.

        All indices are spin-orbital indices. The factor 1/4 is appropriate to
        antisymmetrized two-electron integrals. A general FermionOperator may
        instead be supplied through the main constructor.
        """
        h = np.asarray(one_body, dtype=complex)
        if h.ndim != 2 or h.shape[0] != h.shape[1] or h.shape[0] % 2 or not np.all(np.isfinite(h)):
            raise ValueError("one_body must be a finite even-dimensional square matrix")
        op = FermionOperator.scalar(constant)
        for a, b in zip(*np.nonzero(h), strict=True):
            op += h[a, b]*create(int(a))*annihilate(int(b))
        if interaction is not None:
            interaction = np.asarray(interaction, dtype=complex)
            if interaction.shape != (h.shape[0],)*4 or not np.all(np.isfinite(interaction)):
                raise ValueError("two-electron integral tensor has wrong shape or non-finite entries")
            for a, b, c, d in zip(*np.nonzero(interaction), strict=True):
                op += interaction[a, b, c, d]/4*create(int(a))*create(int(b))*annihilate(int(d))*annihilate(int(c))
        return cls(h.shape[0]//2, op)


def _linear_map(u, v, start):
    result = []
    for row in range(u.shape[0]):
        terms = {}
        for j in np.flatnonzero(abs(u[row]) > 1e-15):
            terms[-(start+int(j)+1)] = u[row, j]
        for j in np.flatnonzero(abs(v[row]) > 1e-15):
            terms[start+int(j)+1] = v[row, j]
        result.append(FermionOperator({(key,): value for key, value in terms.items()}))
    return result


def _coupled_bath_map(baths, phases, direct):
    """Diagonalize the spin-conserving quadratic environment in Nambu space."""
    sizes = [b.levels for b in baths]
    offsets = np.r_[0, np.cumsum(sizes)]
    n = int(offsets[-1])
    xi = np.concatenate([b.xi for b in baths])
    gaps = np.concatenate([np.full(b.levels, b.delta) for b in baths])
    hu, hd = np.diag(xi).astype(complex), np.diag(xi).astype(complex)
    for (l, m), matrix in direct.items():
        matrix = np.asarray(matrix, complex)
        if matrix[0, 1] != 0 or matrix[1, 0] != 0:
            raise ValueError("coupled-bath vacuum currently requires spin-conserving interlead hopping")
        left = slice(offsets[l], offsets[l+1])
        right = slice(offsets[m], offsets[m+1])
        contact = np.sqrt(baths[l].weights[:, None]*baths[m].weights[None, :])
        for spin, h in enumerate((hu, hd)):
            block = np.exp(.5j*(phases[l]-phases[m]))*matrix[spin, spin]*contact
            h[left, right] += block
            h[right, left] += block.conj().T
    nambu = np.block([[hu, -np.diag(gaps)], [-np.diag(gaps), -hd.T]])
    eigenvalues, eigenvectors = np.linalg.eigh(nambu)
    positive, negative = np.flatnonzero(eigenvalues > 0), np.flatnonzero(eigenvalues < 0)[::-1]
    if len(positive) != n or len(negative) != n or np.min(abs(eigenvalues)) < 1e-12*min(gaps):
        raise ValueError("coupled-bath vacuum must be unpolarized and have a nonzero excitation gap")
    u, v = np.zeros((2*n, 2*n), complex), np.zeros((2*n, 2*n), complex)
    u[0::2, 0::2], v[0::2, 1::2] = eigenvectors[:n, positive], eigenvectors[:n, negative]
    v[1::2, 0::2], u[1::2, 1::2] = eigenvectors[n:, positive].conj(), eigenvectors[n:, negative].conj()
    energies = np.column_stack((eigenvalues[positive], -eigenvalues[negative])).ravel()
    vacuum = float(np.trace(hd).real + eigenvalues[negative].sum())
    return u, v, energies.tolist(), vacuum


def make_model(impurity: Impurity, baths, tunneling, phases=None, direct=None,
               phase_velocities=None, eta_basis=False, compress_pairs=False,
               bath_reference="isolated", coefficient_tolerance=1e-14,
               direct_derivatives=None) -> Hamiltonian:
    """Construct a multi-orbital impurity coupled to one or more BCS reservoirs.

    ``tunneling[l]`` is a (2, 2*n_orbitals) matrix T_l in C_l† T_l d + h.c.,
    with C_l=sum_i sqrt(w_li)c_li. ``direct[(l,m)]`` multiplies C_l† W_lm C_m
    for l<m. Pair potentials are real; phase factors exp(i phi_l/2) are inserted
    in T_l and exp(i(phi_l-phi_m)/2) in W_lm.

    The optional eta/paired-mode reduction is restricted to the symmetric,
    particle-hole-symmetric single-dot junction. It preserves bath QP number.
    Decoupled modes are kept in their vacuum when compress_pairs=True; this
    option targets bound impurity states, not arbitrary continuum excitations.

    ``direct_derivatives[name][(l,m)]`` specifies a physical contact-hopping
    observable with the same convention as ``direct``. It is transformed even
    when the corresponding Hamiltonian parameter is zero. Discarded-mode
    vacuum contractions are retained, including in a compressed bath.
    """
    baths = tuple(baths)
    if not baths or not all(isinstance(b, DiscreteBath) for b in baths):
        raise ValueError("provide at least one DiscreteBath")
    nimp, leads = 2*impurity.orbitals, len(baths)
    phases = np.zeros(leads) if phases is None else np.asarray(phases, float)
    velocity = np.zeros(leads) if phase_velocities is None else np.asarray(phase_velocities, float)
    if phases.shape != (leads,) or velocity.shape != (leads,) or not np.all(np.isfinite(np.r_[phases, velocity])):
        raise ValueError("one finite phase and phase derivative per reservoir are required")
    tunneling = tuple(np.asarray(t, complex) for t in tunneling)
    if len(tunneling) != leads or any(t.shape != (2, nimp) or not np.all(np.isfinite(t)) for t in tunneling):
        raise ValueError("each tunneling matrix must have shape (2, 2*impurity.orbitals)")
    direct = {} if direct is None else dict(direct)
    for (l, m), w in direct.items():
        if not 0 <= l < m < leads or np.asarray(w).shape != (2, 2) or not np.all(np.isfinite(w)):
            raise ValueError("direct hopping requires l<m and finite 2-by-2 spin matrices")
    direct_derivatives = {} if direct_derivatives is None else dict(direct_derivatives)
    for name, matrices in direct_derivatives.items():
        if not isinstance(name, str) or name in impurity.observables or name == "phase_derivative":
            raise ValueError("contact observable names must be distinct from built-in observables")
        for (l, m), w in matrices.items():
            if not 0 <= l < m < leads or np.asarray(w).shape != (2, 2) or not np.all(np.isfinite(w)):
                raise ValueError("contact observables require l<m and finite 2-by-2 spin matrices")
    if compress_pairs and not eta_basis:
        raise ValueError("paired-mode compression requires the verified eta representation")
    if bath_reference not in ("isolated", "coupled"):
        raise ValueError("bath_reference must be isolated or coupled")
    if not np.isfinite(coefficient_tolerance) or coefficient_tolerance < 0:
        raise ValueError("coefficient_tolerance must be nonnegative and finite")
    if eta_basis and bath_reference == "coupled":
        raise ValueError("eta charge representation is implemented for isolated-reservoir QPs")
    if eta_basis:
        valid = (leads == 2 and impurity.orbitals == 1 and
                 impurity.metadata.get("kind") == "anderson" and impurity.metadata.get("detuning") == 0
                 and baths[0].paired and baths[1].paired and baths[0].delta == baths[1].delta
                 and np.array_equal(baths[0].xi, baths[1].xi)
                 and np.array_equal(baths[0].weights, baths[1].weights)
                 and np.allclose(tunneling[0], tunneling[1], rtol=0, atol=1e-14)
                 and np.allclose(tunneling[0], np.eye(2)*tunneling[0][0, 0], rtol=0, atol=1e-14)
                 and abs(tunneling[0][0, 0].imag) < 1e-14)
        w = np.asarray(direct.get((0, 1), np.zeros((2, 2))))
        valid = valid and np.allclose(w, np.diag([1j*w[0, 0].imag, -1j*w[0, 0].imag]), rtol=0, atol=1e-14)
        if not valid:
            raise ValueError("eta basis requires matched symmetric contacts, symmetric gate/bands, and pure W_S")
        if compress_pairs and np.max(abs(w)) > 0:
            raise ValueError("finite W_S couples the otherwise decoupled paired modes")

    sizes = [2*b.levels for b in baths]
    offsets = np.r_[0, np.cumsum(sizes)]
    total = int(offsets[-1])
    u, v = np.zeros((total, total), complex), np.zeros((total, total), complex)
    energies = []
    for l, bath in enumerate(baths):
        for i, (xi, e) in enumerate(zip(bath.xi, bath.energies, strict=True)):
            a = int(offsets[l]) + 2*i
            cu, cv = np.sqrt((1+xi/e)/2), np.sqrt((1-xi/e)/2)
            u[a, a] = u[a+1, a+1] = cu
            v[a, a+1], v[a+1, a] = cv, -cv
            energies.extend([e, e])
    spins = [1, -1]*(total//2)
    isolated_vacuum = float(sum(np.sum(b.xi-b.energies) for b in baths))
    vacuum_shift = 0.
    if bath_reference == "coupled":
        u, v, energies, coupled_vacuum = _coupled_bath_map(baths, phases, direct)
        vacuum_shift = coupled_vacuum-isolated_vacuum
    def contact_covariance(hole_map):
        rows = np.asarray([sum(np.sqrt(weight)*hole_map[int(offsets[l])+2*i+s]
                               for i, weight in enumerate(bath.weights))
                           for l, bath in enumerate(baths) for s in range(2)])
        return rows.conj() @ rows.T
    full_contact_covariance = contact_covariance(v) if direct_derivatives else None
    eta = None
    impurity_maps = [annihilate(i) for i in range(nimp)]
    decoupled = []
    if eta_basis:
        bath, levels = baths[0], baths[0].levels
        rotation = np.zeros((total, total), complex)
        for i in range(levels):
            for s in range(2):
                l, r = 2*i+s, 2*levels+2*(levels-1-i)+s
                plus, minus = 2*i+s, 2*levels+2*i+s
                rotation[l, plus] = rotation[l, minus] = 1/np.sqrt(2)
                rotation[r, plus], rotation[r, minus] = 1/np.sqrt(2), -1/np.sqrt(2)
        u, v = u @ rotation, v @ rotation.conj()
        energies = np.tile(np.repeat(bath.energies, 2), 2).tolist()
        eta = [-(1-2*s)*r for r in (1, -1) for i in range(levels) for s in range(2)]
        common_phase = np.exp(-.25j*(phases[0]+phases[1]))
        impurity_maps = [common_phase*(annihilate(0)+create(1))/np.sqrt(2),
                         common_phase*(annihilate(1)-create(0))/np.sqrt(2)]
        if compress_pairs:
            phi = phases[1]-phases[0]
            cu = np.sqrt((1+bath.xi/bath.energies)/2)
            cv = np.sqrt((1-bath.xi/bath.energies)/2)
            columns, kept_energies, kept_spins, kept_eta = [], [], [], []
            for channel, r in enumerate((1, -1)):
                f = np.exp(-.25j*phi)*cu + r*np.exp(.25j*phi)*cv
                # r=+ uses s, r=- uses a; same coefficients for both spins.
                for i in range((levels+1)//2):
                    indices = sorted(set((i, levels-1-i)))
                    for s in range(2):
                        old = [channel*2*levels+2*j+s for j in indices]
                        coefficients = f[indices]*np.sqrt(bath.weights[indices])
                        norm = np.linalg.norm(coefficients)
                        if norm < 1e-14:
                            decoupled.extend([float(bath.energies[i])]*len(indices))
                            continue
                        column = np.zeros(total, complex)
                        column[old] = coefficients/norm
                        columns.append(column)
                        kept_energies.append(float(bath.energies[i]))
                        kept_spins.append(1-2*s)
                        kept_eta.append(-(1-2*s)*r)
                        decoupled.extend([float(bath.energies[i])]*(len(indices)-1))
            rotation = np.column_stack(columns)
            u, v = u @ rotation, v @ rotation.conj()
            energies, spins, eta = kept_energies, kept_spins, kept_eta
        eta = [1, -1] + eta

    electron_maps = _linear_map(u, v, nimp)
    contacts = []
    for l, bath in enumerate(baths):
        contacts.append([sum(np.sqrt(weight)*electron_maps[int(offsets[l])+2*i+s]
                             for i, weight in enumerate(bath.weights)) for s in range(2)])
    operator = impurity.operator.substitute(impurity_maps) + vacuum_shift
    for i, e in enumerate(energies):
        operator += e*number(nimp+i)
    derivative = FermionOperator()
    for l, hopping in enumerate(tunneling):
        for s, a in zip(*np.nonzero(hopping), strict=True):
            term = np.exp(.5j*phases[l])*hopping[s, a]*contacts[l][s].dagger()*impurity_maps[a]
            operator += hermitian_pair(term)
            derivative += hermitian_pair(.5j*velocity[l]*term)
    for (l, m), hopping in direct.items():
        hopping = np.asarray(hopping, complex)
        for s, t in zip(*np.nonzero(hopping), strict=True):
            term = np.exp(.5j*(phases[l]-phases[m]))*hopping[s, t]*contacts[l][s].dagger()*contacts[m][t]
            if bath_reference == "isolated":
                operator += hermitian_pair(term)
            derivative += hermitian_pair(.5j*(velocity[l]-velocity[m])*term)
    # Observable coefficients need not have the Hamiltonian's units or scale.
    observables = {key: value.substitute(impurity_maps) for key, value in impurity.observables.items()}
    if np.any(velocity):
        observables["phase_derivative"] = derivative
    if direct_derivatives:
        # Projecting each electron factor alone loses the dark-vacuum part of
        # <C_l^dagger C_m>. Add that contraction after normal ordering the
        # retained-mode product; all other dark terms have zero expectation.
        missing_covariance = full_contact_covariance-contact_covariance(v)
        for name, matrices in direct_derivatives.items():
            observable = FermionOperator()
            for (l, m), hopping in matrices.items():
                hopping = np.asarray(hopping, complex)
                for s, t in zip(*np.nonzero(hopping), strict=True):
                    product = contacts[l][s].dagger()*contacts[m][t]
                    product += missing_covariance[2*l+s, 2*m+t]
                    observable += hermitian_pair(np.exp(.5j*(phases[l]-phases[m]))*hopping[s, t]*product)
            observables[name] = observable
    metadata = dict(geometry="single-reservoir" if leads == 1 else f"{leads}-reservoir",
                    impurity=impurity.metadata, bath=[b.record() for b in baths], phases=phases.tolist(),
                    representation="eta-charge" if eta_basis else "physical-impurity",
                    paired_mode_compression=bool(compress_pairs), decoupled_qp_energies=decoupled,
                    energy_reference="isolated BCS reservoirs subtracted",
                    bcs_vacuum_energy=isolated_vacuum,
                    qp_vacuum="isolated BCS reservoirs" if bath_reference == "isolated" else "coupled quadratic reservoirs",
                    bath_reference=bath_reference, bath_vacuum_shift=vacuum_shift,
                    tunneling=[[[[z.real, z.imag] for z in row] for row in t] for t in tunneling],
                    direct=[dict(leads=[l, m], matrix=[[[complex(z).real, complex(z).imag] for z in row]
                                                     for row in np.asarray(w)]) for (l, m), w in direct.items()])
    cleaned = operator.cleaned(coefficient_tolerance)
    metadata.update(coefficient_tolerance=coefficient_tolerance,
                    coefficient_cleanup_bound=float(sum(abs(v) for v in (operator-cleaned).terms.values())),
                    coordinates=dict(kind="electron-to-quasiparticle", version=1,
                                     restriction="dark-vacuum" if compress_pairs else "none",
                                     electron_annihilators=[m.to_records() for m in impurity_maps+electron_maps]))
    return Hamiltonian(cleaned, nimp, (1, -1)*impurity.orbitals + tuple(spins),
                       None if eta is None else tuple(eta), observables, metadata)


def reference_model(bath, u=2., gamma=.4, phi=3*np.pi/5, rho_ws=0., rho_wn=0.,
                    geometry="junction", detuning=0., field=0., symmetry=True, compress=True,
                    bath_reference="isolated", coefficient_tolerance=1e-14):
    """Reference single-dot model; Gamma is the *total* hybridization.

    The dimensional interlead amplitude is W_S=(rho W_S)/rho. It acquires
    sqrt(w_i w_j), including fitted surrogate residues, upon discretization.
    """
    if geometry not in ("single", "junction"):
        raise ValueError("geometry must be single or junction")
    if not np.all(np.isfinite([u, gamma, phi, rho_ws, rho_wn, detuning, field])) or gamma < 0:
        raise ValueError("model parameters must be finite and Gamma nonnegative")
    if geometry == "single" and (rho_ws != 0 or rho_wn != 0):
        raise ValueError("interlead hopping requires two reservoirs")
    leads = 1 if geometry == "single" else 2
    t = np.sqrt(gamma/(np.pi*bath.rho*leads))*np.eye(2)
    eta = bool(symmetry and leads == 2 and bath.paired and detuning == 0 and rho_wn == 0 and bath_reference == "isolated")
    direct = {(0, 1): np.diag([(rho_wn+1j*rho_ws)/bath.rho, (rho_wn-1j*rho_ws)/bath.rho])} if leads == 2 else {}
    response = {"rho_ws_derivative": {(0, 1): np.diag([1j, -1j])/bath.rho}} if leads == 2 else {}
    result = make_model(Impurity.anderson(u, detuning, field), [bath]*leads, [t]*leads,
                        phases=[-phi/2, phi/2] if leads == 2 else [0.], direct=direct,
                        phase_velocities=[-.5, .5] if leads == 2 else None,
                         eta_basis=eta, compress_pairs=bool(eta and compress and rho_ws == 0),
                         bath_reference=bath_reference, coefficient_tolerance=coefficient_tolerance,
                         direct_derivatives=response)
    result.metadata.update(u=u, gamma=gamma, phi=phi if leads == 2 else 0., rho_ws=rho_ws,
                           rho_wn=rho_wn, detuning=detuning, field=field)
    return result
