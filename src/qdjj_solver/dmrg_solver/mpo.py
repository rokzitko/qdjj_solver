"""Exact finite-fermion MPO compilation in declared one-particle coordinates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from tenpy.linalg import np_conserved as npc
from tenpy.models.lattice import TrivialLattice
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.networks.mps import MPSEnvironment
from tenpy.networks.mpo import MPOEnvironment
from tenpy.networks.site import FermionSite, Site, SpinHalfFermionSite

from ..common.io import fingerprint
from ..common.problem import DEFAULT_SECTOR, Hamiltonian, Sector


def _merge_parallel_rows(mpo, max_dense_entries=1_000_000):
    """Merge exactly proportional MPO continuation rows, right to left.

    A star with its impurity in the middle otherwise retains one redundant
    virtual channel per earlier bath mode. Only bitwise-reconstructible rows
    with identical normalized entries AND identical charges are merged. There
    is no singular-value/coefficient threshold. Identity channels are protected.
    Large dense intermediates are skipped, affecting performance only.
    """
    labels = ["wL", "wR", "p", "p*"]
    original_dimensions = mpo.chi
    for i in range(mpo.L-1, 0, -1):
        tensor, previous = mpo.get_W(i).transpose(labels), mpo.get_W(i-1).transpose(labels)
        if max(np.prod(tensor.shape), np.prod(previous.shape)) > max_dense_entries:
            continue
        dense = tensor.to_ndarray()
        old_leg = tensor.get_leg("wL")
        charges = old_leg.to_qflat()
        protected = {mpo.IdL[i], mpo.IdR[i]}-{None}
        rows, representatives, assignments, factors, seen = [], [], [], [], {}
        for j, row in enumerate(dense):
            nonzero = np.flatnonzero(row)
            if not len(nonzero) and j not in protected:
                assignments.append(None)
                factors.append(0.)
                continue
            factor, normalized = 1., row
            key = None
            if j not in protected:
                factor = row.flat[nonzero[0]]
                normalized = row/factor
                # Complex division z/z need not round to exactly one. Primitive
                # fermionic rows consist of exact signed copies of a strength;
                # recognize these copies without a floating-point tolerance.
                for unit in (1., -1., 1j, -1j):
                    if np.iscomplexobj(normalized) or np.isreal(unit):
                        normalized[row == unit*factor] = unit
                # Canonicalize signed zeros only for hashing. Require exact
                # reconstruction before replacing any original tensor entries.
                normalized.real[normalized.real == 0.] = 0.
                if np.iscomplexobj(normalized):
                    normalized.imag[normalized.imag == 0.] = 0.
                if np.array_equal(row, factor*normalized):
                    key = (tuple(charges[j]), normalized.tobytes())
                else:
                    factor, normalized = 1., row
            group = seen.get(key) if key is not None else None
            if group is None:
                group = len(rows)
                rows.append(normalized)
                representatives.append(j)
                if key is not None:
                    seen[key] = group
            assignments.append(group)
            factors.append(factor)
        if len(rows) >= len(dense):
            continue
        transfer = np.zeros((len(dense), len(rows)), dtype=dense.dtype)
        for j, group in enumerate(assignments):
            if group is not None:
                transfer[j, group] = factors[j]
        new_leg = npc.LegCharge.from_qflat(old_leg.chinfo, charges[representatives], qconj=old_leg.qconj)
        new_tensor = npc.Array.from_ndarray(np.asarray(rows),
                                             [new_leg]+tensor.legs[1:], labels=labels,
                                             qtotal=tensor.qtotal, cutoff=0.)
        previous_dense = np.tensordot(previous.to_ndarray(), transfer, axes=[1, 0]).transpose(0, 3, 1, 2)
        new_previous = npc.Array.from_ndarray(previous_dense,
                                               [previous.legs[0], new_leg.conj()]+previous.legs[2:],
                                               labels=labels, qtotal=previous.qtotal, cutoff=0.)
        mpo.set_W(i, new_tensor)
        mpo.set_W(i-1, new_previous)
        if mpo.IdL[i] is not None:
            mpo.IdL[i] = assignments[mpo.IdL[i]]
        if mpo.IdR[i] is not None:
            mpo.IdR[i] = assignments[mpo.IdR[i]]
    mpo.test_sanity()
    mpo.qdjj_original_dimensions = original_dimensions
    return mpo


def _site(modes, hamiltonian, sector, chinfo):
    """One or two fermionic modes, with charges derived from physical mode labels."""
    if len(modes) == 2:
        raw = SpinHalfFermionSite(cons_N=None, cons_Sz=None)
        labels = ["empty", "up", "down", "full"]
        occupations = [(0, 0), (1, 0), (0, 1), (1, 1)]
        names = ("JW", "Cu", "Cdu", "Cd", "Cdd", "Nu", "Nd", "Ntot")
        fermions = ("Cu", "Cdu", "Cd", "Cdd")
    else:
        raw = FermionSite(conserve=None)
        labels = ["empty", "full"]
        occupations = [(0,), (1,)]
        names = ("JW", "C", "Cd", "N")
        fermions = ("C", "Cd")
    charges = []
    for occ in occupations:
        charge = [sum(occ) % 2]
        if sector.twice_sz is not None:
            charge.append(sum(n*hamiltonian.spins[m] for n, m in zip(occ, modes, strict=True)))
        if sector.eta is not None:
            charge.append(sum(n*(hamiltonian.eta_labels[m] == -1)
                              for n, m in zip(occ, modes, strict=True)) % 2)
        if sector.particle_number is not None:
            charge.append(sum(occ))
        charges.append(charge)
    leg = npc.LegCharge.from_qflat(chinfo, charges)
    site = Site(leg, labels, sort_charge=True,
                **{name: raw.get_op(name).to_ndarray() for name in names})
    site.need_JW_string.update(fermions)
    site.charge_to_JW_parity = np.array([1]+[0]*(chinfo.qnumber-1))
    return site


@dataclass
class PreparedModel:
    hamiltonian: Hamiltonian
    sector: Sector
    order: tuple[int, ...]
    groups: tuple[tuple[int, ...], ...]
    sites: list
    lattice: TrivialLattice
    mode_operators: dict
    target_charge: tuple[int, ...]
    coordinate_id: str

    def term(self, key, *, lattice=False):
        result = []
        for op in key:
            mode = abs(op)-1
            if mode not in self.mode_operators:
                raise ValueError("operator acts on an undeclared mode")
            site, annihilation, creation = self.mode_operators[mode]
            result.append((creation if op > 0 else annihilation, [0, site] if lattice else site))
        return result

    def explicit_jw_term(self, key):
        """Tensor-product operators, including the left boundary string for odd terms."""
        factors = [[] for _ in self.sites]
        for name, site in self.term(key):
            for left in range(site):
                factors[left].append("JW")
            factors[site].append(name)
        result = []
        for i, names in enumerate(factors):
            if not names:
                continue
            product = " ".join(names)
            if not np.array_equal(self.sites[i].get_op(product).to_ndarray(), np.eye(self.sites[i].dim)):
                result.append((product, i))
        return result

    def compile(self, operator, *, merge_parallel=True):
        """Compile an even, charge-preserving operator without coefficient truncation."""
        if not operator.terms:
            return None
        coupling = CouplingModel(self.lattice)
        for key, coefficient in operator.terms.items():
            if key:
                coupling.add_local_term(coefficient, self.term(key, lattice=True))
            else:
                coupling.add_onsite_term(coefficient, 0, "Id")
        mpo = coupling.calc_H_MPO(tol_zero=0.)
        return _merge_parallel_rows(mpo) if merge_parallel else mpo

    def model(self, operator):
        return MPOModel(self.lattice, self.compile(operator))

    def matrix_element(self, operator, bra, ket):
        """Exact contraction of a general polynomial, including parity-changing terms."""
        env = MPSEnvironment(bra, ket)
        value = 0j
        for key, coefficient in operator.terms.items():
            term = self.explicit_jw_term(key)
            value += coefficient*(env.expectation_value_term(term, autoJW=False)
                                  if term else bra.overlap(ket))
        return value

    def conserves_charges(self, operator):
        for key in operator.terms:
            self.term(key)  # validate all mode indices
            if len(key) % 2:
                return False
            if self.sector.particle_number is not None and sum(1 if op > 0 else -1 for op in key):
                return False
            if self.sector.twice_sz is not None:
                if sum((1 if op > 0 else -1)*self.hamiltonian.spins[abs(op)-1] for op in key):
                    return False
            if self.sector.eta is not None:
                if np.prod([self.hamiltonian.eta_labels[abs(op)-1] for op in key], dtype=int) != 1:
                    return False
        return True

    def matrix_elements(self, operator, bras, kets):
        if self.conserves_charges(operator):
            mpo = self.compile(operator)
            if mpo is None:
                return np.zeros((len(bras), len(kets)), complex)
            return np.array([[MPOEnvironment(bra, mpo, ket).full_contraction(0)
                              for ket in kets] for bra in bras])
        return np.array([[self.matrix_element(operator, bra, ket) for ket in kets] for bra in bras])

    def onsite_costs(self):
        """Diagonal local energies used only to select distinct product-state seeds."""
        costs = [np.zeros(site.dim) for site in self.sites]
        for key, coefficient in self.hamiltonian.operator.terms.items():
            term = self.term(key)
            if term and len({i for _, i in term}) == 1:
                i = term[0][1]
                matrix = self.sites[i].get_op(" ".join(name for name, _ in term)).to_ndarray()
                costs[i] += (coefficient*np.diag(matrix)).real
        return costs


def prepare(hamiltonian, sector=DEFAULT_SECTOR, *, group_size=2, mode_order=None):
    """Choose an MPS layout without modifying physical coordinates or removing modes."""
    hamiltonian.check_sector(sector)
    if not isinstance(group_size, int) or isinstance(group_size, bool) or group_size not in (1, 2):
        raise ValueError("group_size must be one or two fermionic modes per site")
    nmodes = len(hamiltonian.spins)
    if nmodes == 0:
        raise ValueError("DMRG requires at least one fermionic mode")
    order = tuple(range(nmodes)) if mode_order is None else tuple(mode_order)
    if len(order) != nmodes or sorted(order) != list(range(nmodes)):
        raise ValueError("mode_order must be a permutation of all declared modes")
    groups = tuple(order[i:i+group_size] for i in range(0, nmodes, group_size))
    mods, names, target = [2], ["parity"], [sector.parity]
    if sector.twice_sz is not None:
        mods.append(1)
        names.append("twice_sz")
        target.append(sector.twice_sz)
    if sector.eta is not None:
        mods.append(2)
        names.append("eta_parity")
        target.append(int(sector.eta == -1))
    if sector.particle_number is not None:
        mods.append(1)
        names.append("particle_number")
        target.append(sector.particle_number)
    chinfo = npc.ChargeInfo(mods, names)
    sites = [_site(group, hamiltonian, sector, chinfo) for group in groups]
    lattice = TrivialLattice(sites, bc="open", bc_MPS="finite")
    operators = {}
    for i, group in enumerate(groups):
        if len(group) == 1:
            operators[group[0]] = (i, "C", "Cd")
        else:
            operators[group[0]] = (i, "Cu", "Cdu")
            operators[group[1]] = (i, "Cd", "Cdd")
    identity = fingerprint(dict(physical_coordinates=hamiltonian.coordinate_fingerprint(),
                                groups=groups, charges=names))
    return PreparedModel(hamiltonian, sector, order, groups, sites, lattice,
                         operators, tuple(target), identity)


def residual_norm(prepared, psi, energy):
    """Compute ||(H-E)psi|| with an untruncated MPO application and QR sweep.

    This avoids subtracting <H>² from <H²>. No MPS compression or singular-value
    cutoff is applied to the residual; the temporary bonds can be substantially
    larger than those of psi.
    """
    mpo = prepared.compile(prepared.hamiltonian.operator-energy)
    if mpo is None:
        return 0.
    if len(prepared.sites) == 1:
        matrix = mpo.get_W(0).take_slice([mpo.get_IdL(0), mpo.get_IdR(0)], ["wL", "wR"])
        vector = psi.get_B(0, form="B").take_slice([0, 0], ["vL", "vR"])
        return float(npc.tensordot(matrix, vector, axes=["p*", "p"]).norm())
    image = psi.copy()
    mpo.apply_naively(image)
    right = None
    for i in range(image.L):
        tensor = image.get_B(i, form=None)
        if right is not None:
            tensor = npc.tensordot(right, tensor, axes=["vR", "vL"])
        matrix = tensor.combine_legs([["vL", "p"]], qconj=[+1])
        _, right = npc.qr(matrix, inner_labels=["vR", "vL"], cutoff=0.)
    return float(right.norm())
