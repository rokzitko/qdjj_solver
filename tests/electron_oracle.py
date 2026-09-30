"""Independent physical-electron tensor products for small finite reservoirs."""

import numpy as np
from scipy.sparse import csr_matrix, eye, kron


def electron_operators(modes):
    if modes > 12:
        raise ValueError("the independent ED oracle is limited to twelve modes")
    identity = eye(2, format="csr")
    z = csr_matrix(np.diag([1., -1.]))
    lower = csr_matrix([[0., 1.], [0., 0.]])
    result = []
    for i in range(modes):
        matrix = csr_matrix([[1.]])
        for j in range(modes-1, -1, -1):
            matrix = kron(matrix, lower if j == i else (z if j < i else identity), format="csr")
        result.append(matrix)
    return result


def evaluate(operator, annihilators):
    result = csr_matrix(annihilators[0].shape, dtype=complex)
    for key, value in operator.terms.items():
        term = value*eye(result.shape[0], format="csr", dtype=complex)
        for op in key:
            term = term @ (annihilators[op-1].getH() if op > 0 else annihilators[-op-1])
        result += term
    return result


def electron_hamiltonian(impurity, baths, tunneling, phases, direct=None):
    nimp = 2*impurity.orbitals
    modes = nimp+2*sum(b.levels for b in baths)
    cs = electron_operators(modes)
    matrix = evaluate(impurity.operator, cs)
    contacts, start = [], nimp
    for bath in baths:
        contact = [csr_matrix(matrix.shape, dtype=complex) for _ in range(2)]
        for i, (xi, weight) in enumerate(zip(bath.xi, bath.weights, strict=True)):
            up, down = cs[start+2*i:start+2*i+2]
            matrix += xi*(up.getH() @ up+down.getH() @ down)
            pair = -bath.delta*up.getH() @ down.getH()
            matrix += pair+pair.getH()
            contact[0] += np.sqrt(weight)*up
            contact[1] += np.sqrt(weight)*down
        contacts.append(contact)
        start += 2*bath.levels
    for lead, hopping in enumerate(tunneling):
        for spin in range(2):
            for a in range(nimp):
                term = np.exp(.5j*phases[lead])*hopping[spin, a]*contacts[lead][spin].getH() @ cs[a]
                matrix += term+term.getH()
    for (left, right), hopping in (direct or {}).items():
        for s in range(2):
            for t in range(2):
                term = np.exp(.5j*(phases[left]-phases[right]))*hopping[s, t]*contacts[left][s].getH() @ contacts[right][t]
                matrix += term+term.getH()
    reference = sum(np.sum(b.xi-b.energies) for b in baths)
    return matrix-reference*eye(matrix.shape[0], format="csr")


def sector_indices(modes, parity, twice_sz=None):
    result = []
    for state in range(2**modes):
        occupied = [i for i in range(modes) if state & (1 << i)]
        if len(occupied) % 2 == parity and (twice_sz is None or sum(1-2*(i % 2) for i in occupied) == twice_sz):
            result.append(state)
    return np.array(result)
