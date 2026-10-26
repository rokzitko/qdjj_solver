"""Physical impurity-orbital reduced density matrices and local pair weights."""

import numpy as np

from .algebra import FermionOperator, annihilate, create, number


def orbital_density_operators(orbital=0, *, hamiltonian=None, physical=True):
    """Operators O[a,b] whose expectations equal rho[a,b] in a parity sector.

    Return a dictionary keyed by (row, column). The local basis is
    ``[|0>, |up>, |down>, |up down>]``, with ``|up down>=d_up^dagger d_down^dagger|0>``.
    Consequently ``O[a,b]=|b><a|``, not ``|a><b|``. Entries connecting different local
    parities vanish in fixed-total-parity states and are omitted.

    Without a Hamiltonian the operators use physical electron mode indices
    (2*orbital, 2*orbital+1). With one, transform to its canonical coordinates.
    physical=False explicitly interprets those indices as canonical modes.
    Verified eta compression is safe for impurity-only operators: its impurity
    map has no discarded bath factors. Arbitrary compressed-bath RDMs are not
    covered by this function.
    """
    if isinstance(orbital, (bool, np.bool_)) or not isinstance(orbital, (int, np.integer)) or orbital < 0:
        raise ValueError("orbital must be a nonnegative integer")
    up, down = 2*int(orbital), 2*int(orbital)+1
    if hamiltonian is not None and down >= hamiltonian.nimp:
        raise ValueError("orbital is outside the declared impurity modes")
    nu, nd = number(up), number(down)
    pair = annihilate(down)*annihilate(up)
    operators = {(0, 0): (1-nu)*(1-nd), (1, 1): nu*(1-nd),
                 (2, 2): (1-nu)*nd, (3, 3): nu*nd,
                 (0, 3): pair.dagger(), (3, 0): pair,
                 (1, 2): create(down)*annihilate(up), (2, 1): create(up)*annihilate(down)}
    if hamiltonian is None or not physical:
        return operators
    coordinates = hamiltonian.metadata.get("coordinates", {})
    if coordinates.get("restriction") != "dark-vacuum":
        return {key: hamiltonian.physical_operator(op) for key, op in operators.items()}
    records = coordinates.get("electron_annihilators", [])
    maps = [FermionOperator.from_records(record) for record in records[:hamiltonian.nimp]]
    if (not hamiltonian.metadata.get("paired_mode_compression")
            or hamiltonian.metadata.get("representation") != "eta-charge"
            or len(maps) != hamiltonian.nimp
            or any(abs(i) > hamiltonian.nimp for op in maps for word in op.terms for i in word)):
        raise ValueError("compressed physical orbital density requires a verified impurity-only eta map")
    return {key: op.substitute(maps) for key, op in operators.items()}


def orbital_density_properties(density):
    """Occupations and coherent pair probabilities for one normalized 4x4 RDM.

    ``pair_amplitude=<d_down d_up>=rho[3,0]``. plus/minus refer to the fixed-gauge
    normalized states ``(|0> +/- |up down>)/sqrt(2)``. phase_optimized_plus maximizes
    over the relative pair phase; it is different when pair_amplitude is complex.
    The input is validated, not renormalized or projected onto positive matrices.
    """
    rho = np.asarray(density, dtype=complex)
    if rho.shape != (4, 4) or not np.all(np.isfinite(rho)):
        raise ValueError("orbital density must be a finite 4-by-4 matrix")
    if not np.allclose(rho, rho.conj().T, atol=1e-9, rtol=0):
        raise ValueError("orbital density must be Hermitian")
    if not np.isclose(np.trace(rho), 1, atol=1e-9, rtol=0):
        raise ValueError("orbital density must have trace one")
    if np.linalg.eigvalsh(rho).min() < -1e-9:
        raise ValueError("orbital density must be positive semidefinite")
    p0, pu, pd, p2 = rho.diagonal().real
    pair = complex(rho[3, 0])
    even = float(p0+p2)
    return dict(occupation=float(pu+pd+2*p2), empty=float(p0), up=float(pu), down=float(pd),
                single=float(pu+pd), double=float(p2), even_weight=even, pair_amplitude=pair,
                plus=even/2+pair.real, minus=even/2-pair.real,
                phase_optimized_plus=even/2+abs(pair), phase_optimized_minus=even/2-abs(pair),
                purity=float(np.trace(rho @ rho).real))
