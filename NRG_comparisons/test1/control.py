"""Independent, finite-bath controls; no native QP, NRG, or Mathematica calls.

Physical inputs use one common energy unit and total single-reservoir Gamma.
Bath records use the standard DiscreteBath schema, with weights measuring
rho dxi, rho=1/(2*bandwidth). No input bath weights are renormalized.
"""

import math
from numbers import Real

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse.linalg import eigsh


_DENSE_THRESHOLD = 128


def _validate_physical(physical):
    for key in ("gap", "bandwidth", "u", "gamma", "detuning", "field"):
        value = physical[key]
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f"physical.{key} must be a finite real number")
    if physical["gap"] <= 0 or physical["bandwidth"] <= 0:
        raise ValueError("gap and bandwidth must be positive")
    if physical["u"] < 0 or physical["gamma"] < 0:
        raise ValueError("u and gamma must be nonnegative")
    if physical["geometry"] != "single":
        raise ValueError("the electron control requires geometry='single'")


def wilson_bath(physical, Lambda, nmax):
    """Return the normal star equivalent of the finite Z, z=1 Wilson chain.

    The bath contains f_0,...,f_nmax, hence nmax+1 signed spinful levels and
    nmax bonds. In physical energy units its zero-onsite chain has

        t_n = D*(1-Lambda**-1)/log(Lambda) * Lambda**(-n/2)
              * (1-Lambda**(-n-1))
              / sqrt((1-Lambda**(-2*n-1))*(1-Lambda**(-2*n-3))).

    This is the flat-band Z discretization at z=1, not the conventional
    Yoshida prefactor (1+1/Lambda)/2. The contact hopping is sqrt(2*D*Gamma/pi),
    without an A_Lambda correction. Eigenvector first-row squares are the
    normalized contact measure, NOT weights multiplied by 2*D.

    metadata.chain_hoppings and metadata.impurity_hopping retain physical
    energy units. The bound 0 <= nmax <= 100 is separate from the ED limit;
    QP callers can use longer chains. No superconducting modes are discarded.
    """
    _validate_physical(physical)
    if isinstance(Lambda, (bool, np.bool_)) or not isinstance(Lambda, Real) or not math.isfinite(Lambda):
        raise ValueError("Lambda must be a finite real number")
    if Lambda < 1.8:
        raise ValueError("Lambda < 1.8 is not authorized")
    if isinstance(nmax, (bool, np.bool_)) or not isinstance(nmax, (int, np.integer)) or not 0 <= nmax <= 100:
        raise ValueError("nmax must be an integer between 0 and 100")
    nmax, Lambda = int(nmax), float(Lambda)
    n = np.arange(nmax, dtype=float)
    q = 1 / Lambda
    hoppings = (physical["bandwidth"] * (1-q) / math.log(Lambda) * q**(n/2)
                * (1-q**(n+1)) / np.sqrt((1-q**(2*n+1)) * (1-q**(2*n+3))))
    onsite = np.zeros(nmax+1)
    if nmax == 0:
        # SciPy 1.12's STEV wrapper rejects an empty off-diagonal array.
        xi, vectors = onsite, np.ones((1, 1))
    else:
        # QR retains tiny contact components better than MRRR on long chains.
        xi, vectors = eigh_tridiagonal(onsite, hoppings, lapack_driver="stev")
    weights = vectors[0]**2
    if not np.all(weights > 0) or not np.isclose(weights.sum(), 1., rtol=0, atol=1e-13):
        raise ValueError("Wilson contact measure lost positivity or normalization")
    return {
        "xi": xi.tolist(), "weights": weights.tolist(),
        "delta": float(physical["gap"]), "bandwidth": float(physical["bandwidth"]),
        "metadata": {
            "kind": "wilson-chain", "discretization": "Z", "z": 1.,
            "lambda": Lambda, "nmax": nmax,
            "chain_onsite": onsite.tolist(), "chain_hoppings": hoppings.tolist(),
            "impurity_hopping": math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi),
            "coefficient_units": "physical", "measure_error": float(weights.sum()-1),
        },
    }


def electron_ed(physical, bath_record):
    """Solve a single-reservoir electron Hamiltonian with explicit JW products.

    At most five signed bath levels (twelve electron modes) are supported.
    H_d = U/2*(n_d-1)^2 + detuning*(n_d-1) + field*S_d^z, bath pairing is
    -gap*(c_up^dagger*c_down^dagger + h.c.), and each hopping is
    sqrt(weight*Gamma/(pi*rho)). The disconnected BCS vacuum is subtracted.

    Branch energy and residual are in gap units; signed_gap = E_D-E_S in
    those same units. P0/P1/P2 are impurity charge probabilities; moment is
    <S_d^z>, using the +1/2 member of the doublet, not a Pauli expectation.
    These are branch energies, not a search over all ground-state sectors.

    The lowest even Sz=0 and odd Sz=1/2 subspaces must have S^2=0 and 3/4.
    Unexpected spin (including field-induced mixing) raises ValueError rather
    than mislabeling a state. Degenerate minima of the expected spin retain
    their energy but return None for all four basis-dependent observables.
    At Gamma=0 use full dense sectors to avoid missing repeated eigenvalues
    with single-vector Lanczos. A detected sparse degeneracy also falls back
    to a complete dense eigenspace, including its spin verification.
    """
    _validate_physical(physical)
    gap = float(physical["gap"])
    xi = np.asarray(bath_record["xi"], dtype=float)
    weights = np.asarray(bath_record["weights"], dtype=float)
    if xi.ndim != 1 or xi.size == 0 or weights.shape != xi.shape:
        raise ValueError("xi and weights must be nonempty one-dimensional arrays of equal size")
    if not np.all(np.isfinite(xi)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("bath energies must be finite and weights strictly positive")
    if xi.size > 5:
        raise ValueError("electron ED is limited to five signed bath levels (twelve modes)")
    for key, expected in (("delta", gap), ("bandwidth", physical["bandwidth"])):
        value = bath_record[key]
        if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
                or not math.isfinite(value) or not math.isclose(value, expected, rel_tol=1e-13, abs_tol=0)):
            raise ValueError(f"bath {key} must match the physical parameters")

    modes = 2 + 2*xi.size
    identity = eye(2, format="csr")
    z = csr_matrix([[1., 0.], [0., -1.]])
    lower = csr_matrix([[0., 1.], [0., 0.]])
    annihilators = []
    for mode in range(modes):
        operator = csr_matrix([[1.]])
        for site in range(modes-1, -1, -1):
            factor = lower if site == mode else (z if site < mode else identity)
            operator = kron(operator, factor, format="csr")
        annihilators.append(operator)

    identity = eye(2**modes, format="csr")
    up, down = annihilators[:2]
    n_up, n_down = up.T @ up, down.T @ down
    charge = n_up+n_down-identity
    hamiltonian = (physical["u"]/(2*gap) * (charge @ charge)
                   + physical["detuning"]/gap * charge
                   + physical["field"]/(2*gap) * (n_up-n_down))
    spin_raise = up.T @ down
    contact_hopping = math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi)
    for level, (energy, weight) in enumerate(zip(xi, weights, strict=True)):
        bath_up, bath_down = annihilators[2+2*level:4+2*level]
        hamiltonian += energy/gap * (bath_up.T @ bath_up + bath_down.T @ bath_down)
        pair = -bath_up.T @ bath_down.T
        hopping = contact_hopping/gap * math.sqrt(weight) * (bath_up.T @ up + bath_down.T @ down)
        hamiltonian += pair + pair.T + hopping + hopping.T
        spin_raise += bath_up.T @ bath_down
    bath_reference = float(np.sum(xi - np.hypot(xi, gap)))
    hamiltonian -= bath_reference/gap * identity

    occupations = (np.arange(2**modes)[:, None] >> np.arange(modes)) & 1
    parity = occupations.sum(axis=1) % 2
    twice_sz = occupations[:, ::2].sum(axis=1) - occupations[:, 1::2].sum(axis=1)
    branches = {}
    for name, p, sz in (("singlet", 0, 0.), ("doublet", 1, .5)):
        indices = np.flatnonzero((parity == p) & (twice_sz == 2*sz))
        sector = hamiltonian[indices][:, indices].tocsr()
        scale = max(1., float(abs(sector).sum(axis=1).max()))
        tolerance = max(1e-10, 64*np.finfo(float).eps*scale)
        dense = len(indices) <= _DENSE_THRESHOLD or physical["gamma"] == 0
        if dense:
            energies, vectors = eigh(sector.toarray())
        else:
            energies, vectors = eigsh(sector, k=2, which="SA", tol=1e-12,
                                     v0=np.random.default_rng(137).normal(size=len(indices)))
            order = np.argsort(energies)
            energies, vectors = energies[order], vectors[:, order]
            if energies[1]-energies[0] <= tolerance:
                energies, vectors = eigh(sector.toarray())
                dense = True
        count = int(np.count_nonzero(energies-energies[0] <= tolerance))
        vectors = vectors[:, :count]
        residual = float(np.linalg.norm(sector @ vectors - vectors*energies[:count], axis=0).max())
        if residual > max(1e-9, 256*np.finfo(float).eps*scale):
            raise ValueError(f"{name} eigensolver residual is too large: {residual}")
        # S^2 = S^- S^+ + Sz*(Sz+1). Both named branches are highest weights;
        # check the whole minimum subspace, not a convenient degenerate vector.
        raised = spin_raise[:, indices] @ vectors
        spin_excess = float(np.linalg.norm(raised)**2)
        if spin_excess > 1e-8:
            raise ValueError(f"lowest {name} sector is not purely spin {sz}: "
                             f"S^2 exceeds {sz*(sz+1)} (possibly a mixed-spin degeneracy)")
        row = {
            "energy": float(energies[0]), "residual": residual,
            "spin_squared": float(sz*(sz+1) + spin_excess/count),
            "degeneracy": count, "sector_dimension": len(indices),
            "method": "dense" if dense else "sparse",
            "P0": None, "P1": None, "P2": None, "moment": None,
        }
        if count == 1:
            probability = np.abs(vectors[:, 0])**2
            n_up, n_down = occupations[indices, 0], occupations[indices, 1]
            row.update(
                P0=float(probability @ ((1-n_up)*(1-n_down))),
                P1=float(probability @ (n_up+n_down-2*n_up*n_down)),
                P2=float(probability @ (n_up*n_down)),
                moment=float(probability @ ((n_up-n_down)/2)),
            )
        branches[name] = row
    return {
        "branches": branches,
        "signed_gap": branches["doublet"]["energy"] - branches["singlet"]["energy"],
        "energy_reference": "centered impurity; disconnected BCS vacuum subtracted; in gap units",
        "residual_units": "gap",
    }
