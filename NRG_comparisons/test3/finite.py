"""Bounded finite-bath longitudinal comparison; no campaign is run on import."""

import math

import numpy as np
from threadpoolctl import threadpool_limits

from NRG_comparisons.test3.manual_run.physics import PRIMARY, TOLERANCE, finish_states, state_values


def run_finite(task, physical, bath_record):
    """Run explicitly selected roots/sectors with ED, QP or DMRG.

    ``bath_record`` contains explicit xi/weights/delta/bandwidth, without any
    fitting or weight normalization. Energy and residual outputs use gap units.
    ``qualified`` concerns the three primary roots only, never bath convergence
    or completeness. Selected-sector controls can omit primary branches.
    """
    backend = task["backend"]
    if backend not in ("ed", "qp", "dmrg"):
        raise ValueError("backend must be ed, qp or dmrg")
    if physical["geometry"] != "single":
        raise ValueError("test3 requires single-reservoir geometry")
    gap, bandwidth, u, gamma, detuning, field = (
        float(physical[name]) for name in ("gap", "bandwidth", "u", "gamma", "detuning", "field")
    )
    if not all(math.isfinite(x) for x in (gap, bandwidth, u, gamma, detuning, field)):
        raise ValueError("physical parameters must be finite")
    if gap <= 0 or bandwidth <= 0 or gamma < 0:
        raise ValueError("gap/bandwidth must be positive and gamma nonnegative")
    if bath_record["delta"] != gap or bath_record["bandwidth"] != bandwidth:
        raise ValueError("bath gap/bandwidth differ from physical parameters")
    xi, weights = np.asarray(bath_record["xi"], float), np.asarray(bath_record["weights"], float)
    if (xi.ndim != 1 or not xi.size or weights.shape != xi.shape or
            not np.all(np.isfinite(xi)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0)):
        raise ValueError("explicit bath requires finite xi and positive weights of equal size")
    spins = list(task.get("spins", [0, -1, 1, -2, 2]))
    roots = task.get("roots", 2)
    if not spins or any(type(spin) is not int for spin in spins) or len(set(spins)) != len(spins):
        raise ValueError("spins must be a nonempty list of distinct twice-Sz integers")
    if type(roots) is not int or roots < 1:
        raise ValueError("roots must be a positive integer")
    cutoff = task.get("cutoff")
    if backend != "qp" and cutoff is not None:
        raise ValueError("only QP accepts a cutoff")
    if cutoff is not None and (type(cutoff) is not int or not 0 <= cutoff <= 2*len(xi)):
        raise ValueError("cutoff must lie between zero and the bath mode count")
    states, sectors = [], []
    diagnostics = dict(backend=backend, sectors=sectors, bath_convergence_checked=False,
                       spectral_completeness_established=False, energy_units="gap",
                       residual_scope="projected_cutoff_hamiltonian" if backend == "qp"
                       else "untruncated_finite_hamiltonian")
    with threadpool_limits(limits=1):
        if backend == "ed":
            _electron_ed(spins, roots, xi, weights, gap, bandwidth, u, gamma, detuning, field,
                         states, sectors)
        else:
            _solver(task, backend, spins, roots, cutoff, physical, bath_record, states, sectors)
    coverage = finish_states(states)
    primary = {row["id"]: row["finite_converged"] and not row["degenerate"]
               for row in states if row["id"] in PRIMARY.values()}
    diagnostics["primary_qualified"] = primary
    diagnostics["all_returned_states_converged"] = all(row["finite_converged"] for row in states)
    diagnostics["qualification_scope"] = "three primary root-0 branches; exclusion roots reported separately"
    return dict(states=states, coverage=coverage, diagnostics=diagnostics,
                qualified=set(primary) == set(PRIMARY.values()) and all(primary.values()))


def _electron_ed(spins, roots, xi, weights, gap, bandwidth, u, gamma, detuning, field, states, sectors):
    """Independent electron occupations/Jordan-Wigner signs, at most 256 kets."""
    from scipy.sparse import csr_matrix, diags

    if len(xi) > 3:
        raise ValueError("independent ED is bounded to at most three signed bath levels")
    modes = 2+2*len(xi)
    masks = np.arange(1 << modes)
    occupied = (masks[:, None] >> np.arange(modes)) & 1
    charge = occupied[:, 0]+occupied[:, 1]
    double = occupied[:, 0]*occupied[:, 1]
    md = (occupied[:, 0]-occupied[:, 1])/2
    mb = (occupied[:, 2::2]-occupied[:, 3::2]).sum(axis=1)/2
    diagonal = u/2*(charge-1)**2+detuning*(charge-1)+field*md
    diagonal += (occupied[:, 2::2]+occupied[:, 3::2]) @ xi
    reference = float(np.sum(xi-np.hypot(xi, gap)))
    matrix = diags(diagonal-reference, format="csr")
    annihilators = []
    for mode in range(modes):
        source = masks[occupied[:, mode] == 1]
        signs = 1-2*(occupied[source, :mode].sum(axis=1) % 2)
        annihilators.append(csr_matrix((signs, (source ^ (1 << mode), source)),
                                      shape=matrix.shape, dtype=float))
    for level, weight in enumerate(weights):
        up, down = annihilators[2+2*level:4+2*level]
        pair = -gap*(up.T @ down.T)
        matrix += pair+pair.T
        for spin, bath_c in enumerate((up, down)):
            hopping = math.sqrt(2*bandwidth*gamma*weight/math.pi)*(bath_c.T @ annihilators[spin])
            matrix += hopping+hopping.T
    for spin in spins:
        selected = np.flatnonzero(2*(md+mb) == spin)
        dim = len(selected)
        info = dict(twice_sz=spin, parity=spin % 2, dimension=dim, requested_roots=roots,
                    returned_roots=min(roots, dim), method="electron-space-full-dense",
                    energy_reference=reference, complete_sector_diagonalized=True)
        sectors.append(info)
        if not dim:
            continue
        block = matrix[selected][:, selected].toarray()
        energies, vectors = np.linalg.eigh(block)
        for root in range(min(roots, dim)):
            vector, energy = vectors[:, root], energies[root]
            probabilities = abs(vector)**2
            measured = [float(probabilities @ op[selected]) for op in
                        (charge, double, md, md**2, mb, md*mb, mb**2)]
            values = state_values(*measured[:4], spin/2, m_bath=measured[4],
                                  Czz=measured[5], bath_z2=measured[6])
            # Independent diagonal projectors also check the shared probability map.
            for name, op in (("P0", charge == 0), ("P1", charge == 1), ("P2", double),
                             ("P_up", (charge == 1) & (md > 0)),
                             ("P_down", (charge == 1) & (md < 0))):
                value = float(probabilities @ op[selected])
                if abs(values[name]-value) > TOLERANCE:
                    raise ValueError("electron projector probability mismatch")
                values[name] = value
            values["energy"] = float(energy/gap)
            residual = float(np.linalg.norm(block @ vector-energy*vector)/gap)
            states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                               residual=residual, finite_converged=residual <= 1e-8,
                               degenerate=bool(np.count_nonzero(abs(energies-energy)/gap <= TOLERANCE) > 1)))


def _solver(task, backend, spins, roots, cutoff, physical, bath_record, states, sectors):
    from qdjj_solver import DiscreteBath, Sector, reference_model, solve
    from qdjj_solver.common import number
    from qdjj_solver.common.io import json_value

    gap = physical["gap"]
    bath = DiscreteBath.from_record(bath_record)
    h = reference_model(bath, **{name: physical[name] for name in
                                ("geometry", "u", "gamma", "detuning", "field")},
                        symmetry=False, compress=False)
    dot = (number(0)-number(1))/2
    bath_z = sum((number(i)-number(i+1))/2 for i in range(2, len(h.spins), 2))
    # Form complete electron polynomials before either coordinate change or projection.
    for name, operator in dict(dot_z2=dot*dot, bathSz=bath_z, bathSz2=bath_z*bath_z,
                               dotSz_bathSz=dot*bath_z, totSz=dot+bath_z).items():
        h.observables[name] = h.physical_operator(operator)
    if backend == "qp":
        from qdjj_solver.qp_solver.basis import dimension

        limits = dict(max_dimension=task.get("max_dimension", 1_000_000),
                      max_memory_gib=task.get("max_memory_gib", 8.), max_nnz=task.get("max_nnz", 50_000_000))
        if (type(limits["max_dimension"]) is not int or not 0 < limits["max_dimension"] <= 1_000_000
                or type(limits["max_nnz"]) is not int or not 0 < limits["max_nnz"] <= 50_000_000
                or not 0 < limits["max_memory_gib"] <= 8.):
            raise ValueError("QP resource gate: at most 1M states, 8 GiB, 50M nnz (16 GiB parent)")
        dimensions = [dimension(h.nimp, h.spins, h.bath_modes if cutoff is None else cutoff,
                                Sector(spin % 2, spin)) for spin in spins]
        if max(dimensions) > limits["max_dimension"]:
            raise MemoryError("QP sector dimension exceeds the explicit resource gate")
    else:
        chi = task.get("chi", 32)
        min_sweeps, max_sweeps = task.get("min_sweeps", 8), task.get("max_sweeps", 40)
        if chi not in (32, 64, 128) or not 8 <= min_sweeps <= max_sweeps <= 40:
            raise ValueError("DMRG pilot requires chi 32/64/128 and 8 <= min_sweeps <= max_sweeps <= 40")
        orbitals = len(h.spins)//2
        dimensions = [math.comb(2*orbitals, orbitals+spin) if abs(spin) <= orbitals else 0
                      for spin in spins]
    for spin, dim in zip(spins, dimensions, strict=True):
        count = min(roots, dim)
        info = dict(twice_sz=spin, parity=spin % 2, dimension=dim, requested_roots=roots,
                    returned_roots=count)
        sectors.append(info)
        if not dim:
            continue
        if backend == "qp":
            method = task.get("method", "dense" if dim <= 128 else "sparse")
            if method not in ("dense", "sparse"):
                raise ValueError("QP pilot method must be dense or sparse")
            # Tiny controls retain the whole degenerate subspace. Sparse runs use
            # one lookahead root, but do not claim to exhaust hidden degeneracies.
            if dim <= count+2:
                method = "dense"
            solved_roots = dim if method == "dense" else min(count+1, dim)
            options = dict(method=method, eigenpairs=solved_roots, residual_tolerance=1e-8*gap,
                           tolerance=1e-12, threads=1, seed=task.get("seed", 1729), **limits)
        else:
            solved_roots = count
            options = dict(eigenpairs=count, chi_max=chi, min_sweeps=min_sweeps, max_sweeps=max_sweeps,
                           energy_tolerance=1e-12*gap, residual_tolerance=1e-7*gap,
                           lanczos_probability_tolerance=1e-22, seed_trials=2,
                           calculate_residuals=True, require_convergence=False, threads=1,
                           # TeNPy's two-site sweep requires more than two sites.
                           group_size=1 if len(h.spins) == 4 else 2,
                           seed=task.get("seed", 1729))
        result = solve(h, backend=backend, cutoff=cutoff, sector=Sector(spin % 2, spin), options=options)
        info.update(solved_roots=solved_roots, solver_metadata=json_value(result.metadata),
                    timings=json_value(result.timings),
                    complete_sector_diagonalized=backend == "qp" and solved_roots == dim)
        for root in range(count):
            obs = {name: values[root] for name, values in result.observables.items()}
            if abs(obs["totSz"]-spin/2) > TOLERANCE:
                raise ValueError("direct total Sz does not match the requested sector")
            values = state_values(obs["impurity_charge"], obs["double_occupancy_0"],
                                  obs["impurity_spin_z"], obs["dot_z2"], spin/2,
                                  m_bath=obs["bathSz"], Czz=obs["dotSz_bathSz"], bath_z2=obs["bathSz2"])
            energy = float(result.energies[root]/gap)
            values["energy"] = energy
            residual = None if result.residuals is None else float(result.residuals[root]/gap)
            converged = residual is not None and math.isfinite(residual) and residual <= (1e-8 if backend == "qp" else 1e-7)
            if backend == "dmrg":
                root_info = result.metadata["roots"][root]
                # A failed exclusion root must not disqualify an accurate primary.
                converged &= bool(root_info["sweep_energy_converged"] and
                                  root_info.get("lower_root_overlap", 0.) <= 1e-8 and
                                  abs(result.metadata["gram"][root, root]-1) <= 1e-8 and
                                  np.linalg.norm(result.metadata["gram"][root, :root]) <= 1e-8 and
                                  np.all(result.energies[:root] <= result.energies[root]+1e-12*gap) and
                                  np.all(result.energies[root:] >= result.energies[root]-1e-12*gap))
            states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                               residual=residual, finite_converged=bool(converged),
                               degenerate=bool(np.count_nonzero(abs(result.energies/gap-energy) <= TOLERANCE) > 1)))
