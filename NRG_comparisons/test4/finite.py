"""Bounded physical L/R ASQ junction controls; no calculation on import."""

import math

import numpy as np
from threadpoolctl import threadpool_limits

from NRG_comparisons.test4.physics import (
    PRIMARY, TOLERANCE, detuning, finish_states, paper_offset, state_values, validate_physical,
)


def run_finite(task, physical, bath_record):
    """Return states, coverage, diagnostics, qualified for ED/QP/DMRG.

    One explicit bath record is reused identically in both physical leads,
    without normalizing its weights. Roots may be overridden by string-keyed
    ``roots_by_spin``. Qualification applies to the three primary roots only,
    not bath/cutoff convergence, spectral completeness, or a subgap count.
    """
    backend = task["backend"]
    if backend not in ("ed", "qp", "dmrg"):
        raise ValueError("backend must be ed, qp or dmrg")
    physical = validate_physical(physical)
    gap, bandwidth = physical["gap"], physical["bandwidth"]
    if bath_record["delta"] != gap or bath_record["bandwidth"] != bandwidth:
        raise ValueError("bath gap/bandwidth differ from physical parameters")
    xi, weights = np.asarray(bath_record["xi"], float), np.asarray(bath_record["weights"], float)
    if (xi.ndim != 1 or not xi.size or weights.shape != xi.shape or
            not np.all(np.isfinite(xi)) or not np.all(np.isfinite(weights)) or np.any(weights <= 0)):
        raise ValueError("explicit bath requires finite xi and positive weights of equal size")
    spins = list(task.get("spins", [0, -1, 1, -2, 2]))
    if not spins or any(type(spin) is not int for spin in spins) or len(set(spins)) != len(spins):
        raise ValueError("spins must be a nonempty list of distinct twice-Sz integers")
    roots = task.get("roots", 2)
    by_spin = task.get("roots_by_spin", {})
    if type(roots) is not int or roots < 1:
        raise ValueError("roots must be a positive integer")
    if (not isinstance(by_spin, dict) or any(key not in {str(s) for s in spins} for key in by_spin)
            or any(type(value) is not int or value < 1 for value in by_spin.values())):
        raise ValueError("roots_by_spin requires selected spin string keys and positive integer roots")
    counts = {spin: by_spin.get(str(spin), roots) for spin in spins}
    cutoff = task.get("cutoff")
    if backend != "qp" and cutoff is not None:
        raise ValueError("only QP accepts a cutoff")
    if cutoff is not None and (type(cutoff) is not int or not 0 <= cutoff <= 4*len(xi)):
        raise ValueError("cutoff must lie between zero and the combined L/R bath mode count")
    states, sectors = [], []
    diagnostics = dict(backend=backend, sectors=sectors, requested_spins=spins,
                       bath_levels_per_lead=len(xi), bath_leads=2, bath_modes_total=4*len(xi),
                       identical_bath_record_for_both_leads=True, bath_weight_sum_per_lead=float(weights.sum()),
                       bath_convergence_checked=False, spectral_completeness_established=False,
                       energy_units="gap", energy_reference="centered impurity; both isolated BCS vacua subtracted",
                       paper_energy_offset=paper_offset(physical)/gap, detuning=detuning(physical),
                       current_units="I/(2e*Delta/hbar)", current_is_real=True,
                       residual_scope="projected_cutoff_hamiltonian" if backend == "qp"
                       else "untruncated_finite_hamiltonian")
    with threadpool_limits(limits=1):
        if backend == "ed":
            _electron_ed(spins, counts, xi, weights, physical, states, sectors)
        else:
            _solver(task, backend, spins, counts, cutoff, physical, bath_record, states, sectors)
    coverage = finish_states(states)
    primary = {row["id"]: row["finite_converged"] and not row["degenerate"]
               for row in states if row["id"] in PRIMARY.values()}
    diagnostics["primary_qualified"] = primary
    diagnostics["primary_unique"] = {row["id"]: not row["degenerate"] for row in states
                                     if row["id"] in PRIMARY.values()}
    diagnostics["primary_residual_qualified"] = {
        row["id"]: row["residual"] is not None and row["residual"] <= (1e-7 if backend == "dmrg" else 1e-8)
        for row in states if row["id"] in PRIMARY.values()}
    diagnostics["all_returned_states_converged"] = all(row["finite_converged"] for row in states)
    diagnostics["qualification_scope"] = "three primary root-0 branches; extra roots reported separately"
    return dict(states=states, coverage=coverage, diagnostics=diagnostics,
                qualified=set(primary) == set(PRIMARY.values()) and all(primary.values()))


def _electron_ed(spins, counts, xi, weights, physical, states, sectors):
    """Independent complex Jordan-Wigner construction, at most 1024 electron kets."""
    from scipy.sparse import csr_matrix, diags

    if len(xi) > 2:
        raise ValueError("independent ED is bounded to two signed bath levels per lead")
    gap, bandwidth, u, gamma, nu, field, phi, wn, ws = (
        physical[name] for name in ("gap", "bandwidth", "u", "gamma", "nu", "field", "phi", "wn", "ws")
    )
    modes = 2+4*len(xi)
    masks = np.arange(1 << modes)
    occupied = (masks[:, None] >> np.arange(modes)) & 1
    charge = occupied[:, 0]+occupied[:, 1]
    double = occupied[:, 0]*occupied[:, 1]
    md = (occupied[:, 0]-occupied[:, 1])/2
    mb = (occupied[:, 2::2]-occupied[:, 3::2]).sum(axis=1)/2
    # Literal paper charging energy, independent of the solver's centered form.
    diagonal = u/2*(charge-nu)*(charge-nu)+field*md
    diagonal += (occupied[:, 2::2]+occupied[:, 3::2]) @ np.tile(xi, 2)
    reference = float(2*np.sum(xi-np.hypot(xi, gap)))
    matrix = diags(diagonal-reference, format="csr", dtype=complex)
    derivative = csr_matrix(matrix.shape, dtype=complex)
    annihilators = []
    for mode in range(modes):
        source = masks[occupied[:, mode] == 1]
        signs = 1-2*(occupied[source, :mode].sum(axis=1) % 2)
        annihilators.append(csr_matrix((signs, (source ^ (1 << mode), source)),
                                      shape=matrix.shape, dtype=complex))
    for lead, direction in enumerate((-1, 1)):
        for level, weight in enumerate(weights):
            start = 2+2*(lead*len(xi)+level)
            up, down = annihilators[start:start+2]
            pair = -gap*(up.getH() @ down.getH())
            matrix += pair+pair.getH()
            for spin, bath_c in enumerate((up, down)):
                hopping = (math.sqrt(bandwidth*gamma/math.pi)*math.sqrt(weight)
                           * np.exp(direction*1j*phi/4)*(bath_c.getH() @ annihilators[spin]))
                matrix += hopping+hopping.getH()
                differentiated = direction*1j/4*hopping
                derivative += differentiated+differentiated.getH()
    for i, wi in enumerate(weights):
        for j, wj in enumerate(weights):
            for spin, sigma in enumerate((1, -1)):
                left = annihilators[2+2*i+spin]
                right = annihilators[2+2*len(xi)+2*j+spin]
                bridge = math.sqrt(wi*wj)*(wn+1j*sigma*ws)*np.exp(-1j*phi/2)*(left.getH() @ right)
                matrix += bridge+bridge.getH()
                differentiated = -1j/2*bridge
                derivative += differentiated+differentiated.getH()
    for spin in spins:
        selected = np.flatnonzero(2*(md+mb) == spin)
        dim, roots = len(selected), counts[spin]
        info = dict(twice_sz=spin, parity=spin % 2, dimension=dim, requested_roots=roots,
                    returned_roots=min(roots, dim), solved_roots=dim, method="electron-space-full-dense",
                    energy_reference=reference, complete_sector_diagonalized=True)
        sectors.append(info)
        if not dim:
            continue
        block = matrix[selected][:, selected].toarray()
        current = derivative[selected][:, selected]
        energies, vectors = np.linalg.eigh(block)
        for root in range(min(roots, dim)):
            vector, energy = vectors[:, root], energies[root]
            probabilities = abs(vector)**2
            measured = [float(probabilities @ op[selected]) for op in
                        (charge, double, md, md**2, mb, md*mb, mb**2)]
            values = state_values(*measured[:4], spin/2, m_bath=measured[4],
                                  Czz=measured[5], bath_z2=measured[6])
            for name, op in (("P0", charge == 0), ("P1", charge == 1), ("P2", double),
                             ("P_up", (charge == 1) & (md > 0)),
                             ("P_down", (charge == 1) & (md < 0))):
                value = float(probabilities @ op[selected])
                if abs(values[name]-value) > TOLERANCE:
                    raise ValueError("electron projector probability mismatch")
                values[name] = value
            values["energy"] = float((energy-paper_offset(physical))/gap)
            values["energy_paper"] = values["energy"]+paper_offset(physical)/gap
            # finish_states checks that this Hermitian scalar has no imaginary part.
            values["current"] = np.vdot(vector, current @ vector)/gap
            residual = float(np.linalg.norm(block @ vector-energy*vector)/gap)
            states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                               residual=residual, finite_converged=residual <= 1e-8,
                               degenerate=bool(np.count_nonzero(abs(energies-energy)/gap <= TOLERANCE) > 1)))


def _solver(task, backend, spins, counts, cutoff, physical, bath_record, states, sectors):
    from qdjj_solver import DiscreteBath, Sector, reference_model, solve
    from qdjj_solver.common import number
    from qdjj_solver.common.io import json_value

    gap = physical["gap"]
    bath = DiscreteBath.from_record(bath_record)
    # Leave the solver Hamiltonian centered; the paper constant is output-only.
    h = reference_model(bath, **{name: physical[name] for name in
                                ("geometry", "u", "gamma", "field", "phi")},
                        detuning=detuning(physical), rho_wn=physical["wn"]/(2*physical["bandwidth"]),
                        rho_ws=physical["ws"]/(2*physical["bandwidth"]), symmetry=False, compress=False)
    h.observables.pop("rho_ws_derivative", None)
    dot = (number(0)-number(1))/2
    bath_z = sum((number(i)-number(i+1))/2 for i in range(2, len(h.spins), 2))
    # Complete electron polynomials include both leads, before projection.
    for name, operator in dict(dot_z2=dot*dot, bathSz=bath_z, bathSz2=bath_z*bath_z,
                               dotSz_bathSz=dot*bath_z, totSz=dot+bath_z).items():
        h.observables[name] = h.physical_operator(operator)
    if backend == "qp":
        from qdjj_solver.qp_solver.basis import dimension

        limits = dict(max_dimension=task.get("max_dimension", 1_000_000),
                      max_memory_gib=task.get("max_memory_gib", 8.), max_nnz=task.get("max_nnz", 50_000_000))
        if (type(limits["max_dimension"]) is not int or not 0 < limits["max_dimension"] <= 1_000_000
                or type(limits["max_nnz"]) is not int or not 0 < limits["max_nnz"] <= 50_000_000
                or type(limits["max_memory_gib"]) not in (int, float) or not 0 < limits["max_memory_gib"] <= 8.):
            raise ValueError("QP resource gate: at most 1M states, 8 GiB, 50M nnz")
        dimensions = [dimension(h.nimp, h.spins, h.bath_modes if cutoff is None else cutoff,
                                Sector(spin % 2, spin)) for spin in spins]
        if max(dimensions) > limits["max_dimension"]:
            raise MemoryError("QP sector dimension exceeds the explicit resource gate")
        if task.get("method", "sparse") not in ("dense", "sparse"):
            raise ValueError("QP method must be dense or sparse")
        if any(dim > 512 and (task.get("method") == "dense" or dim <= counts[spin]+2)
               for spin, dim in zip(spins, dimensions, strict=True)):
            raise MemoryError("QP dense control cap is 512 states; request fewer sparse roots")
    else:
        chi = task.get("chi", 32)
        min_sweeps, max_sweeps = task.get("min_sweeps", 8), task.get("max_sweeps", 40)
        if (type(chi) is not int or chi not in (32, 64, 128)
                or type(min_sweeps) is not int or type(max_sweeps) is not int
                or not 8 <= min_sweeps <= max_sweeps <= 40):
            raise ValueError("DMRG control requires chi 32/64/128 and 8 <= min_sweeps <= max_sweeps <= 40")
        orbitals = len(h.spins)//2
        dimensions = [math.comb(2*orbitals, orbitals+spin) if abs(spin) <= orbitals else 0
                      for spin in spins]
    for spin, dim in zip(spins, dimensions, strict=True):
        count = min(counts[spin], dim)
        info = dict(twice_sz=spin, parity=spin % 2, dimension=dim, requested_roots=counts[spin],
                    returned_roots=count)
        sectors.append(info)
        if not dim:
            continue
        if backend == "qp":
            method = task.get("method", "dense" if dim <= 512 else "sparse")
            if dim <= count+2:
                method = "dense"
            # Dense controls see every degeneracy; sparse lookahead does not
            # establish completeness of a repeated eigenspace.
            solved_roots = dim if method == "dense" else count+1
            options = dict(method=method, eigenpairs=solved_roots, residual_tolerance=1e-8*gap,
                           tolerance=1e-12, threads=1, seed=task.get("seed", 1729), **limits)
        else:
            solved_roots = count
            options = dict(eigenpairs=count, chi_max=chi, min_sweeps=min_sweeps, max_sweeps=max_sweeps,
                           energy_tolerance=1e-12*gap, residual_tolerance=1e-7*gap,
                           lanczos_probability_tolerance=1e-22, seed_trials=2,
                           calculate_residuals=True, require_convergence=False, threads=1,
                           group_size=2, seed=task.get("seed", 1729))
        result = solve(h, backend=backend, cutoff=cutoff, sector=Sector(spin % 2, spin), options=options)
        info.update(solved_roots=solved_roots, solver_metadata=json_value(result.metadata),
                    timings=json_value(result.timings),
                    complete_sector_diagonalized=backend == "qp" and solved_roots == dim,
                    degeneracy_scope="complete sector" if backend == "qp" and solved_roots == dim
                    else "computed roots only; hidden degeneracies not excluded")
        for root in range(count):
            obs = {name: values[root] for name, values in result.observables.items()}
            if abs(obs["totSz"]-spin/2) > TOLERANCE:
                raise ValueError("direct total Sz does not match the requested sector")
            values = state_values(obs["impurity_charge"], obs["double_occupancy_0"],
                                  obs["impurity_spin_z"], obs["dot_z2"], spin/2,
                                  m_bath=obs["bathSz"], Czz=obs["dotSz_bathSz"], bath_z2=obs["bathSz2"])
            energy = float(result.energies[root]/gap)
            values.update(energy=energy, energy_paper=energy+paper_offset(physical)/gap,
                          current=obs["phase_derivative"]/gap)
            residual = None if result.residuals is None else float(result.residuals[root]/gap)
            converged = residual is not None and math.isfinite(residual) and residual <= (1e-8 if backend == "qp" else 1e-7)
            if backend == "dmrg":
                root_info = result.metadata["roots"][root]
                converged &= bool(root_info["sweep_energy_converged"] and
                                  root_info.get("lower_root_overlap", 0.) <= 1e-8 and
                                  abs(result.metadata["gram"][root, root]-1) <= 1e-8 and
                                  np.linalg.norm(result.metadata["gram"][root, :root]) <= 1e-8 and
                                  np.all(result.energies[:root] <= result.energies[root]+1e-12*gap) and
                                  np.all(result.energies[root:] >= result.energies[root]-1e-12*gap))
            states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                               residual=residual, finite_converged=bool(converged),
                               degenerate=bool(np.count_nonzero(abs(result.energies/gap-energy) <= TOLERANCE) > 1)))
