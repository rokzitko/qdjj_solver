"""Finite-bath DMRG with sequential root targeting and physical residual checks."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from copy import deepcopy
import logging
from time import perf_counter
from typing import ClassVar

import numpy as np
import tenpy
from tenpy.algorithms import dmrg
from tenpy.linalg import np_conserved as npc
from tenpy.networks.mps import MPS
from threadpoolctl import threadpool_limits

from ..common.problem import DEFAULT_SECTOR, Hamiltonian
from .mpo import prepare, residual_norm


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SolverOptions:
    eigenpairs: int = 1
    chi_max: int = 128
    max_sweeps: int = 60
    min_sweeps: int = 10
    energy_tolerance: float = 1e-11
    entropy_tolerance: float = 1e-8
    residual_tolerance: float = 1e-7
    orthogonality_tolerance: float = 1e-8
    svd_min: float = 1e-14
    lanczos_maxiter: int = 100
    seed_trials: int = 4
    seed: int = 1729
    seed_randomization_steps: int = 4
    excitation_operator: str | None = None
    mixer: bool = True
    mixer_amplitude: float = 1e-5
    mixer_sweeps: int = 6
    group_size: int = 2
    mode_order: tuple[int, ...] | None = None
    observables: tuple[str, ...] | None = None
    calculate_residuals: bool = True
    require_convergence: bool = False
    threads: int = 1
    lanczos_probability_tolerance: float = 1e-14

    def __post_init__(self):
        for name in ("eigenpairs", "chi_max", "max_sweeps", "min_sweeps", "lanczos_maxiter",
                     "seed_trials", "mixer_sweeps", "threads"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("seed", "seed_randomization_steps"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.excitation_operator is not None and not isinstance(self.excitation_operator, str):
            raise ValueError("excitation_operator must name a model observable")
        for name in ("energy_tolerance", "entropy_tolerance", "residual_tolerance",
                     "orthogonality_tolerance", "mixer_amplitude", "lanczos_probability_tolerance"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not np.isfinite(self.svd_min) or self.svd_min < 0:
            raise ValueError("svd_min must be finite and nonnegative")
        if self.min_sweeps > self.max_sweeps:
            raise ValueError("min_sweeps exceeds max_sweeps")
        if self.mixer and self.mixer_sweeps >= self.min_sweeps:
            raise ValueError("min_sweeps must exceed mixer_sweeps")
        if (not isinstance(self.group_size, int) or isinstance(self.group_size, bool)
                or self.group_size not in (1, 2)):
            raise ValueError("group_size must be one or two")
        if self.require_convergence and not self.calculate_residuals:
            raise ValueError("require_convergence needs physical residual checks")


@dataclass
class Eigenstates:
    backend: ClassVar[str] = "dmrg"
    energies: np.ndarray
    states: list
    residuals: np.ndarray | None
    observables: dict[str, np.ndarray]
    timings: dict[str, float]
    metadata: dict
    hamiltonian: Hamiltonian
    prepared: object

    def _check_coordinates(self, other):
        if self.prepared.coordinate_id != other.prepared.coordinate_id:
            raise ValueError("states use different coordinates, layouts, or conserved charges")

    def overlaps(self, other=None):
        other = self if other is None else other
        self._check_coordinates(other)
        return np.array([[bra.overlap(ket) for ket in other.states] for bra in self.states])

    def matrix_elements(self, operator, ket=None):
        """All <self_i|operator|ket_j> in common prepared coordinates."""
        ket = self if ket is None else ket
        self._check_coordinates(ket)
        return self.prepared.matrix_elements(operator, self.states, ket.states)


def _product_seeds(prepared, count):
    """K lowest local-cost product states in the requested charge sector, via DP."""
    chinfo = prepared.sites[0].leg.chinfo
    paths = {tuple(chinfo.make_valid()): [(0., ())]}
    for site, costs in zip(prepared.sites, prepared.onsite_costs(), strict=True):
        next_paths = {}
        for total, candidates in paths.items():
            for state, charge in enumerate(site.leg.to_qflat()):
                new_charge = tuple(chinfo.make_valid(np.asarray(total)+charge))
                bucket = next_paths.setdefault(new_charge, [])
                bucket.extend((cost+costs[state], path+(state,)) for cost, path in candidates)
        paths = {charge: sorted(candidates)[:count] for charge, candidates in next_paths.items()}
    return [path for _, path in paths.get(prepared.target_charge, [])]


def _finite_stats(stats):
    return {name: [float(v) if np.isfinite(v) else None for v in values]
            for name, values in stats.items()}


def _randomize_seed(psi, rng, steps, chi):
    """Apply locally generated charge-block Haar gates without global RNG state."""
    for layer in range(steps):
        for i in range(layer % 2, psi.L-1, 2):
            left, right = psi.sites[i:i+2]
            charges = left.leg.to_qflat()[:, None, :]+right.leg.to_qflat()[None, :, :]
            charges = left.leg.chinfo.make_valid(charges.reshape(left.dim*right.dim, -1))
            blocks = {}
            for index, charge in enumerate(charges):
                blocks.setdefault(tuple(charge), []).append(index)
            gate = np.zeros((left.dim*right.dim,)*2, complex)
            for indices in blocks.values():
                shape = (len(indices),)*2
                matrix = rng.normal(size=shape)+1j*rng.normal(size=shape)
                q, r = np.linalg.qr(matrix)
                diagonal = np.diag(r)
                q *= (diagonal/np.maximum(abs(diagonal), np.finfo(float).tiny)).conj()
                gate[np.ix_(indices, indices)] = q
            operator = npc.Array.from_ndarray(gate.reshape(left.dim, right.dim, left.dim, right.dim),
                                               [left.leg, right.leg, left.leg.conj(), right.leg.conj()],
                                               labels=["p0", "p1", "p0*", "p1*"])
            psi.apply_local_op(i, operator, unitary=True, renormalize=True)
        psi.compress_svd(dict(chi_max=chi, svd_min=1e-14))
    # A one-site seed is already canonical; TeNPy's QR/SVD sweep requires L > 1.
    if psi.L > 1:
        psi.canonical_form(renormalize=True)
    psi.norm = 1.


def _operator_seed(prepared, source, operator):
    """Normalized fluctuation O|psi>-<O>|psi>, or None if it annihilates psi."""
    mean = prepared.matrix_element(operator, source, source)
    fluctuation = operator-mean
    scale = float(sum(abs(c) for c in fluctuation.terms.values()))
    if scale == 0.:
        return None
    mpo = prepared.compile(fluctuation/scale)
    candidate = source.copy()
    mpo.apply_naively(candidate)
    norm_squared = float(candidate.overlap(candidate, ignore_form=True).real)
    if norm_squared <= 1e-24:
        return None
    candidate.canonical_form(renormalize=True)
    candidate.norm = 1.
    return candidate


def _single_site(prepared, eigenpairs):
    site = prepared.sites[0]
    selected = np.flatnonzero(np.all(site.leg.to_qflat() == prepared.target_charge, axis=1))
    if eigenpairs > len(selected):
        raise ValueError("more eigenpairs requested than states in the sector")
    matrix = np.zeros((site.dim, site.dim), complex)
    for key, coefficient in prepared.hamiltonian.operator.terms.items():
        term = prepared.term(key)
        matrix += coefficient*site.get_op(" ".join(name for name, _ in term) if term else "Id").to_ndarray()
    energies, vectors = np.linalg.eigh(matrix[np.ix_(selected, selected)])
    states = []
    for vector in vectors[:, :eigenpairs].T:
        local = np.zeros(site.dim, complex)
        local[selected] = vector
        states.append(MPS.from_product_state([site], [local], bc="finite", dtype=complex,
                                            permute=False, unit_cell_width=1))
    return energies[:eigenpairs], states


DEFAULT_OPTIONS = SolverOptions()


def solve(hamiltonian, cutoff=None, sector=DEFAULT_SECTOR, options=DEFAULT_OPTIONS, initial=None):
    """Lowest states in the complete finite-mode space defined by ``hamiltonian``.

    ``cutoff`` must be None: a QP cutoff and a Schmidt-rank limit are different
    approximations. ``initial`` may be an earlier Eigenstates in identical
    coordinates/layout for bond-dimension or parameter continuation.

    Convergence here concerns the finite Hamiltonian, not bath discretization.
    A small residual identifies an eigenstate; independent seeds and root/bath
    ladders are still needed to establish that the desired roots were found.
    """
    if cutoff is not None:
        raise ValueError("DMRG has no QP cutoff; set cutoff=None and control chi_max")
    start = perf_counter()
    prepared = prepare(hamiltonian, sector, group_size=options.group_size, mode_order=options.mode_order)
    if initial is not None:
        if not isinstance(initial, Eigenstates) or prepared.coordinate_id != initial.prepared.coordinate_id:
            raise ValueError("initial states have incompatible coordinates or MPS layout")
        if sector != initial.prepared.sector:
            raise ValueError("initial states belong to a different sector")
    names = list(hamiltonian.observables) if options.observables is None else list(options.observables)
    if set(names)-hamiltonian.observables.keys():
        raise ValueError("unknown requested observable")
    excitation = None
    if options.excitation_operator is not None:
        if options.excitation_operator not in hamiltonian.observables:
            raise ValueError("excitation_operator must name a model observable")
        excitation = hamiltonian.observables[options.excitation_operator]
        if not prepared.conserves_charges(excitation):
            raise ValueError("excitation_operator must preserve the target sector charges")
    times = {"preparation": perf_counter()-start}
    states, energies, diagnostics = [], [], []
    bound = float(sum(abs(v) for v in hamiltonian.operator.terms.values()))
    shift = 0.
    with threadpool_limits(limits=options.threads):
        physical_mpo = prepared.compile(hamiltonian.operator)
        seeds = _product_seeds(prepared, options.eigenpairs*(options.seed_trials+1)+4)
        if len(seeds) < options.eigenpairs:
            raise ValueError("more eigenpairs requested than states in the sector")
        if physical_mpo is not None:
            # Min--max: the kth eigenvalue is bounded above by the largest
            # Ritz value in ANY k-dimensional trial space. Distinct product
            # seeds are exactly orthonormal, so this much tighter upper bound
            # avoids an extensive search shift from high-energy empty modes.
            trial = [MPS.from_product_state(prepared.sites, seed, bc="finite", dtype=complex,
                                            permute=False, unit_cell_width=1)
                     for seed in seeds[:options.eigenpairs]]
            trial_h = prepared.matrix_elements(hamiltonian.operator, trial, trial)
            upper = float(np.linalg.eigvalsh((trial_h+trial_h.conj().T)/2)[-1])
            shift = upper+max(1., .01*abs(upper))
        search_model = prepared.model(hamiltonian.operator-shift) if physical_mpo is not None else None
        solve_start = perf_counter()
        if len(prepared.sites) == 1:
            energies, states = _single_site(prepared, options.eigenpairs)
            diagnostics = [dict(sweep_converged=True, sweep_energy_converged=True, sweeps=0, chi_actual=1,
                                method="single-site-exact", seed_trials=[]) for _ in states]
        elif physical_mpo is None:
            states = [MPS.from_product_state(prepared.sites, seed, bc="finite", dtype=complex,
                                              permute=False, unit_cell_width=1)
                      for seed in seeds[:options.eigenpairs]]
            energies = [0.]*len(states)
            diagnostics = [dict(sweep_converged=True, sweep_energy_converged=True, sweeps=0, chi_actual=1,
                                method="zero-operator", seed_trials=[]) for _ in states]
        else:
            engine_options = dict(
                # TrivialLattice stores this finite chain as one bookkeeping
                # cell; its cell size is not a physical cylinder circumference.
                max_N_sites_per_ring=len(prepared.sites),
                # A deliberately coarse chi point is still useful in a ladder.
                # Retain TeNPy's warning and enforce our physical residual/Gram
                # checks through require_convergence rather than this heuristic.
                max_trunc_err=None,
                mixer=options.mixer,
                mixer_params=dict(amplitude=options.mixer_amplitude, decay=2., disable_after=options.mixer_sweeps),
                max_sweeps=options.max_sweeps, min_sweeps=options.min_sweeps, N_sweeps_check=2,
                # TeNPy divides by max(E, 1), not max(abs(E), 1). Negative
                # search energies make this an absolute per-sweep tolerance.
                max_E_err=options.energy_tolerance, max_S_err=options.entropy_tolerance,
                trunc_params=dict(chi_max=options.chi_max, svd_min=options.svd_min, trunc_cut=None),
                diag_method="lanczos", P_tol_to_trunc=None, E_tol_to_trunc=None,
                lanczos_params=dict(N_min=2, N_max=options.lanczos_maxiter, P_tol=options.lanczos_probability_tolerance,
                                    E_tol=options.energy_tolerance/(10*max(1., abs(shift)))))
            for root in range(options.eigenpairs):
                candidates = []
                for attempt in range(options.seed_trials):
                    if attempt == 0 and initial is not None and root < len(initial.states):
                        psi = initial.states[root].copy()
                        seed_description = "continuation"
                    else:
                        index = (root+attempt*options.eigenpairs) % len(seeds)
                        psi = MPS.from_product_state(prepared.sites, seeds[index], bc="finite", dtype=complex,
                                                      permute=False, unit_cell_width=1)
                        seed_description = list(seeds[index])
                        if attempt == 0 and root > 0 and excitation is not None:
                            candidate = _operator_seed(prepared, states[root-1], excitation)
                            if candidate is not None:
                                psi = candidate
                                seed_description = dict(operator=options.excitation_operator, source_root=root-1)
                        if attempt > 0 and options.seed_randomization_steps:
                            # Charge-preserving random gates expose components
                            # in additional, undeclared conserved subspaces;
                            # product seeds alone can lock onto a higher root.
                            rng = np.random.default_rng(np.random.SeedSequence([options.seed, root, attempt]))
                            _randomize_seed(psi, rng, options.seed_randomization_steps, min(8, options.chi_max))
                            seed_description = dict(product=seed_description,
                                                    random_seed=[options.seed, root, attempt],
                                                    randomization_steps=options.seed_randomization_steps)
                    engine = dmrg.TwoSiteDMRGEngine(psi, search_model, deepcopy(engine_options),
                                                    orthogonal_to=states)
                    _, psi = engine.run()
                    psi.canonical_form(renormalize=True)
                    psi.norm = 1.
                    energy = float(physical_mpo.expectation_value(psi).real)
                    last_energies = engine.sweep_stats["E"][-3:]
                    stable_energy = (len(last_energies) >= 3 and
                                     np.ptp(last_energies) <= 4*options.energy_tolerance)
                    info = dict(sweep_converged=bool(engine.is_converged()), sweeps=int(engine.sweeps),
                                sweep_energy_converged=bool(stable_energy),
                                chi_actual=int(max(psi.chi)), seed=seed_description,
                                lower_root_overlap=max((abs(psi.overlap(state)) for state in states), default=0.),
                                norm_error=float(np.linalg.norm(psi.norm_test())),
                                sweep_statistics=_finite_stats(engine.sweep_stats))
                    candidates.append((energy, psi, info))
                    logger.info("root=%d trial=%d E=%.14g sweeps=%d chi=%d lower-overlap=%.2e elapsed=%.1fs",
                                root, attempt, energy, engine.sweeps, info["chi_actual"],
                                info["lower_root_overlap"], perf_counter()-solve_start)
                eligible = [i for i, (_, _, d) in enumerate(candidates)
                            if d["lower_root_overlap"] <= options.orthogonality_tolerance]
                best = min(eligible or range(len(candidates)), key=lambda i: candidates[i][0])
                energy, psi, info = candidates[best]
                info["seed_trials"] = [dict(energy=e, sweep_converged=d["sweep_converged"],
                                            lower_root_overlap=d["lower_root_overlap"],
                                            seed=d["seed"]) for e, _, d in candidates]
                states.append(psi)
                energies.append(energy)
                diagnostics.append(info)
        times["diagonalization"] = perf_counter()-solve_start
        measure_start = perf_counter()
        energies = np.asarray(energies, dtype=float)
        residuals = (np.array([residual_norm(prepared, psi, energy)
                               for psi, energy in zip(states, energies, strict=True)])
                     if options.calculate_residuals else None)
        gram = np.array([[bra.overlap(ket) for ket in states] for bra in states])
        gram_error = float(np.linalg.norm(gram-np.eye(len(states))))
        observables = {}
        for name in names:
            op = hamiltonian.observables[name]
            observables[name] = np.real_if_close(np.diag(prepared.matrix_elements(op, states, states)))
        times["observables_and_residuals"] = perf_counter()-measure_start
    times["total"] = perf_counter()-start
    ordered = bool(np.all(np.diff(energies) >= -options.energy_tolerance))
    # Entropy can rotate within a degenerate eigenspace even when the physical
    # state residual and sweep energies have converged. Preserve TeNPy's full
    # stop flag separately rather than equating that flag with eigenaccuracy.
    converged = bool(all(d["sweep_energy_converged"] for d in diagnostics) and ordered
                     and gram_error <= options.orthogonality_tolerance
                     and residuals is not None and np.all(residuals <= options.residual_tolerance))
    metadata = dict(backend="dmrg", tenpy=tenpy.__version__, options=asdict(options), sector=asdict(sector),
                    qp_cutoff=None, bath_modes=hamiltonian.bath_modes, impurity_modes=hamiltonian.nimp,
                    mode_order=list(prepared.order), groups=prepared.groups, coordinate_id=prepared.coordinate_id,
                    energy_reference=hamiltonian.metadata.get("energy_reference", "as specified in Hamiltonian"),
                    restriction=hamiltonian.metadata.get("coordinates", {}).get("restriction", "as specified"),
                    search_energy_shift=shift, operator_norm_bound=bound,
                    mpo_dimensions=None if physical_mpo is None else physical_mpo.chi,
                    mpo_original_dimensions=None if physical_mpo is None else physical_mpo.qdjj_original_dimensions,
                    search_shift_bound="min-max product-state Ritz bound" if physical_mpo is not None else None,
                    roots=diagnostics, gram=gram, gram_error=gram_error, roots_energy_ordered=ordered,
                    residual_method="untruncated-(H-E)-MPS-QR" if residuals is not None else None,
                    finite_problem_converged=converged, bath_convergence_checked=False)
    result = Eigenstates(energies, states, residuals, observables, times, metadata, hamiltonian, prepared)
    if options.require_convergence and not converged:
        raise RuntimeError(f"finite-problem convergence failed: residuals={residuals}, "
                           f"Gram error={gram_error:.3e}, ordered={ordered}")
    return result
