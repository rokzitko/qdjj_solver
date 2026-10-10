"""Checkpoint one test4 sector prefix, without assigning global excitations.

The parent owns dependency resolution, aggregation, and the 4-CPU/32-GiB/time
guard. This helper uses one numerical thread and never launches other workers.
Only a returned solve can be checkpointed; an interrupted solve has no partial
Eigenstates object in the current core API.
"""

from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import logging
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from qdjj_solver import DiscreteBath, Sector, load_result, reference_model, save_result, solve
from qdjj_solver.common import number
from qdjj_solver.common.io import complex_array, json_value
from qdjj_solver.common.problem import Hamiltonian

from NRG_comparisons.test4.physics import (
    TOLERANCE, detuning, finish_states, paper_offset, state_values, validate_physical,
)


def build_model(physical, bath):
    """The frozen finite.py Hamiltonian and observables, with unchanged weights.

    ``bath`` is an explicit record or DiscreteBath, reused in both physical leads.
    Gamma is total hybridization; the impurity and both BCS references stay centered.
    """
    physical = validate_physical(physical)
    bath = bath if isinstance(bath, DiscreteBath) else DiscreteBath.from_record(bath)
    if bath.delta != physical["gap"] or bath.bandwidth != physical["bandwidth"]:
        raise ValueError("bath gap/bandwidth differ from physical parameters")
    h = reference_model(bath, **{name: physical[name] for name in
                                ("geometry", "u", "gamma", "field", "phi")},
                        detuning=detuning(physical), rho_wn=physical["wn"]/(2*physical["bandwidth"]),
                        rho_ws=physical["ws"]/(2*physical["bandwidth"]), symmetry=False, compress=False)
    h.observables.pop("rho_ws_derivative", None)
    dot = (number(0)-number(1))/2
    bath_z = sum((number(i)-number(i+1))/2 for i in range(2, len(h.spins), 2))
    for name, operator in dict(dot_z2=dot*dot, bathSz=bath_z, bathSz2=bath_z*bath_z,
                               dotSz_bathSz=dot*bath_z, totSz=dot+bath_z).items():
        h.observables[name] = h.physical_operator(operator)
    return h


def rows_from_result(result, physical):
    """Decode a trusted LIVE or LOADED DMRG result into the 14 shared values.

    Physical parameters are checked against the entire Hamiltonian, including
    observables and both baths. Nondegeneracy means only no tie in this prefix;
    a single root does not establish uniqueness or lowest-root completeness.
    """
    physical = validate_physical(physical)
    h, metadata = result.hamiltonian, result.metadata
    if build_model(physical, h.metadata["bath"][0]).fingerprint() != h.fingerprint():
        raise ValueError("result Hamiltonian differs from test4 physical parameters/observables")
    if result.backend != "dmrg" or metadata.get("backend") != "dmrg":
        raise ValueError("rows require a DMRG result")
    spin = metadata["sector"]["twice_sz"]
    if type(spin) is not int or metadata["sector"] != asdict(Sector(spin % 2, spin)):
        raise ValueError("result must have a fixed twice-Sz sector")
    energies = np.asarray(result.energies)
    count = len(energies)
    if (energies.shape != (count,) or count not in (1, 2)
            or len(metadata["roots"]) != count or metadata["options"]["eigenpairs"] != count):
        raise ValueError("result must contain a complete one- or two-root prefix")
    gram = complex_array(metadata["gram"])
    if gram.shape != (count, count):
        raise ValueError("Gram shape differs from the root prefix")
    gap = physical["gap"]
    gram_ok = bool(np.linalg.norm(gram-np.eye(count)) <= 1e-8)
    ordered = bool(np.all(np.diff(energies) >= -1e-12*gap))
    states = []
    for root, energy in enumerate(energies):
        obs = {name: values[root] for name, values in result.observables.items()}
        if not abs(obs["totSz"]-spin/2) <= TOLERANCE:
            raise ValueError("direct total Sz does not match the requested sector")
        values = state_values(obs["impurity_charge"], obs["double_occupancy_0"],
                              obs["impurity_spin_z"], obs["dot_z2"], spin/2,
                              m_bath=obs["bathSz"], Czz=obs["dotSz_bathSz"], bath_z2=obs["bathSz2"])
        values.update(energy=energy/gap, energy_paper=(energy+paper_offset(physical))/gap,
                      current=obs["phase_derivative"]/gap)
        residual = None if result.residuals is None else float(result.residuals[root]/gap)
        info = metadata["roots"][root]
        norm_error = info.get("norm_error", math.inf)
        converged = bool(residual is not None and math.isfinite(residual) and 0 <= residual <= 1e-7
                         and info["sweep_energy_converged"] and gram_ok and ordered
                         and 0 <= info.get("lower_root_overlap", math.inf) <= 1e-8
                         and 0 <= norm_error <= 1e-8)
        states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                           residual=residual, finite_converged=converged,
                           degenerate=bool(np.count_nonzero(abs(energies-energy)/gap <= TOLERANCE) > 1)))
    # Reuse all shared scalar/identity checks, but discard single-sector excitation labels.
    finish_states(states)
    for row in states:
        del row["excitation"], row["subgap"]
    return states


def _sha256(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"checkpoint requires an ordinary file: {path}")
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def run_state(task, physical, bath_record, directory, initial=None) -> dict:
    """Solve and immediately save an intact sector prefix to eigenstates.json.

    Required task fields: backend='dmrg', twice_sz (int), roots (1 or 2), chi
    (32/64/128/256). Defaults: seed_trials=2, min_sweeps=8, max_sweeps=40,
    seed=1729, layout='centered'. Other parent bookkeeping fields are not solver
    options. The existing directory may contain only request.json, launch.json,
    worker.log. Checkpoint/result outputs are never overwritten.

    ``initial`` is None or {'path': absolute JSON path, 'sha256': byte hash} of
    a trusted checkpoint from this workflow with a compatible root-0 prefix.
    Two roots always mean genuine sequential optimization, including root 0.
    Status, request hash, resource monitoring and global excitations belong to
    the parent, not the immutable checkpoint or this per-sector report.
    """
    start = perf_counter()
    if not isinstance(task, dict) or task.get("backend") != "dmrg":
        raise ValueError("backend must be dmrg")
    spin, roots, chi = (task.get(name) for name in ("twice_sz", "roots", "chi"))
    if type(spin) is not int:
        raise ValueError("twice_sz must be an integer")
    if type(roots) is not int or roots not in (1, 2):
        raise ValueError("roots must be 1 or 2 (a prefix starting at root 0)")
    if type(chi) is not int or chi not in (32, 64, 128, 256):
        raise ValueError("chi must be 32, 64, 128 or 256")
    trials, seed = task.get("seed_trials", 2), task.get("seed", 1729)
    if type(trials) is not int or trials < 1:
        raise ValueError("seed_trials must be a positive integer")
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    minimum, maximum = task.get("min_sweeps", 8), task.get("max_sweeps", 40)
    if (type(minimum) is not int or type(maximum) is not int
            or not 8 <= minimum <= maximum <= 40):
        raise ValueError("require 8 <= min_sweeps <= max_sweeps <= 40")
    if task.get("layout", "centered") != "centered":
        raise ValueError("layout must be centered (L QP, impurity pair, R QP)")
    if task.get("cutoff") is not None or any(name in task for name in
            ("mode_order", "group_size", "threads", "orthogonal_to", "eigenpairs", "chi_max")):
        raise ValueError("unsupported solver override; use the declared task fields")
    physical = validate_physical(physical)
    bath_record = deepcopy(json_value(bath_record))
    directory = Path(directory)
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("directory must already exist and be an ordinary directory")
    if any(path.name not in {"request.json", "launch.json", "worker.log"}
           or path.is_symlink() or not path.is_file() for path in directory.iterdir()):
        raise ValueError("existing checkpoint or unrelated results; use a fresh directory")

    with threadpool_limits(limits=1):
        h = build_model(physical, bath_record)
        levels = len(bath_record["xi"])
        orbitals = 1+2*levels
        dimension = math.comb(2*orbitals, orbitals+spin) if abs(spin) <= orbitals else 0
        if roots > dimension:
            raise ValueError("requested root prefix exceeds the sector dimension")
        sector = Sector(spin % 2, spin)
        # Only an MPS permutation: the bath remains the original two-lead QP star.
        order = list(range(2, 2+2*levels))+[0, 1]+list(range(2+2*levels, len(h.spins)))
        options = dict(eigenpairs=roots, chi_max=chi, min_sweeps=minimum, max_sweeps=maximum,
                       energy_tolerance=1e-12*physical["gap"], residual_tolerance=1e-7*physical["gap"],
                       lanczos_probability_tolerance=1e-22, seed_trials=trials, seed=seed,
                       calculate_residuals=True, require_convergence=False, threads=1,
                       group_size=2, mode_order=order)
        h_hash = h.fingerprint()
        workflow = dict(format_version=1, physical=physical, bath=bath_record, layout="centered",
                        root_prefix=list(range(roots)), task=deepcopy(json_value(task)),
                        initial=deepcopy(initial))
        previous = None
        load_start = perf_counter()
        if initial is not None:
            if (not isinstance(initial, dict) or set(initial) != {"path", "sha256"}
                    or not isinstance(initial["path"], str) or not Path(initial["path"]).is_absolute()
                    or not isinstance(initial["sha256"], str) or len(initial["sha256"]) != 64
                    or any(c not in "0123456789abcdef" for c in initial["sha256"])):
                raise ValueError("initial requires an absolute checkpoint path and SHA256 byte hash")
            path = Path(initial["path"])
            if _sha256(path) != initial["sha256"]:
                raise ValueError("initial checkpoint JSON checksum mismatch")
            record = json.loads(path.read_text(encoding="utf-8"))
            if ((record.get("format"), record.get("format_version"), record.get("backend"))
                    != ("qdjj-eigenstates", 1, "dmrg")):
                raise ValueError("initial must be a DMRG checkpoint")
            if (record["hamiltonian_sha256"] != h_hash
                    or Hamiltonian.from_record(record["hamiltonian"]).fingerprint() != h_hash):
                raise ValueError("initial full Hamiltonian fingerprint differs (physical/bath/observables)")
            metadata = record["calculation"]
            old_options = metadata["options"]
            count = old_options["eigenpairs"]
            old_workflow = metadata.get("test4_dmrg_states", {})
            if (type(count) is not int or count not in (1, 2) or count > roots
                    or len(record["energies"]) != count or len(metadata["roots"]) != count
                    or old_workflow.get("format_version") != 1
                    or old_workflow.get("root_prefix") != list(range(count))):
                raise ValueError("initial must retain its complete root-0 prefix; no slicing/relabeling")
            if (old_workflow.get("physical") != physical or old_workflow.get("bath") != bath_record
                    or old_workflow.get("layout") != "centered"):
                raise ValueError("initial frozen physical/bath/layout metadata differ")
            if metadata.get("backend") != "dmrg" or metadata["sector"] != asdict(sector):
                raise ValueError("initial sector differs")
            from qdjj_solver.dmrg_solver.mpo import prepare

            prepared = prepare(h, sector, group_size=2, mode_order=order)
            if (metadata.get("qp_cutoff") is not None or old_options["group_size"] != 2
                    or old_options["mode_order"] != order or metadata["mode_order"] != order
                    or metadata["groups"] != json_value(prepared.groups)
                    or metadata["coordinate_id"] != prepared.coordinate_id):
                raise ValueError("initial mode order/grouping/coordinate layout differs")
            # The core verifies the companion hash and its binding to this entire JSON record.
            previous = load_result(path)
            if (_sha256(path) != initial["sha256"] or previous.hamiltonian.fingerprint() != h_hash
                    or previous.prepared.sector != sector
                    or previous.prepared.coordinate_id != prepared.coordinate_id
                    or len(previous.states) != count):
                raise ValueError("initial checkpoint changed or loaded inconsistent states")
        load_seconds = perf_counter()-load_start

        logger = logging.getLogger("qdjj_solver.dmrg_solver.solver")
        old_level = logger.level
        handler = None
        if not logger.hasHandlers():
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s"))
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        try:
            result = solve(h, backend="dmrg", cutoff=None, sector=sector, options=options, initial=previous)
            result.metadata["test4_dmrg_states"] = workflow
            checkpoint_path = directory/"eigenstates.json"
            save_start = perf_counter()
            save_result(result, checkpoint_path, save_states=True)
            save_seconds = perf_counter()-save_start
        finally:
            logger.setLevel(old_level)
            if handler is not None:
                logger.removeHandler(handler)
                handler.close()

        # Nothing below may rewrite the checkpoint, including failed row conversion.
        record = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        artifact = record["state_artifact"]
        artifact_hash = _sha256(directory/artifact["file"])
        if artifact_hash != artifact["sha256"]:
            raise ValueError("saved MPS checkpoint checksum mismatch")
        checkpoint = dict(path=checkpoint_path.name, sha256=_sha256(checkpoint_path),
                          state_artifact=dict(file=artifact["file"], sha256=artifact_hash),
                          hamiltonian_sha256=h_hash, coordinate_id=result.metadata["coordinate_id"],
                          twice_sz=spin, roots=roots, chi=chi)
        states = rows_from_result(result, physical)
        selection = []
        for root, info in enumerate(result.metadata["roots"]):
            trial_data = info["seed_trials"]
            eligible = [trial["energy"]/physical["gap"] for trial in trial_data
                        if trial["lower_root_overlap"] <= 1e-8]
            spread = float(np.ptp(eligible)) if eligible else None
            selection.append(dict(root=root, selected_seed=json_value(info["seed"]),
                                  trial_count=len(trial_data), eligible_trials=len(eligible),
                                  eligible_energy_spread_over_gap=spread,
                                  independent_trial_agreement=len(eligible) >= 2 and spread <= TOLERANCE))
        diagnostics = dict(
            backend="dmrg", physical=physical, bath=bath_record, task=deepcopy(json_value(task)), chi=chi,
            layout="centered", layout_description="L QP modes, impurity pair, R QP modes; star, not chain",
            hamiltonian_sha256=h_hash, coordinate_id=result.metadata["coordinate_id"],
            sector=asdict(sector), options=json_value(result.metadata["options"]),
            mode_order=order, group_size=2, dimension=dimension,
            requested_roots=roots, returned_roots=len(states),
            coverage_scope="single-sector root-0 prefix; not global coverage",
            qualification_scope="all returned roots; finite residual/sweep/Gram/order/norm checks and no observed ties",
            degeneracy_scope="computed roots only; hidden degeneracies not excluded",
            spectral_completeness_established=False, lowest_root_completeness_established=False,
            bath_convergence_checked=False, chi_convergence_checked=False,
            full_qp_comparison_checked=False,
            bath_levels_per_lead=levels, bath_leads=2, bath_modes_total=4*levels,
            identical_bath_record_for_both_leads=True,
            bath_weight_sum_per_lead=float(np.sum(bath_record["weights"])), weights_renormalized=False,
            energy_units="gap", energy_reference="centered impurity; both isolated BCS vacua subtracted",
            paper_energy_offset=paper_offset(physical)/physical["gap"], detuning=detuning(physical),
            current_units="I/(2e*Delta/hbar)", current_is_real=True,
            residual_scope="untruncated_finite_hamiltonian", residual_tolerance_over_gap=1e-7,
            norm_error=[info.get("norm_error") for info in result.metadata["roots"]],
            norm_error_tolerance=1e-8, root_selection_seed_independence=selection,
            trial_data_scope="all core trial summaries; full sweep statistics only for selected candidates",
            solver_metadata=json_value(result.metadata), timings=json_value(result.timings),
            checkpoint_seconds=save_seconds, initial_load_seconds=load_seconds,
            helper_seconds=perf_counter()-start,
            all_returned_states_converged=all(row["finite_converged"] for row in states))
    return dict(states=states, checkpoint=checkpoint, diagnostics=diagnostics,
                qualified=all(row["finite_converged"] and not row["degenerate"] for row in states))
