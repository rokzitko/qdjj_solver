"""Version-1 finite worker for the approved test4 continuum campaign.

The parent owns the 72-hour campaign, 8-hour task, 64-GiB process watchdog,
4-CPU allocation, serial scheduling, exceptions and final receipts. This worker
uses one numerical thread. Finite residuals never certify bath convergence.
No legacy limits, Hamiltonians or workers are modified.
"""

from copy import deepcopy
from dataclasses import asdict
import json
import logging
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from qdjj_solver import Sector, load_result, save_result, solve
from qdjj_solver.common.io import json_value
from qdjj_solver.common.problem import Hamiltonian

from NRG_comparisons.test4.dmrg_states import _sha256, build_model, rows_from_result
from NRG_comparisons.test4.physics import (
    PRIMARY, TOLERANCE, detuning, finish_states, paper_offset, state_values, validate_physical,
)


FORMAT_VERSION = 1


def run_task(task, physical, bath_record, directory, initial=None) -> dict:
    """Return the frozen finite (QP) or dmrg_states (DMRG) report schema.

    QP requires method='sparse', cutoff=4/6/8/10 and the five-sector inventory
    0:3, -1:2, +1:2, -2:1, +2:1, plus one lookahead per sector. Optional limits
    may only lower the 10M-state/48-GiB/500M-nnz caps. No QP vectors are saved.

    DMRG requires twice_sz=0/-1/+1, roots=1/2, chi=128/256/512/1024. Defaults
    are layout='centered', min_sweeps=8, max_sweeps=40, seed_trials=2, seed=1729;
    layout, trials and seed are fixed. Save the entire returned prefix BEFORE
    conversion, even when unconverged. No single-sector global reference is set.

    The ordinary directory must exist and contain only ordinary request.json,
    launch.json and worker.log files. Outputs are never reused or overwritten.
    DMRG initial is None or {'path': absolute JSON path, 'sha256': byte hash};
    only intact same-Hamiltonian/bath/sector/layout prefixes may be continued.
    Parent bookkeeping fields are retained, not forwarded as solver options.
    """
    start = perf_counter()
    if not isinstance(task, dict) or task.get("backend") not in ("qp", "dmrg"):
        raise ValueError("backend must be qp or dmrg")
    backend = task["backend"]
    seed = task.get("seed", 1729)
    if type(seed) is not int or seed != 1729:
        raise ValueError("seed must be 1729")
    if any(name in task for name in (
            "options", "solver", "sector", "mode_order", "group_size", "threads", "orthogonal_to",
            "eigenpairs", "chi_max", "initial", "tolerance", "residual_tolerance", "ncv", "maxiter")):
        raise ValueError("unsupported solver override; use the declared task fields")
    if backend == "qp" and initial is not None:
        raise ValueError("QP does not accept initial checkpoints")
    if backend == "dmrg":
        spin, roots, chi = (task.get(name) for name in ("twice_sz", "roots", "chi"))
        if type(spin) is not int or spin not in (0, -1, 1):
            raise ValueError("twice_sz must be 0, -1 or +1")
        if type(roots) is not int or roots not in (1, 2):
            raise ValueError("roots must be 1 or 2 (a prefix starting at root 0)")
        if type(chi) is not int or chi not in (128, 256, 512, 1024):
            raise ValueError("chi must be 128, 256, 512 or 1024")
        trials = task.get("seed_trials", 2)
        if type(trials) is not int or trials != 2:
            raise ValueError("seed_trials must be 2")
        minimum, maximum = task.get("min_sweeps", 8), task.get("max_sweeps", 40)
        if (type(minimum) is not int or type(maximum) is not int
                or not 8 <= minimum <= maximum <= 40):
            raise ValueError("require 8 <= min_sweeps <= max_sweeps <= 40")
        if task.get("layout", "centered") != "centered":
            raise ValueError("layout must be centered (L QP, impurity pair, R QP)")
        if task.get("cutoff") is not None:
            raise ValueError("DMRG requires cutoff=None")
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
        if backend == "qp":
            return _run_qp(task, physical, bath_record, directory, h)
        levels = len(bath_record["xi"])
        orbitals = 1+2*levels
        dimension = math.comb(2*orbitals, orbitals+spin)
        if roots > dimension:
            raise ValueError("requested root prefix exceeds the sector dimension")
        sector = Sector(spin % 2, spin)
        order = list(range(2, 2+2*levels))+[0, 1]+list(range(2+2*levels, len(h.spins)))
        options = dict(eigenpairs=roots, chi_max=chi, min_sweeps=minimum, max_sweeps=maximum,
                       energy_tolerance=1e-12*physical["gap"], residual_tolerance=1e-7*physical["gap"],
                       lanczos_probability_tolerance=1e-22, seed_trials=trials, seed=seed,
                       calculate_residuals=True, require_convergence=False, threads=1,
                       group_size=2, mode_order=order)
        h_hash = h.fingerprint()
        workflow = dict(format_version=1, continuum_finite_version=FORMAT_VERSION,
                        physical=physical, bath=bath_record, layout="centered",
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
            artifact = record["state_artifact"]
            filename = artifact["file"]
            if not isinstance(filename, str) or not filename or Path(filename).name != filename:
                raise ValueError("MPS checkpoint filename must be local to the result directory")
            companion = path.parent/filename
            if _sha256(companion) != artifact["sha256"]:
                raise ValueError("MPS checkpoint checksum mismatch")
            # Core loading also verifies the MPS binding to the ENTIRE scalar JSON.
            previous = load_result(path)
            if (_sha256(path) != initial["sha256"] or _sha256(companion) != artifact["sha256"]
                    or previous.hamiltonian.fingerprint() != h_hash
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

        # Conversion failures must leave the immutable checkpoint available to the parent.
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
            backend="dmrg", continuum_finite_version=FORMAT_VERSION,
            physical=physical, bath=bath_record, task=deepcopy(json_value(task)), chi=chi,
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


def _run_qp(task, physical, bath_record, directory, h):
    from qdjj_solver.qp_solver.basis import dimension

    cutoff = task.get("cutoff")
    if type(cutoff) is not int or cutoff not in (4, 6, 8, 10) or cutoff > h.bath_modes:
        raise ValueError("QP cutoff must be 4, 6, 8 or 10, not exceeding the bath mode count")
    if task.get("method") != "sparse":
        raise ValueError("QP requires explicit method='sparse'")
    spins = task.get("spins", [0, -1, 1, -2, 2])
    if (not isinstance(spins, list) or any(type(spin) is not int for spin in spins)
            or spins != [0, -1, 1, -2, 2]):
        raise ValueError("QP requires spins=[0, -1, 1, -2, 2]")
    roots = task.get("roots", 2)
    by_spin = task.get("roots_by_spin", {"0": 3, "-2": 1, "2": 1})
    if (type(roots) is not int or roots < 1 or not isinstance(by_spin, dict)
            or any(key not in {str(s) for s in spins} for key in by_spin)
            or any(type(value) is not int or value < 1 for value in by_spin.values())):
        raise ValueError("roots_by_spin requires selected spin string keys and positive integer roots")
    counts = {spin: by_spin.get(str(spin), roots) for spin in spins}
    if counts != {0: 3, -1: 2, 1: 2, -2: 1, 2: 1}:
        raise ValueError("QP requires the complete root inventory 0:3, -1:2, +1:2, -2:1, +2:1")
    limits = dict(max_dimension=task.get("max_dimension", 10_000_000),
                  max_memory_gib=task.get("max_memory_gib", 48.), max_nnz=task.get("max_nnz", 500_000_000))
    if (type(limits["max_dimension"]) is not int or not 0 < limits["max_dimension"] <= 10_000_000
            or type(limits["max_nnz"]) is not int or not 0 < limits["max_nnz"] <= 500_000_000
            or type(limits["max_memory_gib"]) not in (int, float) or not 0 < limits["max_memory_gib"] <= 48.):
        raise ValueError("QP resource gate: at most 10M states, 48 GiB, 500M nnz")
    dimensions = [dimension(h.nimp, h.spins, cutoff, Sector(spin % 2, spin)) for spin in spins]
    if max(dimensions) > limits["max_dimension"]:
        raise MemoryError("QP sector dimension exceeds the explicit resource gate")
    if any(dim <= counts[spin]+2 for spin, dim in zip(spins, dimensions, strict=True)):
        raise ValueError("QP requires room for the full inventory and sparse lookahead; no dense fallback")
    gap = physical["gap"]
    states, sectors = [], []
    diagnostics = dict(
        backend="qp", continuum_finite_version=FORMAT_VERSION, sectors=sectors, requested_spins=list(spins),
        bath_levels_per_lead=len(bath_record["xi"]), bath_leads=2, bath_modes_total=h.bath_modes,
        identical_bath_record_for_both_leads=True,
        bath_weight_sum_per_lead=float(np.sum(bath_record["weights"])), weights_renormalized=False,
        bath_convergence_checked=False, spectral_completeness_established=False,
        energy_units="gap", energy_reference="centered impurity; both isolated BCS vacua subtracted",
        paper_energy_offset=paper_offset(physical)/gap, detuning=detuning(physical),
        current_units="I/(2e*Delta/hbar)", current_is_real=True,
        residual_scope="projected_cutoff_hamiltonian", qp_weights_scope="all solved roots, including lookahead")
    for spin, dim in zip(spins, dimensions, strict=True):
        count = counts[spin]
        options = dict(method="sparse", eigenpairs=count+1, residual_tolerance=1e-8*gap,
                       tolerance=1e-12, ncv=40, threads=1, seed=1729, **limits)
        result = solve(h, backend="qp", cutoff=cutoff, sector=Sector(spin % 2, spin), options=options)
        if len(result.energies) != count+1:
            raise ValueError("QP solve did not return the complete requested prefix and lookahead")
        result.metadata["test4_continuum_finite"] = dict(
            format_version=FORMAT_VERSION, task=deepcopy(json_value(task)), physical=physical, bath=bath_record,
            requested_roots=count, root_prefix=list(range(count+1)),
            coverage_scope="single-sector scalar receipt; not a complete task or global reference")
        receipt_path = directory/f"qp-sz{spin:+d}.json"
        save_start = perf_counter()
        save_result(result, receipt_path, save_states=False)
        save_seconds = perf_counter()-save_start
        info = dict(twice_sz=spin, parity=spin % 2, dimension=dim, requested_roots=count,
                    returned_roots=count, solved_roots=count+1, solver_metadata=json_value(result.metadata),
                    timings=json_value(result.timings), complete_sector_diagonalized=False,
                    degeneracy_scope="computed roots only; hidden degeneracies not excluded",
                    qp_weights=json_value(result.qp_weights), scalar_receipt=dict(
                        path=receipt_path.name, sha256=_sha256(receipt_path)), checkpoint_seconds=save_seconds)
        sectors.append(info)
        for root in range(count):
            obs = {name: values[root] for name, values in result.observables.items()}
            if not abs(obs["totSz"]-spin/2) <= TOLERANCE:
                raise ValueError("direct total Sz does not match the requested sector")
            values = state_values(obs["impurity_charge"], obs["double_occupancy_0"],
                                  obs["impurity_spin_z"], obs["dot_z2"], spin/2,
                                  m_bath=obs["bathSz"], Czz=obs["dotSz_bathSz"], bath_z2=obs["bathSz2"])
            energy = float(result.energies[root]/gap)
            values.update(energy=energy, energy_paper=energy+paper_offset(physical)/gap,
                          current=obs["phase_derivative"]/gap)
            residual = None if result.residuals is None else float(result.residuals[root]/gap)
            converged = residual is not None and math.isfinite(residual) and 0 <= residual <= 1e-8
            states.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                               residual=residual, finite_converged=bool(converged),
                               degenerate=bool(np.count_nonzero(abs(result.energies/gap-energy) <= TOLERANCE) > 1)))
        # Do not retain giant vectors while allocating the next sector's solve.
        del result
    coverage = finish_states(states)
    primary = {row["id"]: row["finite_converged"] and not row["degenerate"]
               for row in states if row["id"] in PRIMARY.values()}
    diagnostics["primary_qualified"] = primary
    diagnostics["primary_unique"] = {row["id"]: not row["degenerate"] for row in states
                                     if row["id"] in PRIMARY.values()}
    diagnostics["primary_residual_qualified"] = {
        row["id"]: row["residual"] is not None and row["residual"] <= 1e-8
        for row in states if row["id"] in PRIMARY.values()}
    diagnostics["all_returned_states_converged"] = all(row["finite_converged"] for row in states)
    diagnostics["qualification_scope"] = "three primary root-0 branches; extra roots reported separately"
    return dict(states=states, coverage=coverage, diagnostics=diagnostics,
                qualified=set(primary) == set(PRIMARY.values()) and all(primary.values()))
