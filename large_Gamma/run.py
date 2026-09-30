"""Restartable, resource-bounded large-Gamma scan. Run from the source checkout.

Each point/bath/backend runs in a fresh process. Completed QP cutoffs and DMRG
bond/sector steps are published atomically. MPS continuation is used inside a
live worker; a restarted worker uses fresh seeds at the first unfinished step.
Only scalar measurements are written, including in the working directory.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import logging
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic, sleep

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
from threadpoolctl import threadpool_limits

from qdjj_solver import DiscreteBath, chain_expansion, cosh_grid, fit_surrogate, solve
from qdjj_solver.common.io import fingerprint, numerical_environment, source_provenance, write_json
from large_Gamma.study import (DEFAULT_INPUT, OBSERVABLES, SECTORS, empirical_resolution,
                               file_hash, make_hamiltonian, points, quadratic_continuum,
                               quadratic_finite, read_json, refinement_changes, resource_estimate,
                               sector_key, validate_config, values_from_sectors)

THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")


def peak_rss():
    try:
        import resource
        value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(value if sys.platform == "darwin" else 1024*value)
    except ImportError:
        return None


def live_rss(pid):
    if os.name != "posix":
        return None
    try:
        return 1024*int(subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(pid)], text=True, stderr=subprocess.DEVNULL).strip())
    except (ValueError, subprocess.CalledProcessError):
        return None


def prepare_archive(directory, config):
    directory = Path(directory)
    identity = dict(format="qdjj-large-Gamma", format_version=1, input=config,
                    implementation=source_provenance(),
                    scripts={name: file_hash(ROOT/"large_Gamma"/name) for name in ("run.py", "study.py")},
                    environment=numerical_environment())
    path = directory/"manifest.json"
    if path.exists():
        previous = read_json(path)
        if previous != identity:
            raise ValueError("resume requires identical inputs, numerical code, and environment; use a new run directory")
    else:
        write_json(path, identity)
    return fingerprint(identity)


def prepare_bath(directory, family, levels, config, bandwidth=None):
    bandwidth = config["model"]["bandwidth"] if bandwidth is None else bandwidth
    spec = {k: v for k, v in family.items() if k != "levels"}
    spec.update(levels=levels, delta=1., bandwidth=bandwidth)
    path = Path(directory)/"baths"/f"{fingerprint(spec)}.json"
    if path.exists():
        record = read_json(path)
        if record["spec"] != spec or record["bath_key"] != fingerprint(record["bath"]):
            raise ValueError("bath cache identity mismatch")
        return record
    kind = spec["kind"]
    options = {k: v for k, v in spec.items() if k not in ("kind", "levels")}
    if kind == "cosh-grid":
        bath = cosh_grid(levels//2, **options)
    elif kind == "surrogate":
        bath = fit_surrogate(levels, **options)
    else:
        bath = chain_expansion(levels, **options)
    bath_record = bath.record()
    bath_record["metadata"].pop("fit_seconds", None)
    record = dict(spec=spec, bath=bath_record, bath_key=fingerprint(bath_record))
    write_json(path, record)
    return record


def task_record(point, bath, kind, variant="standard"):
    return dict(point=point, bath_key=bath["bath_key"], family=bath["spec"]["kind"],
                levels=bath["spec"]["levels"], bandwidth=bath["spec"]["bandwidth"],
                kind=kind, variant=variant)


def result_path(directory, task):
    return Path(directory)/"tasks"/f"{fingerprint(task)}.json"


def load_task(path, identity=None, task=None):
    result = read_json(path)
    if (result["records_sha256"] != fingerprint(result["records"])
            or identity is not None and result["identity"] != identity
            or task is not None and result["task"] != task):
        raise ValueError("task identity or measurement checksum mismatch")
    return result


def state_record(state, sector, seconds, chi=None):
    obs = {k: complex(state.observables[k][0]) for k in OBSERVABLES}
    if max(abs(v.imag) for v in obs.values()) > 1e-9:
        raise ValueError("a Hermitian observable has a complex expectation")
    meta = state.metadata
    record = dict(sector=asdict(sector), energy=float(state.energies[0]),
                  observables={k: v.real for k, v in obs.items()},
                  residual=float(state.residuals[0]), seconds=seconds,
                  peak_rss_bytes=peak_rss(), hamiltonian_sha256=state.hamiltonian.fingerprint())
    if state.backend == "qp":
        record.update(projected_converged=True, dimension=meta["dimension"],
                      qp_weights=state.qp_weights[0].tolist(),
                      solver={k: meta[k] for k in ("method", "nnz", "estimated_solver_bytes", "real_arithmetic")})
    else:
        root = meta["roots"][0]
        record.update(finite_problem_converged=meta["finite_problem_converged"], chi_max=chi,
                      chi_actual=root["chi_actual"], sweeps=root["sweeps"],
                      sweep_energy_converged=root["sweep_energy_converged"],
                      trials=root["seed_trials"], gram_error=meta["gram_error"])
    return record


def execute_task(request, destination):
    config, task, identity = request["config"], request["task"], request["identity"]
    bath = DiscreteBath.from_record(request["bath"])
    if fingerprint(request["bath"]) != task["bath_key"]:
        raise ValueError("request refers to different bath coefficients")
    path = Path(destination)
    records = load_task(path, identity, task)["records"] if path.exists() else []
    start = monotonic()

    def publish(complete=False):
        write_json(path, dict(identity=identity, task=task, records=records,
                              records_sha256=fingerprint(records), complete=complete,
                              last_worker_seconds=monotonic()-start, peak_rss_bytes=peak_rss()))

    def existing(label):
        return next((r for r in records if r["label"] == label), None)

    def append(row):
        records.append(row)
        publish()

    logging.getLogger("tenpy").setLevel(logging.ERROR)
    h = make_hamiltonian(bath, task["point"], config,
                         compress=False if task["variant"] == "uncompressed" else None)
    options = dict(config["qp_options"] if task["kind"] == "qp" else config["dmrg_options"])
    if task["variant"] == "tight":
        if task["kind"] == "qp":
            options.update(tolerance=options["tolerance"]*.1, ncv=2*options["ncv"],
                           residual_tolerance=options["residual_tolerance"]*.1)
        else:
            options.update(max_sweeps=2*options["max_sweeps"],
                           energy_tolerance=options["energy_tolerance"]*.1,
                           seed=options["seed"]+1009, lanczos_maxiter=2*options["lanczos_maxiter"])
    if task["variant"] == "centered" and task["kind"] == "dmrg":
        middle = 2+2*((len(h.spins)-2)//4)
        options["mode_order"] = list(range(2, middle))+[0, 1]+list(range(middle, len(h.spins)))

    if task["point"]["u"] == 0 and task["kind"] == "qp" and existing("quadratic") is None:
        append(dict(label="quadratic", method="BdG", status="ok", accepted=True,
                    values=quadratic_finite(bath, task["point"]),
                    continuum=quadratic_continuum(task["point"], bath.bandwidth)))

    if task["kind"] == "qp":
        cutoffs = list(config["cutoffs"])
        if max(resource_estimate(h, None, s, options)["dimension"] for s in SECTORS) <= config["resources"]["ed_max_dimension"]:
            cutoffs.append(None)
        for cutoff in cutoffs:
            label = "ED" if cutoff is None else f"{cutoff}QP"
            if existing(label) is not None:
                continue
            estimates = [resource_estimate(h, cutoff, s, options) for s in SECTORS]
            if any(not estimate["allowed"] for estimate in estimates):
                append(dict(label=label, method=label, status="resource_limited", estimates=estimates))
                continue
            started = monotonic()
            try:
                states = []
                effective_cutoff = h.bath_modes if cutoff is None else min(cutoff, h.bath_modes)
                for sector in SECTORS:
                    tick = monotonic()
                    state = solve(h, cutoff=effective_cutoff, sector=sector, options=options)
                    states.append(state_record(state, sector, monotonic()-tick))
                    del state
                append(dict(label=label, method=label, status="ok", accepted=True,
                            effective_cutoff=effective_cutoff, full_space=effective_cutoff == h.bath_modes,
                            values=values_from_sectors(states), sectors=states,
                            options=options, seconds=monotonic()-started,
                            dark_threshold=float(min(h.metadata["decoupled_qp_energies"], default=np.inf))
                            if h.metadata["decoupled_qp_energies"] else None))
            except (ValueError, RuntimeError, MemoryError) as error:
                append(dict(label=label, method=label, status="failed", exception=type(error).__name__,
                            message=str(error), seconds=monotonic()-started))
    else:
        previous_states = {}
        history = [r for r in records if r.get("method") == "DMRG" and r["status"] == "ok"]
        for chi in config["chi_values"]:
            label = f"chi{chi}"
            old = existing(label)
            if old is not None:
                if old.get("bond_converged"):
                    break
                continue
            states = []
            started = monotonic()
            failed = False
            for sector in SECTORS:
                key = f"{label}_{sector_key(sector)}"
                old = existing(key)
                if old is not None and old["status"] == "ok":
                    states.append(old["state"])
                    continue
                tick = monotonic()
                try:
                    state = solve(h, backend="dmrg", sector=sector, options=dict(options, chi_max=chi),
                                  initial=previous_states.get(sector_key(sector)))
                    measured = state_record(state, sector, monotonic()-tick, chi)
                    previous_states[sector_key(sector)] = state
                    states.append(measured)
                    append(dict(label=key, status="ok", state=measured))
                except (ValueError, RuntimeError, MemoryError) as error:
                    append(dict(label=key, status="failed", exception=type(error).__name__, message=str(error)))
                    failed = True
                    break
            if failed:
                append(dict(label=label, method="DMRG", status="failed", chi=chi))
                break
            row = dict(label=label, method="DMRG", status="ok", chi=chi, options=options,
                       accepted=all(s["finite_problem_converged"] for s in states),
                       values=values_from_sectors(states), sectors=states, seconds=monotonic()-started,
                       dark_threshold=min(h.metadata["decoupled_qp_energies"], default=None))
            history.append(row)
            changes = refinement_changes(history, config["convergence"]["stable_steps"])
            row["bond_changes"] = changes
            row["bond_converged"] = bool(changes is not None
                and all(r["accepted"] for r in history[-config["convergence"]["stable_steps"]-1:])
                and max(changes.values())*config["convergence"]["safety_factor"] <= config["convergence"]["target"])
            append(row)
            if row["bond_converged"]:
                break
    publish(complete=True)


def execute_queue(jobs, directory, identity, config, *, workers, deadline=None):
    """Bound concurrency by reservations, then monitor actual RSS on POSIX."""
    directory = Path(directory)
    pending, running = list(jobs), []
    limits = config["resources"]
    env = dict(os.environ, **{name: "1" for name in THREAD_VARIABLES}, PYTHONDONTWRITEBYTECODE="1")
    completed = 0
    try:
        while pending or running:
            if deadline is not None and monotonic() >= deadline:
                pending.clear()
                if not running:
                    break
            while pending and len(running) < workers:
                task, bath, reservation = pending[0]
                if running and sum(r["reservation"] for r in running)+reservation > limits["total_memory_gib"]:
                    break
                pending.pop(0)
                path = result_path(directory, task)
                if path.exists() and load_task(path, identity, task)["complete"]:
                    continue
                request = directory/"requests"/path.name
                write_json(request, dict(identity=identity, task=task, bath=bath["bath"], config=config))
                log = directory/"logs"/f"{path.stem}.log"
                log.parent.mkdir(parents=True, exist_ok=True)
                handle = log.open("w", encoding="utf-8")
                process = subprocess.Popen([sys.executable, "-B", str(Path(__file__).resolve()),
                                            "--worker", str(request), "--result", str(path)],
                                           env=env, stdout=handle, stderr=subprocess.STDOUT)
                running.append(dict(process=process, path=path, handle=handle, reservation=reservation,
                                    started=monotonic(), task=task))
                print(f"Started {task['kind']} {task['point']} {task['family']} L={task['levels']} {task['variant']}", flush=True)
            rss = {r["process"].pid: live_rss(r["process"].pid) or 0 for r in running}
            total_rss = sum(rss.values())
            for item in list(running):
                process = item["process"]
                reason = None
                if rss[process.pid] > limits["worker_memory_gib"]*2**30:
                    reason = "worker RSS limit"
                elif total_rss > limits["total_memory_gib"]*2**30 and rss[process.pid] == max(rss.values()):
                    reason = "aggregate RSS limit"
                elif monotonic()-item["started"] > limits["task_timeout_hours"]*3600:
                    reason = "task time limit"
                if reason:
                    process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait()
                code = process.poll()
                if code is None:
                    continue
                item["handle"].close()
                path = item["path"]
                if reason or code != 0:
                    result = load_task(path, identity, item["task"]) if path.exists() else dict(
                        identity=identity, task=item["task"], records=[])
                    result.update(complete=True, interruption=reason or f"worker exit {code}",
                                  records_sha256=fingerprint(result["records"]))
                    write_json(path, result)
                running.remove(item)
                completed += 1
                print(f"Finished {completed}: {path.stem} ({reason or code}); queued={len(pending)}", flush=True)
            if running:
                sleep(limits["poll_seconds"])
    finally:
        for item in running:
            item["process"].terminate()
            try:
                item["process"].wait(timeout=10)
            except subprocess.TimeoutExpired:
                item["process"].kill()
                item["process"].wait()
            item["handle"].close()


def needs_refinement(directory, point, family, config, results=None):
    histories = {name: [] for name in [*(f"{q}QP" for q in config["cutoffs"]), "DMRG"]}
    if results is None:
        results = [load_task(path) for path in (Path(directory)/"tasks").glob("*.json")]
    for result in results:
        task = result["task"]
        if (task["point"] != point or task["family"] != family["kind"]
                or task["variant"] != "standard" or task["bandwidth"] != config["model"]["bandwidth"]):
            continue
        for name in histories:
            rows = [r for r in result["records"] if r.get("method") == name and r["status"] == "ok"]
            if rows:
                row = dict(rows[-1], levels=task["levels"])
                if name == "DMRG":
                    row["accepted"] = row.get("bond_converged", False)
                histories[name].append(row)
    for history in histories.values():
        history.sort(key=lambda row: row["levels"])
        resolution = empirical_resolution(history, config)
        if resolution is None or max(resolution.values()) > config["convergence"]["target"]:
            return True
    return False


def run(config, directory, *, profile="full", backend="all", workers=None,
        max_levels=None, wall_hours=None, dry_run=False, controls=False, selected_points=None):
    validate_config(config)
    selected = points(config) if selected_points is None else selected_points
    families = config["baths"][:1]
    if profile == "pilot":
        pilot = read_json(ROOT/"large_Gamma/input/pilot.json")
        selected = pilot["points"]
        max_levels = min(max_levels or max(pilot["levels"]), max(pilot["levels"]))
    if controls:
        selected = [p for p in selected if p["u"] in config["controls"]["u_values"]
                    and p["gamma"] in config["controls"]["gamma_values"]]
        families = config["baths"]
    directory = Path(directory)
    identity = None if dry_run else prepare_archive(directory, config)
    deadline = None if wall_hours is None else monotonic()+3600*wall_hours
    kinds = ("qp", "dmrg") if backend == "all" else (backend,)
    for family in families:
        for levels in family["levels"]:
            if max_levels is not None and levels > max_levels:
                continue
            if deadline is not None and monotonic() >= deadline:
                return
            if dry_run:
                # Cosh dimensions depend on mode/sector counts, not numerical quadrature coefficients.
                bath = cosh_grid((levels+1)//2, bandwidth=config["model"]["bandwidth"])
                h = make_hamiltonian(bath, selected[0], config)
                counts = {q: max(resource_estimate(h, q, s, config["qp_options"])["dimension"] for s in SECTORS)
                          for q in config["cutoffs"]}
                print(f"{family['kind']} L={levels}: {len(selected)} points, max eta-sector QP dimensions {counts}")
                continue
            try:
                bath = prepare_bath(directory, family, levels, config)
            except (ValueError, RuntimeError) as error:
                write_json(directory/"bath_failures"/f"{family['kind']}_{levels}.json",
                           dict(family=family, levels=levels, error=str(error)))
                continue
            h = make_hamiltonian(DiscreteBath.from_record(bath["bath"]), selected[0], config)
            estimates = [resource_estimate(h, q, s, config["qp_options"])
                         for q in config["cutoffs"] for s in SECTORS]
            qp_reservation = min(config["resources"]["worker_memory_gib"],
                                 max([1.]+[2+2*r["workspace_bytes"]/2**30 for r in estimates if r["allowed"]]))
            jobs = []
            previous = [load_task(path) for path in (directory/"tasks").glob("*.json")]
            for point in selected:
                if levels != family["levels"][0] and not needs_refinement(directory, point, family, config, previous):
                    continue
                for kind in kinds:
                    task = task_record(point, bath, kind)
                    reservation = qp_reservation if kind == "qp" else config["resources"]["dmrg_reservation_gib"]
                    jobs.append((task, bath, reservation))
            execute_queue(jobs, directory, identity, config, workers=workers or config["resources"]["jobs"], deadline=deadline)
            from large_Gamma.analyze import export
            export(directory, directory/"summary", plots=False)
    if controls and not dry_run:
        run_controls(config, directory, selected, kinds, identity,
                     workers or config["resources"]["jobs"], deadline, max_levels)


def run_controls(config, directory, selected, kinds, identity, workers, deadline, max_levels):
    """Independent representation, optimization, derivative and bandwidth checks."""
    settings = config["controls"]
    primary = config["baths"][0]
    jobs = []
    # Small uncompressed systems are tractable in electron space and in full ED.
    small = prepare_bath(directory, primary, settings["uncompressed_levels"], config)
    for point in selected:
        for kind in kinds:
            for variant in ("standard", "uncompressed"):
                jobs.append((task_record(point, small, kind, variant), small, 1.))
    check_level = min(settings["optimization_levels"], max_levels or settings["optimization_levels"])
    check = prepare_bath(directory, primary, check_level, config)
    for point in selected:
        if point["gamma"] not in settings["optimization_gammas"]:
            continue
        for kind in kinds:
            variants = ["standard", "tight"]
            if kind == "dmrg":
                variants.append("centered")
            for variant in variants:
                jobs.append((task_record(point, check, kind, variant), check,
                             config["resources"]["dmrg_reservation_gib"] if kind == "dmrg" else 2.))
            for sign, label in ((-1, "phase_minus"), (1, "phase_plus")):
                shifted = dict(point, phi=point["phi"]+sign*settings["phase_step"])
                task = dict(task_record(shifted, check, kind, label), anchor=point)
                jobs.append((task, check, config["resources"]["dmrg_reservation_gib"] if kind == "dmrg" else 2.))
    execute_queue(jobs, directory, identity, config, workers=workers, deadline=deadline)
    for bandwidth in settings["bandwidths"]:
        for levels in settings["bandwidth_levels"]:
            if max_levels is not None and levels > max_levels:
                continue
            if deadline is not None and monotonic() >= deadline:
                return
            bath = prepare_bath(directory, primary, levels, config, bandwidth)
            h = make_hamiltonian(DiscreteBath.from_record(bath["bath"]), selected[0], config)
            estimates = [resource_estimate(h, q, s, config["qp_options"])
                         for q in config["cutoffs"] for s in SECTORS]
            reservation = min(config["resources"]["worker_memory_gib"],
                              max([1.]+[2+2*r["workspace_bytes"]/2**30 for r in estimates if r["allowed"]]))
            jobs = [(task_record(point, bath, kind), bath,
                     config["resources"]["dmrg_reservation_gib"] if kind == "dmrg" else reservation)
                    for point in selected if point["gamma"] == max(settings["gamma_values"]) for kind in kinds]
            execute_queue(jobs, directory, identity, config, workers=workers, deadline=deadline)
    from large_Gamma.analyze import export
    export(directory, Path(directory)/"summary", plots=False)


def refine(config, directory, *, backend, workers, max_levels, wall_hours):
    from large_Gamma.analyze import collect, compare, refinement_points
    deadline = None if wall_hours is None else monotonic()+3600*wall_hours
    for _ in range(config["refinement"]["max_rounds"]):
        _, _, rows, _, _ = collect(directory)
        comparisons, _ = compare(rows, config)
        selected = refinement_points(comparisons, config)
        existing_points = {fingerprint(row["point"]) for row in rows}
        selected = [point for point in selected if fingerprint(point) not in existing_points]
        if not selected or deadline is not None and monotonic() >= deadline:
            break
        run(config, directory, backend=backend, workers=workers, max_levels=max_levels,
            wall_hours=None if deadline is None else (deadline-monotonic())/3600,
            selected_points=selected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=ROOT/"runs/large_Gamma/study")
    parser.add_argument("--profile", choices=("pilot", "full"), default="full")
    parser.add_argument("--backend", choices=("all", "qp", "dmrg"), default="all")
    parser.add_argument("--jobs", type=int)
    parser.add_argument("--max-levels", type=int)
    parser.add_argument("--wall-hours", type=float)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--controls", action="store_true")
    parser.add_argument("--refine", action="store_true", help="add adaptive Gamma points after the base scan")
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--result", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.wall_hours is not None and args.wall_hours <= 0:
        parser.error("--wall-hours must be positive")
    with threadpool_limits(limits=1):
        if args.worker:
            execute_task(read_json(args.worker), args.result)
        else:
            run(read_json(args.input), args.output, profile=args.profile, backend=args.backend,
                workers=args.jobs, max_levels=args.max_levels, wall_hours=args.wall_hours,
                dry_run=args.dry_run, controls=args.controls)
            if args.refine and not args.dry_run:
                refine(read_json(args.input), args.output, backend=args.backend, workers=args.jobs,
                       max_levels=args.max_levels, wall_hours=args.wall_hours)


if __name__ == "__main__":
    main()
