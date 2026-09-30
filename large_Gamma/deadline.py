"""Hard-deadline continuation with exact-reference-validated DMRG observables.

Completed numerical evidence is imported with its original provenance. The
revised accuracy assessment is separate from the original solver diagnostics.
Every child is in a process group and is terminated at the computation cutoff;
publication and report generation have their own bounded time allowance.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from copy import deepcopy
from datetime import datetime, timezone
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic, sleep, time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from threadpoolctl import threadpool_limits

from qdjj_solver import DiscreteBath, Sector, solve
from qdjj_solver.common.io import fingerprint, numerical_environment, source_provenance, write_json
from large_Gamma import analyze, run, study

DEFAULT_PROFILE = ROOT/"large_Gamma/input/deadline.json"


def timestamp(value):
    parsed = datetime.fromisoformat(value)
    if parsed.tzinfo is None:
        raise ValueError("deadline timestamps must include a timezone")
    return parsed.timestamp()


def configuration(profile):
    if profile["format_version"] != 1 or timestamp(profile["compute_stop"]) >= timestamp(profile["report_deadline"]):
        raise ValueError("invalid deadline profile")
    if (min(profile["finite_observable_target"], profile["bath_target"], profile["jobs"]) <= 0
            or profile["exact_reference_floor"] > profile["finite_observable_target"]):
        raise ValueError("invalid deadline accuracy or worker count")
    config = deepcopy(study.read_json(ROOT/profile["base_input"]))
    config["chi_values"] = profile["chi_values"]
    config["dmrg_options"].update(profile["dmrg_options"],
                                   residual_tolerance=profile["unreferenced_residual_tolerance"])
    config["convergence"].update(target=profile["bath_target"],
                                  finite_target=profile["finite_observable_target"],
                                  exact_reference_floor=profile["exact_reference_floor"])
    config["resources"].update(jobs=profile["jobs"], total_memory_gib=profile["total_memory_gib"])
    study.validate_config(config)
    return config


def references(rows):
    selected = {}
    for row in rows:
        if row["method"] in ("ED", "BdG") and row["accepted"]:
            key = (row["point_id"], row["bath_key"], row["variant"])
            if key not in selected or row["method"] == "ED":
                selected[key] = row
    return selected


def assess(row, oracle, profile, history=()):
    """Certify listed observables, preserving finite-state residual diagnostics."""
    result = deepcopy(row)
    target = profile["finite_observable_target"]
    if oracle is not None:
        errors = {k: abs(row["values"][k]-oracle["values"][k]) for k in study.VALUES}
        errors = {k: v+profile["exact_reference_floor"] for k, v in errors.items()}
        result.update(observable_converged=max(errors.values()) <= target,
                      accuracy_evidence="exact same-bath observable comparison",
                      observable_errors=errors, method_resolution=errors,
                      reference_method=oracle["method"],
                      reference_record_sha256=fingerprint(oracle))
    elif row.get("bond_converged", False):
        resolution = {k: max(profile["exact_reference_floor"], 4*(row.get("bond_changes") or {}).get(k, 0.))
                      for k in study.VALUES}
        result.update(observable_converged=max(resolution.values()) <= target,
                      accuracy_evidence="retained stricter bond convergence", method_resolution=resolution)
    else:
        changes = study.refinement_changes([*history, row], stable_steps=2)
        residual_ok = all(s["residual"] <= profile["unreferenced_residual_tolerance"]
                          and s["sweep_energy_converged"] for s in row["sectors"])
        resolution = None if changes is None else {k: 4*v for k, v in changes.items()}
        passed = bool(residual_ok and resolution is not None and max(resolution.values()) <= target)
        result.update(observable_converged=passed,
                      accuracy_evidence="two stable bond steps and finest-state residual",
                      bond_changes=changes, method_resolution=resolution)
    return result


def prepare(directory, profile):
    directory = Path(directory)
    config = configuration(profile)
    donor = ROOT/profile["source_run"]
    original, baths, rows, _, _ = analyze.collect(donor)
    expected = study.read_json(ROOT/profile["base_input"])
    if original["input"] != expected:
        raise ValueError("source run and the declared base input disagree")
    if (original["implementation"] != source_provenance()
            or original["scripts"] != {name: study.file_hash(ROOT/"large_Gamma"/name) for name in ("run.py", "study.py")}):
        raise ValueError("source numerical implementation changed; cannot import it as the same model")
    identity_data = dict(format="qdjj-large-Gamma-deadline", format_version=1,
                         input=config, deadline_profile=profile,
                         implementation=source_provenance(), environment=numerical_environment(),
                         scripts={name: study.file_hash(ROOT/"large_Gamma"/name)
                                  for name in ("deadline.py", "run.py", "study.py")},
                         imported_source=dict(identity=fingerprint(original), manifest=original))
    manifest_path = directory/"manifest.json"
    if manifest_path.exists() and study.read_json(manifest_path) != identity_data:
        raise ValueError("deadline resume requires identical numerical settings and code")
    write_json(manifest_path, identity_data)
    identity = fingerprint(identity_data)
    for record in baths.values():
        write_json(directory/"baths"/f"{fingerprint(record['spec'])}.json", record)
    oracles = references(rows)
    for path in sorted((donor/"tasks").glob("*.json")):
        old = run.load_task(path, fingerprint(original))
        target = run.result_path(directory, old["task"])
        if target.exists():
            run.load_task(target, identity, old["task"])
            continue
        records = deepcopy(old["records"])
        oracle = oracles.get((study.point_key(old["task"]["point"]), old["task"]["bath_key"], old["task"]["variant"]))
        for i, row in enumerate(records):
            if row.get("method") == "DMRG" and row["status"] == "ok":
                records[i] = assess(row, oracle, profile)
                records[i]["original_record_sha256"] = fingerprint(row)
        complete = old["complete"]
        if old["task"]["kind"] == "dmrg":
            complete = any(r.get("observable_converged", False) for r in records)
        write_json(target, dict(identity=identity, task=old["task"], records=records,
                                records_sha256=fingerprint(records), complete=complete,
                                imported_from=dict(identity=old["identity"], task_id=path.stem,
                                                   records_sha256=old["records_sha256"])))
    return config, identity


def selected_sectors(oracle, floor):
    if oracle is None or not oracle.get("sectors"):
        return study.SECTORS
    result = []
    for parity in (0, 1):
        choices = [s for s in oracle["sectors"] if s["sector"]["parity"] == parity]
        lowest = min(s["energy"] for s in choices)
        # Retain the complete eta-degenerate subspace when necessary.
        result.extend(Sector(**s["sector"]) for s in choices if s["energy"]-lowest <= floor)
    return tuple(result)


def worker(request, destination):
    config, profile, task = request["config"], request["profile"], request["task"]
    if request["numerical_scripts"] != {name: study.file_hash(ROOT/"large_Gamma"/name)
                                       for name in ("deadline.py", "run.py", "study.py")}:
        raise ValueError("numerical code changed after this job was scheduled")
    bath = DiscreteBath.from_record(request["bath"])
    if fingerprint(request["bath"]) != task["bath_key"]:
        raise ValueError("wrong bath in deadline job")
    h = study.make_hamiltonian(bath, task["point"], config)
    path = Path(destination)
    old = run.load_task(path, request["identity"], task) if path.exists() else {}
    records = old.get("records", [])
    started = monotonic()

    def publish(complete=False):
        write_json(path, dict(identity=request["identity"], task=task, records=records,
                              records_sha256=fingerprint(records), complete=complete,
                              last_worker_seconds=monotonic()-started, peak_rss_bytes=run.peak_rss(),
                              **({"imported_from": old["imported_from"]} if "imported_from" in old else {})))

    def append(row):
        records.append(row)
        publish()

    logging.getLogger("tenpy").setLevel(logging.ERROR)
    oracle = request.get("oracle")
    if oracle is not None:
        if oracle["bath_key"] != task["bath_key"] or oracle["point"] != task["point"]:
            raise ValueError("reference and calculation have different physics")
        if any(s["hamiltonian_sha256"] != h.fingerprint() for s in oracle.get("sectors", [])):
            raise ValueError("reference uses a different Hamiltonian")
    if task["kind"] in ("ed", "quadratic", "qp"):
        if task["point"]["u"] == 0 and not any(r.get("method") == "BdG" for r in records):
            append(dict(label="quadratic", method="BdG", status="ok", accepted=True,
                        values=study.quadratic_finite(bath, task["point"]),
                        continuum=study.quadratic_continuum(task["point"], bath.bandwidth)))
        methods = []
        if task["kind"] == "ed" and task["point"]["u"] != 0:
            methods.append(("ED", None))
        if task.get("include_qp", False) or task["kind"] == "qp":
            methods += [(f"{q}QP", q) for q in config["cutoffs"]]
        for method, cutoff in methods:
            if any(r.get("method") == method and r["status"] == "ok" for r in records):
                continue
            options = dict(config["qp_options"])
            if method == "ED":
                options.update(max_dimension=6000000, max_memory_gib=16., ncv=24)
            tick = monotonic()
            try:
                estimates = [study.resource_estimate(h, cutoff, s, options) for s in study.SECTORS]
                if any(not r["allowed"] for r in estimates):
                    raise MemoryError("deadline ED/QP preflight exceeded")
                states = []
                effective = h.bath_modes if cutoff is None else min(cutoff, h.bath_modes)
                for sector in study.SECTORS:
                    before = monotonic()
                    result = solve(h, cutoff=effective, sector=sector, options=options)
                    states.append(run.state_record(result, sector, monotonic()-before))
                    del result
                append(dict(label=method, method=method, status="ok", accepted=True,
                            values=study.values_from_sectors(states), sectors=states, options=options,
                            full_space=effective == h.bath_modes, effective_cutoff=effective,
                            seconds=monotonic()-tick, dark_threshold=min(h.metadata["decoupled_qp_energies"], default=None)))
            except (ValueError, RuntimeError, MemoryError) as error:
                append(dict(label=method, method=method, status="failed", message=str(error), seconds=monotonic()-tick))
    else:
        sectors = selected_sectors(oracle, profile["exact_reference_floor"])
        previous, history = {}, []
        options = dict(config["dmrg_options"])
        for chi in config["chi_values"]:
            if any(r.get("observable_converged", False) for r in records):
                break
            tick, states = monotonic(), []
            try:
                for sector in sectors:
                    before = monotonic()
                    state = solve(h, backend="dmrg", sector=sector, options=dict(options, chi_max=chi),
                                  initial=previous.get(study.sector_key(sector)))
                    previous[study.sector_key(sector)] = state
                    record = run.state_record(state, sector, monotonic()-before, chi)
                    states.append(record)
                    append(dict(label=f"deadline_chi{chi}_{study.sector_key(sector)}", status="ok", state=record))
                row = dict(label=f"deadline_chi{chi}", method="DMRG", status="ok", chi=chi,
                           values=study.values_from_sectors(states), sectors=states, options=options,
                           accepted=all(s["finite_problem_converged"] for s in states), bond_converged=False,
                           seconds=monotonic()-tick,
                           sector_selection="exact reference minima" if len(sectors) < 4 else "both eta signs",
                           dark_threshold=min(h.metadata["decoupled_qp_energies"], default=None))
                row = assess(row, oracle, profile, history)
                history.append(row)
                append(row)
            except (ValueError, RuntimeError, MemoryError) as error:
                append(dict(label=f"deadline_chi{chi}", method="DMRG", status="failed", chi=chi, message=str(error)))
                break
    publish(complete=True)


def ordered_points(config, profile):
    return [dict(u=float(u), gamma=float(g), phi=float(phi)) for g in profile["gamma_priority"]
            if g in config["gamma_values"] for u in config["u_values"] for phi in config["phases"]]


def jobs(directory, config, profile):
    """Prioritize complete coverage, then high-Gamma and threshold-region refinements."""
    primary = config["baths"][0]
    result = []

    def add(point, levels, kind, priority, *, family=primary, bandwidth=100., include_qp=False):
        try:
            bath = run.prepare_bath(directory, family, levels, config, bandwidth)
        except (ValueError, RuntimeError) as error:
            write_json(Path(directory)/"bath_failures"/f"{family['kind']}_{levels}_{bandwidth}.json",
                       dict(family=family, levels=levels, bandwidth=bandwidth, error=str(error)))
            return
        task = run.task_record(point, bath, kind)
        if include_qp:
            task["include_qp"] = True
        if kind == "ed":
            timeout = profile["task_seconds"]["ed12" if levels >= 12 else "ed10"]
            memory = 12. if levels >= 12 and point["u"] != 0 else 3.
        elif kind == "dmrg":
            timeout = profile["task_seconds"]["dmrg8" if levels <= 8 else "dmrg_refined"]
            memory = 2. if levels <= 8 else 6.
        else:
            timeout, memory = profile["task_seconds"]["qp"], 3.
        result.append(dict(task=task, bath=bath, priority=priority, memory_gib=memory, timeout=timeout))

    grid = ordered_points(config, profile)
    for point in grid:
        add(point, 8, "dmrg", 0)
        add(point, 10, "ed", 1, include_qp=True)
        if point["gamma"] in profile["dmrg16_gammas"]:
            add(point, 16, "dmrg", 2)
        add(point, 12, "ed", 3 if point["gamma"] in profile["refined_gammas"] else 6)
        if point["gamma"] in profile["refined_gammas"]:
            add(point, 12, "dmrg", 4)
        if point["u"] == 0:
            for levels in profile["quadratic_levels"]:
                add(point, levels, "quadratic", 1)
        if point["u"] in profile["control_u"] and point["gamma"] in profile["control_gammas"]:
            for family in config["baths"][1:]:
                for levels in profile["surrogate_levels"]:
                    add(point, levels, "ed", 5, family=family, include_qp=True)
            if point["gamma"] == 10.:
                for bandwidth in profile["bandwidth_controls"]:
                    add(point, 10, "ed", 5, bandwidth=bandwidth, include_qp=True)
    # Additional Gamma samples improve the finite-bath breakdown brackets cheaply.
    _, _, rows, _, _ = analyze.collect(directory)
    comparisons, _ = analyze.compare(rows, config)
    groups = defaultdict(list)
    for point in analyze.refinement_points(comparisons, config):
        groups[(point["u"], point["phi"])].append(point)
    refined = []
    for index in range(max(map(len, groups.values()), default=0)):
        for u in config["u_values"]:
            for phi in config["phases"]:
                group = groups[(u, phi)]
                if index < len(group):
                    refined.append(group[index])
    for point in refined[:profile["max_refinement_points"]]:
        add(point, 8, "ed", 7, include_qp=True)
    return sorted(result, key=lambda item: item["priority"])


def terminate(process, grace=3.):
    if process.poll() is not None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
    except ProcessLookupError:
        process.wait()
        return
    try:
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


def worker_command(request, destination):
    return [sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(request), "--result", str(destination)]


def schedule(queue, directory, config, identity, profile, *, on_snapshot=None):
    directory = Path(directory)
    pending, active = list(queue), []
    stop_at = timestamp(profile["compute_stop"])
    last_publish, completed, stopped = monotonic(), 0, 0
    env = dict(os.environ, **{k: "1" for k in run.THREAD_VARIABLES}, PYTHONDONTWRITEBYTECODE="1")
    numerical_scripts = {name: study.file_hash(ROOT/"large_Gamma"/name) for name in ("deadline.py", "run.py", "study.py")}
    _, _, rows, _, _ = analyze.collect(directory)
    oracles = references(rows)

    def stop_job(item, reason):
        terminate(item["process"], profile["termination_grace_seconds"])
        path = item["path"]
        record = run.load_task(path, identity, item["job"]["task"]) if path.exists() else dict(
            identity=identity, task=item["job"]["task"], records=[])
        record.update(complete=True, interruption=reason, records_sha256=fingerprint(record["records"]))
        write_json(path, record)

    try:
        while pending or active:
            remaining = stop_at-time()
            if remaining <= 0:
                for item in active:
                    stop_job(item, "hard computation deadline")
                    item["log"].close()
                stopped += len(active)
                active.clear()
                break
            while pending and len(active) < profile["jobs"]:
                used = sum(item["job"]["memory_gib"] for item in active)
                index = next((i for i, job in enumerate(pending)
                              if used+job["memory_gib"] <= profile["total_memory_gib"]), None)
                if index is None:
                    break
                job = pending.pop(index)
                task = job["task"]
                path = run.result_path(directory, task)
                if path.exists() and run.load_task(path, identity, task)["complete"]:
                    continue
                oracle = oracles.get((study.point_key(task["point"]), task["bath_key"], task["variant"]))
                if task["kind"] in ("ed", "quadratic") and oracle is not None and not task.get("include_qp"):
                    continue
                request = directory/"requests"/path.name
                write_json(request, dict(config=config, profile=profile, identity=identity, task=task,
                                         bath=job["bath"]["bath"], oracle=oracle, numerical_scripts=numerical_scripts))
                log_path = directory/"logs"/f"{path.stem}.log"
                log_path.parent.mkdir(parents=True, exist_ok=True)
                log = log_path.open("w", encoding="utf-8")
                process = subprocess.Popen(worker_command(request, path), env=env, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                active.append(dict(job=job, path=path, process=process, log=log, started=monotonic()))
                print(f"Started {task['kind']} L={task['levels']} {task['point']}", flush=True)
            rss = {item["process"].pid: run.live_rss(item["process"].pid) or 0 for item in active}
            for item in list(active):
                process = item["process"]
                reason = None
                if monotonic()-item["started"] >= item["job"]["timeout"]:
                    reason = "deadline-profile task time limit"
                elif rss[process.pid] > profile["worker_memory_gib"]*2**30:
                    reason = "worker memory limit"
                elif sum(rss.values()) > profile["total_memory_gib"]*2**30 and rss[process.pid] == max(rss.values()):
                    reason = "aggregate memory limit"
                if reason:
                    stop_job(item, reason)
                    stopped += 1
                code = process.poll()
                if code is None:
                    continue
                if code != 0 and reason is None:
                    stop_job(item, f"worker exit {code}")
                    stopped += 1
                item["log"].close()
                active.remove(item)
                completed += 1
                if item["path"].exists():
                    result = run.load_task(item["path"], identity)
                    for row in result["records"]:
                        if row.get("method") in ("ED", "BdG") and row["status"] == "ok":
                            enriched = dict(row, **result["task"], point_id=study.point_key(result["task"]["point"]))
                            oracles[(enriched["point_id"], enriched["bath_key"], enriched["variant"])] = enriched
                print(f"Finished {completed}; pending={len(pending)}; active={len(active)}", flush=True)
            write_json(directory/"progress.json", dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                       pending=len(pending), active=[dict(pid=i["process"].pid, task=i["job"]["task"]) for i in active],
                       completed_new_tasks=completed, stopped_tasks=stopped,
                       compute_stop=profile["compute_stop"], report_deadline=profile["report_deadline"]))
            if on_snapshot and monotonic()-last_publish >= profile["publish_interval_seconds"]:
                # Publishing runs with a bounded timeout shorter than the time to cutoff.
                allowance = min(profile["publication_timeout_seconds"], max(0., stop_at-time()-5))
                if allowance > 10:
                    on_snapshot(allowance)
                last_publish = monotonic()
            if active:
                sleep(min(profile["poll_seconds"], max(0., stop_at-time())))
    finally:
        for item in active:
            stop_job(item, "deadline controller interrupted")
            item["log"].close()
    return dict(completed_new_tasks=completed, stopped_tasks=stopped, unscheduled_tasks=len(pending))


def publish(directory, output, *, final=False):
    summary = analyze.export(directory, output, plots=True)
    from large_Gamma.deadline_report import write_report
    write_report(directory, output, summary, final=final)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, default=DEFAULT_PROFILE)
    parser.add_argument("--publish-only", action="store_true")
    parser.add_argument("--final", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        if args.worker:
            worker(study.read_json(args.worker), args.result)
            return
        profile = study.read_json(args.profile)
        directory, output = ROOT/profile["output"], ROOT/profile["publish"]
        if args.publish_only:
            publish(directory, output, final=args.final)
            return
        if time() >= timestamp(profile["compute_stop"]):
            parser.error("the computation deadline has passed; use --publish-only --final")
        config, identity = prepare(directory, profile)
        publish(directory, output)
        queue = jobs(directory, config, profile)
        status = dict(pid=os.getpid(), state="running", started_utc=datetime.now(timezone.utc).isoformat(),
                      compute_stop=profile["compute_stop"], report_deadline=profile["report_deadline"],
                      input_sha256=study.file_hash(args.profile), total_planned_tasks=len(queue))
        write_json(directory/"campaign.json", status)

        def snapshot(timeout):
            try:
                subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()), "--profile", str(args.profile),
                                "--publish-only"], timeout=timeout, check=True, start_new_session=True)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                print(f"Intermediate publication deferred: {error}", flush=True)

        def stop(signum, frame):
            raise KeyboardInterrupt

        signal.signal(signal.SIGTERM, stop)
        signal.signal(signal.SIGINT, stop)
        try:
            status.update(schedule(queue, directory, config, identity, profile, on_snapshot=snapshot))
            status["state"] = "finalizing"
        except KeyboardInterrupt:
            status["state"] = "interrupted_finalizing"
        except (OSError, ValueError, RuntimeError) as error:
            status.update(state="failed_finalizing", error=str(error))
        write_json(directory/"campaign.json", status)
        # All numerical children have stopped before this phase starts.
        before = time()
        try:
            allowance = min(300., timestamp(profile["report_deadline"])-time()-60.)
            if allowance <= 0:
                raise RuntimeError("no time left to publish before the report deadline")
            subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()), "--profile", str(args.profile),
                            "--publish-only", "--final"], timeout=allowance, check=True, start_new_session=True)
            analyze.verify_archive(output)
            status.update(state="final_report_published", finished_utc=datetime.now(timezone.utc).isoformat(),
                          publication_seconds=time()-before)
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            status.update(state="publication_failed", error=str(error))
        write_json(directory/"campaign.json", status)
        print(status, flush=True)


if __name__ == "__main__":
    main()
