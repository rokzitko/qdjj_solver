#!/usr/bin/env python3
"""Budgeted reproduction of the published finite-bath endpoints from fresh states."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import signal
import sys
import time


HERE = Path(__file__).resolve().parent
if str(HERE.parents[1]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[1]))
from NRG_comparisons.test1 import run, targeted_nrg


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def preflight(task):
    """Count the declared uncompressed basis before constructing any Hamiltonian."""
    from qdjj_solver.common.problem import Sector
    from qdjj_solver.qp_solver.basis import dimension

    if task["backend"] == "dmrg":
        return dict(backend="dmrg", bond=task["chi"], exact_rank_ceiling=4 ** ((task["bath"]["levels"] + 1) // 2))
    spins = (1, -1) * (1 + task["bath"]["levels"])
    counts = {name: dimension(2, spins, task["cutoff"], Sector(p, p))
              for name, p in (("singlet", 0), ("doublet", 1))}
    if max(counts.values()) > task["max_dimension"]:
        raise MemoryError(f"focused task exceeds its declared dimension guard: {counts}")
    return dict(backend="qp", dimensions=counts, max_dimension=task["max_dimension"],
                memory_guard_gib=96, nnz_guard=500_000_000)


def compare(a, b, target):
    result = {}
    for observable in sorted(a.keys() & b.keys()):
        if a[observable] is None or b[observable] is None:
            continue
        change = a[observable] - b[observable]
        threshold = target["absolute"] + target["relative"] * abs(a[observable])
        result[observable] = dict(change=change, threshold=threshold,
                                  ratio=abs(change) / threshold, stable_change=abs(change) <= threshold)
    return result


def verify_sources(identity):
    for entry in identity["sources"].values():
        if sha(HERE / entry["path"]) != entry["sha256"]:
            raise ValueError(f"focused source changed: {entry['path']}")


def validate_bath(physical, specification, record):
    from qdjj_solver.common.baths import DiscreteBath, cosh_grid
    bath = DiscreteBath.from_record(record)
    if bath.levels != specification["levels"] or bath.delta != physical["gap"] or bath.bandwidth != physical["bandwidth"]:
        raise ValueError("focused bath dimensions/scales differ from the declared recipe")
    if specification["kind"] == "cosh":
        expected = cosh_grid(bath.levels // 2, physical["gap"], physical["bandwidth"]).record()
        import numpy as np
        if not (np.allclose(record["xi"], expected["xi"], rtol=1e-13, atol=1e-14)
                and np.allclose(record["weights"], expected["weights"], rtol=1e-13, atol=1e-14)):
            raise ValueError("focused cosh coefficients do not match the declared recipe")
    elif specification["kind"] == "surrogate":
        expected = dict(kind="surrogate", levels=bath.levels, frequency_min=1e-3 * bath.delta,
                        frequency_cutoff=specification["frequency_cutoff"], relative_weight=specification["relative_weight"],
                        starts=4, seed=1729, frequency_points=1000, max_nfev=specification["max_nfev"])
        if any(bath.metadata.get(k) != v for k, v in expected.items()):
            raise ValueError("focused surrogate settings do not match the declared recipe")
    else:
        raise ValueError("unsupported focused bath family")


def check_comparisons(identity):
    verify_sources(identity)
    references = {name: json.loads((HERE / entry["path"]).read_text())
                  for name, entry in identity["sources"].items() if name.startswith("comparison:")}
    for name, row in references.items():
        if row.get("status") != "completed" or row.get("physical") != identity["manifest"]["physical"]:
            raise ValueError(f"focused comparison physical identity/status mismatch: {name}")
        json.dumps(row, allow_nan=False)
    return references


def summarize(campaign, identity):
    references = check_comparisons(identity)
    records, failures = [], []
    for entry in identity["specification"]["tasks"]:
        task = entry["calculation"]
        taskdir = campaign / "tasks" / (task["backend"] + "-" + run.digest(task)[:16])
        results = sorted(taskdir.glob("attempt-*/result.json"))
        completed = [p for p in results if json.loads(p.read_text()).get("status") == "completed"]
        for path in results:
            if path not in completed:
                failures.append(dict(path=str(path.relative_to(campaign)), result=json.loads(path.read_text())))
        if completed:
            row = json.loads(completed[-1].read_text())
            if row["task"] != task or row["physical"] != identity["manifest"]["physical"]:
                raise ValueError("focused result identity mismatch")
            records.append(dict(record=row, source=str(completed[-1].relative_to(campaign)), sha256=sha(completed[-1])))
    comparisons = []
    candidates = [(r["source"], r["record"]) for r in records]
    previous = [(name, row) for name, row in references.items()]
    for name, row in candidates:
        target = identity["manifest"]["targets"][row["task"]["backend"]]
        for other_name, other in previous:
            task, before = row["task"], other["task"]
            matching = (task["backend"] == "qp" and before["backend"] == "qp" and
                        (task["cutoff"] == before["cutoff"] or task["bath"] == before["bath"]))
            matching |= task["backend"] == "dmrg" and before["backend"] in ("qp", "dmrg") and task["bath"] == before["bath"]
            if matching:
                comparisons.append(dict(records=[name, other_name], new_finite_converged=row["finite_converged"],
                                        old_finite_converged=other["finite_converged"],
                                        changes=compare(row["values"], other["values"], target)))
        previous.append((name, row))
    result = dict(format_version=1, manifest=identity, records=records, failures=failures, comparisons=comparisons,
                  complete=len(records) == len(identity["specification"]["tasks"]),
                  continuum_converged=False,
                  scope="Focused differences only. Missing cutoff/bond/bath refinements are not certified by agreement.")
    run.write_json(campaign / "summary.json", result)
    lines = ["# Focused QP/DMRG Follow-Up", "", result["scope"], "",
             f"Completed: {len(records)}/{len(identity['specification']['tasks'])}; failed attempts: {len(failures)}.", "",
             "| Backend | Bath | Levels | Cutoff/Bond | Signed gap | Finite checks |",
             "|---|---|---:|---:|---:|---|"]
    for entry in records:
        row, task = entry["record"], entry["record"]["task"]
        lines.append(f"| {task['backend']} | {task['bath']['kind']} | {task['bath']['levels']} | "
                     f"{task.get('cutoff', task.get('chi'))} | {row['values']['signed_gap']} | {row['finite_converged']} |")
    (campaign / "REPORT.md").write_text("\n".join(lines) + "\n")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "prepare", "run", "report"))
    parser.add_argument("--campaign", default="finite-reproduction")
    parser.add_argument("--group", choices=("window", "bath", "dmrg"))
    parser.add_argument("--max-tasks", type=int)
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign) or (args.max_tasks is not None and args.max_tasks < 1):
        raise ValueError("invalid campaign name or task limit")
    manifest, _ = run.load_inputs()
    spec = json.loads((HERE / "input/qp-dmrg-focused.json").read_text())
    if spec["format_version"] != 1 or any(isinstance(h, bool) or not math.isfinite(h) or not 0 < h <= 4 for h in spec["group_hours"].values()):
        raise ValueError("invalid focused group budget")
    counts = [preflight(e["calculation"]) for e in spec["tasks"]]
    if args.command == "plan":
        print(json.dumps(dict(manifest=manifest, specification=spec, preflight=counts), indent=2))
        return
    runs = HERE / "runs"
    runs.mkdir(exist_ok=True)
    with run.lock(runs / ".campaign.lock"):
        campaign = runs / args.campaign
        campaign.mkdir(exist_ok=True)
        path = campaign / "manifest.json"
        if args.command == "report":
            summarize(campaign, json.loads(path.read_text()))
            return
        baths = json.loads((HERE / "input/final-baths.json").read_text())
        for entry in spec["tasks"]:
            validate_bath(manifest["physical"], entry["calculation"]["bath"], baths[entry["bath_input"]])
        identity = dict(manifest=manifest, specification=spec, preflight=counts, sources={}, baths=baths,
                         provenance=run.provenance(["qp", "dmrg"]))
        if path.exists() and json.loads(path.read_text()) != identity:
            raise ValueError("focused inputs/implementation changed; use a new campaign")
        run.write_json(path, identity)
        if args.command == "prepare":
            print(f"Prepared {campaign}: {counts}")
            return
        count = 0
        for entry in spec["tasks"]:
            group, task = entry["group"], entry["calculation"]
            if args.group is not None and group != args.group:
                continue
            if args.max_tasks is not None and count >= args.max_tasks:
                break
            taskdir = campaign / "tasks" / (task["backend"] + "-" + run.digest(task)[:16])
            attempts = sorted(taskdir.glob("attempt-*"))
            if any((p / "result.json").exists() and json.loads((p / "result.json").read_text()).get("status") == "completed" for p in attempts):
                continue
            if attempts and not args.retry_failed:
                raise RuntimeError("focused failure requires explicit review/retry")
            global_remaining = manifest["limits"]["campaign_hours"] * 3600 - targeted_nrg.spent_seconds()
            group_spent = 0.
            for prior in campaign.glob("tasks/*/attempt-*/result.json"):
                request = json.loads((prior.parent / "request.json").read_text())
                if request["group"] == group:
                    group_spent += json.loads(prior.read_text())["monitored_seconds"]
            remaining = min(global_remaining, spec["group_hours"][group] * 3600 - group_spent)
            if remaining <= 0:
                continue
            directory = taskdir / f"attempt-{len(attempts) + 1:03d}"
            directory.mkdir(parents=True)
            run.write_json(taskdir / "task.json", task)
            request = dict(group=group, task=task, physical=manifest["physical"])
            run.write_json(directory / "bath.json", identity["baths"][entry["bath_input"]])
            request["bath_sha256"] = sha(directory / "bath.json")
            run.write_json(directory / "request.json", request)
            print(f"FOCUSED {group}: {task}", flush=True)
            started = time.monotonic()
            result = run.execute(directory, manifest["limits"], remaining)
            bath_path = directory / "bath.json"
            validate_bath(manifest["physical"], task["bath"], json.loads(bath_path.read_text()))
            if sha(bath_path) != request["bath_sha256"]:
                raise ValueError("worker did not replay the retained bath")
            print(f"  {result['status']} ({time.monotonic() - started:.1f}s): {result.get('error', result.get('values'))}", flush=True)
            count += 1
            summarize(campaign, identity)
            if result["status"] != "completed":
                raise RuntimeError("focused calculation requires review")
        result = summarize(campaign, identity)
        print(f"{campaign / 'REPORT.md'}: {len(result['records'])} completed", flush=True)


if __name__ == "__main__":
    def interrupted(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupted)
    main()
