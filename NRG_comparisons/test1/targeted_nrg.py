#!/usr/bin/env python3
"""Nested-twist NRG refinement, from scratch or with an optional existing baseline."""

from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re
import signal
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from NRG_comparisons.test1 import nrg_dataset, run
from NRG_comparisons.test1.nrg import render_param


def plan_for(stage, specification):
    if specification["format_version"] != 1:
        raise ValueError("unknown targeted NRG specification")
    settings = specification["stages"][stage]
    anchors = deepcopy(specification["anchors"])
    for anchor, nz in zip(anchors, settings["nz"], strict=True):
        if anchor["lambda"] < 1.8:
            raise ValueError("Lambda < 1.8 is not authorized")
        n = anchor["nmax"]
        anchor.update(nz=nz, checks=dict(
            keep=[2000, 4000, 8000] if anchor["lambda"] in settings["tighter_keep_lambda"] else settings["keep_checks"],
            keepenergy=[10, 12, 14], keepmin=[300, 500, 1000], nmax=[n - 8, n - 4, n]))
        anchor["tail_checks"] = [n + 4, n + 8] if settings["tail_checks"] and nz > 8 else []
    return dict(policy="targeted-policy-v1", stage=stage, anchors=anchors)


def tasks_for(plan):
    """Reference anchors first; auxiliary probes never change reference settings."""
    from NRG_comparisons.test1.twist_analysis import assess
    # Validate the uncertainty plan; main also checks every literal input deck.
    assess([], plan, {"nrg": {"absolute": 1e-6, "relative": 1e-3}})
    anchors, probes = {}, {}
    for a in plan["anchors"]:
        for k in range(1, a["nz"] + 1):
            numerical = dict(z=k / a["nz"], untruncated=False,
                             **{key: a[key] for key in ("lambda", "nmax", "keep", "keepenergy", "keepmin")})
            task = dict(backend="nrg", clean=False, numerical=numerical, units="gap")
            anchors[run.digest(task)] = task
            for axis, values in a["checks"].items():
                for value in values:
                    if value != a[axis]:
                        probe = task | {"numerical": numerical | {axis: value}}
                        probes[run.digest(probe)] = probe
            if k == 1:
                for nmax in a["tail_checks"]:
                    probe = task | {"numerical": numerical | {"nmax": nmax}}
                    probes[run.digest(probe)] = probe
            if k == a["nz"]:
                unit_check = task | {"units": "bandwidth"}
                probes[run.digest(unit_check)] = unit_check
    controls = [t for t in run.make_tasks("controls", {"wilson_nmax": [2, 4]}, {})
                if t["backend"] == "nrg"]
    return controls + list(anchors.values()) + [t for key, t in probes.items() if key not in anchors], set(anchors)


def implementation():
    result = run.provenance(["nrg"], nrg_only=True)
    result["workflow"].update({name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ("targeted_nrg.py", "twist_analysis.py", "nrg_dataset.py")})
    result["postprocessor"] = nrg_dataset.fingerprint()
    return result


def spent_seconds():
    """Use the existing shared ledger; charge unfinished reservations conservatively."""
    spent = 0.

    def charge(seconds):
        if isinstance(seconds, bool) or not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds < 0:
            raise ValueError("ledger durations must be finite nonnegative numbers")
        return seconds

    for attempt in (HERE / "runs").glob("*/tasks/*/attempt-*"):
        result_path, journal_path = attempt / "result.json", attempt / "launch.json"
        result = json.loads(result_path.read_text()) if result_path.exists() else {}
        if "monitored_seconds" in result:
            spent += charge(result["monitored_seconds"])
        elif journal_path.exists():
            journal = json.loads(journal_path.read_text())
            seconds = charge(journal["timeout_seconds"])
            if journal["worker_pid"] and run.process_memory(journal["worker_pid"]):
                raise RuntimeError(f"unfinished worker session exists: {attempt}")
            spent += seconds
            result.update(status="failed", error="interrupted or unmonitored attempt",
                          monitored_seconds=seconds, accounting="conservative reservation")
            run.write_json(result_path, result)
    return spent


def report(campaign, identity, snapshot):
    from NRG_comparisons.test1.twist_analysis import assess

    records = list(snapshot["records"])
    sources = [dict(kind="baseline", source=s) for s in snapshot["sources"]]
    known = {run.digest(r["task"]): i for i, r in enumerate(records)}
    inventory = [dict(task=r["task"], status="completed", record=i) for i, r in enumerate(records)]
    failures = []
    for task in identity["tasks"]:
        key = run.digest(task)
        if key in known:
            continue
        taskdir = campaign / "tasks" / ("nrg-" + key[:16])
        attempts = sorted(taskdir.glob("attempt-*/result.json"))
        completed = [p for p in attempts if json.loads(p.read_text()).get("status") == "completed"]
        for path in attempts:
            if path not in completed:
                failures.append(dict(path=str(path.relative_to(campaign)), result=json.loads(path.read_text())))
        index = None
        if completed:
            path = completed[-1]
            row = json.loads(path.read_text())
            if row["task"] != task or row["physical"] != identity["manifest"]["physical"]:
                raise ValueError("targeted result identity mismatch")
            index = len(records)
            records.append(row)
            sources.append(dict(kind="new", path=str(path.relative_to(campaign)),
                                sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        inventory.append(dict(task=task, record=index,
                              status="completed" if index is not None else "failed" if attempts else "unstarted"))
    result = assess(records, identity["plan"], identity["manifest"]["targets"], inventory)
    result.update(record_sources=sources, failed_attempts=failures, baseline_records=len(snapshot["records"]),
                  new_completed=len(records) - len(snapshot["records"]),
                  analysis_sha256=hashlib.sha256((HERE / "twist_analysis.py").read_bytes()).hexdigest())
    run.write_json(campaign / "summary.json", result)
    lines = ["# Targeted NRG Refinement", "", "Empirical uncertainty estimates, not rigorous bounds.", "",
             f"Stage: `{identity['plan']['stage']}`. Reused: {len(snapshot['records'])}; new completed: {result['new_completed']}; failed attempts: {len(failures)}.", "",
             "| Observable | Extrapolated candidate | Empirical uncertainty | Target | Status |",
             "|---|---:|---:|---:|---|"]
    for name, row in result["nrg_extrapolation"].items():
        lines.append(f"| {name} | {row['value']} | {row['empirical_uncertainty']} | {row['threshold']} | {row['status']} |")
    (campaign / "REPORT.md").write_text("\n".join(lines) + "\n")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("plan", "prepare", "run", "report"))
    parser.add_argument("--stage", choices=("gap16", "energy32", "retention2"), default="retention2")
    parser.add_argument("--campaign", default="nrg-reproduction")
    parser.add_argument("--source", help="optional existing campaign to reuse; otherwise start from scratch")
    parser.add_argument("--mode", choices=("all", "anchors", "probes"), default="all")
    parser.add_argument("--lambda", dest="lambda_value", type=float)
    parser.add_argument("--max-tasks", type=int)
    parser.add_argument("--retry-failed", action="store_true")
    args = parser.parse_args(argv)
    for name in (args.campaign, args.source):
        if name is not None and not re.fullmatch(r"[A-Za-z0-9_-]+", name):
            raise ValueError("campaign/source must be simple directory names")
    if args.campaign == args.source or (args.max_tasks is not None and args.max_tasks < 1):
        raise ValueError("invalid source campaign or task limit")
    manifest, _ = run.load_inputs()
    manifest = dict(manifest, limits=dict(manifest["limits"]))
    for key, cap in (("cpus", 15), ("memory_gib", 126), ("campaign_hours", 47.5)):
        manifest["limits"][key] = min(manifest["limits"][key], cap)
    specification = json.loads((HERE / "input/nrg-targeted.json").read_text())
    plan = plan_for(args.stage, specification)
    tasks, anchors = tasks_for(plan)
    for task in tasks:
        render_param(manifest["physical"], task["numerical"], task["units"], task["clean"])
    if args.command == "plan":
        print(json.dumps(dict(manifest=manifest, plan=plan, tasks=tasks), indent=2))
        return
    campaign = HERE / "runs" / args.campaign
    (HERE / "runs").mkdir(exist_ok=True)
    with run.lock(HERE / "runs/.campaign.lock"):
        campaign.mkdir(exist_ok=True)
        path = campaign / "manifest.json"
        snapshot = dict(format_version=1, records=[], sources=[], verification_bytes=0)
        if path.exists():
            identity = json.loads(path.read_text())
            if "baseline" in identity:
                snapshot = nrg_dataset.load_snapshot(campaign / identity["baseline"]["path"], identity["baseline"]["sha256"])
            if args.command == "report":
                report(campaign, identity, snapshot)
                return
            if any(identity[key] != value for key, value in dict(manifest=manifest, plan=plan, tasks=tasks,
                                                               provenance=implementation(), source_campaign=args.source).items()):
                raise ValueError("targeted inputs/implementation changed; use a new campaign")
        else:
            if args.command == "report":
                raise FileNotFoundError("targeted campaign is not prepared")
            limits = manifest["limits"]
            if spent_seconds() >= limits["campaign_hours"] * 3600:
                raise RuntimeError("shared campaign budget exhausted before baseline verification")
            attempts = sorted((campaign / "tasks/audit").glob("attempt-*"))
            if attempts and not args.retry_failed:
                raise RuntimeError("failed baseline verification requires explicit review and retry-failed")
            audit = campaign / "tasks/audit" / f"attempt-{len(attempts) + 1:03d}"
            audit.mkdir(parents=True)
            snapshot_path = campaign / "baseline.json"
            if snapshot_path.exists():
                snapshot_path = audit / "baseline.json"
            started = time.monotonic()
            try:
                if args.source is not None:
                    snapshot = nrg_dataset.freeze(HERE / "runs" / args.source, manifest["physical"], snapshot_path)
            except BaseException as error:
                run.write_json(audit / "result.json", dict(status="failed", error=str(error), monitored_seconds=time.monotonic() - started))
                raise
            run.write_json(audit / "result.json", dict(status="completed", monitored_seconds=time.monotonic() - started,
                                                       verification_bytes=snapshot["verification_bytes"]))
            identity = dict(manifest=manifest, policy="targeted-policy-v1", plan=plan, tasks=tasks,
                            provenance=implementation(), source_campaign=args.source)
            if args.source is not None:
                identity["baseline"] = dict(path=str(snapshot_path.relative_to(campaign)),
                                            sha256=hashlib.sha256(snapshot_path.read_bytes()).hexdigest())
            run.write_json(path, identity)
        reused = {run.digest(r["task"]) for r in snapshot["records"]}
        if args.command == "prepare":
            report(campaign, identity, snapshot)
            print(f"Prepared {args.campaign}: {len(reused)} audited records, {sum(run.digest(t) not in reused for t in tasks)} new tasks")
            return
        count = 0
        for task in tasks:
            key = run.digest(task)
            if (key in reused or (args.mode == "anchors" and key not in anchors)
                    or (args.mode == "probes" and key in anchors)
                    or (args.lambda_value is not None and task["numerical"]["lambda"] != args.lambda_value)):
                continue
            if args.max_tasks is not None and count >= args.max_tasks:
                break
            taskdir = campaign / "tasks" / ("nrg-" + key[:16])
            attempts = sorted(taskdir.glob("attempt-*"))
            if any((p / "result.json").exists() and json.loads((p / "result.json").read_text()).get("status") == "completed" for p in attempts):
                continue
            if attempts and not args.retry_failed:
                raise RuntimeError(f"unfinished/failed targeted task requires explicit review: {taskdir}")
            remaining = manifest["limits"]["campaign_hours"] * 3600 - spent_seconds()
            if remaining <= 0:
                break
            directory = taskdir / f"attempt-{len(attempts) + 1:03d}"
            directory.mkdir(parents=True)
            run.write_json(taskdir / "task.json", task)
            run.write_json(directory / "request.json", dict(task=task, physical=manifest["physical"], nrg_refinement=True))
            (directory / "param").write_text(render_param(manifest["physical"], task["numerical"], task["units"], task["clean"]))
            print(f"TARGETED NRG {taskdir.name}: {json.dumps(task, sort_keys=True)}", flush=True)
            result = run.execute(directory, manifest["limits"], remaining)
            count += 1
            report(campaign, identity, snapshot)
            print(f"  {result['status']}: {result.get('error', result.get('values'))}", flush=True)
            if result["status"] != "completed":
                raise RuntimeError("targeted NRG paused on failure; review before retry")
        assessment = report(campaign, identity, snapshot)
        print(f"{campaign / 'REPORT.md'}: {assessment['new_completed']} new records", flush=True)


if __name__ == "__main__":
    def interrupted(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupted)
    main()
