"""Run and publish successive study stages under a multi-day scheduling budget.

Numerical work is delegated to run.py's bounded worker controller. This process
can wait for specified earlier controllers before starting, avoiding duplicate
jobs. Campaign status and logs live in the ignored working directory.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import signal
import subprocess
import sys
from time import monotonic, sleep

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qdjj_solver.common.io import write_json
from large_Gamma.study import DEFAULT_INPUT, file_hash, read_json


def stages():
    return [
        ("pilot-L8", ["--profile", "pilot", "--max-levels", "8"]),
        ("grid-L8", ["--max-levels", "8"]),
        ("controls-L8", ["--controls", "--max-levels", "8"]),
        ("pilot-L16", ["--profile", "pilot", "--max-levels", "16"]),
        ("grid-L16", ["--max-levels", "16"]),
        ("controls-L16", ["--controls", "--max-levels", "16"]),
        ("grid-L24", ["--max-levels", "24"]),
        ("grid-L32", ["--max-levels", "32"]),
        ("controls-L32", ["--controls", "--max-levels", "32"]),
        ("grid-L64", ["--max-levels", "64"]),
        ("Gamma-refinement", ["--max-levels", "64", "--refine"]),
    ]


def active_study_process(pid):
    if os.name == "posix":
        try:
            command = subprocess.check_output(["ps", "-o", "command=", "-p", str(pid)],
                                              text=True, stderr=subprocess.DEVNULL)
            return "large_Gamma" in command and "python" in command.lower()
        except subprocess.CalledProcessError:
            return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=ROOT/"runs/large_Gamma/study")
    parser.add_argument("--publish", type=Path, default=ROOT/"large_Gamma/output")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--max-hours", type=float, default=168.)
    parser.add_argument("--wait-for", type=int, nargs="*", default=[])
    args = parser.parse_args()
    if args.max_hours <= 0 or args.jobs < 1:
        parser.error("positive scheduling budget and worker count required")
    directory = args.output
    directory.mkdir(parents=True, exist_ok=True)
    status = dict(pid=os.getpid(), started_utc=datetime.now(timezone.utc).isoformat(),
                  max_hours=args.max_hours, wait_for=args.wait_for,
                  campaign_sha256=file_hash(__file__), completed_stages=[], state="waiting")
    status_path = directory/"campaign.json"
    write_json(status_path, status)
    child = None

    def stop(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, stop)
    try:
        while any(active_study_process(pid) for pid in args.wait_for):
            sleep(10)
        started = monotonic()
        for label, extra in stages():
            remaining = args.max_hours-(monotonic()-started)/3600
            if remaining <= 0:
                status["state"] = "scheduling_budget_exhausted"
                break
            command = [sys.executable, "-B", str(ROOT/"large_Gamma/run.py"),
                       "--input", str(args.input), "--output", str(directory),
                       "--jobs", str(args.jobs), "--wall-hours", str(remaining), *extra]
            status.update(state="running", stage=label, command=command)
            write_json(status_path, status)
            print(f"Starting campaign stage {label}; scheduling budget {remaining:.2f} h", flush=True)
            child = subprocess.Popen(command)
            code = child.wait()
            child = None
            # Publish all durable evidence even if the controller was interrupted.
            subprocess.run([sys.executable, "-B", str(ROOT/"large_Gamma/analyze.py"), str(directory),
                            "--output", str(args.publish), "--plots"], check=True)
            if code:
                status.update(state="controller_failed", exit_code=code)
                break
            status["completed_stages"].append(label)
        else:
            summary = read_json(args.publish/"summary.json")
            status["state"] = "converged" if summary["complete"] else "ladders_finished_with_unresolved_targets"
    except KeyboardInterrupt:
        if child is not None:
            child.send_signal(signal.SIGINT)
            child.wait()
        status["state"] = "interrupted"
    except (OSError, subprocess.CalledProcessError, ValueError) as error:
        status.update(state="campaign_failed", exception=type(error).__name__, message=str(error))
    finally:
        status["updated_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(status_path, status)
    print(status["state"], flush=True)


if __name__ == "__main__":
    main()
