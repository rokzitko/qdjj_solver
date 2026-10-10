#!/usr/bin/env python3
"""Explicit, sequential test3 controls and cost pilots; never a production scan."""

import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from NRG_comparisons.test1.run import digest, lock, process_memory, write_json
from NRG_comparisons.test3.manual_run.nrg import extract_result, render_param
from NRG_comparisons.test3.manual_run.physics import OBSERVABLES, PRIMARY


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def tasks(stage):
    """Frozen pilot inventory. No automatic retries or production expansion."""
    if stage == "controls":
        result = []
        for label, nmax, changes in (("atomic", 1, {"gamma": 0.}), ("chain", 2, {}),
                                     ("zero", 2, {"field": 0.}), ("reverse", 2, {"field": -.03})):
            result.append(dict(id=f"{label}-nrg", backend="nrg", physical_changes=changes,
                               numerical=dict(nmax=nmax, keep=512, keepmin=0, keepenergy=-1.,
                                              untruncated=True, z=1., **{"lambda": 2.})))
            for backend in ("ed", "qp", "dmrg") if label in {"atomic", "chain"} else ("ed",):
                result.append(dict(id=f"{label}-{backend}", backend=backend, physical_changes=changes,
                                   bath_source=f"{label}-nrg", chi=32, cutoff=None, roots=2,
                                   spins=list(range(-nmax-2, nmax+3))))
        result.append(dict(id="chain-bandwidth", backend="nrg", units="bandwidth",
                           numerical=dict(nmax=2, keep=512, keepmin=0, keepenergy=-1., untruncated=True,
                                          z=1., **{"lambda": 2.})))
        for step in (2e-4, 1e-4):
            for sign in (-1, 1):
                result.append(dict(id=f"hf-{step:g}-{sign:+d}", backend="ed", bath_source="chain-nrg",
                                   physical_changes={"field": .03+sign*step}, roots=2, spins=[0, -1, 1]))
        return result
    if stage != "pilot":
        raise ValueError("stage must be controls or pilot")
    result = [dict(id=f"cosh8-qp-q{cutoff if cutoff is not None else 'full'}", backend="qp",
                   bath="cosh8", cutoff=cutoff, roots=2) for cutoff in (2, 4, 6, None)]
    result += [dict(id=f"cosh8-dmrg-chi{chi}", backend="dmrg", bath="cosh8", chi=chi, roots=2)
               for chi in (32, 64, 128)]
    result += [dict(id=f"nrg-k{keep}", backend="nrg", numerical=dict(
        nmax=42, keep=keep, keepmin=keep//5, keepenergy=14., z=1., **{"lambda": 2.}))
        for keep in (500, 1000, 2000)]
    return result


def chain_bath(result):
    """Use the actual generated normal chain, not an independently rounded grid."""
    import numpy as np
    from scipy.linalg import eigh_tridiagonal

    chain = result["diagnostics"]["wilson_chain"]
    xi, vectors = eigh_tridiagonal(chain["zeta_over_gap"], chain["xi_over_gap"][:-1], lapack_driver="stev")
    weights = vectors[0]**2
    if np.any(weights <= 0) or abs(weights.sum()-1) > 1e-12:
        raise ValueError("generated chain lost its contact measure")
    return dict(xi=xi.tolist(), weights=weights.tolist(), delta=1., bandwidth=100.,
                metadata={"kind": "generated-nrg-chain", "chain": chain, "weights_renormalized": False})


def check_controls(results):
    """Absolute control bounds are finite-chain checks, not continuum targets."""
    checks = []

    def rows(name):
        return {(r["twice_sz"], r["root"]): r for r in results[name]["states"]}

    for label, maximum in (("atomic", 3), ("chain", 4), ("zero", 4), ("reverse", 4)):
        expected = {(spin, root) for spin in range(-maximum, maximum+1)
                    for root in range(1 if abs(spin) == maximum else 2)}
        missing = expected - rows(f"{label}-ed").keys()
        checks.append(dict(name=f"{label}-ed declared root inventory", passed=not missing,
                           missing=[list(key) for key in sorted(missing)]))

    def compare(new, old):
        a, b = rows(new), rows(old)
        checks.append(dict(name=f"{new} covers all {old} control roots", passed=b.keys() <= a.keys(),
                           missing=[list(key) for key in sorted(b.keys()-a.keys())]))
        errors = []
        for key in a.keys() & b.keys():
            names = ("energy",) if a[key]["degenerate"] or b[key]["degenerate"] else OBSERVABLES
            errors.extend(abs(a[key]["values"][name]-b[key]["values"][name]) for name in names)
        if not all((spin, 0) in a.keys() & b.keys() for spin in PRIMARY):
            raise ValueError("control lacks a primary state")
        error = max(errors)
        checks.append(dict(name=f"{new} vs {old}", maximum_error=error, tolerance=1e-8, passed=error < 1e-8))

    for label in ("atomic", "chain"):
        for backend in ("nrg", "qp", "dmrg"):
            compare(f"{label}-{backend}", f"{label}-ed")
    for label in ("zero", "reverse"):
        compare(f"{label}-nrg", f"{label}-ed")
    compare("chain-bandwidth", "chain-nrg")
    atomic = rows("atomic-nrg")
    error = max(abs(atomic[key]["values"]["energy"]-energy)
                for key, energy in (((0, 0), .5), ((0, 1), 1.5), ((-1, 0), -.015), ((1, 0), .015)))
    checks.append(dict(name="atomic energies and Zeeman sign", maximum_error=error, tolerance=1e-9, passed=error < 1e-9))
    positive, negative = rows("chain-nrg"), rows("reverse-nrg")
    errors = []
    for spin in PRIMARY:
        a, b = positive[spin, 0]["values"], negative[-spin, 0]["values"]
        for name in ("energy", "P0", "P1", "P2", "dot_z2", "bath_z2", "Czz", "Czz_connected"):
            errors.append(abs(a[name]-b[name]))
        errors.extend(abs(a[name]+b[name]) for name in ("m_dot", "m_bath"))
    error = max(errors)
    checks.append(dict(name="NRG spin reversal", maximum_error=error, tolerance=1e-8, passed=error < 1e-8))
    for step in (2e-4, 1e-4):
        minus, plus, measured = rows(f"hf-{step:g}--1"), rows(f"hf-{step:g}-+1"), rows("chain-qp")
        error = max(abs((plus[spin, 0]["values"]["energy"]-minus[spin, 0]["values"]["energy"])/(2*step)
                        - measured[spin, 0]["values"]["m_dot"]) for spin in PRIMARY)
        checks.append(dict(name=f"Hellmann-Feynman h={step:g}", maximum_error=error,
                           tolerance=1e-7, passed=error < 1e-7))
    for label in ("atomic", "chain", "zero", "reverse"):
        for backend in ("ed", "qp", "dmrg") if label in {"atomic", "chain"} else ("ed",):
            checks.append(dict(name=f"{label}-{backend} primary residual/uniqueness qualification",
                               passed=results[f"{label}-{backend}"]["qualified"]))
            checks.append(dict(name=f"{label}-{backend} all control-root residuals and orthogonality",
                               passed=results[f"{label}-{backend}"]["diagnostics"]["all_returned_states_converged"]))
    return dict(checks=checks, passed=all(row["passed"] for row in checks))


def worker(directory):
    """Foreground child owning its process group, including licensed descendants."""
    request = json.loads((directory / "request.json").read_text())
    parent = json.loads((directory / "launch.json").read_text())["supervisor_pid"]

    def terminate(signum, frame):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        os.killpg(os.getpgrp(), signal.SIGTERM)
        # This fallback must survive loss of the supervising parent itself.
        time.sleep(2)
        os.killpg(os.getpgrp(), signal.SIGKILL)

    if os.getpgrp() != os.getpid():
        raise RuntimeError("worker must own its process group")
    signal.signal(signal.SIGTERM, terminate)
    if ctypes.CDLL(None, use_errno=True).prctl(1, signal.SIGTERM, 0, 0, 0) != 0 or os.getppid() != parent:
        terminate(signal.SIGTERM, None)
    task, physical = request["task"], request["physical"]
    try:
        if task["backend"] == "nrg":
            deck = render_param(physical, task["numerical"], task.get("units", "gap"))
            if (directory / "param").read_text() != deck:
                raise ValueError("prepared NRG input changed")
            tools = request["tools"]
            scratch = directory / "scratch"
            scratch.mkdir()
            environment = dict(os.environ, NRG_WORKDIR=str(scratch), TMPDIR=str(scratch),
                               MKL_THREADING_LAYER="INTEL", OMP_NUM_THREADS="4", MKL_NUM_THREADS="4",
                               OPENBLAS_NUM_THREADS="1", MKL_DYNAMIC="FALSE", OMP_DYNAMIC="FALSE")
            with lock(Path("/tmp") / f"qdjj-mathematica-{os.getuid()}.lock"):
                active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "MathKernel|WolframKernel|Mathematica"],
                                        capture_output=True)
                if active.returncode != 1:
                    raise RuntimeError("Mathematica already active or process check failed; no process was stopped")
                with (directory / "nrginit.log").open("w") as stream:
                    subprocess.run([tools["nrginit"], tools["math"]], cwd=directory, env=environment,
                                   stdout=stream, stderr=subprocess.STDOUT, check=True)
                if not (directory / "data").is_file() or "Success!" not in (directory / "nrginit.log").read_text():
                    raise RuntimeError("nrginit did not produce data and Success!")
            with (directory / "nrg.log").open("w") as stream:
                subprocess.run([tools["nrg"]], cwd=directory, env=environment,
                               stdout=stream, stderr=subprocess.STDOUT, check=True)
            if not (directory / "DONE").is_file():
                raise RuntimeError("NRG did not produce DONE")
            result = extract_result(directory)
            result["artifacts"] = {name: sha256(directory / name) for name in
                                   ("param", "data", "raw.h5", "DONE", "nrg.log", "nrginit.log")}
        else:
            from NRG_comparisons.test3.finite import run_finite
            result = run_finite(task, physical, request["bath"])
        write_json(directory / "result.json", dict(result, status="completed", request_sha256=digest(request)))
    except Exception as error:
        write_json(directory / "result.json", dict(status="failed", error=f"{type(error).__name__}: {error}"))
        raise


def execute(directory, limits, remaining):
    timeout = min(limits["task_seconds"], remaining)
    if timeout <= 0:
        raise RuntimeError("pilot campaign budget exhausted")

    def constrain():
        amount = int(limits["memory_gib"] * 1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (amount, amount))
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:limits["cpus"]])

    write_json(directory / "launch.json", dict(supervisor_pid=os.getpid(), timeout_seconds=timeout))
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONDONTWRITEBYTECODE="1")
    start, peak, failure = time.monotonic(), 0, None
    with (directory / "worker.log").open("w") as log:
        child = subprocess.Popen([sys.executable, "-B", str(HERE / "run.py"), "worker", str(directory)],
                                 stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, env=environment,
                                 start_new_session=True, preexec_fn=constrain)
        try:
            while child.poll() is None:
                peak = max(peak, process_memory(child.pid))
                if peak > limits["memory_gib"]*1024**3 or time.monotonic()-start >= timeout:
                    failure = "memory/time resource limit reached"
                    break
                time.sleep(.25)
        finally:
            for signum in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(child.pid, signum)
                except ProcessLookupError:
                    pass
                try:
                    child.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    pass
            child.wait()
    path = directory / "result.json"
    result = json.loads(path.read_text()) if path.exists() else {}
    if failure or child.returncode or not result:
        result.update(status="failed", error=failure or result.get("error", f"worker exit {child.returncode}"))
    result.update(monitored_seconds=time.monotonic()-start, peak_rss_gib=peak/1024**3)
    write_json(path, result)
    return result


def campaign(stage, directory):
    directory = directory.resolve()
    if not directory.is_relative_to(HERE / "runs"):
        raise ValueError("fresh calculation directory must be inside test3/runs")
    inputs = json.loads((HERE / "input/physical.json").read_text())
    limits = inputs["limits"]
    if limits != dict(cpus=4, memory_gib=16, task_seconds=1200, campaign_seconds=7200, mathematica_instances=1):
        raise ValueError("pilot resource ceilings must not be changed without a new reviewed profile")
    directory.mkdir(parents=True, exist_ok=True)
    with lock(directory / ".campaign.lock"):
        tools = {name: shutil.which(os.environ.get(variable, name)) for name, variable in
                 (("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH"))}
        if any(path is None for path in tools.values()):
            raise RuntimeError("load NRG/Mathematica environment or set NRGINIT, NRG, MATH")
        active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "nrg|nrginit|MathKernel|WolframKernel|Mathematica"],
                                capture_output=True)
        if active.returncode != 1:
            raise RuntimeError("an external numerical task is active, or process check failed; nothing stopped")
        identity = dict(inputs=inputs, tasks={s: tasks(s) for s in ("controls", "pilot")},
                        tools={k: dict(path=v, sha256=sha256(v)) for k, v in tools.items()},
                        source_hashes={str(p.relative_to(ROOT)): sha256(p) for p in
                                       (HERE / "run.py", HERE / "finite.py", HERE / "manual_run/nrg.py",
                                        HERE / "manual_run/physics.py", ROOT / "NRG_comparisons/test1/run.py")},
                        baths_sha256=sha256(HERE / "input/baths.json"))
        initializer = Path(tools["nrginit"]).resolve().parent
        sources = list(initializer.glob("*.m")) + list((initializer.parent / "lib64").glob("libnrgljubljana*.so.*"))
        sources += [p for p in (Path.home() / "nrg/custommodels.m", Path.home() / "nrg/customoperators.m") if p.is_file()]
        identity["external_sources"] = {str(p.resolve()): sha256(p) for p in sorted(sources)}
        manifest_path = directory / "campaign.json"
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            if manifest["identity"] != identity:
                raise ValueError("campaign code/inputs/tools changed; preserve it and use a fresh campaign")
        else:
            from qdjj_solver.common.io import numerical_environment, source_provenance
            # A fresh directory after a code fix does not renew this pilot's budget.
            starts = [json.loads(p.read_text())["started_unix"] for p in directory.parent.glob("*/campaign.json")]
            manifest = dict(identity=identity, started_unix=min([time.time(), *starts]), environment=numerical_environment(),
                            solver=source_provenance())
            write_json(manifest_path, manifest)
        results = {}
        for task in tasks("controls") + (tasks("pilot") if stage == "pilot" else []):
            if task["id"] == "cosh8-qp-q2":
                gate = check_controls(results)
                write_json(directory / "controls.json", gate)
                if not gate["passed"]:
                    raise RuntimeError("finite-chain controls failed; pilot escalation blocked")
            physical = inputs["physical"] | task.get("physical_changes", {})
            request = dict(task=task, physical=physical, tools=tools)
            if task["backend"] != "nrg":
                request["bath"] = (chain_bath(results[task["bath_source"]]) if "bath_source" in task
                                   else json.loads((HERE / "input/baths.json").read_text())[task["bath"]])
            work = directory / task["id"]
            if work.exists():
                if json.loads((work / "request.json").read_text()) != request:
                    raise ValueError("prepared request changed")
                result = json.loads((work / "result.json").read_text())
                if result["status"] != "completed" or result["request_sha256"] != digest(request):
                    raise ValueError(f"{task['id']} previously failed; no automatic retry")
                for name, value in result.get("artifacts", {}).items():
                    if sha256(work / name) != value:
                        raise ValueError("raw NRG artifact changed")
            else:
                work.mkdir()
                write_json(work / "request.json", request)
                if task["backend"] == "nrg":
                    (work / "param").write_text(render_param(physical, task["numerical"], task.get("units", "gap")))
                print(f"Starting {task['id']}", flush=True)
                result = execute(work, limits, limits["campaign_seconds"]-(time.time()-manifest["started_unix"]))
                print(f"{task['id']}: {result['status']}, {result['monitored_seconds']:.1f}s, "
                      f"{result['peak_rss_gib']:.3f} GiB", flush=True)
                if result["status"] != "completed":
                    raise RuntimeError(f"{task['id']} failed: {result['error']}; see {work}/worker.log")
            results[task["id"]] = result
        gate = check_controls(results)
        write_json(directory / "controls.json", gate)
        print(f"Controls passed: {gate['passed']}. Production has not been launched.", flush=True)
        return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("controls", "pilot", "worker"))
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    if args.stage == "worker":
        worker(args.directory.resolve())
    else:
        campaign(args.stage, args.directory)


if __name__ == "__main__":
    main()
