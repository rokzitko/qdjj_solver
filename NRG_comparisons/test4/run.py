#!/usr/bin/env python3
"""Resource-bounded ASQ model validation and two-channel cost pilots, Lambda >= 4."""

import argparse
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
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
from NRG_comparisons.test4.nrg import extract_result, render_param
from NRG_comparisons.test4.physics import OBSERVABLES, PRIMARY, PHYSICAL, paper_offset


MODELS = ("model.m", "modeloperators.m")
LIMITS = dict(cpus=4, memory_gib=16, task_seconds=1200, campaign_seconds=7200, mathematica_instances=1)


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def tasks(stage):
    if stage == "controls":
        result = []
        cases = (("atomic", {"gamma": 0.}, ("ed", "qp")),
                 ("junction", {}, ("ed", "qp", "dmrg")),
                 ("normal", {"wn": 20.}, ("ed",)),
                 ("spin", {"ws": 10.}, ("ed",)),
                 ("bridge", {"wn": 20., "ws": 10.}, ("ed", "qp")),
                 ("direct", {"gamma": 0., "wn": 20., "ws": 10.}, ("ed",)),
                 ("free", {"u": 0.}, ("ed", "qp")),
                 ("reverse", {"wn": 20., "ws": -10., "field": -.03}, ("ed",)),
                 ("phase-reverse", {"wn": 20., "ws": -10., "phi": -PHYSICAL["phi"]}, ("ed",)))
        for label, changes, backends in cases:
            result.append(dict(id=f"{label}-nrg", backend="nrg", physical_changes=changes,
                               numerical=dict(nmax=1, keep=1024, keepmin=0, keepenergy=-1.,
                                              untruncated=True, z=1., **{"lambda": 4.})))
            for backend in backends:
                result.append(dict(id=f"{label}-{backend}", backend=backend, physical_changes=changes,
                                   bath_source=f"{label}-nrg", chi=32, cutoff=None, roots=2,
                                   roots_by_spin={"0": 3}, spins=list(range(-5, 6))))
        for label, extra in (("junction-bandwidth", {"units": "bandwidth"}),
                             ("junction-sds", {"representation": "sds"})):
            result.append(dict(id=label, backend="nrg", **extra, numerical=dict(
                nmax=1, keep=1024, keepmin=0, keepenergy=-1., untruncated=True, z=1., **{"lambda": 4.})))
        for label, changes in (("junction", {}), ("bridge", {"wn": 20., "ws": 10.}),
                               ("direct", {"gamma": 0., "wn": 20., "ws": 10.})):
            for step in (2e-4, 1e-4):
                for sign in (-1, 1):
                    result.append(dict(id=f"hf-{label}-{step:g}-{sign:+d}", backend="ed",
                                       bath_source=f"{label}-nrg", roots=2, spins=[0, -1, 1],
                                       physical_changes=changes | {"phi": PHYSICAL["phi"]+sign*step}))
        return result
    if stage != "pilot":
        raise ValueError("stage must be controls or pilot; production is not authorized")
    result = [dict(id=f"nrg-k{keep}", backend="nrg", numerical=dict(
        nmax=20, keep=keep, keepmin=keep//5, keepenergy=14., z=1., **{"lambda": 4.}))
        for keep in (200, 500, 1000)]
    result += [dict(id="nrg-length22", backend="nrg", numerical=dict(
        nmax=22, keep=500, keepmin=100, keepenergy=14., z=1., **{"lambda": 4.}))]
    result += [dict(id="cosh2-qp-full", backend="qp", bath="cosh2", cutoff=None, roots=3)]
    result += [dict(id=f"cosh4-qp-q{q if q is not None else 'full'}", backend="qp", bath="cosh4",
                    cutoff=q, roots=2, roots_by_spin={"0": 3, "-2": 1, "2": 1}) for q in (2, 4, 6, None)]
    result += [dict(id=f"cosh4-dmrg-chi{chi}", backend="dmrg", bath="cosh4", chi=chi,
                    roots=2, roots_by_spin={"0": 3, "-2": 1, "2": 1}) for chi in (32, 64, 128)]
    return result


def chain_bath(result):
    import numpy as np
    from scipy.linalg import eigh_tridiagonal

    chains = result["diagnostics"]["wilson_chains"]
    for key in ("xi_over_gap", "zeta_over_gap"):
        if not np.allclose(chains[0][key], chains[1][key], rtol=0, atol=1e-12):
            raise ValueError("finite reference requires identical normal chains in both leads")
    xi, vectors = eigh_tridiagonal(chains[0]["zeta_over_gap"], chains[0]["xi_over_gap"][:-1], lapack_driver="stev")
    weights = vectors[0]**2
    if np.any(weights <= 0) or abs(weights.sum()-1) > 1e-12:
        raise ValueError("generated chain lost the normalized contact measure")
    return dict(xi=xi.tolist(), weights=weights.tolist(), delta=1., bandwidth=result["physical"]["bandwidth"],
                metadata=dict(kind="generated-nrg-chain", wilson_chains=chains, weights_renormalized=False,
                              levels_per_lead=len(xi), identical_physical_leads=True))


def check_controls(results):
    checks = []
    inventory = tasks("controls")

    def rows(taskid):
        return {(r["twice_sz"], r["root"]): r for r in results[taskid]["states"]}

    def check(name, error, tolerance):
        checks.append(dict(name=name, maximum_error=error, tolerance=tolerance, passed=error < tolerance))

    def compare(new, old):
        a, b = rows(new), rows(old)
        missing = b.keys()-a.keys()
        checks.append(dict(name=f"{new} covers {old}", passed=not missing, missing=[list(k) for k in sorted(missing)]))
        errors = []
        for key in a.keys() & b.keys():
            observables = ("energy", "energy_paper") if a[key]["degenerate"] or b[key]["degenerate"] else OBSERVABLES
            errors.extend(abs(a[key]["values"][name]-b[key]["values"][name]) for name in observables)
        check(f"{new} vs {old}", max(errors, default=math.inf), 1e-8)

    for task in inventory:
        if task["backend"] == "nrg":
            continue
        expected = {(spin, root) for spin in task["spins"] for root in range(min(
            task.get("roots_by_spin", {}).get(str(spin), task["roots"]), math.comb(10, 5+spin)))}
        missing = expected-rows(task["id"]).keys()
        checks.append(dict(name=f"{task['id']} declared inventory", passed=not missing,
                           missing=[list(k) for k in sorted(missing)]))
        checks.append(dict(name=f"{task['id']} all-root residual/orthogonality",
                           passed=results[task["id"]]["diagnostics"]["all_returned_states_converged"]))
        if task["id"].startswith("hf-"):
            continue  # Different phases: compare their derivative below, not their energies.
        if task["backend"] == "ed":
            compare(task["bath_source"], task["id"])
        else:
            compare(task["id"], task["id"].rsplit("-", 1)[0]+"-ed")
    for name in ("junction-bandwidth", "junction-sds"):
        compare(name, "junction-nrg")
    atomic = rows("atomic-nrg")
    errors = []
    for key, energy in (((0, 0), .5), ((0, 1), 1.5), ((-1, 0), -.015), ((1, 0), .015)):
        value = atomic[key]["values"]
        errors += [abs(value["energy"]-energy), abs(value["energy_paper"]-energy-paper_offset(PHYSICAL)), abs(value["current"])]
    check("atomic energies, paper constant, Zeeman sign and zero current", max(errors), 1e-9)
    for label in ("junction", "bridge", "direct"):
        for step in (2e-4, 1e-4):
            minus, plus = rows(f"hf-{label}-{step:g}--1"), rows(f"hf-{label}-{step:g}-+1")
            measured = rows(f"{label}-nrg")
            error = max(abs((plus[spin, 0]["values"]["energy"]-minus[spin, 0]["values"]["energy"])/(2*step)
                            - measured[spin, 0]["values"]["current"]) for spin in PRIMARY)
            check(f"{label} current Hellmann-Feynman h={step:g}", error, 1e-7)
    direct_current = max(abs(r["values"]["current"]) for r in rows("direct-nrg").values())
    checks.append(dict(name="bridge carries current with Gamma=0", measured=direct_current, passed=direct_current > 1e-7))
    for target, flipspin in (("reverse", True), ("phase-reverse", False)):
        original, changed = rows("bridge-nrg"), rows(target+"-nrg")
        errors = []
        for spin in PRIMARY:
            a, b = original[spin, 0]["values"], changed[-spin if flipspin else spin, 0]["values"]
            for name in ("energy", "energy_paper", "P0", "P1", "P2", "Czz", "Czz_connected", "dot_z2", "bath_z2"):
                errors.append(abs(a[name]-b[name]))
            errors.extend(abs(a[name]+b[name] if flipspin else a[name]-b[name]) for name in ("m_dot", "m_bath"))
            errors.append(abs(a["current"]-b["current"] if flipspin else a["current"]+b["current"]))
        check(f"{target} symmetry", max(errors), 1e-8)
    return dict(checks=checks, passed=all(row["passed"] for row in checks))


def worker(directory):
    request = json.loads((directory / "request.json").read_text())
    parent = json.loads((directory / "launch.json").read_text())["supervisor_pid"]

    def terminate(signum, frame):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        os.killpg(os.getpgrp(), signal.SIGTERM)
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
            deck = render_param(physical, task["numerical"], task.get("units", "gap"), task.get("representation", "paper"))
            if (directory / "param").read_text() != deck or any(sha256(directory / name) != value for name, value in request["models"].items()):
                raise ValueError("prepared model or NRG input changed")
            scratch = directory / "scratch"
            scratch.mkdir()
            env = dict(os.environ, NRG_WORKDIR=str(scratch), TMPDIR=str(scratch), MKL_THREADING_LAYER="INTEL",
                       OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="1", MKL_DYNAMIC="FALSE", OMP_DYNAMIC="FALSE")
            with lock(Path("/tmp") / f"qdjj-mathematica-{os.getuid()}.lock"):
                active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "MathKernel|WolframKernel|Mathematica"], capture_output=True)
                if active.returncode != 1:
                    raise RuntimeError("Mathematica active or process check failed; no unrelated process stopped")
                with (directory / "nrginit.log").open("w") as log:
                    subprocess.run([request["tools"]["nrginit"], request["tools"]["math"]], cwd=directory, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
                if not (directory / "data").is_file() or "Success!" not in (directory / "nrginit.log").read_text():
                    raise RuntimeError("initializer did not produce data and Success!")
            with (directory / "nrg.log").open("w") as log:
                subprocess.run([request["tools"]["nrg"]], cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            if not (directory / "DONE").is_file():
                raise RuntimeError("NRG did not produce DONE")
            result = extract_result(directory)
            result["artifacts"] = {name: sha256(directory / name) for name in (*MODELS, "param", "data", "raw.h5", "DONE", "nrg.log", "nrginit.log")}
        else:
            from NRG_comparisons.test4.finite import run_finite
            result = run_finite(task, physical, request["bath"])
        write_json(directory / "result.json", dict(result, status="completed", request_sha256=digest(request)))
    except Exception as error:
        write_json(directory / "result.json", dict(status="failed", error=f"{type(error).__name__}: {error}"))
        raise


def execute(directory, remaining):
    import resource

    timeout = min(LIMITS["task_seconds"], remaining)
    if timeout <= 0:
        raise RuntimeError("test4 pilot budget exhausted")

    def constrain():
        amount = LIMITS["memory_gib"] * 1024**3
        resource.setrlimit(resource.RLIMIT_AS, (amount, amount))
        os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:LIMITS["cpus"]])

    write_json(directory / "launch.json", dict(supervisor_pid=os.getpid(), timeout_seconds=timeout))
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    start, peak, failure = time.monotonic(), 0, None
    with (directory / "worker.log").open("w") as log:
        child = subprocess.Popen([sys.executable, "-B", str(HERE / "run.py"), "worker", str(directory)], cwd=ROOT,
                                 env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True, preexec_fn=constrain)
        try:
            while child.poll() is None:
                peak = max(peak, process_memory(child.pid))
                if peak > LIMITS["memory_gib"]*1024**3 or time.monotonic()-start >= timeout:
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
        raise ValueError("calculation directory must be inside test4/runs")
    inputs = json.loads((HERE / "input/physical.json").read_text())
    if inputs["limits"] != LIMITS or inputs["minimum_lambda"] != 4:
        raise ValueError("reviewed pilot limits and minimum Lambda must not change")
    directory.mkdir(parents=True, exist_ok=True)
    # Unlike per-campaign locks, this also prevents concurrent sibling campaigns.
    with lock(HERE / "runs/.session.lock"):
        tools = {name: shutil.which(os.environ.get(key, name)) for name, key in (("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH"))}
        if any(path is None for path in tools.values()):
            raise RuntimeError("load NRG/Mathematica or set NRGINIT, NRG, MATH")
        active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "nrg|nrginit|MathKernel|WolframKernel|Mathematica"], capture_output=True)
        if active.returncode != 1:
            raise RuntimeError("external calculation active or process check failed; nothing stopped")
        sources = [HERE / name for name in ("run.py", "finite.py", "nrg.py", "physics.py", *MODELS)]
        sources.append(ROOT / "NRG_comparisons/test1/run.py")
        identity = dict(inputs=inputs, tasks={s: tasks(s) for s in ("controls", "pilot")},
                        source_hashes={str(p.relative_to(ROOT)): sha256(p) for p in sources},
                        models={name: sha256(HERE / name) for name in MODELS},
                        tools={k: dict(path=v, sha256=sha256(v)) for k, v in tools.items()},
                        baths_sha256=sha256(HERE / "input/baths.json"))
        initializer = Path(tools["nrginit"]).resolve().parent
        external = list(initializer.glob("*.m")) + list((initializer.parent / "lib64").glob("libnrgljubljana*.so.*"))
        external += [p for p in (Path.home() / "nrg/custommodels.m", Path.home() / "nrg/customoperators.m") if p.is_file()]
        identity["external_sources"] = {str(p.resolve()): sha256(p) for p in sorted(external)}
        path = directory / "campaign.json"
        if path.exists():
            manifest = json.loads(path.read_text())
            if manifest["identity"] != identity:
                raise ValueError("code/inputs/tools changed; preserve this attempt and use a fresh campaign")
        else:
            from qdjj_solver.common.io import numerical_environment, source_provenance
            starts = [json.loads(p.read_text())["started_unix"] for p in directory.parent.glob("*/campaign.json")]
            manifest = dict(identity=identity, started_unix=min([time.time(), *starts]),
                            environment=numerical_environment(), solver=source_provenance())
            write_json(path, manifest)
        results = {}
        for task in tasks("controls") + (tasks("pilot") if stage == "pilot" else []):
            if task["id"] == tasks("pilot")[0]["id"]:
                gate = check_controls(results)
                write_json(directory / "controls.json", gate)
                if not gate["passed"]:
                    raise RuntimeError("model controls failed; cost pilots blocked")
            physical = inputs["physical"] | task.get("physical_changes", {})
            request = dict(task=task, physical=physical, tools=tools)
            if task["backend"] == "nrg":
                request["models"] = identity["models"]
            else:
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
                        raise ValueError("raw artifact or model snapshot changed")
            else:
                work.mkdir()
                write_json(work / "request.json", request)
                if task["backend"] == "nrg":
                    for name in MODELS:
                        shutil.copyfile(HERE / name, work / name)
                    (work / "param").write_text(render_param(physical, task["numerical"], task.get("units", "gap"), task.get("representation", "paper")))
                print(f"Starting {task['id']}", flush=True)
                result = execute(work, LIMITS["campaign_seconds"]-(time.time()-manifest["started_unix"]))
                print(f"{task['id']}: {result['status']}, {result['monitored_seconds']:.1f}s, {result['peak_rss_gib']:.3f} GiB", flush=True)
                if result["status"] != "completed":
                    raise RuntimeError(f"{task['id']} failed: {result['error']}; see {work}/worker.log")
            results[task["id"]] = result
        gate = check_controls(results)
        write_json(directory / "controls.json", gate)
        print(f"Controls passed: {gate['passed']}. Production not launched.", flush=True)
        return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare-inputs", "controls", "pilot", "worker"))
    parser.add_argument("directory", nargs="?", type=Path)
    args = parser.parse_args()
    if args.stage == "prepare-inputs":
        from qdjj_solver import cosh_grid
        target = HERE / "input/baths.json"
        if target.exists():
            raise ValueError("explicit baths already exist; refusing to overwrite")
        write_json(target, {f"cosh{n}": cosh_grid(n//2, delta=1., bandwidth=100.).record() for n in (2, 4)})
    elif args.directory is None:
        parser.error("directory is required")
    elif args.stage == "worker":
        worker(args.directory.resolve())
    else:
        campaign(args.stage, args.directory)


if __name__ == "__main__":
    main()
