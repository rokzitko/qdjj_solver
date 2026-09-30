#!/usr/bin/env python3
"""Bounded, resumable static comparison; external NRG installation is read-only."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def load_inputs():
    manifest = json.loads((HERE / "input/physical.json").read_text())
    profiles = json.loads((HERE / "input/profiles.json").read_text())
    physical, limits = manifest["physical"], manifest["limits"]
    if manifest["format_version"] != 1 or physical["geometry"] != "single":
        raise ValueError("this runner supports only the version-1 single-reservoir manifest")
    if physical["detuning"] != 0 or physical["field"] != 0:
        raise ValueError("test1 requires zero detuning and field; later cases are separate")
    for key in ("gap", "bandwidth", "u", "gamma", "detuning", "field"):
        if isinstance(physical[key], bool) or not isinstance(physical[key], (int, float)) or not math.isfinite(physical[key]):
            raise ValueError(f"physical {key} must be finite and real")
    if min(physical["gap"], physical["bandwidth"]) <= 0 or min(physical["u"], physical["gamma"]) < 0:
        raise ValueError("positive energy scales and nonnegative U/Gamma required")
    if type(limits["cpus"]) is not int:
        raise ValueError("cpus must be an integer")
    for key, maximum in (("cpus", 16), ("memory_gib", 128), ("campaign_hours", 48), ("task_hours", 48)):
        if isinstance(limits[key], bool) or not 0 < limits[key] <= maximum:
            raise ValueError(f"{key} exceeds the authorized budget")
    if not math.isfinite(limits["minimum_lambda"]) or limits["minimum_lambda"] < 1.8 or limits["mathematica_instances"] != 1:
        raise ValueError("Lambda >= 1.8 and one Mathematica instance are mandatory")
    return manifest, profiles


def make_tasks(profile, specification, physical):
    """Explicit independent axes, not a Cartesian product of all NRG cutoffs."""
    tasks = []
    if profile in ("fulltest1", "nrg-refinement"):
        tasks = make_tasks("controls", {"wilson_nmax": specification.get("wilson_nmax", [2, 4])}, physical)
        if profile == "nrg-refinement":
            tasks = [t for t in tasks if t["backend"] == "nrg"]
        for entry in specification["bath_ladders"] if profile == "fulltest1" else []:
            bath = {k: entry[k] for k in ("kind", "levels")}
            if bath["kind"] == "surrogate":
                bath.update(frequency_cutoff=100 * physical["gap"], relative_weight=1.0)
            tasks.extend(dict(backend="qp", bath=bath, cutoff=q,
                              max_dimension=specification["qp_max_dimension"]) for q in entry["cutoffs"])
            tasks.extend(dict(backend="dmrg", bath=bath, chi=chi, save_states=True,
                              lanczos_probability_tolerance=specification["lanczos_probability_tolerance"])
                         for chi in entry["bonds"])
        for levels in specification["window_levels"] if profile == "fulltest1" else []:
            for window in specification["window_cutoffs_over_gap"]:
                bath = dict(kind="surrogate", levels=levels, frequency_cutoff=window * physical["gap"], relative_weight=1.0)
                tasks.extend(dict(backend="dmrg", bath=bath, chi=chi, save_states=True,
                                  lanczos_probability_tolerance=specification["lanczos_probability_tolerance"])
                             for chi in specification["window_bonds"])
        anchors, probes, remainder, unit_checks = [], [], [], []
        for lam in specification["lambdas"]:
            if not math.isfinite(lam) or lam < 1.8:
                raise ValueError("Lambda < 1.8 is not authorized")
            a = (1 - 1 / lam) / math.log(lam)
            target = specification["terminal_scale_over_gap"] * physical["gap"]
            length = 1 + 2 * (1 - min(specification["z_values"]) + math.log(physical["bandwidth"] * a / target) / math.log(lam))
            base = 2 * math.ceil(length / 2)
            for z in specification["z_values"]:
                anchor = dict(z=z, nmax=base + max(specification["chain_extensions"]),
                              keep=max(specification["keeps"]), keepenergy=max(specification["energy_cutoffs"]),
                              keepmin=max(specification["keepmins"]), untruncated=False, **{"lambda": lam})
                anchors.append(dict(backend="nrg", numerical=anchor, units="gap", clean=False))
                settings = [anchor | {"keep": k} for k in specification["keeps"][:-1]]
                settings += [anchor | {"keepenergy": e} for e in specification["energy_cutoffs"][:-1]]
                settings += [anchor | {"keepmin": k} for k in specification["keepmins"][:-1]]
                settings += [anchor | {"nmax": base + n} for n in specification["chain_extensions"][:-1]]
                selected = probes if z in (min(specification["z_values"]), 1.0) else remainder
                selected.extend(dict(backend="nrg", numerical=n, units="gap", clean=False) for n in settings)
                if z == 1:
                    unit_checks.append(dict(backend="nrg", numerical=anchor, units="bandwidth", clean=False))
        if profile == "nrg-refinement":
            # Finish anchors and their unit checks before any truncation axes.
            ordered = []
            for lam in specification["lambdas"]:
                ordered.extend(t for t in anchors + unit_checks if t["numerical"]["lambda"] == lam)
            refinements = [t for lam in specification["lambdas"] for t in probes + remainder if t["numerical"]["lambda"] == lam]
            return tasks + ordered + refinements
        return tasks + anchors + unit_checks + probes + remainder
    if profile == "controls":
        for nmax in specification["wilson_nmax"]:
            bath = dict(kind="wilson", levels=nmax + 1, nmax=nmax, **{"lambda": 2.0})
            tasks.extend([dict(backend="ed", bath=bath), dict(backend="qp", bath=bath, cutoff=None),
                          dict(backend="dmrg", bath=bath, chi=128)])
            numerical = dict(z=1.0, nmax=nmax, keep=4096, keepenergy=-1, keepmin=0,
                             untruncated=True, **{"lambda": 2.0})
            for units in ("gap", "bandwidth"):
                for clean in (False, True):
                    tasks.append(dict(backend="nrg", bath=bath, numerical=numerical, units=units, clean=clean))
        return tasks
    for kind in ("cosh", "surrogate"):
        for levels in specification[f"{kind}_levels"]:
            bath = dict(kind=kind, levels=levels)
            if kind == "surrogate":
                # Broad physical fitting window: absolute energies are targets too.
                bath.update(frequency_cutoff=100 * physical["gap"], relative_weight=1.0)
            for cutoff in specification["cutoffs"]:
                tasks.append(dict(backend="qp", bath=bath, cutoff=cutoff))
            if levels <= 8:
                tasks.append(dict(backend="qp", bath=bath, cutoff=None))
            for chi in specification["bonds"]:
                task = dict(backend="dmrg", bath=bath, chi=chi)
                if "lanczos_probability_tolerance" in specification:
                    task["lanczos_probability_tolerance"] = specification["lanczos_probability_tolerance"]
                tasks.append(task)
    for lam in specification["lambdas"]:
        if not math.isfinite(lam) or lam < 1.8:
            raise ValueError("Lambda < 1.8 is not authorized")
        a = (1 - 1 / lam) / math.log(lam)
        target = specification["terminal_scale_over_gap"] * physical["gap"]
        # Same length for every twist, set by the slowest-decaying z in the grid.
        length = 1 + 2 * (1 - min(specification["z_values"]) + math.log(physical["bandwidth"] * a / target) / math.log(lam))
        base = 2 * math.ceil(length / 2)
        largest = base + max(specification["chain_extensions"])
        for z in specification["z_values"]:
            anchor = dict(z=z, nmax=largest, keep=max(specification["keeps"]),
                          keepenergy=max(specification["energy_cutoffs"]), keepmin=specification["keepmin"],
                          untruncated=False, **{"lambda": lam})
            settings = [anchor | {"keep": k} for k in specification["keeps"]]
            settings += [anchor | {"keepenergy": e} for e in specification["energy_cutoffs"]]
            settings += [anchor | {"nmax": base + n} for n in specification["chain_extensions"]]
            unique = {digest(n): n for n in settings}
            for numerical in unique.values():
                tasks.append(dict(backend="nrg", numerical=numerical, units="gap", clean=False))
            if z == 1:
                tasks.append(dict(backend="nrg", numerical=anchor, units="bandwidth", clean=False))
    return tasks


def construct_bath(physical, specification):
    from qdjj_solver.common.baths import cosh_grid, fit_surrogate
    if specification["kind"] == "wilson":
        from NRG_comparisons.test1.control import wilson_bath
        return wilson_bath(physical, specification["lambda"], specification["nmax"])
    if specification["kind"] == "cosh":
        if specification["levels"] % 2:
            raise ValueError("cosh levels must be even")
        bath = cosh_grid(specification["levels"] // 2, delta=physical["gap"], bandwidth=physical["bandwidth"])
    elif specification["kind"] == "surrogate":
        extra = {name: specification[name] for name in ("max_nfev",) if name in specification}
        bath = fit_surrogate(specification["levels"], delta=physical["gap"], bandwidth=physical["bandwidth"],
                             frequency_cutoff=specification["frequency_cutoff"], frequency_min=1e-3 * physical["gap"],
                             relative_weight=specification["relative_weight"], starts=4, seed=1729, **extra)
    else:
        raise ValueError("unknown bath kind")
    return bath.record()


def flatten(branches):
    values = {f"{branch}.{name}": state[name] for branch, state in branches.items()
              for name in ("energy", "P0", "P1", "P2", "moment") if name in state}
    values["signed_gap"] = branches["doublet"]["energy"] - branches["singlet"]["energy"]
    return values


def run_internal(task, physical, directory):
    import numpy as np
    from qdjj_solver.common.baths import DiscreteBath
    from qdjj_solver.common.models import reference_model
    from qdjj_solver.common.problem import Sector
    from qdjj_solver.common.io import numerical_environment, source_provenance

    payload = json.loads((directory / "request.json").read_text()) if (directory / "request.json").exists() else {}
    expected_bath = payload.get("bath_sha256")
    if not (directory / "bath.json").exists():
        cache = directory.parents[2] / ("bath-" + digest(task["bath"])[:16] + ".json")
        if not cache.exists():
            if expected_bath is not None:
                raise ValueError("pinned bath is missing; refusing to refit")
            write_json(cache, construct_bath(physical, task["bath"]))
        shutil.copyfile(cache, directory / "bath.json")
    if expected_bath is not None and hashlib.sha256((directory / "bath.json").read_bytes()).hexdigest() != expected_bath:
        raise ValueError("pinned bath changed after task preparation")
    bath_record = json.loads((directory / "bath.json").read_text())
    bath = DiscreteBath.from_record(bath_record)
    if task["backend"] == "ed":
        from NRG_comparisons.test1.control import electron_ed
        result = electron_ed(physical, bath_record)
        return dict(values=flatten(result["branches"]), finite_converged=True, diagnostics=result,
                    environment=numerical_environment(), implementation=source_provenance())
    model = reference_model(bath, u=physical["u"], gamma=physical["gamma"], detuning=physical["detuning"],
                            field=physical["field"], geometry="single", phi=0, symmetry=False, compress=False)
    # Validate branch labels using the full physical total-spin operator.
    from qdjj_solver.common.algebra import create, annihilate, number
    modes = 2 + 2 * bath.levels
    sz = sum((number(i) - number(i + 1)) / 2 for i in range(0, modes, 2))
    raising = sum(create(i) * annihilate(i + 1) for i in range(0, modes, 2))
    model.observables["total_spin_squared"] = model.physical_operator(raising.dagger() * raising + sz * sz + sz)
    if task["backend"] == "qp":
        from qdjj_solver.qp_solver.solver import solve, SolverOptions
        options = SolverOptions(method="sparse", eigenpairs=2, tolerance=1e-12, residual_tolerance=1e-8,
                                max_dimension=task.get("max_dimension", 8_000_000), max_memory_gib=96, max_nnz=500_000_000,
                                ncv=30, threads=1)
    else:
        from qdjj_solver.dmrg_solver.solver import solve, SolverOptions
        options = SolverOptions(chi_max=task["chi"], eigenpairs=2,
                                seed_trials=3, max_sweeps=80, min_sweeps=12,
                                residual_tolerance=1e-7, energy_tolerance=1e-12, threads=1,
                                lanczos_probability_tolerance=task.get("lanczos_probability_tolerance", 1e-14))
    branches, diagnostics = {}, {}
    finite = True
    for branch, parity, spin in (("singlet", 0, 0), ("doublet", 1, 1)):
        kwargs = {"cutoff": task["cutoff"]} if task["backend"] == "qp" else {}
        if task["backend"] == "dmrg" and (directory / "request.json").exists():
            initial = json.loads((directory / "request.json").read_text()).get("initial", {}).get(branch)
            if initial is not None:
                from qdjj_solver.dmrg_solver.io import load_result
                path = Path(initial["path"])
                if hashlib.sha256(path.read_bytes()).hexdigest() != initial["sha256"]:
                    raise ValueError("continuation checkpoint JSON changed after task preparation")
                previous = load_result(path)
                if previous.hamiltonian.fingerprint() != model.fingerprint():
                    raise ValueError("benchmark continuation requires the identical finite Hamiltonian")
                kwargs["initial"] = previous
        state = solve(model, sector=Sector(parity, spin), options=options, **kwargs)
        if task["backend"] == "dmrg" and task.get("save_states"):
            from qdjj_solver.dmrg_solver.io import save_result
            save_result(state, directory / f"{branch}-state.json", save_states=True)
        obs = {name: float(np.real(values[0])) for name, values in state.observables.items()}
        charge, double = obs["impurity_charge"], obs["double_occupancy_0"]
        branches[branch] = dict(energy=float(state.energies[0] / physical["gap"]),
                                P0=1 - charge + double, P1=charge - 2 * double, P2=double,
                                moment=obs["impurity_spin_z"])
        spin_error = abs(obs["total_spin_squared"] - (0 if spin == 0 else .75))
        converged = bool(np.max(state.residuals) <= options.residual_tolerance and spin_error <= 1e-6)
        degenerate = len(state.energies) > 1 and abs(state.energies[1] - state.energies[0]) <= 1e-9 * physical["gap"]
        if degenerate:
            for name in ("P0", "P1", "P2", "moment"):
                branches[branch][name] = None
            converged = False
        if task["backend"] == "dmrg":
            converged &= state.metadata["finite_problem_converged"]
        finite &= converged
        diagnostics[branch] = dict(residual_over_gap=float(state.residuals[0] / physical["gap"]),
                                   residuals_over_gap=(state.residuals / physical["gap"]).tolist(),
                                   lowest_roots_over_gap=(state.energies / physical["gap"]).tolist(),
                                   degenerate_minimum=bool(degenerate),
                                   total_spin_squared=obs["total_spin_squared"], finite_converged=bool(converged),
                                   solver_metadata=state.metadata, timings=state.timings)
        if task["backend"] == "qp":
            diagnostics[branch]["qp_weights"] = state.qp_weights[0].tolist()
            diagnostics[branch]["dimension"] = state.basis.dimension
    return dict(values=flatten(branches), finite_converged=bool(finite), diagnostics=diagnostics,
                environment=numerical_environment(), implementation=source_provenance())


@contextmanager
def lock(path):
    import fcntl
    with Path(path).open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"another process owns {path}") from None
        yield


def executable(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"required executable not found: {name}")
    return str(Path(path).resolve())


def run_nrg(task, physical, directory):
    from NRG_comparisons.test1.nrg import render_param

    deck = render_param(physical, task["numerical"], task["units"], task["clean"])
    if (directory / "param").read_text() != deck:
        raise ValueError("generated param has changed since task preparation")
    initializer, solver = executable("nrginit"), executable("nrg")
    kernel = executable("math")
    environment = dict(os.environ, MKL_THREADING_LAYER="INTEL", MKL_NUM_THREADS="4", OMP_NUM_THREADS="4",
                       OPENBLAS_NUM_THREADS="1", MKL_DYNAMIC="FALSE")
    scratch = directory / "scratch"
    scratch.mkdir(exist_ok=True)
    environment.update(NRG_WORKDIR=str(scratch), TMPDIR=str(scratch))
    # All our kernels share this lock; pre-existing non-cooperating kernels are
    # detected, never killed. No settings or files in the NRG install are changed.
    with lock(Path("/tmp") / f"qdjj-mathematica-{os.getuid()}.lock"):
        active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "MathKernel|WolframKernel|Mathematica"],
                                capture_output=True, text=True)
        if active.returncode != 1:
            raise RuntimeError("Mathematica is already running, or kernel detection failed; retry after it exits")
        with (directory / "nrginit.log").open("w") as log:
            subprocess.run([initializer, kernel], cwd=directory, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        if not (directory / "data").is_file() or "Success!" not in (directory / "nrginit.log").read_text():
            raise RuntimeError("initializer did not produce fresh data and Success!; exit status alone is insufficient")
    with (directory / "nrg.log").open("w") as log:
        subprocess.run([solver], cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    if not (directory / "DONE").is_file():
        raise RuntimeError("NRG did not produce a fresh DONE marker")
    result = postprocess_nrg(task, physical, directory)
    artifacts = {}
    for name in ("param", "data", "nrginit.log", "nrg.log", "raw.h5"):
        with (directory / name).open("rb") as handle:
            artifacts[name] = hashlib.file_digest(handle, "sha256").hexdigest()
    return dict(result, artifacts=artifacts)


def postprocess_nrg(task, physical, directory):
    """Read an existing finite NRG flow without writes or external executables."""
    import numpy as np
    from NRG_comparisons.test1.nrg import extract_result

    result = extract_result(directory, physical, task["numerical"], task["units"])
    chain = result["wilson_chain"]
    if chain is None:
        raise ValueError("absolute energies require the generated Wilson chain")
    diagonal = np.array(chain["zeta_over_gap"])
    hopping = np.array(chain["xi_over_gap"][:-1])
    if np.max(abs(np.array(chain["delta_over_gap"]) - 1)) > 1e-12 or np.max(abs(np.array(chain["kappa_over_gap"]))) > 1e-12:
        raise ValueError("analytic vacuum subtraction requires uniform unit gap and no anomalous hopping")
    xi = np.linalg.eigvalsh(np.diag(diagonal) + np.diag(hopping, 1) + np.diag(hopping, -1))
    # NRG bath onsite terms are centered. There is no sum(xi) contribution.
    vacuum = -float(np.hypot(xi, 1).sum())
    result["analytic_bath_vacuum_over_gap"] = vacuum
    result["bath_subtracted_ground_energy"] = result["raw_total_energy_over_gap"] - vacuum
    branches = {}
    nondegenerate = True
    for name, branch in result["branches"].items():
        branches[name] = {key: branch[key] for key in ("P0", "P1", "P2", "moment") if key in branch}
        branches[name]["energy"] = branch["raw_total_energy_over_gap"] - vacuum - physical["detuning"] / physical["gap"]
        if branch.get("degenerate_minima", 1) > 1 and not task["clean"]:
            nondegenerate = False
            for observable in ("P0", "P1", "P2", "moment"):
                branches[name][observable] = None
    if task["clean"]:
        values = {"clean.vacuum_error": result["bath_subtracted_ground_energy"]}
    else:
        values = flatten(branches)
    return dict(values=values, finite_converged=nondegenerate, diagnostics=result,
                finite_convergence_meaning="completed finite NRG flow, not truncation convergence")


def worker(directory):
    directory = Path(directory).resolve()
    payload = json.loads((directory / "request.json").read_text())
    task, physical = payload["task"], payload["physical"]
    start = time.monotonic()
    launch = directory / "launch.json"
    if launch.exists():
        # Linux delivers SIGTERM if our supervisor dies, including SIGKILL.
        # Kill our owned session as well, so no licensed kernel survives us.
        import ctypes

        def terminate(signum, frame):
            signal.signal(signal.SIGTERM, signal.SIG_IGN)
            write_json(directory / "result.json", dict(task=task, physical=physical, status="failed",
                       error="worker or supervisor terminated", elapsed_seconds=time.monotonic() - start))
            os.killpg(os.getpgrp(), signal.SIGTERM)
            time.sleep(1)
            os.killpg(os.getpgrp(), signal.SIGKILL)

        if os.getpgrp() != os.getpid():
            raise RuntimeError("supervised worker must own its process group")
        signal.signal(signal.SIGTERM, terminate)
        if ctypes.CDLL(None, use_errno=True).prctl(1, signal.SIGTERM, 0, 0, 0) != 0:
            raise RuntimeError("cannot install parent-death protection")
        if os.getppid() != json.loads(launch.read_text())["supervisor_pid"]:
            terminate(signal.SIGTERM, None)
    try:
        if payload.get("nrg_refinement"):
            if task["backend"] != "nrg":
                raise ValueError("NRG-only requests cannot run internal solvers")
            wait_for_nrg_slot()
        result = (run_nrg if task["backend"] == "nrg" else run_internal)(task, physical, directory)
        result.update(task=task, physical=physical, status="completed", elapsed_seconds=time.monotonic() - start)
        # Solver metadata may contain numpy scalars/arrays. Normalize explicitly.
        from qdjj_solver.common.io import json_value
        normalized = json_value(result)
        write_json(directory / "result.json", normalized)
    except Exception as error:
        write_json(directory / "result.json", dict(task=task, physical=physical, status="failed",
                   error=f"{type(error).__name__}: {error}", elapsed_seconds=time.monotonic() - start))
        raise


def process_memory(group):
    """Aggregate Linux RSS for a worker process group, including nrginit children."""
    total = 0
    page = os.sysconf("SC_PAGE_SIZE")
    for path in Path("/proc").iterdir():
        if not path.name.isdigit():
            continue
        try:
            fields = (path / "stat").read_text().rsplit(")", 1)[1].split()
            if int(fields[2]) == group:
                total += int(fields[21]) * page
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return total


def execute(directory, limits, remaining_seconds):
    import resource
    if remaining_seconds <= 0:
        result = dict(status="failed", error="campaign wall-time limit exhausted", monitored_seconds=0)
        write_json(directory / "result.json", result)
        return result
    timeout = min(limits["task_hours"] * 3600, remaining_seconds)
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONDONTWRITEBYTECODE="1")

    def constrain():
        # Address-space limit complements the aggregate-RSS watchdog. The single
        # foreground worker and its descendants never share a budget with another.
        amount = int(limits["memory_gib"] * 1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (amount, amount))
        allowed = sorted(os.sched_getaffinity(0))[:int(limits["cpus"])]
        os.sched_setaffinity(0, allowed)

    start, peak, failure, interrupted = time.monotonic(), 0, None, False
    journal = dict(supervisor_pid=os.getpid(), worker_pid=None, started_unix=time.time(), timeout_seconds=timeout)
    write_json(directory / "launch.json", journal)
    with (directory / "worker.log").open("w") as log:
        process = subprocess.Popen([sys.executable, "-B", str(HERE / "run.py"), "worker", str(directory)],
                                   cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True, preexec_fn=constrain)
        journal["worker_pid"] = process.pid
        write_json(directory / "launch.json", journal)
        try:
            while process.poll() is None:
                peak = max(peak, process_memory(process.pid))
                if peak > limits["memory_gib"] * 1024**3:
                    failure = "aggregate RSS limit exceeded"
                    break
                if time.monotonic() - start > timeout:
                    failure = "task or campaign wall-time limit exceeded"
                    break
                time.sleep(.25)
        except KeyboardInterrupt:
            failure, interrupted = "supervisor interrupted", True
        finally:
            # Clean descendants even if the launcher exits before its child.
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
    result_path = directory / "result.json"
    result = json.loads(result_path.read_text()) if result_path.exists() else {}
    if failure or process.returncode != 0 or not result:
        result.update(status="failed", error=failure or result.get("error", f"worker exit {process.returncode}"))
    result.update(peak_rss_gib=peak / 1024**3, monitored_seconds=time.monotonic() - start)
    write_json(result_path, result)
    if interrupted:
        raise KeyboardInterrupt
    return result


def wait_for_nrg_slot():
    """Wait, never kill another owned kernel/test; the worker watchdog bounds this."""
    while True:
        active = subprocess.run(["pgrep", "-u", str(os.getuid()), "-x", "ctest|MathKernel|WolframKernel|Mathematica"],
                                capture_output=True, text=True)
        if active.returncode == 1:
            return
        if active.returncode != 0:
            raise RuntimeError("ctest/Mathematica process detection failed")
        print(f"Waiting for owned ctest/Mathematica processes: {active.stdout.strip()}", flush=True)
        time.sleep(5)


def nrg_runtime(provenance):
    """Comparable NRG runtime identity, including legacy mixed-campaign records."""
    environment = provenance["environment"]
    return dict(environment={k: environment[k] for k in ("python", "numpy", "system", "release", "machine")} |
                {"packages": {k: environment["packages"][k] for k in ("numpy", "h5py")}},
                nrg_sha256=provenance["workflow"]["nrg.py"],
                **{k: provenance[k] for k in ("executables", "initializer_sources", "custom_modules")})


def provenance(backends, *, nrg_only=False):
    from qdjj_solver.common.io import numerical_environment, source_provenance
    if nrg_only and set(backends) != {"nrg"}:
        raise ValueError("NRG-only provenance requires exclusively NRG tasks")
    paths = [HERE / name for name in ("run.py", "nrg.py", "analysis.py", "full_analysis.py")] if nrg_only else HERE.glob("*.py")
    files = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    result = dict(workflow=files, environment=numerical_environment())
    if not nrg_only:
        result["solver"] = source_provenance()
    if "nrg" in backends:
        binaries = {name: executable(name) for name in ("nrginit", "nrg", "math")}
        result["executables"] = {name: dict(path=path, sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
                                 for name, path in binaries.items()}
        folder = Path(binaries["nrginit"]).parent
        result["initializer_sources"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.glob("*.m")}
        # nrginit also loads these optional research definitions; hash, don't copy.
        result["custom_modules"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in (Path.home() / "nrg").glob("custom*.m")}
    if nrg_only:
        result["environment"] = nrg_runtime(result)["environment"]
        path = Path(sys.executable).resolve()
        result["python_executable"] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return result


def checked_nrg_reuse(campaign, identity, source=None):
    """Freeze/verify references, never copy attempts or add their elapsed time twice."""
    from NRG_comparisons.test1.nrg import render_param

    runs = (HERE / "runs").resolve()

    def inside(path):
        path = path.resolve()
        if not path.is_relative_to(runs):
            raise ValueError("NRG reuse paths must stay within this case/runs")
        return path

    def sha(path):
        with inside(path).open("rb") as handle:
            return hashlib.file_digest(handle, "sha256").hexdigest()

    tasks = {digest(t): t for t in identity["tasks"]}
    if source is not None:
        if not re.fullmatch(r"[A-Za-z0-9_-]+", source) or source == campaign.name:
            raise ValueError("reuse source must be a different campaign name")
        manifest_path = inside(runs / source / "manifest.json")
        archived = json.loads(manifest_path.read_text())
        references = []
        for taskdir in sorted((manifest_path.parent / "tasks").iterdir()):
            completed = [p for p in sorted(taskdir.glob("attempt-*/result.json"))
                         if json.loads(inside(p).read_text()).get("status") == "completed"]
            if not completed:
                continue
            path = inside(completed[-1])
            result = json.loads(path.read_text())
            if result["task"]["backend"] == "nrg":
                references.append(dict(path=os.path.relpath(path, campaign), sha256=sha(path), task_sha256=digest(result["task"])))
        identity["nrg_reuse"] = dict(manifest_path=os.path.relpath(manifest_path, campaign),
                                    manifest_sha256=sha(manifest_path), provenance=archived["provenance"], results=references)
    reuse = identity.get("nrg_reuse")
    if reuse is None:
        return [], []
    if any(Path(p).is_absolute() for p in [reuse["manifest_path"], *(r["path"] for r in reuse["results"])]):
        raise ValueError("NRG reuse paths must be relative within this case/runs")
    manifest_path = inside(campaign / reuse["manifest_path"])
    if sha(manifest_path) != reuse["manifest_sha256"]:
        raise ValueError("NRG reuse source manifest changed")
    archived = json.loads(manifest_path.read_text())
    physical = identity["manifest"]["physical"]
    if any(archived["manifest"][key] != identity["manifest"][key] for key in ("format_version", "energy_unit", "physical")):
        raise ValueError("NRG reuse physical identity or units/version mismatch")
    if archived["provenance"] != reuse["provenance"] or nrg_runtime(archived["provenance"]) != nrg_runtime(identity["provenance"]):
        raise ValueError("NRG reuse runtime/renderer provenance mismatch")
    records, sources, seen = [], [], set()
    for reference in reuse["results"]:
        path = inside(campaign / reference["path"])
        if not path.is_relative_to(manifest_path.parent / "tasks") or sha(path) != reference["sha256"]:
            raise ValueError("NRG reused result path or SHA changed")
        result = json.loads(path.read_text())
        key = reference["task_sha256"]
        task = result["task"]
        if (key in seen or key not in tasks or task != tasks[key] or task["backend"] != "nrg"
                or task not in archived["tasks"] or result.get("status") != "completed"
                or result["physical"] != physical or "monitored_seconds" not in result):
            raise ValueError("NRG reuse task/physical/completion identity mismatch")
        request = json.loads(inside(path.parent / "request.json").read_text())
        if request["task"] != task or request["physical"] != physical:
            raise ValueError("NRG reuse request mismatch")
        deck = render_param(physical, task["numerical"], task["units"], task["clean"])
        if inside(path.parent / "param").read_bytes() != deck.encode():
            raise ValueError("NRG reuse rendered param mismatch")
        artifacts = result["artifacts"]
        if not {"param", "data", "nrginit.log", "nrg.log", "raw.h5"} <= artifacts.keys():
            raise ValueError("NRG reuse missing artifact hashes")
        for name, checksum in artifacts.items():
            if Path(name).name != name or sha(path.parent / name) != checksum:
                raise ValueError(f"NRG reuse artifact changed: {name}")
        seen.add(key)
        records.append(result)
        sources.append(os.path.relpath(path, campaign))
    return records, sources


def summarize(campaign, targets):
    from NRG_comparisons.test1.analysis import analyze
    records, sources = [], []
    identity = json.loads((campaign / "manifest.json").read_text())
    archived = identity["manifest"]
    if targets != archived["targets"]:
        raise ValueError("report targets must match the archived campaign manifest")
    records, sources = checked_nrg_reuse(campaign, identity) if "nrg_reuse" in identity else ([], [])
    reused_count = len(records)
    for taskdir in sorted((campaign / "tasks").iterdir()):
        attempts = sorted(taskdir.glob("attempt-*/result.json"))
        completed = [p for p in attempts if json.loads(p.read_text()).get("status") == "completed"]
        if completed:
            path = completed[-1]
            if any(r["task"] == json.loads(path.read_text())["task"] for r in records):
                raise ValueError("duplicate local/reused task result")
            records.append(json.loads(path.read_text()))
            if records[-1]["physical"] != archived["physical"]:
                raise ValueError("result physical identity does not match campaign manifest")
            sources.append(str(path.relative_to(campaign)))
    if identity.get("policy") == "fulltest1-v1":
        from NRG_comparisons.test1.full_analysis import assess
        by_task = {digest(r["task"]): i for i, r in enumerate(records)}
        inventory = []
        for task in identity["tasks"]:
            identifier = digest(task)
            index = by_task.get(identifier)
            attempted = list((campaign / "tasks" / (task["backend"] + "-" + identifier[:16])).glob("attempt-*/result.json"))
            inventory.append(dict(task=task, record=index,
                                  status="completed" if index is not None else "failed" if attempted else "unstarted"))
        report = assess(records, targets, inventory)
        if identity.get("profile") == "nrg-refinement":
            report["inventory"] = inventory
    else:
        report = analyze(records, targets)
    report["record_sources"] = sources
    report["records"] = records
    if identity.get("profile") == "nrg-refinement":
        report["reused_references"] = [dict(record=i, **identity["nrg_reuse"]["results"][i])
                                       for i in range(reused_count)]
    report["failed_attempts"] = [str(p.relative_to(campaign)) for p in campaign.glob("tasks/*/attempt-*/result.json")
                                 if json.loads(p.read_text()).get("status") != "completed"]
    write_json(campaign / "summary.json", report)
    lines = ["# Test1 Campaign", "", "Empirical evidence only; finite-bath agreement is not continuum convergence.", "",
             f"Completed records: {len(records)}. Failed attempts: {len(report['failed_attempts'])}.", "",
             "| Backend | Observable | Candidate | Status |", "|---|---|---:|---|"]
    for backend, quantities in report["final_reference"].items():
        for name, value in quantities.items():
            lines.append(f"| {backend} | {name} | {value['value']} | {value['status']} |")
    lines += ["", "See `summary.json` for inputs, residuals, refinement differences, unit checks, and exploratory fits.", ""]
    (campaign / "REPORT.md").write_text("\n".join(lines))
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("worker", help=argparse.SUPPRESS).add_argument("directory")
    for name in ("plan", "prepare", "run", "report", "export"):
        command = commands.add_parser(name)
        command.add_argument("--profile", choices=("controls", "pilot", "convergence", "fulltest1", "nrg-refinement"), default="controls")
        command.add_argument("--campaign", default="controls")
        command.add_argument("--backend", choices=("all", "qp", "dmrg", "ed", "nrg"), default="all")
        command.add_argument("--max-tasks", type=int, default=None)
        command.add_argument("--retry-failed", action="store_true")
        command.add_argument("--task", help="run only this exact backend-hash task identifier")
        command.add_argument("--levels", type=int, help="select one finite-bath level count")
        command.add_argument("--kind", choices=("cosh", "surrogate", "wilson"))
        command.add_argument("--lambda", dest="lambda_value", type=float)
        command.add_argument("--chi", type=int)
        command.add_argument("--reuse-from", help="freeze checked NRG result references from this case's campaign")
    args = parser.parse_args(argv)
    if args.command == "worker":
        worker(args.directory)
        return
    manifest, profiles = load_inputs()
    if args.reuse_from and args.profile != "nrg-refinement":
        raise ValueError("reuse-from is only supported by the explicit NRG-only profile")
    if args.profile == "nrg-refinement":
        manifest = dict(manifest, limits=dict(manifest["limits"]))
        for key, cap in (("cpus", 15), ("memory_gib", 126), ("campaign_hours", 47.5)):
            manifest["limits"][key] = min(manifest["limits"][key], cap)
    all_tasks = make_tasks(args.profile, profiles[args.profile], manifest["physical"])
    tasks = all_tasks
    if args.backend != "all":
        tasks = [task for task in tasks if task["backend"] == args.backend]
    for expected, field, container in ((args.levels, "levels", "bath"), (args.kind, "kind", "bath"),
                                       (args.lambda_value, "lambda", "numerical"), (args.chi, "chi", None)):
        if expected is not None:
            tasks = [t for t in tasks if (t.get(container, {}) if container else t).get(field) == expected]
    if args.task is not None:
        tasks = [task for task in tasks if task["backend"] + "-" + digest(task)[:16] == args.task]
        if not tasks:
            raise ValueError("task identifier does not match the selected profile/backend")
    if args.command == "plan":
        print(json.dumps(dict(manifest=manifest, tasks=tasks), indent=2))
        return
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.campaign):
        raise ValueError("campaign must be a simple directory name")
    if args.max_tasks is not None and args.max_tasks < 1:
        raise ValueError("max-tasks must be positive")
    runs = HERE / "runs"
    runs.mkdir(exist_ok=True)
    with lock(runs / ".campaign.lock"):
        campaign = runs / args.campaign
        if args.command in ("report", "export"):
            archived = json.loads((campaign / "manifest.json").read_text())
            report = summarize(campaign, archived["manifest"]["targets"])
            if args.command == "export":
                output = HERE / "output"
                output.mkdir(exist_ok=True)
                compact = dict(manifest=archived, uncertainty_kind="empirical_not_rigorous",
                               prerequisites_passed=report["prerequisites_passed"],
                               same_bath_checks=report["same_bath_checks"], unit_checks=report["unit_checks"],
                               final_reference=report["final_reference"],
                               records=[{k: r[k] for k in ("task", "physical", "values", "finite_converged", "peak_rss_gib",
                                                          "monitored_seconds", "artifacts") if k in r}
                                        for r in report["records"]], source_campaign=args.campaign)
                compact["analysis_sha256"] = hashlib.sha256(Path(__file__).with_name("analysis.py").read_bytes()).hexdigest()
                for saved, original in zip(compact["records"], report["records"], strict=True):
                    diagnostics = original.get("diagnostics", {})
                    saved["diagnostics"] = {
                        branch: {key: value for key, value in diagnostics.get(branch, {}).items()
                                 if key in ("residual_over_gap", "lowest_roots_over_gap", "degenerate_minimum",
                                            "residuals_over_gap", "total_spin_squared", "finite_converged", "dimension", "qp_weights",
                                            "native_residuals_over_gap", "residual_evaluator_difference", "exact_schmidt_tail_weights",
                                            "overlap_squared_with_qp", "observable_differences", "qp_energies", "qp_residuals",
                                            "norm_errors", "solve_seconds")}
                        for branch in ("singlet", "doublet") if branch in diagnostics}
                    if original["task"].get("study") == "residual":
                        for branch, values in saved["diagnostics"].items():
                            values["trials"] = [{k: trial[k] for k in ("energy", "native_residual", "vector_norm", "sweeps",
                                                                     "effective_lanczos", "last_local_updates")}
                                                for trial in diagnostics[branch]["trial_diagnostics"]]
                            values["solver_options"] = diagnostics[branch]["solver_metadata"]["options"]
                    if original["task"]["backend"] == "nrg":
                        saved["diagnostics"] = {k: diagnostics[k] for k in
                            ("wilson_chain", "analytic_bath_vacuum_over_gap", "bath_subtracted_ground_energy", "shell_convergence")
                            if k in diagnostics}
                compact["baths"] = {p.name: json.loads(p.read_text()) for p in campaign.glob("bath-*.json")}
                for source in report["record_sources"]:
                    path = (campaign / source).parent / "bath.json"
                    if path.exists():
                        bath = json.loads(path.read_text())
                        if digest(bath) not in {digest(b) for b in compact["baths"].values()}:
                            compact["baths"][f"bath-{digest(bath)[:16]}.json"] = bath
                compact["failures"] = [dict(path=p, result=json.loads((campaign / p).read_text())) for p in report["failed_attempts"]]
                compact["refinements"] = report["refinements"]
                compact["z_grids"] = report["z_grids"]
                compact["nrg_extrapolation"] = report["nrg_extrapolation"]
                write_json(output / f"{args.campaign}.json", compact)
            return
        backends = [t["backend"] for t in tasks]
        implementation = provenance(backends, nrg_only=True) if args.profile == "nrg-refinement" else provenance(backends)
        identity = dict(manifest=manifest, profile=args.profile, profile_settings=profiles[args.profile],
                        provenance=implementation, tasks=all_tasks)
        if args.profile in ("fulltest1", "nrg-refinement"):
            identity["policy"] = "fulltest1-v1"
        campaign.mkdir(exist_ok=True)
        identity_path = campaign / "manifest.json"
        reused = []
        if args.profile == "nrg-refinement":
            if identity_path.exists():
                prior = json.loads(identity_path.read_text())
                if "nrg_reuse" in prior:
                    identity["nrg_reuse"] = prior["nrg_reuse"]
                if args.reuse_from and prior.get("nrg_reuse", {}).get("manifest_path") != f"../{args.reuse_from}/manifest.json":
                    raise ValueError("frozen NRG reuse source cannot change")
                reused, _ = checked_nrg_reuse(campaign, identity)
            else:
                reused, _ = checked_nrg_reuse(campaign, identity, args.reuse_from)
        if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
            raise ValueError("campaign input or implementation changed; use a new campaign name, never mix results")
        write_json(identity_path, identity)
        (campaign / "tasks").mkdir(exist_ok=True)
        spent = 0
        for attempt in runs.glob("*/tasks/*/attempt-*"):
            result_path, journal_path = attempt / "result.json", attempt / "launch.json"
            result = json.loads(result_path.read_text()) if result_path.exists() else {}
            journal = json.loads(journal_path.read_text()) if journal_path.exists() else {}
            if "monitored_seconds" in result:
                spent += result["monitored_seconds"]
            elif journal:
                if journal["worker_pid"] and process_memory(journal["worker_pid"]):
                    raise RuntimeError(f"unfinished worker session still exists for {attempt}; do not overlap jobs")
                # Unknown-duration interrupted attempts consume their full reserved
                # allowance, not zero. Preserve them and require an explicit retry.
                spent += journal["timeout_seconds"]
                result.update(status="failed", error="interrupted or unmonitored attempt",
                              monitored_seconds=journal["timeout_seconds"], accounting="conservative reservation")
                write_json(result_path, result)
        started, count = time.monotonic(), 0
        for task in tasks:
            if any(r["task"] == task for r in reused):
                continue
            remaining = manifest["limits"]["campaign_hours"] * 3600 - spent - (time.monotonic() - started)
            if remaining <= 0 or (args.max_tasks is not None and count >= args.max_tasks):
                break
            taskdir = campaign / "tasks" / (task["backend"] + "-" + digest(task)[:16])
            taskdir.mkdir(exist_ok=True)
            write_json(taskdir / "task.json", task)
            attempts = sorted(taskdir.glob("attempt-*"))
            completed = any((p / "result.json").exists() and json.loads((p / "result.json").read_text()).get("status") == "completed"
                            for p in attempts)
            pending = [p for p in attempts if not (p / "result.json").exists() and not (p / "launch.json").exists()]
            if args.profile == "nrg-refinement" and attempts and not completed and not pending and not args.retry_failed:
                failures = [json.loads((p / "result.json").read_text()) for p in attempts]
                if any(r.get("error") not in ("aggregate RSS limit exceeded", "task or campaign wall-time limit exceeded") for r in failures):
                    raise RuntimeError(f"NRG queue paused on failed task: {taskdir.name}; explicit review/retry required")
            if completed or (attempts and not pending and not args.retry_failed):
                continue
            if pending:
                directory = pending[-1]
                if args.command == "prepare":
                    continue
            else:
                directory = taskdir / f"attempt-{len(attempts) + 1:03d}"
                directory.mkdir()
            request_path = directory / "request.json"
            if request_path.exists():
                request = json.loads(request_path.read_text())
                if request["task"] != task or request["physical"] != manifest["physical"]:
                    raise ValueError("prepared request no longer matches its frozen task")
            else:
                request = dict(task=task, physical=manifest["physical"])
                if args.profile == "nrg-refinement":
                    request["nrg_refinement"] = True
                if "bath_input" in identity:
                    request["bath_record"] = identity["bath_input"]
                if task["backend"] == "dmrg" and task.get("save_states"):
                    predecessors = []
                    for path in campaign.glob("tasks/dmrg-*/attempt-*/result.json"):
                        prior = json.loads(path.read_text())
                        settings = prior.get("task", {})
                        if (prior.get("status") == "completed" and settings.get("chi", 0) < task["chi"]
                                and {k: v for k, v in settings.items() if k != "chi"} == {k: v for k, v in task.items() if k != "chi"}
                                and all((path.parent / f"{b}-state.json").is_file() for b in ("singlet", "doublet"))):
                            predecessors.append((settings["chi"], path.parent))
                    if predecessors:
                        parent = max(predecessors, key=lambda row: row[0])[1]
                        request["initial"] = {b: dict(path=str(parent / f"{b}-state.json"),
                                                      sha256=hashlib.sha256((parent / f"{b}-state.json").read_bytes()).hexdigest())
                                              for b in ("singlet", "doublet")}
                write_json(request_path, request)
            if task["backend"] == "nrg" and not (directory / "param").exists():
                from NRG_comparisons.test1.nrg import render_param
                (directory / "param").write_text(render_param(manifest["physical"], task["numerical"], task["units"], task["clean"]))
            if args.command == "prepare":
                count += 1
                continue
            print(f"Running {taskdir.name}: {json.dumps(task, sort_keys=True)}", flush=True)
            remaining = manifest["limits"]["campaign_hours"] * 3600 - spent - (time.monotonic() - started)
            if remaining <= 0:
                break
            result = execute(directory, manifest["limits"], remaining)
            print(f"  {result['status']}: {result.get('error', result.get('values', {}))}", flush=True)
            count += 1
            summarize(campaign, manifest["targets"])
            if args.profile == "nrg-refinement" and result["status"] != "completed":
                if result.get("error") not in ("aggregate RSS limit exceeded", "task or campaign wall-time limit exceeded"):
                    raise RuntimeError(f"NRG queue paused on failed task: {taskdir.name}; {result.get('error')}")
        if args.command == "run":
            report = summarize(campaign, manifest["targets"])
            print(f"Report: {campaign / 'REPORT.md'} ({len(report['records'])} completed records)")


if __name__ == "__main__":
    def interrupted(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupted)
    main()
