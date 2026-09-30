"""Independent electron-chain MPS references for the bath documentation study.

Uses the existing Lanczos contact-measure transformation with the unified MPS
solver, so reference calculations have physical residuals rather than a legacy
variance subtraction. This source-checkout utility is not a new core bath kind.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import logging
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from qdjj_solver import (DiscreteBath, Hamiltonian, Impurity, Sector, annihilate,
                         create, hermitian_pair, number, solve)
from qdjj_solver.common.io import write_json
from qdjj_solver.dmrg_solver.legacy import chain_coefficients

try:
    from tools import benchmark_baths as bench
except ModuleNotFoundError:
    import benchmark_baths as bench


def chain_hamiltonian(bath, u, gamma, geometry="single", compress=False):
    """Same finite single-reservoir Anderson model, in physical chain electrons."""
    if geometry != "single" or compress or not np.isfinite(gamma) or gamma < 0:
        raise ValueError("this independent reference requires an uncompressed single reservoir")
    onsite, hopping, norm, vacuum = chain_coefficients(bath)
    if len(onsite) != bath.levels:
        raise ValueError("full-spectrum chain reference requires distinct normal levels")
    impurity = Impurity.anderson(u)
    operator = impurity.operator-vacuum
    tunnel = np.sqrt(gamma/(np.pi*bath.rho))*norm
    for j, energy in enumerate(onsite):
        up, down = 2+2*j, 3+2*j
        operator += energy*(number(up)+number(down))
        operator -= hermitian_pair(bath.delta*create(up)*create(down))
        for spin in (0, 1):
            previous = spin if j == 0 else 2*j+spin
            value = tunnel if j == 0 else hopping[j-1]
            operator += hermitian_pair(value*create(previous)*annihilate(up+spin))
    return Hamiltonian(operator, 2, (1, -1)*(bath.levels+1), observables=impurity.observables,
                       metadata=dict(bath=[bench.bath_record(bath)], u=u, gamma=gamma,
                                     energy_reference="isolated BCS reservoir subtracted",
                                     bcs_vacuum_energy=vacuum,
                                     coordinates=dict(kind="physical-electron-chain", restriction="none",
                                                      onsite=onsite.tolist(), hopping=hopping.tolist(),
                                                      contact_norm=norm)))


def point(bath, task):
    start = perf_counter()
    h = chain_hamiltonian(bath, **task["model"])
    built = perf_counter()
    options = dict(task["dmrg_options"], **task["solver"].get("options", {}))
    states = [bench.compact_state(solve(h, backend="dmrg", cutoff=None, sector=Sector(p, p), options=options))
              for p in task["sectors"]]
    finished = perf_counter()
    values = {}
    for parity, state in zip(task["sectors"], states, strict=True):
        values[f"energy_{parity}"] = state["energy"]
        values.update({f"{name}_{parity}": value for name, value in state["observables"].items()})
    values["signed_gap"] = values["energy_1"]-values["energy_0"]
    values["kappa"] = 1-2*values["impurity_spin_z_1"]
    times = dict(model=built-start, solve=finished-built, reused_bath_total=finished-start)
    for name in {key for state in states for key in state["timings_seconds"]}:
        times[f"solver_{name}"] = sum(state["timings_seconds"].get(name, 0.) for state in states)
    return dict(states=states, values=values, times=times,
                converged=all(s["finite_problem_converged"] for s in states))


def worker(request, destination):
    task = bench.read_json(request)
    if task["type"] == "bath":
        return bench.worker(request, destination)
    start = perf_counter()
    logging.getLogger("tenpy").setLevel(logging.ERROR)
    try:
        if task["solver"].get("layout") != "electron-chain" or task["solver"].get("backend") != "dmrg":
            raise ValueError("chain refinement requires the declared electron-chain DMRG workflow")
        bath = DiscreteBath.from_record(task["bath"])
        with threadpool_limits(limits=1):
            point(bath, task)
            runs = []
            for _ in range(task["repetitions"]):
                gc.collect()
                runs.append(point(bath, task))
        result = dict(status="ok", restriction="none", values=runs[0]["values"], states=runs[0]["states"],
                      converged=all(r["converged"] for r in runs),
                      repeat_value_spread={k: float(np.ptp([r["values"][k] for r in runs])) for k in runs[0]["values"]},
                      timings_seconds={k: bench.timing_summary([r["times"][k] for r in runs]) for k in runs[0]["times"]})
    except (ValueError, RuntimeError, MemoryError, ImportError) as error:
        result = dict(status="failed", exception=type(error).__name__, message=str(error))
    result.update(task=task, peak_process_rss_bytes=bench.peak_rss_bytes(), worker_seconds=perf_counter()-start)
    write_json(destination, result)


def run(config, output):
    output = Path(output)
    (output/"jobs").mkdir(parents=True, exist_ok=True)
    environment = bench.environment_record()
    identity = dict(format="qdjj-bath-chain-reference", format_version=1, input=config,
                    profile="full", include_dmrg=True,
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    helper_sha256=hashlib.sha256(Path(bench.__file__).read_bytes()).hexdigest(),
                    implementation_sha256=environment["implementation"]["sha256"])
    path = output/"manifest.json"
    if path.exists():
        manifest = bench.read_json(path)
        if any(manifest.get(k) != v for k, v in identity.items()):
            raise ValueError("resume requires identical chain input and implementation")
    else:
        manifest = dict(identity, environment=environment, completed=[], commands=[])
    manifest["commands"].append(sys.argv[1:])
    write_json(path, manifest)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", **{k: "1" for k in bench.THREAD_VARIABLES})

    def execute(task):
        key = bench.task_id(task)
        destination = output/"jobs"/f"{key}.json"
        if key not in manifest["completed"]:
            request = output/"request.json"
            write_json(request, task)
            print(task.get("case", task.get("spec")), task.get("solver", ""), flush=True)
            process = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()), "--worker",
                                      str(request), str(destination)], env=env, capture_output=True, text=True)
            if process.returncode:
                write_json(destination, dict(task=task, status="failed", exception="WorkerExit",
                                             message=f"exit {process.returncode}: {process.stdout[-3000:]} {process.stderr[-3000:]}"))
            manifest["completed"].append(key)
            write_json(path, manifest)
        return bench.read_json(destination)

    bath_tasks, later = bench.make_tasks(config, "full", True)
    baths = {}
    for key, task in bath_tasks.items():
        result = execute(dict(task, bath_key=key))
        if result["status"] == "ok":
            baths[key] = result["bath"]
    for task in later:
        if task["bath_key"] in baths:
            execute(dict(task, bath=baths[task["bath_key"]]))
    bench.export_results(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=bench.ROOT/"docs/benchmarks/baths/chain-input.json")
    parser.add_argument("--output", type=Path, default=bench.ROOT/"results/bath-chain-reference")
    parser.add_argument("--worker", nargs=2, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
    else:
        config = bench.read_json(args.input)
        if config.get("format_version") != 2:
            parser.error("chain reference input must have format_version=2")
        run(config, args.output)


if __name__ == "__main__":
    main()
