"""Reproducible reservoir accuracy/cost study; see docs/bath_representations.md.

Each task runs serially in a fresh process. Timed repetitions follow one warm-up;
imports, reference integrals and output serialization are outside timed regions.
Exact baths, individual timings, failures and process high-water RSS are retained.
The quadratic continuum oracle is deliberately independent of the Fock/MPS solvers.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import logging
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre
from threadpoolctl import threadpool_info, threadpool_limits

from qdjj_solver import DiscreteBath, Sector, cosh_grid, fit_surrogate, reference_model, solve
from qdjj_solver.common.baths import discrete_g, hybridization_g
from qdjj_solver.common.io import fingerprint, numerical_environment, source_provenance, write_json

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "docs/benchmarks/baths/input.json"
THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def bath_spec(kind, levels, bandwidth, **fit_options):
    return dict(kind=kind, levels=levels, bandwidth=bandwidth, delta=1., **fit_options)


def construct_bath(spec):
    options = dict(spec)
    kind, levels = options.pop("kind"), options.pop("levels")
    if kind == "surrogate":
        return fit_surrogate(levels, **options)
    if kind == "cosh-grid":
        if levels % 2:
            raise ValueError("cosh-grid requires an even signed-level count")
        return cosh_grid(levels//2, **options)
    if kind == "linear-gl":
        # Independent control, not an additional core bath-kind string.
        x, w = roots_legendre(levels)
        return DiscreteBath(options["bandwidth"]*x, w/2, **options,
                            metadata=dict(kind="direct-xi-gauss-legendre"))
    raise ValueError(f"unknown benchmark bath family: {kind}")


def bath_record(bath):
    record = bath.record()
    record["metadata"] = dict(record["metadata"])
    record["metadata"].pop("fit_seconds", None)
    return record


def timing_summary(samples):
    values = np.asarray(samples, dtype=float)
    if values.ndim != 1 or not len(values) or not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("timing samples must be finite nonnegative values")
    return dict(samples=values.tolist(), median=float(np.median(values)),
                q25=float(np.quantile(values, .25)), q75=float(np.quantile(values, .75)),
                minimum=float(values.min()), maximum=float(values.max()))


def kernel_metrics(bath, minimum=1e-3, cutoff=100.):
    """Evaluate on a shifted/dense held-out grid, including zero and the tail."""
    inside = np.geomspace(minimum, cutoff, 2003)
    outside = np.r_[0., np.geomspace(minimum*1e-3, minimum*.99, 509),
                    np.geomspace(cutoff*1.01, max(10*bath.bandwidth, 100*cutoff), 997)]
    result = {}
    for label, omega in (("in", inside), ("out", outside)):
        exact = hybridization_g(omega, bath.delta, bath.bandwidth)
        error = discrete_g(bath, omega)-exact
        result.update({f"{label}_max_abs": float(np.max(abs(error))),
                       f"{label}_max_rel": float(np.max(abs(error/exact))),
                       f"{label}_rms_abs": float(np.sqrt(np.mean(error**2)))})
    for power in (0, 2, 4):
        # Dimensionless moments; continuum integral of (xi/D)^p is 1/(p+1).
        result[f"moment_{power}_error"] = float(
            np.dot(bath.weights, (bath.xi/bath.bandwidth)**power)-1/(power+1))
    result["minimum_qp_energy"] = float(bath.energies.min())
    weak_exact = 2/np.pi*quad(lambda x: 1/(5+np.hypot(x, bath.delta))**2,
                              0, bath.bandwidth, epsabs=1e-13, epsrel=1e-13)[0]
    weak_finite = float(np.sum(bath.weights/(np.pi*bath.rho)/(5+bath.energies)**2))
    result.update(weak_spin_slope=weak_finite, weak_spin_slope_continuum=weak_exact,
                  weak_spin_slope_error=abs(weak_finite-weak_exact))
    return result


def quadratic_observables(*, bandwidth, gamma, phi, detuning=0., delta=1., bath=None):
    """Finite-band U=0 subgap root, zero-T current and charge.

    Integrating out the flat BCS leads gives
      A(w) = w^2 (1+Gamma*g)^2 + eps^2 + (Gamma*Delta*cos(phi/2)*g)^2.
    dE/dphi = Gamma^2 Delta^2 sin(phi)/(2*pi) integral g^2/A dw;
    n = 1 - 2*eps/pi integral 1/A dw.
    This includes the continuum contribution to the current. A supplied symmetric
    bath replaces the continuum g by its finite positive-pole sum.
    """
    if gamma <= 0 or delta <= 0 or bandwidth <= 0:
        raise ValueError("positive Gamma, gap and half-bandwidth required")
    if bath is not None and (not bath.paired or bath.delta != delta or bath.bandwidth != bandwidth):
        raise ValueError("quadratic scalar reference requires a matching symmetric bath")
    cosine = np.cos(phi/2)

    def g(omega):
        return (float(hybridization_g(omega, delta, bandwidth)) if bath is None
                else float(discrete_g(bath, omega)))

    def denominator(w, value):
        return w*w*(1+gamma*value)**2 + detuning**2 + (gamma*delta*cosine*value)**2

    current_prefactor = gamma**2*delta**2*np.sin(phi)/(2*np.pi)
    # Avoid the nondifferentiable zero-energy crossing at eps=0, phi=pi.
    if abs(cosine) < 1e-13 and detuning == 0:
        raise ValueError("choose a point away from the zero-energy quadratic crossing")
    current_integral, current_error = quad(
        lambda w: g(w)**2/denominator(w, g(w)), 0, np.inf, epsabs=2e-11, epsrel=2e-11,
        limit=300)
    charge_integral, charge_error = quad(
        lambda w: 1/denominator(w, g(w)), 0, np.inf, epsabs=2e-11, epsrel=2e-11, limit=300)
    edge = delta if bath is None else float(bath.energies.min())

    def retarded_g(energy):
        if bath is None:
            a = np.sqrt((delta-energy)*(delta+energy))
            return 2/np.pi*np.arctan(bandwidth/a)/a
        return float(np.sum(bath.weights/(np.pi*bath.rho) /
                            ((bath.energies-energy)*(bath.energies+energy))))

    def pole_equation(energy):
        value = retarded_g(energy)
        return energy*(1+gamma*value)-np.hypot(detuning, gamma*delta*cosine*value)

    bound = brentq(pole_equation, 0., edge*(1-1e-12), xtol=1e-13)
    return dict(excitation=float(bound), current=float(current_prefactor*current_integral),
                charge=float(1-2*detuning/np.pi*charge_integral),
                current_quadrature_error=float(abs(current_prefactor)*current_error),
                charge_quadrature_error=float(abs(2*detuning/np.pi)*charge_error))


def quadratic_bdg(bath, gamma, phi, detuning=0.):
    """Independent physical-electron BdG diagonalization, for validating the oracle."""
    levels = bath.levels
    n = 1+2*levels
    normal = np.diag(np.r_[detuning, bath.xi, bath.xi]).astype(complex)
    derivative = np.zeros_like(normal)
    pair = -np.diag(np.r_[0., np.full(2*levels, bath.delta)])
    for lead, velocity in enumerate((-.5, .5)):
        sl = slice(1+lead*levels, 1+(lead+1)*levels)
        hopping = np.sqrt(gamma/(2*np.pi*bath.rho)*bath.weights)*np.exp(.5j*velocity*phi)
        normal[sl, 0], normal[0, sl] = hopping, hopping.conj()
        derivative[sl, 0] = .5j*velocity*hopping
        derivative[0, sl] = derivative[sl, 0].conj()
    matrix = np.block([[normal, pair], [pair, -normal.conj()]])
    response = np.block([[derivative, np.zeros_like(pair)],
                         [np.zeros_like(pair), -derivative.conj()]])
    energies, vectors = np.linalg.eigh(matrix)
    occupied = vectors[:, :n]
    return dict(excitation=float(energies[n]),
                energy=float(energies[:n].sum()+2*bath.energies.sum()),
                current=float(np.trace(occupied.conj().T @ response @ occupied).real),
                charge=float(1+np.sum(abs(occupied[0])**2)-np.sum(abs(occupied[n])**2)))


def benchmark_bath(task):
    construct_bath(task["spec"])  # unmeasured warm-up; no fitted-bath cache
    samples, baths = [], []
    for _ in range(task["repetitions"]):
        start = perf_counter()
        bath = construct_bath(task["spec"])
        samples.append(perf_counter()-start)
        baths.append(bath_record(bath))
    if any(record != baths[0] for record in baths[1:]):
        raise RuntimeError("fixed-seed bath coefficients changed between timing repetitions")
    bath = DiscreteBath.from_record(baths[0])
    minimum = task["spec"].get("frequency_min", 1e-3)
    cutoff = task["spec"].get("frequency_cutoff", 100.)
    quadratic = []
    for point in task["quadratic_points"]:
        parameters = dict(bandwidth=bath.bandwidth, delta=bath.delta, **point)
        exact = quadratic_observables(**parameters)
        finite = quadratic_observables(**parameters, bath=bath)
        quadratic.append(dict(parameters=parameters, continuum=exact, finite=finite,
                              errors={key: abs(finite[key]-exact[key])
                                      for key in ("excitation", "current", "charge")}))
    return dict(bath=baths[0], construction_seconds=timing_summary(samples),
                kernel=kernel_metrics(bath, minimum, cutoff), quadratic=quadratic)


def build_model(bath, parameters):
    parameters = dict(parameters)
    kind = parameters.pop("kind", "reference")
    if kind == "reference":
        return reference_model(bath, **parameters)
    if kind == "double-dot":
        # Reuse the reviewed physical model; do not duplicate its conventions.
        if str(ROOT) not in sys.path:
            sys.path.insert(0, str(ROOT))
        from applications.zonda_2023_double_dot.run import model
        return model(bath, **parameters)
    raise ValueError(f"unknown benchmark model: {kind}")


def compact_state(state):
    metadata = state.metadata
    record = dict(energy=float(state.energies[0]), residual=float(state.residuals[0]),
                  observables={k: float(np.real(v[0])) for k, v in state.observables.items()},
                  timings_seconds=state.timings, sector=metadata["sector"],
                  options=metadata["options"], bath_modes=metadata["bath_modes"],
                  hamiltonian_sha256=state.hamiltonian.fingerprint())
    if state.backend == "qp":
        record.update({k: metadata[k] for k in ("dimension", "method", "qp_cutoff", "nnz",
                                               "estimated_solver_bytes", "real_arithmetic")})
        record["qp_weights"] = state.qp_weights[0].tolist()
        record["finite_problem_converged"] = True  # only the projected finite problem
    else:
        record.update({k: metadata[k] for k in ("finite_problem_converged", "mpo_dimensions",
                                               "gram_error", "mode_order")})
        record["roots"] = [{k: root[k] for k in ("chi_actual", "sweeps", "sweep_energy_converged",
                                                 "seed_trials")} for root in metadata["roots"]]
    return record


def calculate_point(bath, task):
    start = perf_counter()
    h = build_model(bath, task["model"])
    built = perf_counter()
    settings = task["solver"]
    backend = settings.get("backend", "qp")
    options = dict(task["qp_options"] if backend == "qp" else task["dmrg_options"])
    options.update(settings.get("options", {}))
    states = [compact_state(solve(h, sector=Sector(p, p), backend=backend,
                                  cutoff=settings.get("cutoff"), options=options))
              for p in task["sectors"]]
    end = perf_counter()
    by_parity = dict(zip(task["sectors"], states, strict=True))
    values = {}
    for p, state in by_parity.items():
        values[f"energy_{p}"] = state["energy"]
        for name, value in state["observables"].items():
            values[f"{name}_{p}"] = value
    if 1 in by_parity:
        values["kappa"] = 1-2*by_parity[1]["observables"]["impurity_spin_z"]
    if 0 in by_parity and 1 in by_parity:
        values["signed_gap"] = by_parity[1]["energy"]-by_parity[0]["energy"]
        parity = int(values["signed_gap"] < 0)
        if "phase_derivative" in by_parity[parity]["observables"]:
            values["ground_current"] = by_parity[parity]["observables"]["phase_derivative"]
    times = dict(model=built-start, solve=end-built, reused_bath_total=end-start)
    for name in {key for state in states for key in state["timings_seconds"]}:
        times[f"solver_{name}"] = sum(state["timings_seconds"].get(name, 0.) for state in states)
    return dict(values=values, states=states, times=times,
                converged=all(s["finite_problem_converged"] for s in states),
                restriction=h.metadata["coordinates"]["restriction"])


def benchmark_solve(task):
    bath = DiscreteBath.from_record(task["bath"])
    calculate_point(bath, task)  # imports, allocation and backend initialization warm-up
    runs = []
    for _ in range(task["repetitions"]):
        gc.collect()
        runs.append(calculate_point(bath, task))
    return dict(values=runs[0]["values"], states=runs[0]["states"],
                restriction=runs[0]["restriction"],
                converged=all(run["converged"] for run in runs),
                repeat_value_spread={k: float(np.ptp([run["values"][k] for run in runs]))
                                     for k in runs[0]["values"]},
                timings_seconds={k: timing_summary([run["times"][k] for run in runs])
                                 for k in runs[0]["times"]})


def benchmark_legacy(task):
    from qdjj_solver.dmrg_solver.legacy import dmrg_chain_reference, dmrg_reference
    bath = DiscreteBath.from_record(task["bath"])
    function = dmrg_chain_reference if task["layout"] == "chain" else dmrg_reference
    parameters = dict(task["parameters"])
    function(bath, **parameters)
    samples, records = [], []
    for _ in range(task["repetitions"]):
        gc.collect()
        start = perf_counter()
        records.append(function(bath, **parameters))
        samples.append(perf_counter()-start)
    record = records[0]
    return dict(values={k: record[k] for k in ("energy", "impurity_spin_z", "phase_derivative")},
                diagnostics={k: record[k] for k in ("chi_actual", "sweeps", "variance", "norm_error",
                                                   "sweep_converged", "last_energy_change")},
                seconds=timing_summary(samples),
                repeat_energy_spread=float(np.ptp([r["energy"] for r in records])))


def peak_rss_bytes():
    try:
        import resource
    except ImportError:  # resource is not available on Windows
        return None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else 1024*peak)


def worker(request, destination):
    task = read_json(request)
    logging.getLogger("tenpy").setLevel(logging.ERROR)
    start = perf_counter()
    try:
        with threadpool_limits(limits=1):
            dispatch = {"bath": benchmark_bath, "solve": benchmark_solve, "legacy": benchmark_legacy}
            result = dict(status="ok", **dispatch[task["type"]](task))
    except (ValueError, RuntimeError, MemoryError, ImportError) as error:
        result = dict(status="failed", exception=type(error).__name__, message=str(error))
    result.update(task=task, peak_process_rss_bytes=peak_rss_bytes(),
                  worker_seconds=perf_counter()-start)
    write_json(destination, result)


def task_id(task):
    return fingerprint(task)[:20]


def make_tasks(config, profile, include_dmrg):
    """Expand the explicit input matrix; all fit settings enter bath/task identity."""
    smoke = profile == "smoke"
    bath_tasks, later = {}, []
    fitting = config["fitting"]

    def register(kind, levels, bandwidth, overrides=None):
        spec = bath_spec(kind, levels, bandwidth,
                         **(dict(fitting, **(overrides or {})) if kind == "surrogate" else {}))
        key = task_id(spec)
        bath_tasks[key] = dict(type="bath", spec=spec,
                               repetitions=1 if smoke else config["bath_repetitions"],
                               quadratic_points=config["quadratic_points"][:1] if smoke
                               else config["quadratic_points"])
        return key

    for bandwidth in ([10.] if smoke else config["bandwidths"]):
        for family, sizes in config["kernel_levels"].items():
            for levels in (sizes[:1] if smoke else sizes):
                register(family, levels, bandwidth)
    if not smoke:
        for tuning in config["fit_sensitivity"]:
            register("surrogate", tuning["levels"], tuning["bandwidth"], tuning["options"])
    for study in config["studies"]:
        if smoke and study["name"] != "knight":
            continue
        cases = study["cases"][:1] if smoke else study["cases"]
        for case in cases:
            for family in (study["families"][:2] if smoke else study["families"]):
                for levels in (study["levels"][:1] if smoke else study["levels"]):
                    bath_key = register(family, levels, study["bandwidth"], study.get("fit_options"))
                    for settings in study["solvers"]:
                        if settings.get("backend", "qp") == "dmrg" and not include_dmrg:
                            continue
                        later.append(dict(type="solve", study=study["name"], case=case["label"],
                                          model=case["model"], sectors=study["sectors"],
                                          bath_key=bath_key, solver=settings,
                                          qp_options=config["qp_options"], dmrg_options=config["dmrg_options"],
                                          repetitions=1 if smoke else config["solve_repetitions"]))
    if not smoke and include_dmrg:
        for levels in config["legacy_levels"]:
            for family in ("cosh-grid", "surrogate"):
                key = register(family, levels, 10.)
                for chi in (32, 64):
                    for layout in ("chain", "star"):
                        later.append(dict(type="legacy", bath_key=key, layout=layout,
                                          repetitions=config["solve_repetitions"],
                                          parameters=dict(u=10., gamma=1., phi=.8, geometry="single",
                                                          twice_sz=1, chi=chi, max_sweeps=60,
                                                          energy_tolerance=1e-11, threads=1,
                                                          calculate_variance=True)))
    return bath_tasks, later


def environment_record():
    # Import before hashing so an editable installation's native extension is included.
    from qdjj_solver.qp_solver import _core
    assert _core is not None

    def command(*args):
        try:
            return subprocess.check_output(args, cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    return dict(numerical=numerical_environment(), implementation=source_provenance(),
                processor=platform.processor(), cpu_count=os.cpu_count(),
                cpu_brand=command("sysctl", "-n", "machdep.cpu.brand_string") if sys.platform == "darwin" else None,
                physical_memory_bytes=command("sysctl", "-n", "hw.memsize") if sys.platform == "darwin" else None,
                threadpools=threadpool_info(), threads={key: os.environ.get(key) for key in THREAD_VARIABLES},
                git_revision=command("git", "rev-parse", "HEAD"), git_status=command("git", "status", "--short"))


def write_csv(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def export_results(directory, destination=None):
    """Archive only completed tasks listed in the manifest; no solver/fit calls."""
    directory = Path(directory)
    destination = directory if destination is None else Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    manifest = read_json(directory/"manifest.json")
    results = [read_json(directory/"jobs"/f"{key}.json") for key in manifest["completed"]]
    bath_results = {r["task"]["bath_key"]: r for r in results
                    if r["task"]["type"] == "bath" and r["status"] == "ok"}
    kernel_rows, quadratic_rows, solve_rows, legacy_rows, failures = [], [], [], [], []
    for result in results:
        task = result["task"]
        if result["status"] != "ok":
            failures.append(dict(task=task, exception=result["exception"], message=result["message"]))
            continue
        bath_result = bath_results[task["bath_key"]]
        spec = bath_result["task"]["spec"]
        identity = dict(bath_key=task["bath_key"], family=spec["kind"], levels=spec["levels"],
                        bandwidth=spec["bandwidth"])
        construction = bath_result["construction_seconds"]["median"]
        if task["type"] == "bath":
            kernel_rows.append(dict(**identity, **result["kernel"], construction_seconds=construction,
                                    construction_q25=result["construction_seconds"]["q25"],
                                    construction_q75=result["construction_seconds"]["q75"],
                                    **{k: spec.get(k, "") for k in ("frequency_cutoff", "frequency_min",
                                                                   "relative_weight", "starts", "seed",
                                                                   "frequency_points")}))
            for point in result["quadratic"]:
                quadratic_rows.append(dict(identity | point["parameters"],
                                           **{f"{k}_error": v for k, v in point["errors"].items()},
                                           **{f"{k}_continuum": v for k, v in point["continuum"].items()},
                                           **{f"{k}_finite": v for k, v in point["finite"].items()},
                                           construction_seconds=construction))
        elif task["type"] == "solve":
            times = result["timings_seconds"]
            row = dict(**identity, study=task["study"], case=task["case"],
                       backend=task["solver"].get("backend", "qp"),
                       cutoff=task["solver"].get("cutoff"),
                       chi=task["solver"].get("options", {}).get("chi_max"),
                       restriction=result["restriction"], converged=result["converged"],
                       **result["values"], max_residual=max(s["residual"] for s in result["states"]),
                       max_dimension=max((s.get("dimension", 0) for s in result["states"])),
                       construction_seconds=construction,
                       cold_total_seconds=construction+times["reused_bath_total"]["median"],
                       peak_process_rss_bytes=result["peak_process_rss_bytes"])
            row.update({f"{k}_seconds": v["median"] for k, v in times.items()})
            row.update(reused_q25=times["reused_bath_total"]["q25"],
                       reused_q75=times["reused_bath_total"]["q75"])
            solve_rows.append(row)
        else:
            legacy_rows.append(dict(**identity, layout=task["layout"], chi=task["parameters"]["chi"],
                                    **result["values"], **result["diagnostics"],
                                    seconds=result["seconds"]["median"],
                                    seconds_q25=result["seconds"]["q25"], seconds_q75=result["seconds"]["q75"],
                                    peak_process_rss_bytes=result["peak_process_rss_bytes"]))
    for name, rows in (("kernel", kernel_rows), ("quadratic", quadratic_rows),
                       ("many_body", solve_rows), ("layouts", legacy_rows)):
        write_csv(destination/f"{name}.csv", rows)
    # Exact input and per-task evidence, including timing samples and solver settings.
    write_json(destination/"input.json", manifest["input"])
    write_json(destination/"manifest.json", manifest)
    write_json(destination/"failures.json", failures)
    write_json(destination/"baths.json", {k: r["bath"] for k, r in bath_results.items()})
    compact = []
    for result in results:
        record = dict(result)
        record["task"] = {k: v for k, v in result["task"].items() if k != "bath"}
        record.pop("bath", None)  # already stored once in baths.json
        compact.append(record)
    write_json(destination/"measurements.json", compact)
    print(f"Exported {len(results)} tasks ({len(failures)} failures) to {destination}", flush=True)


def run(config, profile, output, include_dmrg=False, stage="all"):
    output = Path(output)
    (output/"jobs").mkdir(parents=True, exist_ok=True)
    bath_tasks, later = make_tasks(config, profile, include_dmrg)
    environment = environment_record()
    identity = dict(format="qdjj-bath-benchmark", format_version=1, input=config,
                    profile=profile, include_dmrg=include_dmrg,
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    implementation_sha256=environment["implementation"]["sha256"])
    manifest_path = output/"manifest.json"
    if manifest_path.exists():
        manifest = read_json(manifest_path)
        if any(manifest.get(key) != value for key, value in identity.items()):
            raise ValueError("resume requires identical input, runner, profile and solver implementation")
    else:
        manifest = dict(identity, environment=environment, completed=[], commands=[])
    manifest["commands"].append(sys.argv[1:])
    write_json(manifest_path, manifest)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", **{k: "1" for k in THREAD_VARIABLES})

    def execute(task):
        key = task_id(task)
        destination = output/"jobs"/f"{key}.json"
        if key not in manifest["completed"]:
            request = output/"request.json"
            write_json(request, task)
            description = task.get("case", task.get("layout", task.get("spec")))
            print(f"{task['type']} {description} {task.get('solver', '')}", flush=True)
            process = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()),
                                      "--worker", str(request), str(destination)], env=env,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if process.returncode:
                # Preserve externally killed or unexpected failures, not just convergence failures.
                write_json(destination, dict(task=task, status="failed", exception="WorkerExit",
                                             message=f"exit {process.returncode}: {process.stdout[-6000:]}"))
            manifest["completed"].append(key)
            write_json(manifest_path, manifest)
        result = read_json(destination)
        if result["status"] != "ok":
            print(f"  {result['exception']}: {result['message']}", flush=True)
        return result

    baths = {}
    for key, task in bath_tasks.items():
        result = execute(dict(task, bath_key=key))
        if result["status"] == "ok":
            baths[key] = result["bath"]
    if stage != "bath":
        for task in later:
            if stage != "all" and ((stage == "solve") != (task["type"] == "solve")):
                continue
            if task["bath_key"] in baths:
                execute(dict(task, bath=baths[task["bath_key"]]))
    export_results(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--output", type=Path, default=ROOT/"results/bath-benchmarks")
    parser.add_argument("--dmrg", action="store_true", help="include optional MPS and legacy layout jobs")
    parser.add_argument("--stage", choices=("all", "bath", "solve", "legacy"), default="all")
    parser.add_argument("--export", type=Path, help="export completed --output records here without running solvers")
    parser.add_argument("--worker", nargs=2, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(*args.worker)
    elif args.export:
        export_results(args.output, args.export)
    else:
        config = read_json(args.input)
        if config.get("format_version") != 1:
            parser.error("unsupported benchmark input version")
        if min(config["bath_repetitions"], config["solve_repetitions"]) < 1:
            parser.error("positive repetition counts required")
        run(config, args.profile, args.output, args.dmrg, args.stage)


if __name__ == "__main__":
    main()
