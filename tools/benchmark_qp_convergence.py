"""Reproduce the fixed-bath doublet cutoff study in docs/qp_convergence.md.

The check profile repeats Q=1..4; full repeats Q=1..6. --dmrg adds the
three bath/bond reference checks. Completed points are saved atomically and
reused only for identical inputs and code. --dry-run reports sector dimensions.
The optional one-time paper importer is not needed to use the shipped archive.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import gc
import hashlib
import json
import logging
from pathlib import Path
import subprocess
import sys

import numpy as np
from threadpoolctl import threadpool_limits

from qdjj_solver import DiscreteBath, Sector, cosh_grid, reference_model, solve
from qdjj_solver.common.io import fingerprint, numerical_environment, source_provenance, write_json
from qdjj_solver.qp_solver import dimension

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from examples.qp_convergence import doublet_values, values_from_state

DEFAULT_INPUT = ROOT/"docs/benchmarks/qp_convergence/input.json"
OBSERVABLES = ("impurity_charge", "double_occupancy_0", "impurity_spin_z")
VALUES = ("energy", "q_d", "P_0", "P_1", "P_2", "C_d_bath")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_config(config):
    # The inferred scalar spin correlation requires this physical state/model.
    if config.get("format_version") != 1:
        raise ValueError("unsupported study input version")
    model, sector = config["model"], Sector(**config["sector"])
    if (model["field"] != 0 or model["rho_ws"] != 0 or model["rho_wn"] != 0
            or model["detuning"] != 0 or model["bath_reference"] != "isolated"
            or not model["compress"] or not model["symmetry"] or sector != Sector(1, 1, 1)):
        raise ValueError("this study follows the symmetric, zero-field eta=+1 doublet")
    if (any(n < 2 or n % 2 for n in config["levels"])
            or any(q < 1 for q in config["cutoffs"])):
        raise ValueError("positive cutoffs and even signed bath-level counts are required")


def task_list(config, profile, include_dmrg):
    validate_config(config)
    tasks = [dict(id=f"qp_N{n}_Q{q}", levels=n, backend="qp", cutoff=q, chi_max=None,
                  options=config["qp_options"])
             for n in config["levels"] for q in config["cutoffs"]
             if profile == "full" or q <= 4]
    if include_dmrg:
        tasks += [dict(id=f"dmrg_N{r['levels']}_chi{r['chi_max']}", backend="dmrg", cutoff=None,
                       **r, options=dict(config["dmrg_options"], chi_max=r["chi_max"]))
                  for r in config["references"]]
    return tasks


def validate_values(values, weights=None):
    if set(values) != set(VALUES) or not all(np.isfinite(v) for v in values.values()):
        raise ValueError("missing or non-finite doublet observable")
    probabilities = [values[f"P_{i}"] for i in range(3)]
    if (min(probabilities) < -2e-9 or max(probabilities) > 1+2e-9
            or abs(sum(probabilities)-1) > 2e-9 or abs(values["P_0"]-values["P_2"]) > 2e-9):
        raise ValueError("invalid particle-hole-symmetric charge probabilities")
    # For a spin-half doublet, -P1/3 <= q_d <= P1 (bath spin zero or one).
    if not -values["P_1"]/3-2e-9 <= values["q_d"] <= values["P_1"]+2e-9:
        raise ValueError("spin values are incompatible with the assumed doublet")
    if weights is not None and (min(weights) < -2e-9 or abs(sum(weights)-1) > 2e-9):
        raise ValueError("QP weights are not normalized probabilities")


def save_archive(output, manifest, baths, records):
    output = Path(output)
    write_json(output/"baths.json", baths)
    write_json(output/"measurements.json", records)
    manifest = dict(manifest, completed=[r["id"] for r in records],
                    measurements_sha256=fingerprint(records), baths_sha256=fingerprint(baths))
    write_json(output/"manifest.json", manifest)


def load_archive(directory):
    directory = Path(directory)
    manifest = read_json(directory/"manifest.json")
    baths, records = read_json(directory/"baths.json"), read_json(directory/"measurements.json")
    if (manifest["completed"] != [r["id"] for r in records]
            or len({r["id"] for r in records}) != len(records)
            or manifest["measurements_sha256"] != fingerprint(records)
            or manifest["baths_sha256"] != fingerprint(baths)):
        raise ValueError("archive manifest and evidence disagree")
    for key, bath in baths.items():
        if fingerprint(bath) != key:
            raise ValueError("bath checksum mismatch")
    for row in records:
        bath = DiscreteBath.from_record(baths[row["bath_key"]])
        if bath.levels != row["levels"]:
            raise ValueError("record refers to the wrong bath")
        derived = doublet_values(row["values"]["energy"], *(row["observables"][k] for k in OBSERVABLES))
        if derived != row["values"]:
            raise ValueError("derived observables disagree with measured values")
        validate_values(derived, row["qp_weights"])
    return manifest, baths, records


def export_archive(source, destination):
    manifest, baths, records = load_archive(source)
    write_json(Path(destination)/"input.json", manifest["input"])
    save_archive(destination, manifest, baths, records)
    print(f"Exported {len(records)} checked records to {destination}")


def run(config, profile, output, include_dmrg=False, dry_run=False):
    tasks = task_list(config, profile, include_dmrg)
    bath_objects = {n: cosh_grid(n//2, **config["bath"]) for n in {t["levels"] for t in tasks}}
    models = {n: reference_model(b, **config["model"]) for n, b in bath_objects.items()}
    sector = Sector(**config["sector"])
    if dry_run:
        for task in tasks:
            if task["backend"] == "qp":
                h = models[task["levels"]]
                d = dimension(h.nimp, h.spins, task["cutoff"], sector, h.eta_labels)
                print(f"{task['id']}: {d:,} states; one real vector = {8*d/1024**2:.2f} MiB")
            else:
                print(f"{task['id']}: unrestricted MPS, chi_max={task['chi_max']}")
        print("Vector sizes exclude basis, sparse matrix and Krylov work arrays; solver limits remain active.")
        return
    output = Path(output)
    identity = dict(format="qdjj-qp-convergence", format_version=1, origin="current",
                    input=config, profile=profile, include_dmrg=include_dmrg,
                    implementation=source_provenance(),
                    scripts_sha256={"runner": file_hash(__file__),
                                    "example": file_hash(ROOT/"examples/qp_convergence.py")})
    if (output/"manifest.json").exists():
        manifest, baths, records = load_archive(output)
        if any(manifest.get(k) != v for k, v in identity.items()):
            raise ValueError("resume requires identical input, profile, runner and solver implementation")
    else:
        manifest = dict(identity, environment=numerical_environment())
        baths, records = {}, []
    write_json(output/"input.json", config)
    for task in tasks:
        if task["id"] in {r["id"] for r in records}:
            continue
        bath, h = bath_objects[task["levels"]], models[task["levels"]]
        bath_record = bath.record()
        key = fingerprint(bath_record)
        baths[key] = bath_record
        print(f"Solving {task['id']}", flush=True)
        # States are released after each point; the solver preflights storage.
        with threadpool_limits(limits=1):
            state = solve(h, sector=sector, cutoff=task["cutoff"], backend=task["backend"],
                          options=task["options"])
        values = values_from_state(state)
        weights = state.qp_weights[0].tolist() if task["backend"] == "qp" else None
        validate_values(values, weights)
        diagnostics = {k: state.metadata[k] for k in (
            "finite_problem_converged", "sweep_converged", "variance", "basis_bytes",
            "estimated_solver_bytes", "real_arithmetic", "nnz", "roots", "gram_error",
            "roots_energy_ordered", "residual_method") if k in state.metadata}
        records.append(dict(task, origin="current", bath_key=key, model=config["model"],
                            sector=asdict(sector), options=state.metadata["options"],
                            values=values, qp_weights=weights,
                            observables={k: float(np.real(state.observables[k][0])) for k in OBSERVABLES},
                            residual=float(state.residuals[0]), dimension=state.metadata.get("dimension"),
                            residual_kind="projected" if task["backend"] == "qp" else "full finite Hamiltonian",
                            timings_seconds=state.timings, diagnostics=diagnostics))
        save_archive(output, manifest, baths, records)
        print(f"  E={values['energy']:.12f}, q_d={values['q_d']:.10f}", flush=True)
        del state
        gc.collect()


def import_paper(source, output, config):
    """Extract compact measured evidence; retain checksums of all original files.

    This is a one-time migration of reviewed records, not a runtime dependency.
    Original bath coefficients are retained separately even when they differ
    only by quadrature-library roundoff. No DMRG residual is invented from a
    variance obtained by subtracting nearly equal extensive expectations.
    """
    validate_config(config)
    source, output = Path(source), Path(output)
    if (output/"manifest.json").exists():
        raise ValueError("import into a new archive directory")
    baths, records = {}, []

    def source_info(path):
        return dict(path=f"qp_solver/benchmarks/asq_reference/data/{path.name}", sha256=file_hash(path))

    def register(bath):
        key = fingerprint(bath)
        baths[key] = bath
        return key

    for n in config["levels"]:
        for q in config["cutoffs"]:
            path = source/"data"/f"qp-convergence_ws0_N{n}_Q{q}.json"
            if not path.exists():
                path = source/"data"/f"qp_ws0_N{n}_Q{q}.json"
            original = read_json(path)
            states = [r["states"][0] for r in original["repetitions"]]
            values = []
            for state in states:
                meta = state["calculation"]
                for name in ("u", "gamma", "phi", "rho_ws", "rho_wn", "field", "detuning"):
                    if meta[name] != config["model"][name]:
                        raise ValueError(f"wrong paper parameter {name} in {path.name}")
                if (Sector(**meta["sector"]) != Sector(**config["sector"])
                        or meta["qp_cutoff"] != q or meta["bath_reference"] != "isolated"
                        or not meta["paired_mode_compression"] or max(state["residuals"]) > 2e-9):
                    raise ValueError(f"wrong or unconverged paper sector in {path.name}")
                obs = {k: float(state["observables"][k][0]) for k in OBSERVABLES}
                values.append(doublet_values(state["energies"][0], *(obs[k] for k in OBSERVABLES)))
                validate_values(values[-1], state["qp_weights"][0])
            state, meta = states[0], states[0]["calculation"]
            bath = meta["bath"][0]
            if len(bath["xi"]) != n or any(bath[k] != v for k, v in config["bath"].items()):
                raise ValueError("wrong paper bath")
            spread = {k: max(abs(v[k]-values[0][k]) for v in values) for k in VALUES}
            if max(spread.values()) > 2e-9:
                raise ValueError("repeated paper calculations disagree")
            records.append(dict(id=f"qp_N{n}_Q{q}", origin="paper", backend="qp", levels=n,
                                cutoff=q, chi_max=None, bath_key=register(bath), model=config["model"],
                                sector=asdict(Sector(**meta["sector"])), values=values[0],
                                observables={k: float(state["observables"][k][0]) for k in OBSERVABLES},
                                qp_weights=state["qp_weights"][0], options=meta["options"],
                                dimension=meta["dimension"], residual=max(s["residuals"][0] for s in states),
                                residual_kind="projected", repetition_spread=spread,
                                diagnostics={k: meta[k] for k in ("basis_bytes", "estimated_solver_bytes",
                                                                 "real_arithmetic", "nnz")},
                                timings_seconds=state["timings_seconds"],
                                total_seconds=original["repetitions"][0]["seconds"],
                                peak_resident_bytes=original["peak_resident_bytes"],
                                environment=original["environment"], source=source_info(path)))
    for reference in config["references"]:
        n, chi = reference["levels"], reference["chi_max"]
        path = source/"data"/f"reference_ws0_N{n}_chi{chi}.json"
        original = read_json(path)
        state = original["results"][0]
        if (not state["sweep_converged"] or state["rho_ws"] != 0 or state["twice_sz"] != 1
                or not state["paired_mode_compression"] or state["chi_limit"] != chi
                or any(state[k] != config["model"][k] for k in ("u", "gamma", "phi"))):
            raise ValueError("wrong or unconverged paper reference")
        bath = state["bath"]
        if len(bath["xi"]) != n or any(bath[k] != v for k, v in config["bath"].items()):
            raise ValueError("wrong reference bath")
        obs = {k: state[k] for k in OBSERVABLES}
        values = doublet_values(state["energy"], *(obs[k] for k in OBSERVABLES))
        validate_values(values)
        records.append(dict(id=f"dmrg_N{n}_chi{chi}", origin="paper", backend="dmrg", levels=n,
                            cutoff=None, chi_max=chi, bath_key=register(bath), model=config["model"],
                            sector=asdict(Sector(1, 1)), values=values, observables=obs,
                            qp_weights=None, dimension=None, residual=None,
                            residual_kind="not recorded; legacy variance diagnostic only",
                            options={k: state[k] for k in ("chi_limit", "max_sweeps", "energy_tolerance", "threads")},
                            diagnostics={k: state[k] for k in ("chi_actual", "variance", "last_discarded_weight",
                                                             "last_energy_change", "norm_error", "sweep_converged",
                                                             "sweeps", "sweep_energies")},
                            total_seconds=state["seconds"], environment=original["environment"],
                            source=source_info(path)))
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    manifest = dict(format="qdjj-qp-convergence", format_version=1, origin="paper",
                    input=config, source_repository="ASQ_variational", source_revision_context=revision,
                    provenance_note="File checksums identify imported working-tree contents; revision is context.",
                    importer_sha256=file_hash(__file__), example_sha256=file_hash(ROOT/"examples/qp_convergence.py"))
    write_json(output/"input.json", config)
    save_archive(output, manifest, baths, records)
    load_archive(output)
    print(f"Imported {len(records)} compact records and {len(baths)} exact bath records into {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=ROOT/"results/qp-convergence")
    parser.add_argument("--profile", choices=("check", "full"), default="check")
    parser.add_argument("--dmrg", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--import-paper", type=Path, help="one-time import of the asq_reference directory")
    parser.add_argument("--export", type=Path, help="export the checked --output archive without solving")
    args = parser.parse_args()
    logging.basicConfig(level=logging.ERROR)
    config = read_json(args.input)
    if args.export:
        export_archive(args.output, args.export)
    elif args.import_paper:
        import_paper(args.import_paper, args.output, config)
    else:
        run(config, args.profile, args.output, args.dmrg, args.dry_run)


if __name__ == "__main__":
    main()
