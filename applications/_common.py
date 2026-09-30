"""Small shared utilities for the literature reproductions."""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
from io import StringIO
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from qdjj_solver import DiscreteBath, Sector, chain_expansion, cosh_grid, fit_surrogate, save_result, solve
from qdjj_solver.common.io import numerical_environment, source_provenance, write_json


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


@lru_cache(maxsize=32)
def fitted_bath(levels, bandwidth, frequency_cutoff, delta=1.):
    return fit_surrogate(levels, delta=delta, bandwidth=bandwidth,
                         frequency_cutoff=frequency_cutoff, seed=1729, starts=4)


def bath_for(settings, bandwidth, delta=1.):
    """Replay an explicit bath or deterministically construct the requested one."""
    if "bath_record" in settings:
        bath = DiscreteBath.from_record(settings["bath_record"])
        if bath.bandwidth != bandwidth or bath.delta != delta:
            raise ValueError("stored bath does not match the physical energy scales")
        return bath
    if settings.get("bath_kind", "surrogate") == "cosh":
        return cosh_grid(settings["pairs"], delta=delta, bandwidth=bandwidth)
    if settings.get("bath_kind") == "chain-expansion":
        return chain_expansion(settings["levels"], delta=delta, bandwidth=bandwidth,
                               wide_band=settings.get("wide_band", False))
    return fitted_bath(settings["levels"], bandwidth, settings["frequency_cutoff"], delta)


def eigenstate(h, sector, settings, roots=1):
    backend = settings.get("backend", "qp")
    if backend == "qp":
        options = dict(tolerance=1e-12, residual_tolerance=2e-8, threads=1,
                       seed=1729, max_memory_gib=8.)
    else:
        options = dict(chi_max=128, require_convergence=True, seed_trials=2,
                       energy_tolerance=1e-11, residual_tolerance=1e-6,
                       threads=1, seed=1729)
    options.update(settings.get("solver", {}))
    options["eigenpairs"] = roots
    return solve(h, cutoff=settings.get("cutoff"), sector=sector,
                 backend=backend, options=options)


def parity_states(h, settings):
    return [eigenstate(h, Sector(p, p), settings) for p in (0, 1)]


def scalar(result, observable):
    return float(np.real(result.observables[observable][0]))


def state_summary(result):
    """Compact sweep data; complete standard-format records are optional."""
    record = dict(backend=result.backend, energies=result.energies,
                  residuals=result.residuals, observables=result.observables,
                  timings_seconds=result.timings,
                  hamiltonian_sha256=result.hamiltonian.fingerprint(),
                  sector=result.metadata["sector"],
                  options=result.metadata["options"])
    if result.backend == "qp":
        record.update(dimension=result.basis.dimension, qp_cutoff=result.basis.cutoff,
                      qp_weights=result.qp_weights)
    else:
        record["diagnostics"] = {k: v for k, v in result.metadata.items()
                                 if k not in result.hamiltonian.metadata}
    return record


class Run:
    """Archive inputs, baths, implementation identity and all scalar eigensolutions."""

    def __init__(self, case, config, profile, output, raw=False, *, half_filling=True):
        if config["model"]["delta"] != 1.:
            raise ValueError("application energies must be expressed in units of Delta=1")
        if half_filling and config["model"]["detuning"] != 0.:
            raise ValueError("this paper sweep requires particle-hole symmetry (detuning=0)")
        self.case = Path(case)
        self.config = config
        self.profile = profile
        self.output = Path(output)
        self.output.mkdir(parents=True, exist_ok=True)
        self.raw = raw
        self.states = []
        self.baths = {}
        self.started = perf_counter()
        self.script_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in sorted(self.case.glob("*.py"))}
        self.script_hashes["../_common.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def bath(self, settings, *, bandwidth=None, delta=1.):
        bandwidth = self.config["model"]["bandwidth"] if bandwidth is None else bandwidth
        bath = bath_for(settings, bandwidth, delta)
        record = bath.record()
        record["metadata"] = dict(record["metadata"])
        # Fitting time is not part of bath identity.
        record["metadata"].pop("fit_seconds", None)
        key = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
        self.baths[key] = record
        return bath

    def keep(self, label, result):
        self.states.append(dict(label=label, **state_summary(result)))
        if self.raw:
            save_result(result, self.output/"raw"/f"{len(self.states):05d}.json")

    def finish(self, tables, summary):
        import csv
        for name, rows in tables.items():
            if not rows:
                continue
            with (self.output/f"{name}.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        write_json(self.output/"eigenstates.json", self.states)
        write_json(self.output/"summary.json", summary)
        write_json(self.output/"manifest.json", dict(
            format="qdjj-application", format_version=1, case=self.case.name,
            profile=self.profile, input=self.config, baths=self.baths,
            environment=numerical_environment(), implementation=source_provenance(),
            scripts_sha256=self.script_hashes, elapsed_seconds=perf_counter()-self.started,
            state_records=len(self.states)))
        print(json.dumps(summary, indent=2, allow_nan=False))
        print(f"Wrote {self.output}", flush=True)


def arguments(case, description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", type=Path, default=Path(case)/"input"/"parameters.json")
    parser.add_argument("--profile", choices=("quick", "paper", "convergence"), default="quick")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--raw", action="store_true", help="also save standard qdjj-eigenstates JSON")
    args = parser.parse_args()
    args.output = args.output or Path(case)/"runs"/args.profile
    config = read_json(args.input)
    if config.get("format_version") != 1:
        parser.error("unsupported application input version")
    return args, config


def read_csv(path):
    import csv
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def plot_arguments(case):
    parser = argparse.ArgumentParser(description="Plot archived application output (no solver calls).")
    parser.add_argument("--output", type=Path, default=Path(case)/"output"/"paper")
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "qdjj-applications"})
    return args


def save_figure(fig, path):
    svg = StringIO()
    fig.savefig(svg, format="svg", bbox_inches="tight", metadata={"Date": None})
    Path(path).write_text("\n".join(line.rstrip() for line in svg.getvalue().splitlines())+"\n",
                          encoding="utf-8", newline="\n")
    fig.savefig(Path(path).with_suffix(".png"), bbox_inches="tight", dpi=150)
