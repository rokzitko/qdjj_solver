"""Regenerate bath-comparison figures and tables from archived measurements only."""

from __future__ import annotations

import argparse
from io import StringIO
from pathlib import Path

import numpy as np

from qdjj_solver import DiscreteBath
from qdjj_solver.common.baths import discrete_g, hybridization_g
from qdjj_solver.common.io import fingerprint, write_json

try:  # direct script and namespace-package imports used by tests
    from tools.benchmark_baths import make_tasks, read_json, write_csv
except ModuleNotFoundError:
    from benchmark_baths import make_tasks, read_json, write_csv


FAMILIES = ("cosh-grid", "surrogate", "linear-gl")
LABELS = {"cosh-grid": "Cosh Gauss–Legendre", "surrogate": "Surrogate", "linear-gl": "Direct-ξ Gauss–Legendre"}
COLORS = {"cosh-grid": "#2374ab", "surrogate": "#d1495b", "linear-gl": "#37875b"}


def is_default(spec, config):
    return spec["kind"] != "surrogate" or all(spec.get(k) == v for k, v in config["fitting"].items())


def canonical_study(name):
    return name.removesuffix("-refinement").removesuffix("-odd")


def primary_observable(study):
    return "kappa" if study == "knight" else "signed_gap"


def solver_policy(settings, ignore_chi=False):
    """Concrete permutations grow with bath size; compare the declared layout policy."""
    normalized = dict(settings)
    normalized["options"] = dict(settings.get("options", {}))
    normalized["options"].pop("mode_order", None)
    if ignore_chi:
        normalized["options"].pop("chi_max", None)
    return normalized


def compare_interacting(records, specs, config):
    """Use a labelled empirical reference, retaining separate stability indicators.

    Prefer the largest converged cosh-grid DMRG refinement, otherwise the largest
    unrestricted cosh-grid finite problem. Never turn a failed MPS run into a
    reference, or silently call the largest available bath an exact continuum.
    """
    groups = {}
    for record in records:
        task = record["task"]
        if (task["type"] != "solve" or record["status"] != "ok"
                or task["study"] in ("coordinates", "layout-control")):
            continue
        key = canonical_study(task["study"]), task["case"]
        groups.setdefault(key, []).append(record)
    rows, references = [], []
    for (study, case), candidates in sorted(groups.items()):
        observable = primary_observable(study)

        def spec(record):
            return specs[record["task"]["bath_key"]]

        def full(record):
            settings = record["task"]["solver"]
            return settings.get("backend", "qp") == "qp" and settings.get("cutoff") is None

        physics = {fingerprint(dict(model=r["task"].get("model"), sectors=r["task"].get("sectors"),
                                    bandwidth=spec(r)["bandwidth"], delta=spec(r)["delta"])) for r in candidates}
        if len(physics) != 1:
            raise ValueError(f"different physical problems share benchmark case {study}/{case}")
        cosh = [r for r in candidates if spec(r)["kind"] == "cosh-grid"]
        mps = [r for r in cosh if r["task"]["solver"].get("backend") == "dmrg"]
        converged_mps = [r for r in mps if r["converged"]]
        pool = converged_mps or [r for r in cosh if full(r) and r["converged"]]
        if not pool:
            continue
        reference = max(pool, key=lambda r: (spec(r)["levels"],
                                            r["task"]["solver"].get("options", {}).get("chi_max", 0),
                                            -max(s["residual"] for s in r["states"])))
        value = reference["values"][observable]
        size = spec(reference)["levels"]
        solver = reference["task"]["solver"]
        smaller = [r for r in pool if spec(r)["levels"] < size and
                   solver_policy(r["task"]["solver"]) == solver_policy(solver)]
        smaller = max(smaller, key=lambda r: spec(r)["levels"]) if smaller else None
        lower_chi = [r for r in mps if spec(r)["levels"] == size and
                     solver_policy(r["task"]["solver"], True) == solver_policy(solver, True) and
                     r["task"]["solver"].get("options", {}).get("chi_max", 0) <
                     solver.get("options", {}).get("chi_max", 0)]
        lower_chi = max(lower_chi, key=lambda r: r["task"]["solver"]["options"]["chi_max"]) if lower_chi else None
        other = [r for r in candidates if spec(r)["kind"] == "surrogate" and
                 r["converged"] and spec(r)["levels"] == size and
                 solver_policy(r["task"]["solver"]) == solver_policy(solver) and is_default(spec(r), config)]
        differences = dict(
            bath_step=None if smaller is None else abs(value-smaller["values"][observable]),
            bond_step=None if lower_chi is None else abs(value-lower_chi["values"][observable]),
            other_family=None if not other else abs(value-other[0]["values"][observable]))
        label = ("empirical MPS reference" if smaller is not None
                 else "finite MPS reference (bath unresolved)") if converged_mps else "largest finite ED bath"
        references.append(dict(study=study, case=case, observable=observable, value=value,
                               bath_key=reference["task"]["bath_key"], levels=size, solver=solver,
                               label=label,
                               stability_indicators=differences,
                               lower_bond_converged=None if lower_chi is None else lower_chi["converged"],
                               max_residual=max(s["residual"] for s in reference["states"])))
        for record in candidates:
            if not full(record) or not record["converged"] or not is_default(spec(record), config):
                continue
            task = record["task"]
            times = record["timings_seconds"]["reused_bath_total"]
            rows.append(dict(study=study, case=case, family=spec(record)["kind"],
                             levels=spec(record)["levels"], bath_key=task["bath_key"],
                             observable=observable, value=record["values"][observable],
                             reference=value, difference=abs(record["values"][observable]-value),
                             reused_seconds=times["median"], reused_q25=times["q25"], reused_q75=times["q75"],
                             peak_process_rss_bytes=record["peak_process_rss_bytes"],
                             reference_label=references[-1]["label"]))
    return rows, references


def coordinate_comparisons(records, specs):
    selected = [r for r in records if r["status"] == "ok" and r["task"]["type"] == "solve"
                and r["task"]["study"] == "coordinates"]
    exact = {(r["task"]["bath_key"], r["task"]["case"]): r for r in selected
             if r["task"]["solver"].get("backend", "qp") == "qp"
             and r["task"]["solver"].get("cutoff") is None}
    rows = []
    for record in selected:
        task = record["task"]
        reference = exact.get((task["bath_key"], "isolated"))
        if reference is None:
            continue
        rows.append(dict(bath_key=task["bath_key"], family=specs[task["bath_key"]]["kind"],
                         levels=specs[task["bath_key"]]["levels"], coordinates=task["case"],
                         backend=task["solver"].get("backend", "qp"), cutoff=task["solver"].get("cutoff"),
                         chi=task["solver"].get("options", {}).get("chi_max"), converged=record["converged"],
                         max_energy_error=max(abs(record["values"][f"energy_{p}"]-
                                                  reference["values"][f"energy_{p}"]) for p in (0, 1)),
                         gap_error=abs(record["values"]["signed_gap"]-reference["values"]["signed_gap"]),
                         seconds=record["timings_seconds"]["reused_bath_total"]["median"]))
    return rows


def legacy_comparisons(records, specs):
    controls = {}
    for record in records:
        task = record["task"]
        if (record["status"] == "ok" and task["type"] == "solve" and
                task["study"] == "layout-control" and task["solver"].get("backend", "qp") == "qp"
                and task["solver"].get("cutoff") is None and task["sectors"] == [1]):
            model = task["model"]
            controls[task["bath_key"], model["u"], model["gamma"], model["geometry"]] = record
    rows = []
    for record in records:
        task = record["task"]
        if record["status"] != "ok" or task["type"] != "legacy":
            continue
        parameters = task["parameters"]
        if parameters["geometry"] != "single" or parameters["twice_sz"] != 1:
            continue
        reference = controls.get((task["bath_key"], parameters["u"], parameters["gamma"], parameters["geometry"]))
        if reference is None:
            continue
        rows.append(dict(bath_key=task["bath_key"], family=specs[task["bath_key"]]["kind"],
                         levels=specs[task["bath_key"]]["levels"], layout=task["layout"],
                         chi=parameters["chi"], energy=record["values"]["energy"],
                         exact_energy=reference["values"]["energy_1"],
                         energy_error=abs(record["values"]["energy"]-reference["values"]["energy_1"]),
                         spin_error=abs(record["values"]["impurity_spin_z"]-
                                        reference["values"]["impurity_spin_z_1"]),
                         seconds=record["seconds"]["median"], **record["diagnostics"]))
    return rows


def render_table(headers, rows):
    def format_value(value):
        if value is None:
            return "—"
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)
    return ("| " + " | ".join(headers) + " |\n|" + "|".join(["---"]*len(headers)) + "|\n" +
            "\n".join("| " + " | ".join(format_value(v) for v in row) + " |" for row in rows) + "\n")


def save_figure(fig, path):
    stream = StringIO()
    fig.savefig(stream, format="svg", bbox_inches="tight", metadata={"Date": None})
    Path(path).write_text("\n".join(line.rstrip() for line in stream.getvalue().splitlines())+"\n",
                          encoding="utf-8", newline="\n")


def generate(directory, supplements=()):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "qdjj-bath-benchmark",
                                "font.size": 10})
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter

    directory = Path(directory)
    records = read_json(directory/"measurements.json")
    config = read_json(directory/"input.json")
    baths = read_json(directory/"baths.json")
    manifests = [read_json(directory/"manifest.json")]
    for supplement in supplements:
        supplement = Path(supplement)
        manifests.append(read_json(supplement/"manifest.json"))
        if manifests[-1]["implementation_sha256"] != manifests[0]["implementation_sha256"]:
            raise ValueError("supplement uses a different solver implementation")
        for key, bath in read_json(supplement/"baths.json").items():
            if key in baths and baths[key] != bath:
                raise ValueError("supplement changed a bath's coefficients or provenance")
            baths[key] = bath
        existing = {(r["task"]["type"], r["task"]["bath_key"]) for r in records}
        # Retain primary-run fitting timings for shared baths; later fits are
        # independently retained in the supplemental archive, not averaged in.
        records.extend(r for r in read_json(supplement/"measurements.json")
                       if r["task"]["type"] != "bath" or ("bath", r["task"]["bath_key"]) not in existing)
    bath_results = {r["task"]["bath_key"]: r for r in records
                    if r["task"]["type"] == "bath" and r["status"] == "ok"}
    specs = {k: r["task"]["spec"] for k, r in bath_results.items()}
    defaults = [r for r in bath_results.values() if is_default(r["task"]["spec"], config)]
    bands = config["bandwidths"]

    fig, axes = plt.subplots(3, len(bands), figsize=(12, 10), squeeze=False)
    for column, bandwidth in enumerate(bands):
        for family in FAMILIES:
            points = sorted([r for r in defaults if r["task"]["spec"]["bandwidth"] == bandwidth and
                             r["task"]["spec"]["kind"] == family and r["task"]["spec"]["levels"] % 2 == 0],
                            key=lambda r: r["task"]["spec"]["levels"])
            if not points:
                continue
            levels = [r["task"]["spec"]["levels"] for r in points]
            values = [[r["kernel"]["in_max_abs"] for r in points],
                      [max(p["errors"]["excitation"] for p in r["quadratic"]) for r in points],
                      [max(p["errors"]["current"] for p in r["quadratic"]) for r in points]]
            for row, errors in enumerate(values):
                axes[row, column].loglog(levels, np.maximum(errors, 1e-16), "o-", ms=3,
                                         color=COLORS[family], label=LABELS[family])
        axes[0, column].set_title(f"D/Δ = {bandwidth:g}")
        for row in range(3):
            axes[row, column].grid(alpha=.2, which="both")
            axes[row, column].set_xlabel("Signed normal levels per reservoir")
    for row, label in enumerate(("max |g error| in [0.001, 100]", "max U=0 excitation error / Δ",
                                  "max U=0 current error / (2eΔ/ℏ)")):
        axes[row, 0].set_ylabel(label)
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, directory/"accuracy.svg")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    omega = np.geomspace(1e-5, 1e5, 1001)
    for family in FAMILIES:
        selected = [r for r in defaults if r["task"]["spec"]["bandwidth"] == 100 and
                    r["task"]["spec"]["kind"] == family]
        even = sorted([r for r in selected if r["task"]["spec"]["levels"] % 2 == 0],
                      key=lambda r: r["task"]["spec"]["levels"])
        axes[0].loglog([r["task"]["spec"]["levels"] for r in even],
                       [r["construction_seconds"]["median"] for r in even], "o-", ms=3,
                       color=COLORS[family], label=LABELS[family])
        for result in selected:
            if result["task"]["spec"]["levels"] != 8:
                continue
            bath = DiscreteBath.from_record(baths[result["task"]["bath_key"]])
            error = discrete_g(bath, omega)-hybridization_g(omega, bandwidth=100.)
            axes[1].loglog(omega, np.maximum(abs(error), 1e-16), color=COLORS[family])
            axes[2].loglog(omega, np.maximum(abs(error/hybridization_g(omega, bandwidth=100.)), 1e-16),
                           color=COLORS[family])
    axes[0].set(xlabel="Signed levels", ylabel="Construction / fit time (s)", title="D/Δ=100; median of five")
    axes[0].legend(fontsize=8)
    for ax, label in zip(axes[1:], ("Absolute g error", "Relative g error"), strict=True):
        ax.axvspan(.001, 100, alpha=.08, color="black")
        ax.set(xlabel="Imaginary-axis frequency / Δ", ylabel=label, title="D/Δ=100; L=8")
    for ax in axes:
        ax.grid(alpha=.2, which="both")
    fig.tight_layout()
    save_figure(fig, directory/"construction_and_window.svg")
    plt.close(fig)

    rows, references = compare_interacting(records, specs, config)
    for row in rows:
        row["construction_seconds"] = bath_results[row["bath_key"]]["construction_seconds"]["median"]
        row["cold_seconds"] = row["construction_seconds"]+row["reused_seconds"]
    write_csv(directory/"interacting_comparison.csv", rows)
    coordinates = coordinate_comparisons(records, specs)
    write_csv(directory/"coordinate_comparison.csv", coordinates)
    write_csv(directory/"layout_comparison.csv", legacy_comparisons(records, specs))
    refinement_rows = []
    for record in records:
        task = record["task"]
        if (record["status"] != "ok" or task["type"] != "solve" or
                not task["study"].endswith("-refinement")):
            continue
        settings, spec = task["solver"], specs[task["bath_key"]]
        observable = primary_observable(canonical_study(task["study"]))
        refinement_rows.append(dict(
            study=canonical_study(task["study"]), case=task["case"], family=spec["kind"], levels=spec["levels"],
            backend=settings.get("backend", "qp"), layout=settings.get("layout", "declared-star"),
            cutoff=settings.get("cutoff"), chi=settings.get("options", {}).get("chi_max"),
            observable=observable, value=record["values"][observable], converged=record["converged"],
            max_residual=max(s["residual"] for s in record["states"]),
            residual_scope="projected" if settings.get("cutoff") is not None else "declared finite Hamiltonian",
            seconds=record["timings_seconds"]["reused_bath_total"]["median"]))
    write_csv(directory/"refinement.csv", refinement_rows)

    panels = [("knight", "g1-p0.5"), ("spectrum", "g0.5"), ("spectrum", "g20")]
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))
    for column, (study, case) in enumerate(panels):
        for family in FAMILIES:
            selected = sorted([r for r in rows if r["study"] == study and r["case"] == case and
                               r["family"] == family and r["levels"] % 2 == 0], key=lambda r: r["levels"])
            if not selected:
                continue
            errors = np.maximum([r["difference"] for r in selected], 1e-16)
            for index, x in enumerate(([r["levels"] for r in selected],
                                        [r["reused_seconds"] for r in selected],
                                        [r["cold_seconds"] for r in selected])):
                axes[index, column].loglog(x, errors, "o-", ms=4, color=COLORS[family], label=LABELS[family])
                for a, b, row in zip(x, errors, selected, strict=True):
                    axes[index, column].annotate(str(row["levels"]), (a, b), xytext=(3, 4),
                                                 textcoords="offset points", fontsize=7)
        axes[0, column].set_title("Knight shift: Γ/Δ=1, φ=π/2" if study == "knight"
                                  else f"Spectrum: U/Δ=15, Γ/Δ={case[1:]}")
        shown_levels = sorted({r["levels"] for r in rows if r["study"] == study and r["case"] == case
                               and r["levels"] % 2 == 0})
        axes[0, column].set_xticks(shown_levels, [str(level) for level in shown_levels])
        axes[0, column].xaxis.set_minor_formatter(NullFormatter())
        reference = next((r for r in references if r["study"] == study and r["case"] == case), None)
        if reference and reference["stability_indicators"]["bath_step"] is None:
            for ax in axes[:, column]:
                ax.text(.04, .04, "Bath reference unresolved", transform=ax.transAxes,
                        fontsize=8, bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))
        for row, xlabel in enumerate(("Signed levels", "Reused-bath wall time (s)", "Fresh-fit total time (s)")):
            axes[row, column].set_xlabel(xlabel)
            axes[row, column].grid(alpha=.2, which="both")
    for row in range(3):
        axes[row, 0].set_ylabel("Absolute difference from reference")
    axes[0, 0].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, directory/"interacting_cost.svg")
    plt.close(fig)

    points = [r for r in records if r["status"] == "ok" and r["task"]["type"] == "solve"
              and r["task"]["study"] == "knight" and r["task"]["case"] == "g1-p0.5"
              and r["task"]["solver"].get("backend", "qp") == "qp"
              and r["task"]["solver"].get("cutoff") is None
              and is_default(specs[r["task"]["bath_key"]], config)]
    if points:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        bars = [r for r in points if specs[r["task"]["bath_key"]]["kind"] != "linear-gl"
                and specs[r["task"]["bath_key"]]["levels"] >= 4]
        heights = []
        for record in bars:
            timing = record["timings_seconds"]
            parts = ("model", "solver_basis", "solver_hamiltonian", "solver_diagonalization",
                     "solver_observables_and_residuals")
            overhead = np.median([max(0., total-sum(timing[k]["samples"][i] for k in parts))
                                  for i, total in enumerate(timing["reused_bath_total"]["samples"])])
            heights.append([bath_results[record["task"]["bath_key"]]["construction_seconds"]["median"],
                            timing["model"]["median"],
                            timing["solver_basis"]["median"]+timing["solver_hamiltonian"]["median"],
                            timing["solver_diagonalization"]["median"],
                            timing["solver_observables_and_residuals"]["median"], overhead])
        bottom = np.zeros(len(bars))
        for index, label in enumerate(("Bath / fit", "Model", "Basis + assembly", "Diagonalization",
                                        "Observables + residual", "Other / scalar records")):
            values = np.array(heights)[:, index]
            axes[0].bar(np.arange(len(bars)), values, bottom=bottom, label=label)
            bottom += values
        axes[0].set_xticks(np.arange(len(bars)),
                          [("C" if specs[r["task"]["bath_key"]]["kind"] == "cosh-grid" else "S")+
                           str(specs[r["task"]["bath_key"]]["levels"]) for r in bars])
        axes[0].set(yscale="log", ylabel="Stacked component medians (s)",
                    xlabel="C: cosh grid; S: surrogate; number: signed levels", title="Knight checkpoint cost breakdown")
        axes[0].legend(fontsize=7)
        for family in FAMILIES:
            selected = [r for r in points if specs[r["task"]["bath_key"]]["kind"] == family and
                        r["peak_process_rss_bytes"] is not None]
            axes[1].loglog([r["states"][0]["dimension"] for r in selected],
                           [r["peak_process_rss_bytes"]/2**20 for r in selected], "o-",
                           color=COLORS[family], label=LABELS[family])
        axes[1].set(xlabel="Sector Fock-space dimension", ylabel="Peak process RSS (MiB)",
                    title="Same checkpoint; unrestricted retained space")
        axes[1].legend(fontsize=8)
        for ax in axes:
            ax.grid(alpha=.2, axis="y", which="both")
        fig.tight_layout()
        save_figure(fig, directory/"cost_breakdown.svg")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    for ax, case in zip(axes, ("g0.5", "g20"), strict=True):
        inset = ax.inset_axes([.54, .48, .43, .43])
        for family in ("cosh-grid", "surrogate"):
            for parity in (0, 1) if family == "surrogate" else (0,):
                selected = sorted([r for r in rows if r["study"] == "spectrum" and r["case"] == case
                                   and r["family"] == family and r["levels"] % 2 == parity], key=lambda r: r["levels"])
                ax.plot([r["levels"] for r in selected], [r["value"] for r in selected],
                        "o--" if parity else "o-", color=COLORS[family],
                        label=LABELS[family]+(" odd" if parity else " even"))
                inset.plot([r["levels"] for r in selected], [r["value"] for r in selected],
                           "o--" if parity else "o-", color=COLORS[family], ms=3)
        reference = next((r for r in references if r["study"] == "spectrum" and r["case"] == case), None)
        if reference:
            ax.axhline(reference["value"], color="black", lw=.8, label="Reference (see convergence checks)")
            inset.axhline(reference["value"], color="black", lw=.8)
            if reference["stability_indicators"]["bath_step"] is None:
                ax.text(.02, .04, "Bath reference unresolved", transform=ax.transAxes, fontsize=8,
                        bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))
        edge = -1 if case == "g0.5" else 1
        ax.axhline(edge, color="gray", ls=":", label="Physical gap edge")
        inset.axhline(edge, color="gray", ls=":")
        inset.set(xlim=(5.5, 10.5), ylim=(-.94, -.90) if case == "g0.5" else (.96, 1.05))
        inset.tick_params(labelsize=7)
        inset.grid(alpha=.2)
        ax.set(xlabel="Signed levels", ylabel="(Eodd − Eeven) / Δ", title=f"U/Δ=15; Γ/Δ={case[1:]}")
        ax.grid(alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8)
    fig.tight_layout(rect=(0, .16, 1, 1))
    save_figure(fig, directory/"odd_even.svg")
    plt.close(fig)

    # Time-to-tolerance is a measured selection over the tested sizes, not interpolation.
    targets = []
    for study, case in panels:
        reference = next((r for r in references if r["study"] == study and r["case"] == case), None)
        if reference is None:
            continue
        indicators = reference["stability_indicators"]
        supported = all(indicators[k] is not None for k in ("bath_step", "bond_step"))
        spread = max((indicators[k] for k in ("bath_step", "bond_step") if indicators[k] is not None), default=0.)
        for tolerance in (1e-2, 1e-3, 1e-4, 1e-5):
            for family in FAMILIES:
                candidates = [r for r in rows if r["study"] == study and r["case"] == case and
                              r["family"] == family and r["levels"] % 2 == 0 and r["difference"] <= tolerance]
                for cost in ("cold_seconds", "reused_seconds"):
                    best = min(candidates, key=lambda r: r[cost]) if candidates else None
                    targets.append(dict(study=study, case=case, tolerance=tolerance, family=family, cost=cost,
                                        reference_resolved=bool(supported and spread < tolerance/5 and
                                                                reference["max_residual"] < tolerance/10),
                                        levels=None if best is None else best["levels"],
                                        seconds=None if best is None else best[cost],
                                        difference=None if best is None else best["difference"]))
    write_csv(directory/"time_to_accuracy.csv", targets)
    failures = [r for r in records if r["status"] != "ok"]
    unconverged = [r["task"] for r in records if r["status"] == "ok" and
                   r["task"]["type"] == "solve" and not r["converged"]]
    planned = [make_tasks(m["input"], m["profile"], m["include_dmrg"]) for m in manifests]
    planned_bath_keys = {key for batch, _ in planned for key in batch}
    planned_later = [task for _, tasks in planned for task in tasks]
    unavailable = [task for task in planned_later if task["bath_key"] not in bath_results]
    write_json(directory/"summary.json", dict(
        tasks=len(records), failures=len(failures), unconverged=unconverged,
        planned_tasks=len(planned_bath_keys)+len(planned_later), unavailable_bath_tasks=unavailable,
        supplemental_manifests=[str(Path(p)/"manifest.json") for p in supplements],
        worker_seconds_including_warmups=sum(r.get("worker_seconds", 0.) for r in records),
        empirical_references=references, time_to_accuracy=targets,
        interpretation="Stability indicators are empirical differences, not rigorous continuum error bounds."))

    table_rows = []
    for bandwidth in bands:
        for family in FAMILIES:
            record = next((r for r in defaults if r["task"]["spec"]["bandwidth"] == bandwidth and
                           r["task"]["spec"]["kind"] == family and r["task"]["spec"]["levels"] == 8), None)
            if record:
                table_rows.append([bandwidth, family, record["construction_seconds"]["median"],
                                   record["kernel"]["in_max_abs"], record["kernel"]["in_max_rel"],
                                   max(p["errors"]["excitation"] for p in record["quadratic"]),
                                   max(p["errors"]["current"] for p in record["quadratic"])])
    text = "# Bath-accuracy and calculation-time tables\n\nRegenerate with `tools/plot_bath_benchmarks.py`.\n\n"
    text += ("Energies and frequencies are in gap units, with $\\Delta=1$. "
             "Absolute kernel errors and the weak-hybridization Knight-shift slope have units of $1/\\Delta$; "
             "relative errors are dimensionless. "
             "Currents are in units of $2e\\Delta/\\hbar$. "
             "See [the calculation details](https://github.com/rokzitko/qdjj_solver/blob/main/docs/benchmarks/baths/README.md) "
             "for reference uncertainties and timing definitions.\n\n")
    text += "## Equal size: eight signed levels per reservoir\n\n"
    text += render_table(["Half-bandwidth", "Family", "Construction (s)", "Max absolute $g$ error",
                          "Max relative $g$ error", "Max $U=0$ excitation error", "Max $U=0$ current error"], table_rows)
    text += "\nThe $U=0$ maxima cover all three noninteracting parameter points.\n"
    text += "\n## Surrogate sensitivity at eight levels\n\n"
    sensitivity = []
    for record in bath_results.values():
        spec = record["task"]["spec"]
        if spec["kind"] != "surrogate" or spec["levels"] != 8 or spec["bandwidth"] == 10:
            continue
        sensitivity.append(dict(bath_key=record["task"]["bath_key"], **spec,
                                weak_spin_slope_error=record["kernel"]["weak_spin_slope_error"],
                                maximum_excitation_error=max(p["errors"]["excitation"] for p in record["quadratic"]),
                                maximum_current_error=max(p["errors"]["current"] for p in record["quadratic"]),
                                construction_seconds=record["construction_seconds"]["median"]))
    write_csv(directory/"fit_sensitivity.csv", sensitivity)
    text += render_table(["$D$", "Minimum frequency", "Maximum fit frequency", "Weight exponent",
                          "Frequency points", "Initial guesses", "Random seed",
                          "Weak spin slope error", "Max $U=0$ excitation error", "Max $U=0$ current error", "Fit (s)"],
                         [[r[k] for k in ("bandwidth", "frequency_min", "frequency_cutoff", "relative_weight",
                                         "frequency_points", "starts", "seed", "weak_spin_slope_error",
                                         "maximum_excitation_error", "maximum_current_error", "construction_seconds")]
                          for r in sensitivity])
    text += "\n## Interacting reference checks\n\n"
    text += render_table(["Study", "Point", "Reference levels", "Reference type", "Observable", "Value",
                          "Change with bath size", "Change with bond dimension", "Difference with another bath family", "Max residual"],
                         [[r["study"], r["case"], r["levels"], r["label"], r["observable"], r["value"],
                           *r["stability_indicators"].values(), r["max_residual"]] for r in references])
    text += "\n## Solver and bath refinements\n\n"
    text += "QP residuals with a finite cutoff are projected residuals; MPS residuals use the declared finite Hamiltonian.\n\n"
    text += render_table(["Study", "Point", "Bath family", "Levels", "Method", "Representation / ordering",
                          "QP cutoff", "MPS bond dimension", "Value", "Finite-problem checks passed", "Max residual"],
                         [[r[k] for k in ("study", "case", "family", "levels", "backend", "layout", "cutoff", "chi",
                                         "value", "converged", "max_residual")] for r in refinement_rows])
    for study, case in panels:
        if study == "knight":
            text += f"\n## Knight shift: {case}\n\n"
            text += "Reported value: $\\kappa$. The point name refers to the parameters in `input.json`.\n\n"
        else:
            text += f"\n## Parity excitation: {case}\n\n"
            text += ("Reported value: $(E_{\\mathrm{odd}}-E_{\\mathrm{even}})/\\Delta$. "
                     "The point name refers to the parameters in `input.json`.\n\n")
        text += render_table(["Bath family", "Levels", "Value", "Difference from reference",
                              "Time with saved bath (s)", "Time including bath construction (s)", "Peak memory (MiB)"],
                             [[r["family"], r["levels"], r["value"], r["difference"], r["reused_seconds"],
                               r["cold_seconds"], None if r["peak_process_rss_bytes"] is None else
                               r["peak_process_rss_bytes"]/2**20]
                              for r in rows if r["study"] == study and r["case"] == case])
    (directory/"tables.md").write_text(text, encoding="utf-8", newline="\n")
    print(f"Generated figures, tables and assessments in {directory}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="export directory containing measurements.json")
    parser.add_argument("--supplement", type=Path, action="append", default=[],
                        help="include an additional convergence archive (repeatable)")
    args = parser.parse_args()
    generate(args.directory, args.supplement)


if __name__ == "__main__":
    main()
