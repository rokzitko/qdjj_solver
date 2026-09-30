"""Validate, tabulate and plot the self-contained doublet cutoff archive.

Matplotlib is needed only by render(), not by the numerical archive checks.
All reference resolutions are empirical observable-by-observable stability
estimates. In particular, a legacy DMRG variance is not a residual certificate.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from qdjj_solver.common.io import write_json
from tools.benchmark_qp_convergence import VALUES, load_archive, validate_config

LABELS = {"energy": r"$E/\Delta$", "q_d": r"$q_d$", "P_0": r"$P_0$",
          "P_1": r"$P_1$", "P_2": r"$P_2$", "C_d_bath": r"$\langle\mathbf{S}_d\cdot\mathbf{S}_{\rm bath}\rangle$"}


def matching_baths(left, right):
    return (all(left[k] == right[k] for k in ("bandwidth", "delta"))
            and len(left["xi"]) == len(right["xi"])
            and all(np.allclose(left[k], right[k], rtol=2e-14, atol=1e-14) for k in ("xi", "weights")))


def compare_current(directory, records, baths):
    _, fresh_baths, fresh = load_archive(directory)
    archived = {r["id"]: r for r in records}
    differences, reference_differences = [], []
    for row in fresh:
        old = archived[row["id"]]
        if (row["origin"] != "current" or row["model"] != old["model"]
                or not matching_baths(fresh_baths[row["bath_key"]], baths[old["bath_key"]])):
            raise ValueError("verification calculation has different physics or coordinates")
        difference = {k: abs(row["values"][k]-old["values"][k]) for k in VALUES}
        if row["backend"] == "dmrg":
            # The legacy reference targeted both eta sectors; the bound doublet
            # is in eta=+1. Preserve this distinction in the actual records.
            if (max(difference.values()) > 1e-8 or not row["diagnostics"]["finite_problem_converged"]
                    or any(row["sector"][k] != old["sector"][k] for k in ("parity", "twice_sz"))):
                raise ValueError("current DMRG does not reproduce the empirical reference")
            reference_differences.append(dict(id=row["id"], absolute_changes=difference,
                                              full_residual=row["residual"]))
            continue
        if row["sector"] != old["sector"]:
            raise ValueError("verification calculation changed the QP sector")
        weight_difference = float(np.max(abs(np.array(row["qp_weights"])-old["qp_weights"])))
        if max(*difference.values(), weight_difference, row["residual"]) > 2e-9:
            raise ValueError("current solver does not reproduce the archived QP calculation")
        differences.append(dict(id=row["id"], absolute_changes=difference,
                                maximum_qp_weight_change=weight_difference, residual=row["residual"]))
    if not differences:
        raise ValueError("verification archive has no matching QP points")
    return dict(points=differences, reference_points=reference_differences,
                maximum_absolute_changes={k: max(r["absolute_changes"][k] for r in differences) for k in VALUES},
                maximum_qp_weight_change=max(r["maximum_qp_weight_change"] for r in differences))


def analyze(directory, verification=None):
    manifest, baths, records = load_archive(directory)
    config = manifest["input"]
    validate_config(config)
    by_id = {r["id"]: r for r in records}
    levels = sorted(config["levels"])
    if len(levels) != 2 or len(config["references"]) != 3:
        raise ValueError("the report requires two baths and a three-point bath/bond reference check")
    expected = {f"qp_N{n}_Q{q}" for n in levels for q in config["cutoffs"]}
    expected.update(f"dmrg_N{r['levels']}_chi{r['chi_max']}" for r in config["references"])
    if expected != set(by_id):
        raise ValueError("incomplete or unexpected cutoff/reference points in archive")
    for row in records:
        if row["model"] != config["model"]:
            raise ValueError("mixed physical models in cutoff archive")
        n = row["levels"]
        baseline = by_id[f"qp_N{n}_Q{config['cutoffs'][0]}"]
        if not matching_baths(baths[row["bath_key"]], baths[baseline["bath_key"]]):
            raise ValueError("cutoff sequence or reference changed the finite bath")
        if row["backend"] == "qp":
            if row["sector"] != baseline["sector"] or row["residual"] > 2e-9:
                raise ValueError("inconsistent sector or excessive projected residual")
        elif (row["cutoff"] is not None or
              (row["origin"] == "paper" and not row["diagnostics"].get("sweep_converged")) or
              (row["origin"] == "current" and not row["diagnostics"].get("finite_problem_converged"))):
            raise ValueError("an unconverged or QP-truncated state cannot be a reference")
    reference_rows = [by_id[f"dmrg_N{r['levels']}_chi{r['chi_max']}"] for r in config["references"]]
    base, bath_check, reference = reference_rows
    if (base["chi_max"] != bath_check["chi_max"] or base["levels"] != reference["levels"]
            or bath_check["levels"] <= base["levels"] or reference["chi_max"] <= base["chi_max"]):
        raise ValueError("reference checks must vary bath and bond dimension independently")
    stability = {}
    for key in VALUES:
        bath_step = abs(bath_check["values"][key]-base["values"][key])
        bond_step = abs(reference["values"][key]-base["values"][key])
        stability[key] = dict(bath_step=bath_step, bond_step=bond_step,
                              empirical_resolution=max(config["reference_stability_floors"][key],
                                                       4*(bath_step+bond_step)))
    points = []
    previous = dict.fromkeys(levels, np.inf)
    for q in config["cutoffs"]:
        coarse, fine = [by_id[f"qp_N{n}_Q{q}"] for n in levels]
        for row in (coarse, fine):
            if row["values"]["energy"] > previous[row["levels"]]+2e-9:
                raise ValueError("fixed-bath Ritz energies increase with the cutoff")
            previous[row["levels"]] = row["values"]["energy"]
        errors = {k: abs(fine["values"][k]-reference["values"][k]) for k in VALUES}
        grid = {k: abs(fine["values"][k]-coarse["values"][k]) for k in VALUES}
        resolved = {k: errors[k] > stability[k]["empirical_resolution"] for k in VALUES}
        points.append(dict(cutoff=q, **fine["values"], absolute_errors=errors, bath_changes=grid,
                           error_resolved=resolved, dimension=fine["dimension"], residual=fine["residual"],
                           qp_weights=fine["qp_weights"]))
    for key in VALUES:
        stability[key]["maximum_qp_bath_change"] = max(p["bath_changes"][key] for p in points)
        stability[key]["qp_bath_change_below_resolution"] = (
            stability[key]["maximum_qp_bath_change"] < stability[key]["empirical_resolution"]/5)
    summary = dict(qp_levels=levels[-1], bath_check_levels=levels[0],
                   reference_levels=reference["levels"], reference_chi=reference["chi_max"],
                   reference_id=reference["id"], reference_origin=reference["origin"],
                   reference_values=reference["values"], reference_diagnostics=reference["diagnostics"],
                   reference_residual_kind=reference["residual_kind"],
                   resolution_rule="max(stated floor, 4*(bath step + bond step)); empirical, not an error bound",
                   stability=stability, points=points)
    if verification is not None:
        summary["current_solver_verification"] = compare_current(verification, records, baths)
    return summary


def write_tables(output, summary):
    output = Path(output)
    rows = []
    for point in summary["points"]:
        for key in VALUES:
            rows.append(dict(cutoff=point["cutoff"], observable=key, value=point[key],
                             reference=summary["reference_values"][key],
                             absolute_difference=point["absolute_errors"][key],
                             empirical_resolution=summary["stability"][key]["empirical_resolution"],
                             error_resolved=point["error_resolved"][key],
                             bath_change=point["bath_changes"][key],
                             dimension=point["dimension"], projected_residual=point["residual"]))
    with (output/"convergence.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Doublet cutoff convergence", "", "Generated by `tools/plot_qp_convergence.py`.", "",
             "| Cutoff | $E/\\Delta$ | $q_d$ | $P_0$ | $P_1$ | $P_2$ | $C_{d,\\mathrm{bath}}$ |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for row in summary["points"]:
        lines.append(f"| {row['cutoff']} | " + " | ".join(f"{row[k]:.10f}" for k in VALUES) + " |")
    lines += ["| DMRG reference | " + " | ".join(f"{summary['reference_values'][k]:.10f}" for k in VALUES) + " |",
              "", f"QP rows use {summary['qp_levels']} signed normal levels per reservoir; the selected reference uses",
              f"{summary['reference_levels']} levels and bond dimension {summary['reference_chi']}. "
              "Printed digits identify the records, not an accuracy claim.",
              "The correlation is inferred using the SU(2) doublet identity.", "",
              "## Independent bath and bond refinement", "",
              "| Observable | DMRG bath step | DMRG bond step | Empirical resolution | Maximum QP bath step |",
              "|---|---:|---:|---:|---:|"]
    for key, row in summary["stability"].items():
        lines.append(f"| `{key}` | {row['bath_step']:.3e} | {row['bond_step']:.3e} | "
                     f"{row['empirical_resolution']:.1e} | {row['maximum_qp_bath_change']:.3e} |")
    lines += ["", "Energy differences use units of the gap; other differences are dimensionless.",
              "Resolutions are empirical stability estimates, not rigorous error bounds.", ""]
    (output/"table.md").write_text("\n".join(lines), encoding="utf-8")


def render(summary, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "svg.hashsalt": "qdjj-qp-convergence"})
    points = summary["points"]
    cutoffs = [p["cutoff"] for p in points]

    def finish(fig, name):
        fig.tight_layout()
        path = output/f"{name}.svg"
        fig.savefig(path, metadata={"Date": None})
        # Matplotlib leaves spaces at the ends of SVG path-coordinate lines.
        text = "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())
        path.write_text(text+"\n", encoding="utf-8")
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for ax, names in zip(axes.flat, (("energy",), ("q_d",), ("P_0", "P_1", "P_2"), ("C_d_bath",)), strict=True):
        for key, marker in zip(names, ("o", "s", "x"), strict=False):
            line, = ax.plot(cutoffs, [p[key] for p in points], marker=marker, label=LABELS[key])
            ax.axhline(summary["reference_values"][key], color=line.get_color(), linestyle="--", alpha=.6)
        ax.set(xlabel="Maximum bath-QP number Q", xticks=cutoffs)
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle("Fixed-bath cutoff convergence; dashed lines: unrestricted DMRG reference", fontsize=11)
    finish(fig, "observables")

    fig, axes = plt.subplots(2, 3, figsize=(11, 6))
    for ax, key in zip(axes.flat, VALUES, strict=True):
        floor = summary["stability"][key]["empirical_resolution"]
        error = np.array([p["absolute_errors"][key] for p in points])
        resolved = np.array([p["error_resolved"][key] for p in points])
        ax.semilogy(cutoffs, np.maximum(error, floor), color="C0", linewidth=1)
        ax.scatter(np.array(cutoffs)[resolved], error[resolved], color="C0", s=25, zorder=3)
        ax.scatter(np.array(cutoffs)[~resolved], np.full(sum(~resolved), floor),
                   facecolors="white", edgecolors="C0", s=35, zorder=4)
        ax.semilogy(cutoffs, np.maximum([p["bath_changes"][key] for p in points], 1e-12),
                    ":", color=".4", marker="+",
                    label=f"Bath change ({summary['bath_check_levels']} to {summary['qp_levels']})")
        ax.axhspan(1e-12, floor, color=".9", label="Empirical reference resolution")
        ax.set(title=LABELS[key], xlabel="Maximum bath-QP number Q", ylabel="Absolute difference",
               xticks=cutoffs, ylim=(1e-12, max(error.max()*3, floor*100)))
        ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=7, loc="upper right")
    fig.suptitle("Open markers: unresolved differences placed at the reference resolution", fontsize=11)
    finish(fig, "errors")

    fig, ax = plt.subplots(figsize=(7, 4))
    for point in points:
        weights = np.array(point["qp_weights"])
        ax.semilogy(np.arange(len(weights)), weights, "o-", label=f"Q={point['cutoff']}")
    ax.set(xlabel="Bath-QP number n", ylabel=r"Probability $w_n$", xticks=range(max(cutoffs)+1),
           title="QP-number probabilities of the normalized truncated states")
    ax.legend(ncol=3)
    ax.grid(alpha=.2)
    finish(fig, "weights")
    write_json(output/"summary.json", summary)
    write_tables(output, summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--verification", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT/"results/qp-convergence-plots")
    args = parser.parse_args()
    summary = analyze(args.archive, args.verification)
    render(summary, args.output)
    print(f"Wrote tables, checked summary and three figures to {args.output}")


if __name__ == "__main__":
    main()
