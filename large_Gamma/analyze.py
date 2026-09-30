"""Export compact, auditable large-Gamma evidence and observable-specific errors.

Incomplete calculations remain incomplete. In particular, an accurate solution
of a small bath is never promoted to a continuum reference.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import io
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qdjj_solver.common.io import fingerprint, write_json
from large_Gamma.study import VALUES, empirical_resolution, file_hash, point_key, points, read_json, validate_values


def collect(directory):
    from large_Gamma.run import load_task
    directory = Path(directory)
    manifest = read_json(directory/"manifest.json")
    identity = fingerprint(manifest)
    rows, failures, tasks = [], [], []
    baths = {}
    for path in sorted((directory/"bath_failures").glob("*.json")):
        failures.append(dict(reason="bath_construction_failed", **read_json(path)))
    for path in sorted((directory/"baths").glob("*.json")):
        record = read_json(path)
        if fingerprint(record["bath"]) != record["bath_key"]:
            raise ValueError("archived bath checksum mismatch")
        baths[record["bath_key"]] = record
    for path in sorted((directory/"tasks").glob("*.json")):
        result = load_task(path, identity)
        task = result["task"]
        if task["bath_key"] not in baths:
            raise ValueError("task references an absent bath")
        tasks.append(dict(task=task, complete=result["complete"], records_sha256=result["records_sha256"]))
        if result.get("interruption"):
            failures.append(dict(task=task, reason=result["interruption"]))
        for measured in result["records"]:
            if "method" not in measured:
                continue
            if measured["status"] != "ok":
                failures.append(dict(task=task, method=measured["method"], reason=measured["status"],
                                     message=measured.get("message")))
                continue
            validate_values(measured["values"])
            row = dict(measured, **task)
            row["point_id"] = point_key(task["point"])
            if row["method"] == "DMRG":
                row["accepted"] = row.get("observable_converged", row.get("bond_converged", False))
            rows.append(row)
    return manifest, baths, rows, failures, tasks


def index_rows(rows):
    groups = defaultdict(list)
    for row in rows:
        key = (row["point_id"], row["family"], row["bandwidth"], row["variant"], row["method"])
        groups[key].append(row)
    finest = {}
    for key, group in groups.items():
        by_bath = {}
        for row in sorted(group, key=lambda r: (r["levels"], r["accepted"], r.get("chi", 0))):
            by_bath[row["levels"]] = row
        finest[key] = [by_bath[level] for level in sorted(by_bath)]
    return finest


def method_resolution(row, key, config):
    exact_floor = config["convergence"].get("exact_reference_floor", config["convergence"]["target"])
    if row["method"] in ("ED", "BdG", "quadratic continuum"):
        return exact_floor
    if row.get("method_resolution") is not None:
        return max(exact_floor, row["method_resolution"][key])
    return max(config["convergence"].get("finite_target", config["convergence"]["target"]),
               config["convergence"]["safety_factor"]*(row.get("bond_changes") or {}).get(key, 0.))


def bath_resolution(history, config):
    if "finite_target" not in config["convergence"]:
        return empirical_resolution(history, config)
    from large_Gamma.study import refinement_changes
    steps = config["convergence"]["stable_steps"]
    changes = refinement_changes(history, steps)
    recent = history[-steps-1:]
    if changes is None or any(not r["accepted"] for r in recent):
        return None
    # QP residuals certify only the projected problem; cutoff error is measured
    # separately. Retained DMRG states have a per-observable method assessment.
    return {key: max(config["convergence"]["target"], config["convergence"]["safety_factor"]*
                     (changes[key]+max(method_resolution(r, key, config) if r["method"] == "DMRG" else 0.
                                       for r in recent))) for key in VALUES}


def compare(rows, config):
    """Matched-bath differences and separately qualified continuum comparisons."""
    histories = index_rows(rows)
    lookup = {}
    for history in histories.values():
        for row in history:
            lookup[(row["point_id"], row["bath_key"], row["variant"], row["method"])] = row
    comparisons, stability = [], []
    physical_refs = {}
    for key, history in histories.items():
        pid, family, bandwidth, variant, method = key
        if variant != "standard":
            continue
        resolution = bath_resolution(history, config)
        for name in VALUES:
            stability.append(dict(point_id=pid, **history[-1]["point"], family=family,
                                  bandwidth=bandwidth, method=method, levels=history[-1]["levels"],
                                  observable=name, resolution=None if resolution is None else resolution[name],
                                  target_met=resolution is not None and resolution[name] <= config["convergence"]["target"]))
        if method in ("DMRG", "ED") and resolution is not None:
            for name in VALUES:
                if resolution[name] <= config["convergence"]["target"]:
                    current = physical_refs.get((pid, bandwidth, name))
                    if current is None or current[0]["method"] != "quadratic continuum":
                        physical_refs[(pid, bandwidth, name)] = (history[-1], resolution[name])
        if method == "BdG":
            last = history[-1]
            for name in VALUES:
                physical_refs[(pid, bandwidth, name)] = (
                    dict(last, method="quadratic continuum", values=last["continuum"]["values"]),
                    config["convergence"].get("exact_reference_floor", config["convergence"]["target"]))

    def error_row(row, ref, name, scope, resolution, bath_resolved, method_resolution=0.):
        value, reference = row["values"][name], ref["values"][name]
        error = abs(value-reference)
        return dict(point_id=row["point_id"], **row["point"], family=row["family"],
                    bandwidth=row["bandwidth"], method=row["method"], levels=row["levels"],
                    observable=name, scope=scope, value=value, reference=reference,
                    reference_method=ref["method"], reference_levels=ref["levels"],
                    finite_accepted=row["accepted"],
                    absolute_error=error, reference_resolution=resolution,
                    comparison_resolution=resolution+method_resolution,
                    relative_error=error/abs(reference) if abs(reference) > 10*resolution else None,
                    qp_bath_resolved=bath_resolved)

    for key, history in histories.items():
        pid, _, bandwidth, variant, method = key
        if method not in ("2QP", "4QP", "6QP", "DMRG") or variant != "standard":
            continue
        resolution = bath_resolution(history, config)
        for row in history:
            refs = [lookup.get((pid, row["bath_key"], variant, name)) for name in ("BdG", "ED", "DMRG")]
            ref = next((r for r in refs if r is not None and r["accepted"] and r is not row), None)
            if ref is not None:
                for name in VALUES:
                    ref_resolution = method_resolution(ref, name, config)
                    comparisons.append(error_row(row, ref, name, "matched_bath", ref_resolution,
                                                 resolution is not None and resolution[name] <= config["convergence"]["target"]))
        row = history[-1]
        for name in VALUES:
            match = physical_refs.get((pid, bandwidth, name))
            if match is not None:
                ref, ref_resolution = match
                comparisons.append(error_row(row, ref, name, "continuum", ref_resolution,
                                             resolution is not None and resolution[name] <= config["convergence"]["target"],
                                             0. if resolution is None else resolution[name]))
    return comparisons, stability


def threshold_brackets(comparisons, thresholds):
    groups = defaultdict(list)
    for row in comparisons:
        if not row.get("finite_accepted", True):
            continue
        # Continuum threshold claims require the QP bath as well as the reference.
        if row["scope"] == "continuum" and not row["qp_bath_resolved"]:
            continue
        key = tuple(row[k] for k in ("u", "phi", "family", "bandwidth", "method", "observable", "scope"))
        key += (row["levels"] if row["scope"] == "matched_bath" else None,)
        groups[key].append(row)
    output = []
    names = ("u", "phi", "family", "bandwidth", "method", "observable", "scope", "levels")
    for key, rows in sorted(groups.items(), key=lambda item: str(item[0])):
        rows.sort(key=lambda row: row["gamma"])
        for threshold in thresholds:
            labels = ["above" if r["absolute_error"]-r.get("comparison_resolution", r["reference_resolution"]) > threshold else
                      "below" if r["absolute_error"]+r.get("comparison_resolution", r["reference_resolution"]) < threshold else
                      "unresolved" for r in rows]
            intervals = []
            if labels[0] == "above":
                intervals.append((None, rows[0]["gamma"], "already_above_at_first_sample"))
            for i in range(1, len(rows)):
                if labels[i] != labels[i-1]:
                    intervals.append((rows[i-1]["gamma"], rows[i]["gamma"], f"{labels[i-1]}_to_{labels[i]}"))
            if not intervals:
                intervals = [(None, None, f"all_{labels[0]}_on_sampled_points")]
            for lower, upper, status in intervals:
                output.append(dict(zip(names, key, strict=True), threshold=threshold,
                                   gamma_lower=lower, gamma_upper=upper, status=status,
                                   sample_min=rows[0]["gamma"], sample_max=rows[-1]["gamma"],
                                   samples=len(rows)))
    return output


def refinement_points(comparisons, config):
    """Refine certified threshold brackets without presuming monotonic errors."""
    proposed = {}
    for row in threshold_brackets(comparisons, config["convergence"]["thresholds"]):
        lower, upper = row["gamma_lower"], row["gamma_upper"]
        if lower is None or upper is None or row["method"] == "DMRG":
            continue
        width = max(config["refinement"]["absolute_gamma_width"],
                    config["refinement"]["relative_gamma_width"]*(lower+upper)/2)
        if upper-lower > width:
            point = dict(u=row["u"], phi=row["phi"], gamma=float((lower+upper)/2))
            proposed[point_key(point)] = point
    # Also resolve singlet-doublet crossings using reference branch energies.
    gaps = defaultdict(dict)
    for row in comparisons:
        if row["observable"] == "signed_gap" and row["scope"] == "matched_bath":
            key = (row["u"], row["phi"], row["family"], row["bandwidth"], row["levels"])
            gaps[key][row["gamma"]] = row["reference"]
    for key, grid in gaps.items():
        samples = sorted(grid.items())
        for (a, ea), (b, eb) in zip(samples, samples[1:], strict=False):
            width = max(config["refinement"]["absolute_gamma_width"],
                        config["refinement"]["relative_gamma_width"]*(a+b)/2)
            if ea*eb < 0 and b-a > width:
                point = dict(u=key[0], phi=key[1], gamma=float((a+b)/2))
                proposed[point_key(point)] = point
    return list(proposed.values())


def compact_comparisons(comparisons):
    """Wide scalar tables at the finest matched bath; all baths remain in observables.csv."""
    latest = {}
    for row in comparisons:
        key = tuple(row[k] for k in ("point_id", "family", "bandwidth", "method", "observable", "scope"))
        if key not in latest or row["levels"] > latest[key]["levels"]:
            latest[key] = row
    wide = {}
    identifying = ("point_id", "u", "gamma", "phi", "family", "bandwidth", "method", "levels", "scope", "finite_accepted",
                   "reference_method", "reference_levels")
    for row in latest.values():
        key = tuple(row[k] for k in identifying)
        result = wide.setdefault(key, {k: row[k] for k in identifying})
        for source, prefix in (("absolute_error", "error"), ("reference", "reference"),
                               ("comparison_resolution", "resolution"), ("relative_error", "relative_error"),
                               ("qp_bath_resolved", "bath_resolved")):
            result[f"{prefix}_{row['observable']}"] = row[source]
    fields = list(identifying)+[f"{prefix}_{name}" for name in VALUES
                               for prefix in ("error", "reference", "resolution", "relative_error", "bath_resolved")]
    return [{key: row.get(key) for key in fields} for row in wide.values()]


def control_checks(rows, config):
    histories = index_rows(rows)
    lookup = {(r["point_id"], r["bath_key"], r["method"], r["variant"]): r
              for history in histories.values() for r in history}
    checks = []
    for history in histories.values():
        for row in history:
            variant = row["variant"]
            if variant in ("tight", "centered", "uncompressed"):
                other = lookup.get((row["point_id"], row["bath_key"], row["method"], "standard"))
                if other is not None:
                    checks.append(dict(point=row["point"], bath_key=row["bath_key"], method=row["method"],
                                       check=variant, both_finite_accepted=row["accepted"] and other["accepted"],
                                       absolute_changes={k: abs(row["values"][k]-other["values"][k]) for k in VALUES}))
            if variant == "phase_plus":
                anchor = row["anchor"]
                base = lookup.get((point_key(anchor), row["bath_key"], row["method"], "standard"))
                minus_point = dict(anchor, phi=anchor["phi"]-config["controls"]["phase_step"])
                minus = lookup.get((point_key(minus_point), row["bath_key"], row["method"], "phase_minus"))
                if base is not None and minus is not None:
                    differences = {branch: abs((row["values"][f"energy_{branch}"]-minus["values"][f"energy_{branch}"])
                                   /(2*config["controls"]["phase_step"])-base["values"][f"current_{branch}"])
                                   for branch in ("even", "odd")}
                    checks.append(dict(point=anchor, bath_key=row["bath_key"], method=row["method"],
                                       check="phase_derivative", absolute_changes=differences,
                                       both_finite_accepted=all(r["accepted"] for r in (row, minus, base))))
    return checks


def write_csv(path, rows, fields=None, max_bytes=2*1024**2):
    """Bound individual text files; a dataset catalog keeps shards transparent."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else ["status"]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    header = buffer.getvalue()
    chunk, size, paths = [header], len(header.encode("utf-8")), []

    def publish():
        target = path if not paths else path.with_name(f"{path.stem}.part{len(paths)+1:03d}{path.suffix}")
        target.write_text("".join(chunk), encoding="utf-8", newline="")
        paths.append(target.name)

    for row in rows:
        buffer.seek(0)
        buffer.truncate(0)
        writer.writerow(row)
        text = buffer.getvalue()
        length = len(text.encode("utf-8"))
        if length+len(header.encode("utf-8")) > max_bytes:
            raise ValueError("a scalar CSV row exceeds the archive file budget")
        if size+length > max_bytes and len(chunk) > 1:
            publish()
            chunk, size = [header], len(header.encode("utf-8"))
        chunk.append(text)
        size += length
    publish()
    return paths


def export(directory, destination, *, plots=False):
    manifest, baths, rows, failures, tasks = collect(directory)
    config = manifest["input"]
    destination = Path(destination)
    comparisons, stability = compare(rows, config)
    histories = index_rows(rows)
    measurements, bonds, weights = [], [], []
    for history in histories.values():
        for row in history:
            sectors = row.get("sectors", [])
            trial_spreads = [max(t["energy"] for t in s["trials"])-min(t["energy"] for t in s["trials"])
                             for s in sectors if s.get("trials")]
            measured = dict(point_id=row["point_id"], **row["point"], family=row["family"],
                            bandwidth=row["bandwidth"], bath_key=row["bath_key"], levels=row["levels"],
                            variant=row["variant"], method=row["method"], chi=row.get("chi"),
                            finite_accepted=row["accepted"], seconds=row.get("seconds"),
                            accuracy_evidence=row.get("accuracy_evidence"),
                            method_resolution_max=max((row.get("method_resolution") or {}).values(), default=None),
                            state_residual_converged=all(s["finite_problem_converged"] for s in sectors)
                            if row["method"] == "DMRG" else None,
                            max_residual=max((s["residual"] for s in sectors), default=None),
                            max_dimension=max((s.get("dimension", 0) for s in sectors), default=None),
                            max_rss_bytes=max((s.get("peak_rss_bytes") or 0 for s in sectors), default=None),
                            chi_actual=max((s.get("chi_actual", 0) for s in sectors), default=None),
                            max_sweeps=max((s.get("sweeps", 0) for s in sectors), default=None),
                            seed_energy_spread=max(trial_spreads, default=None),
                            dark_threshold=row.get("dark_threshold"), full_space=row.get("full_space"),
                            **{f"energy_p{p}_eta{eta:+d}": next((s["energy"] for s in sectors
                                if s["sector"]["parity"] == p and s["sector"]["eta"] == eta), None)
                               for p in (0, 1) for eta in (1, -1)}, **row["values"])
            gap = row["values"]["signed_gap"]
            measured.update(ground_energy=min(row["values"]["energy_even"], row["values"]["energy_odd"]),
                            ground_parity="unresolved" if abs(gap) <= config["convergence"]["target"]
                            else "odd" if gap < 0 else "even",
                            ground_current=None if abs(gap) <= config["convergence"]["target"] else
                            row["values"]["current_odd" if gap < 0 else "current_even"],
                            continuum_edge_margin=1-abs(gap),
                            dark_edge_margin=None if row.get("dark_threshold") is None else row["dark_threshold"]-abs(gap))
            measurements.append(measured)
            if row is history[-1] and row["method"] in ("2QP", "4QP", "6QP"):
                for parity in (0, 1):
                    state = min((s for s in sectors if s["sector"]["parity"] == parity), key=lambda s: s["energy"])
                    weights.append(dict(point_id=row["point_id"], bath_key=row["bath_key"],
                                        variant=row["variant"], method=row["method"],
                                        parity=parity, eta=state["sector"]["eta"],
                                        **{f"w{q}": state["qp_weights"][q] if q < len(state["qp_weights"]) else 0.
                                           for q in range(7)}))
    for row in rows:
        if row["method"] == "DMRG":
            bonds.append(dict(point_id=row["point_id"], bath_key=row["bath_key"], variant=row["variant"],
                              chi=row["chi"], finite_converged=all(s["finite_problem_converged"] for s in row["sectors"]),
                              bond_converged=row.get("bond_converged", False),
                              observable_converged=row.get("observable_converged", row.get("bond_converged", False)),
                              accuracy_evidence=row.get("accuracy_evidence"),
                              max_residual=max(s["residual"] for s in row["sectors"]),
                              max_observable_step=max((row.get("bond_changes") or {}).values(), default=None),
                              seconds=row["seconds"], **row["values"]))
    # Deterministic output permits meaningful reviews between exports.
    measurements.sort(key=lambda r: (r["phi"], r["u"], r["gamma"], r["family"], r["bandwidth"],
                                     r["variant"], r["levels"], r["method"]))
    compact = compact_comparisons(comparisons)
    tables = {"observables.csv": measurements, "convergence.csv": stability, "comparisons.csv": compact,
              "breakdown.csv": threshold_brackets(comparisons, config["convergence"]["thresholds"]),
              "bond_convergence.csv": bonds, "qp_weights.csv": weights}
    datasets = {name: write_csv(destination/name, data) for name, data in tables.items()}
    previous = destination/"manifest.json"
    if previous.exists():
        old_manifest = read_json(previous)
        for name, old_parts in old_manifest.get("datasets", {}).items():
            for part in set(old_parts)-set(datasets.get(name, [])):
                path = destination/part
                if path.is_file():
                    if file_hash(path) != old_manifest["artifacts"][part]["sha256"]:
                        raise ValueError("a previously generated CSV shard was modified")
                    path.unlink()
    write_json(destination/"baths.json", baths)
    write_json(destination/"failures.json", failures)
    write_json(destination/"validation.json", control_checks(rows, config))
    expected = {point_key(p) for p in points(config)}
    coverage = {}
    for method in ("2QP", "4QP", "6QP", "DMRG"):
        available = {r["point_id"] for r in rows if r["method"] == method and r["variant"] == "standard"
                     and r["bandwidth"] == config["model"]["bandwidth"]}
        accepted = {r["point_id"] for r in rows if r["method"] == method and r["variant"] == "standard"
                    and r["bandwidth"] == config["model"]["bandwidth"] and r["accepted"]}
        # Every observable must be converged within one common bath family/bandwidth.
        resolved_groups = defaultdict(set)
        for r in stability:
            if r["method"] == method and r["target_met"] and r["bandwidth"] == config["model"]["bandwidth"]:
                resolved_groups[(r["point_id"], r["family"])].add(r["observable"])
        fully = {pid for (pid, _), names in resolved_groups.items() if names == set(VALUES)}
        coverage[method] = dict(measured_grid_points=len(expected & available),
                                finite_accepted_grid_points=len(expected & accepted),
                                all_observables_bath_converged_grid_points=len(expected & fully))
    complete = all(c["all_observables_bath_converged_grid_points"] == len(expected) for c in coverage.values())
    summary = dict(format="qdjj-large-Gamma-summary", format_version=1, complete=complete,
                   status="convergence target achieved on base grid" if complete else "incomplete convergence study",
                   expected_base_grid_points=len(expected), measured_rows=len(measurements),
                   coverage=coverage, failures=len(failures), completed_tasks=sum(t["complete"] for t in tasks),
                   partial_tasks=sum(not t["complete"] for t in tasks),
                   empirical_target=config["convergence"]["target"],
                   finite_observable_target=config["convergence"].get("finite_target", config["convergence"]["target"]),
                   resolution_rule="max(target, safety_factor*(max of last two bath changes + bond changes)); empirical",
                   measured_worker_seconds=sum(r.get("seconds", 0.) for r in rows),
                   maximum_resident_bytes=max((r["max_rss_bytes"] or 0 for r in measurements), default=0),
                   next_refinement_points=refinement_points(comparisons, config))
    write_json(destination/"summary.json", summary)
    if plots and measurements:
        render(expand_comparisons(compact), destination/"figures")
    artifacts = {p.relative_to(destination).as_posix(): dict(sha256=file_hash(p), bytes=p.stat().st_size)
                 for p in sorted(destination.rglob("*")) if p.is_file() and p.name != "manifest.json"}
    write_json(destination/"manifest.json", dict(source=manifest, analyzer_sha256=file_hash(__file__),
                                                task_checksums={fingerprint(t["task"]):
                                                    [t["records_sha256"], t["complete"]] for t in tasks},
                                                task_checksum_columns=["records_sha256", "complete"],
                                                artifacts=artifacts, datasets=datasets,
                                                measurements_sha256=fingerprint(measurements)))
    return summary


def expand_comparisons(compact):
    expanded = []
    for row in compact:
        for name in VALUES:
            if row.get(f"error_{name}") in (None, ""):
                continue
            expanded.append(dict(u=float(row["u"]), phi=float(row["phi"]), gamma=float(row["gamma"]),
                                 bandwidth=float(row["bandwidth"]), levels=int(row["levels"]),
                                 scope=row["scope"], family=row["family"], method=row["method"], observable=name,
                                 absolute_error=float(row[f"error_{name}"]),
                                 reference_resolution=float(row[f"resolution_{name}"])))
    return expanded


def verify_archive(directory):
    directory = Path(directory)
    manifest = read_json(directory/"manifest.json")
    for name, record in manifest["artifacts"].items():
        path = directory/name
        if (Path(name).is_absolute() or ".." in Path(name).parts or not path.is_file()
                or path.stat().st_size != record["bytes"] or file_hash(path) != record["sha256"]):
            raise ValueError(f"compact archive checksum mismatch: {name}")
    return expand_comparisons(archive_rows(directory, "comparisons.csv"))


def archive_rows(directory, dataset):
    directory = Path(directory)
    manifest = read_json(directory/"manifest.json")
    rows = []
    for name in manifest.get("datasets", {}).get(dataset, [dataset]):
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("nonrelative archive table")
        with (directory/name).open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def render(comparisons, directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "svg.hashsalt": "large-Gamma-study"})
    for phi in sorted({r["phi"] for r in comparisons}):
        chosen = [r for r in comparisons if r["phi"] == phi and r["bandwidth"] == 100
                  and r["observable"] in ("signed_gap", "q_d", "C_d_bath")
                  and r["method"] in ("2QP", "4QP", "6QP") and r["scope"] == "matched_bath"
                  and r["family"] == "cosh-grid"]
        us = sorted({r["u"] for r in chosen})
        if not us:
            continue
        fig, axes = plt.subplots(len(us), 3, figsize=(11, 2.1*len(us)), squeeze=False)
        for i, u in enumerate(us):
            for j, name in enumerate(("signed_gap", "q_d", "C_d_bath")):
                ax = axes[i, j]
                for method, color in zip(("2QP", "4QP", "6QP"), ("C0", "C1", "C2"), strict=True):
                    group = [r for r in chosen if r["u"] == u and r["observable"] == name and r["method"] == method]
                    levels = sorted({r["levels"] for r in group})
                    level_label = ",".join(map(str, levels))
                    group = sorted(group, key=lambda r: r["gamma"])
                    ax.loglog([r["gamma"] for r in group],
                              [max(r["absolute_error"], r["reference_resolution"]) for r in group],
                              ".-", color=color, label=f"{method}, L={level_label}")
                ax.axhline(1e-4, color=".6", linestyle=":")
                ax.set_title(f"U/Delta={u:g}: {name}")
                ax.set_ylabel("Matched-bath absolute error")
                ax.grid(alpha=.2)
                ax.legend(fontsize=7)
                if i == len(us)-1:
                    ax.set_xlabel("Gamma/Delta")
        fig.suptitle(f"phi={phi:.6g}; finest matched bath; errors clipped at recorded reference resolution")
        fig.tight_layout(rect=(0, 0, 1, .98))
        path = directory/f"cutoff_errors_phi_{phi:.6g}.svg"
        fig.savefig(path, metadata={"Date": None})
        # Matplotlib leaves spaces at the ends of SVG path-coordinate lines.
        text = "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines())
        path.write_text(text+"\n", encoding="utf-8")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT/"large_Gamma/output")
    parser.add_argument("--plots", action="store_true")
    parser.add_argument("--replot", action="store_true", help="verify/replot a compact archive without raw runs or solvers")
    args = parser.parse_args()
    if args.replot:
        comparisons = verify_archive(args.source)
        render(comparisons, args.output/"figures")
        if args.source.resolve() == args.output.resolve():
            # Preserve the numerical analysis provenance when only regenerating figures.
            manifest = read_json(args.source/"manifest.json")
            manifest["plotter_sha256"] = file_hash(__file__)
            for phi in {row["phi"] for row in comparisons}:
                name = f"figures/cutoff_errors_phi_{phi:.6g}.svg"
                path = args.output/name
                if path.is_file():
                    manifest["artifacts"][name] = dict(sha256=file_hash(path), bytes=path.stat().st_size)
            write_json(args.output/"manifest.json", manifest)
        print("Verified compact archive and regenerated figures")
        return
    summary = export(args.source, args.output, plots=args.plots)
    print(summary["status"], summary["coverage"])


if __name__ == "__main__":
    main()
