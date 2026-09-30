"""Offline error budgets and one-to-one matching of the published subgap levels."""

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from qdjj_solver import DiscreteBath
from qdjj_solver.common.baths import discrete_g, hybridization_g
from applications._common import read_csv, read_json, write_json
from applications.analyze import differences

CASE = Path(__file__).resolve().parent


def distinct(values, tolerance):
    result = []
    for value in sorted(values):
        if not result or value-result[-1] > tolerance:
            result.append(value)
    return np.array(result)


def smooth_interpolate(x, rows, key, value, jump):
    """Read a simplified published polyline without interpolating a phase jump."""
    rows = sorted(rows, key=lambda r: float(r[key]))
    xs, ys = np.array([float(r[key]) for r in rows]), np.array([float(r[value]) for r in rows])
    if len(rows) < 2 or not xs[0] <= x <= xs[-1]:
        return None
    right = min(max(1, int(np.searchsorted(xs, x))), len(xs)-1)
    if abs(ys[right]-ys[right-1]) > jump:
        return None
    return float(np.interp(x, xs, ys))


def analyze():
    output = CASE/"output"
    current = read_csv(output/"convergence"/"current.csv")
    shared = read_csv(output/"convergence"/"shared.csv")
    comparison8 = read_csv(output/"paper"/"current_comparison.csv")
    comparison17 = read_csv(output/"paper"/"shared_comparison.csv")
    result = dict(paper=read_json(output/"paper"/"summary.json"), current_convergence={}, shared_convergence={})
    for a, b in (("ChE-W4", "ChE-W6"), ("ChE-W6", "ChE-W8"), ("ChE-F4", "ChE-F6"),
                 ("ChE-F6", "ChE-F8"), ("ChE-W8", "ChE-F8"), ("ChE-F8-q4", "ChE-F8"),
                 ("ChE-F8", "SM8"), ("SM8", "SM10")):
        result["current_convergence"][f"{a}_versus_{b}"] = differences(
            current, a, b, ("u", "phi_over_pi"), ("current", "signed_gap"))
    for a, b in (("ChE-W4", "ChE-W6"), ("ChE-W6", "ChE-W8"), ("ChE-W10", "ChE-W8"),
                 ("ChE-W8-q4", "ChE-W8"), ("ChE-W8-q6", "ChE-W8"), ("ChE-W8", "ChE-F8"),
                 ("ChE-F8", "SM8"), ("SM10", "SM8")):
        result["shared_convergence"][f"{a}_versus_{b}"] = differences(
            shared, a, b, ("u",), ("doublet_singlet_gap", "triplet_singlet_gap", "pairing", "spin_correlation"))
    result["current_reference"] = {}
    for setting in sorted({r["setting"] for r in comparison8}):
        data = [r for r in comparison8 if r["setting"] == setting and r["reference"] == "NRG"]
        same = [r for r in data if r["same_sign"] == "1"]
        result["current_reference"][setting] = dict(
            points=len(data), maximum_error=max(abs(float(r["difference"])) for r in data),
            maximum_same_sign_error=max(abs(float(r["difference"])) for r in same),
            opposite_sign_points=[dict(u=float(r["u"]), phi_over_pi=float(r["phi_over_pi"]))
                                  for r in data if r["same_sign"] == "0"],
            maximum_doublet_error=max(abs(float(r["difference"])) for r in data if float(r["reference_current"]) < -1e-5))
    result["shared_reference"] = {}
    for setting in sorted({r["setting"] for r in comparison17}):
        data = [r for r in comparison17 if r["setting"] == setting]
        metrics = {}
        # Ground-state scalar observables jump at different transition positions.
        # The NRG phase is identified by its pairing/correlation signs.
        for observable in ("pairing", "spin_correlation"):
            selected = [r for r in data if r["observable"] == observable]
            same = []
            for row in selected:
                u = float(row["u"])
                nu = next(float(r["reference_value"]) for r in data if r["observable"] == "pairing" and abs(float(r["u"])-u) < 1e-5)
                corr = next(float(r["reference_value"]) for r in data if r["observable"] == "spin_correlation" and abs(float(r["u"])-u) < 1e-5)
                nrg_spin = .5 if nu > 0 else (0. if corr < 0 else 1.)
                if nrg_spin == float(row["ground_spin"]):
                    same.append(row)
            metrics[observable] = dict(maximum_error=max(abs(float(r["difference"])) for r in selected),
                                        maximum_same_phase_error=max(abs(float(r["difference"])) for r in same),
                                        same_phase_points=len(same), total_points=len(selected))
        result["shared_reference"][setting] = metrics

    # Match distinct displayed energy levels one-to-one: independent nearest
    # neighbours could match several reference levels to one computed state.
    spectrum = read_csv(output/"paper"/"spectrum.csv")
    reference = read_csv(CASE/"reference"/"figure17.csv")
    full_energies = {}
    for state in read_json(output/"paper"/"eigenstates.json"):
        label = state["label"]
        if label["figure"] == 17:
            full_energies.setdefault((label["setting"], label["u"]), []).extend(state["energies"])
    matching = []
    for setting in sorted({r["setting"] for r in spectrum}):
        for u in sorted({float(r["u"]) for r in spectrum if r["setting"] == setting}):
            published = distinct([float(r["value"]) for r in reference if r["method"] == "NRG"
                                  and r["observable"] == "excitation" and abs(float(r["u"])-u) < 3e-5], 2e-5)
            if not len(published):
                continue
            energies = full_energies[setting, u]
            computed = distinct([e-min(energies) for e in energies if e-min(energies) > 2e-5], 2e-5)
            if min(published) < 2e-5:
                computed = np.r_[0., computed]
            i, j = linear_sum_assignment(abs(published[:, None]-computed))
            in_window = computed[computed <= 1.02]
            iw, jw = linear_sum_assignment(abs(published[:, None]-in_window))
            matching.append(dict(setting=setting, u=u, reference_levels=len(published),
                                 computed_levels=len(computed), matched_levels=len(i),
                                 unmatched_levels=len(published)-len(i),
                                 computed_plot_levels=len(in_window), missing_reference_levels=len(published)-len(iw),
                                 maximum_matched_error=float(max(abs(published[i]-computed[j]))),
                                 maximum_plot_window_error=float(max(abs(published[iw]-in_window[jw])))))
    with (output/"paper"/"spectrum_matching.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(matching[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(matching)
    result["spectrum_reference"] = {setting: dict(
        maximum_matched_error=max(r["maximum_matched_error"] for r in matching if r["setting"] == setting),
        maximum_plot_window_error=max(r["maximum_plot_window_error"] for r in matching if r["setting"] == setting),
        unmatched_levels=sum(r["unmatched_levels"] for r in matching if r["setting"] == setting),
        missing_reference_levels=sum(r["missing_reference_levels"] for r in matching if r["setting"] == setting))
        for setting in sorted({r["setting"] for r in matching})}

    result["published_che_curves"] = {}
    scalar_rows = read_csv(output/"paper"/"shared.csv")
    for setting in ("ChE-W2", "ChE-W4", "ChE-W6", "ChE-W8"):
        errors = {name: [] for name in ("excitation", "pairing", "spin_correlation")}
        for row in (r for r in scalar_rows if r["setting"] == setting):
            u = float(row["u"])
            for observable, jump in (("pairing", .04), ("spin_correlation", .2)):
                selected = [r for r in reference if r["method"] == setting and r["observable"] == observable]
                value = smooth_interpolate(u, selected, "u", "value", jump)
                if value is not None:
                    errors[observable].append(abs(float(row[observable])-value))
            levels = distinct([float(r["excitation"]) for r in spectrum if r["setting"] == setting
                               and float(r["u"]) == u and float(r["excitation"]) > 2e-5], 2e-5)
            for branch, level in enumerate(levels[:5]):
                selected = [r for r in reference if r["method"] == setting and r["observable"] == "excitation"
                            and int(r["branch"]) == branch]
                value = smooth_interpolate(u, selected, "u", "value", .2)
                if value is not None:
                    errors["excitation"].append(abs(level-value))
        result["published_che_curves"][setting] = {k: dict(points=len(v), maximum_error=max(v))
                                                  for k, v in errors.items() if v}

    result["bath_error_budget"] = []
    manifest = read_json(output/"convergence"/"manifest.json")
    omega = np.r_[0., np.geomspace(.001, 10., 301)]
    for record in manifest["baths"].values():
        bath = DiscreteBath.from_record(record)
        wide = bath.metadata.get("wide_band", False)
        exact = 1/np.hypot(omega, bath.delta) if wide else hybridization_g(omega, bath.delta, bath.bandwidth)
        result["bath_error_budget"].append(dict(kind=bath.metadata["kind"], levels=bath.levels, wide_band=wide,
                                               maximum_kernel_error=float(max(abs(discrete_g(bath, omega)-exact)))))
    states = read_json(output/"convergence"/"eigenstates.json")
    result["native_solver_cost"] = []
    for figure, key in ((8, "current_settings"), (17, "shared_settings")):
        for settings in manifest["input"]["profiles"]["convergence"][key]:
            selected = [s for s in states if s["label"]["figure"] == figure and s["label"]["setting"] == settings["label"]]
            result["native_solver_cost"].append(dict(figure=figure, setting=settings["label"], sector_solves=len(selected),
                                                     seconds=sum(s["timings_seconds"]["total"] for s in selected)))
    result["runtimes_seconds"] = {p: read_json(output/p/"manifest.json")["elapsed_seconds"] for p in ("quick", "paper", "convergence")}
    if (output/"dmrg"/"summary.json").exists():
        result["independent_dmrg"] = read_json(output/"dmrg"/"summary.json")
    write_json(output/"validation.json", result)
    return result


if __name__ == "__main__":
    import json
    print(json.dumps(analyze(), indent=2))
