"""Offline schema v1: input-record indices, empirical errors, exploratory NRG fits.
Refinements retain endpoint evidence; z_grids average complete grids only.
prerequisites_passed means no available check failed, not that every control
exists. status_changed_count counts gate-induced demotions. No solver launches.
"""

from collections import defaultdict
from copy import deepcopy
from itertools import combinations
import json
import math

import numpy as np


OBSERVABLES = ("signed_gap",) + tuple(
    f"{branch}.{name}" for branch in ("singlet", "doublet")
    for name in ("energy", "P0", "P1", "P2", "moment")
)


def analyze(records, targets):
    """Assess supplied successful payloads; missing evidence stays unresolved.

    Energies are already in gap units with the SAME finite BCS bath subtracted
    and centered impurity. Probabilities/moments are dimensionless. No energy
    shifts, residual-to-observable bounds, or highest-basis exactness are inferred.
    z grids are k/N, k=1..N. Only gap-unit decks enter NRG ladders and fits.
    DMRG identity needs two separated roots per branch, not just a small residual.
    """
    key = lambda obj: json.dumps(obj, sort_keys=True, allow_nan=False)

    def number(value):
        return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) else None

    if any("physical" in r for r in records) and any(
        not isinstance(r.get("physical"), dict) or r["physical"] != records[0].get("physical") for r in records
    ):
        raise ValueError("records must have one identical physical dict with consistent presence")
    report = {"schema_version": 1, "uncertainty_kind": "empirical_not_rigorous",
              "energy_units": "gap", "symmetry_only": ["singlet.moment"],
              "energy_reference": "centered impurity; same finite-chain BCS bath subtracted",
              "excluded": [], "controls": [], "same_bath_checks": [], "unit_checks": [],
              "prerequisites_passed": True, "status_changed_count": 0, "consistency_checks": [], "dmrg_identity_checks": [],
              "refinements": [], "z_grids": [], "final_reference": {},
              "nrg_extrapolation": {}, "cross_solver_discrepancy": {}}
    rows = []
    for index, record in enumerate(records):
        task = record.get("task", {})
        for source in (task.get("bath", {}), task.get("numerical", {})):
            if "lambda" in source and (not math.isfinite(source["lambda"]) or source["lambda"] < 1.8):
                raise ValueError("Lambda must be finite and >= 1.8")
        reason = ("failed or missing payload" if record.get("error") or record.get("success") is False
                  or record.get("status") in ("failed", "error", "timeout") or not record.get("values")
                  else "CLEAN is not an impurity result" if task.get("clean") else None)
        if reason:
            report["excluded"].append({"record": index, "reason": reason})
            continue
        control = (task["backend"] == "ed" or task.get("bath", {}).get("kind") == "wilson"
                   or task.get("numerical", {}).get("untruncated", False))
        if control:
            report["controls"].append(index)
        values = {o: number(record["values"].get(o)) for o in OBSERVABLES}
        identity = "passed"
        if task["backend"] == "dmrg":
            branches = [record.get("diagnostics", {}).get(b, {}) for b in ("singlet", "doublet")]
            roots = [b.get("lowest_roots_over_gap", []) for b in branches]
            identity = "not_checked" if any(len(es) < 2 or any(number(e) is None for e in es[:2]) for es in roots) else (
                "passed" if all(es[1] - es[0] > 1e-9 for es in roots)
                and not any(b.get("degenerate_minimum", False) for b in branches) else "failed")
            report["dmrg_identity_checks"].append({"record": index, "status": identity})
        rows.append({"id": index, "task": task, "values": values, "control": control,
                     "finite": record.get("finite_converged") is True, "identity": identity})

    def threshold(backend, reference):
        target = targets["qp" if backend == "ed" else backend]
        return target["absolute"] + target["relative"] * abs(reference)

    def compare(a, b):
        result = {}
        for o in OBSERVABLES:
            x, y = a["values"][o], b["values"][o]
            delta = number(x - y) if x is not None and y is not None else None
            tolerance = min(threshold(r["task"]["backend"], y) for r in (a, b)) if y is not None else None
            result[o] = {"delta": delta, "absolute_delta": abs(delta) if delta is not None else None,
                         "threshold": tolerance, "passed": bool(a["finite"] and b["finite"]
                         and delta is not None and abs(delta) <= tolerance)}
        return result

    for a, b in combinations(rows, 2):
        ta, tb = a["task"], b["task"]
        ua, ub = ({k: v for k, v in t.items() if k != "units"} for t in (ta, tb))
        unit_check = (ta["backend"] == tb["backend"] == "nrg" and ua == ub
                      and {ta["units"], tb["units"]} == {"gap", "bandwidth"})
        same_bath = (ta["backend"] != tb["backend"] and ta.get("bath") == tb.get("bath")
                     and bool(ta.get("bath")) and all(t["backend"] != "nrg" or
                     (t["units"] == "gap" and t["bath"]["kind"] == "wilson" and t["numerical"]["z"] == 1
                      and all(t["numerical"][k] == t["bath"].get(k) for k in ("lambda", "nmax")))
                     for t in (ta, tb)))
        if unit_check or same_bath:
            report["unit_checks" if unit_check else "same_bath_checks"].append(
                {"records": [a["id"], b["id"]], "observables": compare(a, b)})

    canonical = [r for r in rows if r["task"]["backend"] != "nrg" or r["task"]["units"] == "gap"]
    evidence, measured, windows, final_solver = {}, {}, {}, set()
    for backend, axes in (("qp", ("cutoff", "levels")), ("dmrg", ("chi", "levels")),
                          ("nrg", ("keep", "keepenergy", "keepmin", "nmax"))):
        for axis in axes:
            groups = defaultdict(list)
            for r in canonical:
                if r["task"]["backend"] != backend:
                    continue
                task = deepcopy(r["task"])
                container = task["bath"] if axis == "levels" else task["numerical"] if backend == "nrg" else task
                setting = container.pop(axis)
                if backend == "nrg" and axis == "nmax":
                    for k in ("nmax", "levels"):
                        task.get("bath", {}).pop(k, None)
                order = math.inf if setting is None or (axis == "keepenergy" and setting <= 0) else setting
                groups[key(task)].append((order, setting, r))
            for group in groups.values():
                group.sort(key=lambda entry: entry[0])
                ladder = [entry[2] for entry in group]
                changes = [{"records": [a["id"], b["id"]], "observables": compare(a, b)}
                           for a, b in zip(ladder, ladder[1:], strict=False)]
                unique = len({entry[0] for entry in group}) == len(group)
                if axis in ("cutoff", "chi"):
                    final_solver.add(ladder[-1]["id"])
                endpoints = []
                for end, r in enumerate(ladder):
                    window, adjacent = ladder[max(0, end-2):end+1], changes[max(0, end-2):end]
                    windows[r["id"], axis] = [w["id"] for w in window]
                    finite = all(w["finite"] and w["identity"] == "passed" for w in window)
                    distinct = all(sum(entry[0] == group[j][0] for entry in group) == 1 for j in range(max(0, end-2), end+1))
                    passed = {}
                    for o in OBSERVABLES:
                        deltas = [c["observables"][o]["absolute_delta"] for c in adjacent]
                        if deltas and None not in deltas and finite:
                            measured[r["id"], axis, o] = max(deltas)
                        passed[o] = len(adjacent) == 2 and distinct and finite and all(c["observables"][o]["passed"] for c in adjacent)
                        if passed[o]:
                            evidence[r["id"], axis, o] = max(deltas)
                    endpoints.append({"record": r["id"], "passed": passed})
                report["refinements"].append({"backend": backend, "axis": axis,
                    "records": [r["id"] for r in ladder], "settings": [entry[1] for entry in group],
                    "changes": changes, "passed": passed, "endpoints": endpoints, "duplicate_settings": not unique})

    for backend, axis in (("qp", "cutoff"), ("dmrg", "chi")):
        candidates = [r for r in canonical if r["task"]["backend"] == backend and not r["control"]]
        families = {}
        for r in candidates:
            setting = r["task"][axis]
            rank = (r["task"]["bath"]["levels"], math.inf if setting is None else setting)
            family = r["task"]["bath"]["kind"]
            if family not in families or rank > families[family][0]:
                families[family] = (rank, r)
        selected = [entry[1] for entry in families.values()]
        reference = max(families.values(), key=lambda entry: entry[0])[1] if selected else None
        report["final_reference"][backend] = {}
        for o in OBSERVABLES:
            value = reference["values"][o] if reference else None
            estimates = {a: max((measured.get((r["id"], a, o), math.inf) for r in selected), default=math.inf)
                         for a in (axis, "levels")}
            independent = (max(r["values"][o] for r in selected) - min(r["values"][o] for r in selected)
                           if len(selected) >= 2 and all(r["values"][o] is not None for r in selected) else math.inf)
            estimates["independent_bath"] = independent
            # Bath refinement requires a solver ladder at EACH participating bath.
            bath_rows = [windows[r["id"], "levels"] for r in selected]
            stable = (bool(selected) and all((r["id"], a, o) in evidence for r in selected for a in (axis, "levels"))
                      and all((i, axis, o) in evidence for ids in bath_rows for i in ids))
            uncertainty = estimates[axis] + max(estimates["levels"], independent)
            accepted = (o != "singlet.moment" and value is not None and stable
                        and uncertainty <= threshold(backend, value))
            report["final_reference"][backend][o] = {
                "value": value, "records": [r["id"] for r in selected],
                "reference_record": reference["id"] if reference else None,
                "status": "empirically_converged" if accepted else "unresolved",
                "reason": "two solver/bath refinements and independent families agree" if accepted else
                          "requires two finite-converged solver/bath refinements and independent families; symmetry zeros excluded",
                "threshold": threshold(backend, value) if value is not None else None,
                "estimates": {k: v if math.isfinite(v) else None for k, v in estimates.items()},
                "empirical_uncertainty": uncertainty if math.isfinite(uncertainty) else None}

    groups = defaultdict(list)
    for r in canonical:
        if r["task"]["backend"] == "nrg" and not r["control"]:
            task = deepcopy(r["task"])
            task["numerical"].pop("z")
            groups[key(task)].append(r)
    for task_key, group in groups.items():
        grids = {}
        for size in (4, 8):
            expected = {k / size for k in range(1, size + 1)}
            members = [r for r in group if r["task"]["numerical"]["z"] in expected]
            complete = (len(members) == size and {r["task"]["numerical"]["z"] for r in members} == expected
                        and all(r["finite"] for r in members))
            values = {o: sum(r["values"][o] / size for r in members) if complete
                      and all(r["values"][o] is not None for r in members) else None for o in OBSERVABLES}
            grids[str(size)] = {"complete": complete, "records": [r["id"] for r in members], "values": values}
        task = json.loads(task_key)
        comparisons = compare({"task": task, "finite": grids["4"]["complete"], "values": grids["4"]["values"]},
                              {"task": task, "finite": grids["8"]["complete"], "values": grids["8"]["values"]})
        axes_passed = {o: comparisons[o]["passed"] and all((i, a, o) in evidence
                       for i in grids["8"]["records"] for a in ("keep", "keepenergy", "nmax")) for o in OBSERVABLES}
        report["z_grids"].append({"task": task, "grids": grids, "comparison": comparisons, "all_axes_passed": axes_passed})

    # Highest settings per Lambda, not whichever incomplete/coarse grid looks best.
    points, series = {}, set()
    for grid in report["z_grids"]:
        n = grid["task"]["numerical"]
        fixed = deepcopy(grid["task"])
        for k in ("lambda", "keep", "keepenergy", "keepmin", "nmax"):
            fixed["numerical"].pop(k)
        for k in ("nmax", "levels"):
            fixed.get("bath", {}).pop(k, None)
        series.add(key(fixed))
        rank = (n["nmax"], n["keep"], math.inf if n["keepenergy"] <= 0 else n["keepenergy"], n["keepmin"])
        if n["lambda"] not in points or rank > points[n["lambda"]][0]:
            points[n["lambda"]] = (rank, grid)
    report["final_reference"]["nrg"] = {}
    for o in OBSERVABLES:
        usable = [(lam, g[1]) for lam, g in sorted(points.items()) if g[1]["grids"]["8"]["values"][o] is not None]
        fits, leaveouts = {}, []
        if len(usable) >= 3 and len(series) == 1 and o != "singlet.moment":
            lam = np.array([p[0] for p in usable])
            y = np.array([p[1]["grids"]["8"]["values"][o] for p in usable])
            for form, x in (("log_lambda", np.log(lam)), ("lambda_minus_one", lam - 1)):
                for degree in (1, 2):
                    name = f"{form}_degree_{degree}"
                    fits[name] = number(float(np.polynomial.polynomial.polyfit(x, y, degree)[0]))
                    if len(usable) >= 5:
                        for i in range(len(usable)):
                            value = number(float(np.polynomial.polynomial.polyfit(np.delete(x, i), np.delete(y, i), degree)[0]))
                            leaveouts.append({"form": name, "omitted_lambda": float(lam[i]), "value": value})
        predictions = [*fits.values(), *(fit["value"] for fit in leaveouts)]
        sensitivity = number(max(predictions) - min(predictions)) if predictions and None not in predictions else None
        all_axes = bool(usable) and all(g["all_axes_passed"][o] for _, g in usable)
        missing = [a for a in ("keep", "keepenergy", "nmax") if not usable or any(
            (i, a, o) not in evidence for _, g in usable for i in g["grids"]["8"]["records"])]
        if not usable or any(not g["comparison"][o]["passed"] for _, g in usable):
            missing.append("z4_vs_z8")
        report["nrg_extrapolation"][o] = {
            "status": "exploratory_unvalidated" if fits else "not_attempted", "accepted": False,
            "reason": "far-from-one empirical fits have no universal justification; validation deferred" if fits else
                      "requires >=3 complete finite-converged gap-unit z8 grids of one series; symmetry zeros excluded",
            "lambdas": [p[0] for p in usable], "fits": fits, "leave_one_out": leaveouts,
            "sensitivity": sensitivity, "all_axes_passed": all_axes, "missing_axes": missing,
            "sensitivity_passed": bool(len(usable) >= 5 and sensitivity is not None
                                       and all(sensitivity <= threshold("nrg", v) for v in fits.values()))}
        value = usable[0][1]["grids"]["8"]["values"][o] if usable else None
        ids = usable[0][1]["grids"]["8"]["records"] if usable else []
        estimates = {a: max((measured.get((i, a, o), math.inf) for i in ids), default=math.inf)
                     for a in ("keep", "keepenergy", "nmax")}
        estimates["z4_vs_z8"] = usable[0][1]["comparison"][o]["absolute_delta"] if usable else None
        estimates["lambda_fit_sensitivity"] = sensitivity
        report["final_reference"]["nrg"][o] = {"value": value, "status": "unresolved",
            "reason": "finite-Lambda z8 candidate only; continuum extrapolation unvalidated",
            "records": ids, "reference_record": None,
            "estimates": {k: v if v is not None and math.isfinite(v) else None for k, v in estimates.items()},
            "threshold": threshold("nrg", value) if value is not None else None,
            "empirical_uncertainty": None}
    for a, b in combinations(("qp", "dmrg", "nrg"), 2):
        report["cross_solver_discrepancy"][f"{a}_vs_{b}"] = {
            o: number(abs(x - y)) if (x := report["final_reference"][a][o]["value"]) is not None
            and (y := report["final_reference"][b][o]["value"]) is not None else None for o in OBSERVABLES}
    by_id, inconsistent = {r["id"]: r for r in rows}, set()
    def consistency(kind, ids, o, delta, limit):
        passed = delta is not None and delta <= limit
        report["consistency_checks"].append({"kind": kind, "records": ids, "observable": o,
                                            "absolute_delta": delta, "threshold": limit, "passed": passed})
        if not passed:
            inconsistent.add(o)

    for check in report["unit_checks"] + report["same_bath_checks"]:
        pair = [by_id[i] for i in check["records"]]
        check["prerequisite"] = check in report["unit_checks"] or all(
            r["task"]["backend"] == "ed" or (r["task"]["backend"] == "qp" and r["task"]["cutoff"] is None)
            or (r["task"]["backend"] == "nrg" and r["task"]["numerical"].get("untruncated", False))
            or (r["task"]["backend"] == "dmrg" and r["id"] in final_solver) for r in pair)
        for o, comparison in check["observables"].items():
            if check["prerequisite"]:
                limit = min(1e-8, comparison["threshold"]) if comparison["threshold"] is not None else None
                comparison.update(threshold=limit, passed=bool(comparison["passed"] and comparison["absolute_delta"] <= limit))
                report["prerequisites_passed"] &= comparison["passed"]
            elif {r["task"]["backend"] for r in pair} == {"qp", "dmrg"} and all(
                r["id"] in final_solver and (r["id"], "cutoff" if r["task"]["backend"] == "qp" else "chi", o) in evidence for r in pair
            ):
                error = sum(evidence[r["id"], "cutoff" if r["task"]["backend"] == "qp" else "chi", o] for r in pair)
                limit = max(error, sum(threshold(r["task"]["backend"], r["values"][o]) for r in pair))
                consistency("same_bath", check["records"], o, comparison["absolute_delta"], limit)
    for o in OBSERVABLES:
        refs = [report["final_reference"][b][o] for b in ("qp", "dmrg")]
        if all(r["status"] == "empirically_converged" for r in refs):
            limit = max(sum(r["empirical_uncertainty"] for r in refs), sum(r["threshold"] for r in refs))
            consistency("final_reference", [r["reference_record"] for r in refs], o, report["cross_solver_discrepancy"]["qp_vs_dmrg"][o], limit)
        for r in refs:
            if not report["prerequisites_passed"] or o in inconsistent:
                report["status_changed_count"] += r["status"] == "empirically_converged"
                r.update(status="unresolved", reason="failed control/unit prerequisite" if not report["prerequisites_passed"] else "inconsistent QP/DMRG evidence")
    return report
