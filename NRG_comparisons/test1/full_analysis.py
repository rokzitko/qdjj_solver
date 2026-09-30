"""Offline full-test1 policy, schema v2; empirical evidence, never certification.

Inventory is the frozen plan, not a list of successful runs. Record indices in
the report always address the caller's original payloads. Only save_states and
max_dimension are computational metadata; all other settings define science
groups. No residual is converted into an observable error bound.
"""

from collections import defaultdict
from copy import deepcopy
import json
import math

import numpy as np

from NRG_comparisons.test1.analysis import OBSERVABLES, analyze


LAMBDAS = (1.8, 2.0, 2.5, 3.0, 4.0)
NRG_AXES = ("keep", "keepenergy", "keepmin", "nmax")


def _key(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def _science(task):
    return deepcopy({k: v for k, v in task.items() if k not in ("save_states", "max_dimension")})


def _without(task, axis):
    task = deepcopy(task)
    container = task["numerical"] if task["backend"] == "nrg" else task["bath"] if axis == "levels" else task
    container.pop(axis, None)
    if axis == "nmax":
        for k in ("nmax", "levels"):
            task.get("bath", {}).pop(k, None)
    return task


def _order(setting, axis):
    return math.inf if setting is None or (axis == "keepenergy" and setting <= 0) else setting


def _safe(value):
    if isinstance(value, dict):
        return {k: _safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


def assess(records, targets, inventory):
    """Return a JSON-safe v2 report without modifying inputs or launching solvers.

    A finite-bath endpoint needs two successive refinements, an unrestricted
    calculation, or an independently qualified same-bath counterpart. Test1's
    spinful sites have local dimension four: chi >= 4**floor((levels+1)/2)
    saturates the maximum MPS rank (the worker uses group_size=2).

    Bath errors include uncertainty at BOTH ends of each difference. NRG fits
    propagate point uncertainties using absolute pseudoinverse weights and add
    the full model/range envelope. Acceptance is conditional empirical stability,
    not a rigorous continuum bound or a fit calibrated to QP/DMRG.
    """
    if records and (not isinstance(records[0].get("physical"), dict) or any(
        r.get("physical") != records[0]["physical"] for r in records
    )):
        raise ValueError("records must have one identical physical dict")
    for backend in ("qp", "dmrg", "nrg"):
        if any(isinstance(targets[backend][k], bool) or not math.isfinite(targets[backend][k])
               or targets[backend][k] < 0 for k in ("absolute", "relative")):
            raise ValueError("targets must be finite and nonnegative")
    active, seen = set(), set()
    for entry in inventory:
        task, status, index = entry["task"], entry["status"], entry["record"]
        if status not in ("completed", "failed", "unstarted") or _key(task) in seen:
            raise ValueError("invalid or duplicate inventory entry")
        seen.add(_key(task))
        if status == "completed":
            if type(index) is not int or not 0 <= index < len(records) or records[index]["task"] != task:
                raise ValueError("completed inventory record must match the exact task")
            active.add(index)
        elif index is not None:
            raise ValueError("unfinished inventory entry must have a null record")
    for task in [r["task"] for r in records] + [e["task"] for e in inventory]:
        for source in (task.get("bath", {}), task.get("numerical", {})):
            if "lambda" in source and (not math.isfinite(source["lambda"]) or source["lambda"] < 1.8):
                raise ValueError("Lambda must be finite and >= 1.8")
    normalized = deepcopy(records)
    for i, row in enumerate(normalized):
        row["task"] = _science(row["task"])
        if inventory and i not in active:
            row["error"] = "outside frozen inventory"
    if not inventory:
        active = set(range(len(records)))
    report = analyze(normalized, targets)
    report.update(schema_version=2, policy="fulltest1-v1", final_reference={}, nrg_extrapolation={},
                  consistency_checks=[], status_changed_count=0, missing_gates=[], finite_bath_endpoints=[],
                  bath_series=[], finite_lambda={}, prerequisites_passed=True)
    if not inventory:
        report["missing_gates"].append("frozen_inventory")
    excluded = {e["record"] for e in report["excluded"]}
    rows = {i: r for i, r in enumerate(normalized) if i in active and i not in excluded}
    tasks = [_science(e["task"]) for e in inventory] if inventory else [r["task"] for r in rows.values()]
    identity = {c["record"]: c["status"] == "passed" for c in report["dmrg_identity_checks"]}

    def value(i, o):
        x = rows[i]["values"].get(o)
        return float(x) if isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) else None

    def valid(i, o):
        return (rows[i].get("finite_converged") is True and identity.get(i, True) and value(i, o) is not None)

    def threshold(backend, x):
        t = targets["qp" if backend == "ed" else backend]
        return t["absolute"] + t["relative"] * abs(x)

    def control(t):
        return (t["backend"] == "ed" or t.get("clean", False) or t.get("numerical", {}).get("untruncated", False)
                or t.get("bath", {}).get("kind") == "wilson" or t.get("bath", {}).get("levels", 9) <= 8)

    report["controls"] = [i for i, r in rows.items() if control(r["task"])]
    # Retain legacy diagnostic ladders, but require the last three PLANNED
    # settings, not three surviving successes straddling an unexecuted task.
    planned_axes = defaultdict(set)
    for t in tasks:
        axes = NRG_AXES if t["backend"] == "nrg" else ("cutoff",) if t["backend"] == "qp" else ("chi",) if t["backend"] == "dmrg" else ()
        for axis in axes:
            source = t["numerical"] if t["backend"] == "nrg" else t
            planned_axes[axis, _key(_without(t, axis))].add(_order(source[axis], axis))
    measured, evidence = {}, {}
    for ladder in report["refinements"]:
        axis = ladder["axis"]
        if axis == "levels":
            continue
        ids = ladder["records"]
        for end, i in enumerate(ids):
            window = ids[max(0, end - 2):end + 1]
            settings = [_order(s, axis) for s in ladder["settings"][max(0, end - 2):end + 1]]
            expected = sorted(s for s in planned_axes[axis, _key(_without(rows[i]["task"], axis))] if s <= settings[-1])[-3:]
            for o in OBSERVABLES:
                if len(window) > 1 and all(valid(j, o) for j in window):
                    error = max(abs(value(a, o) - value(b, o)) for a, b in zip(window, window[1:], strict=False))
                    measured[i, axis, o] = error
                    if (len(window) == 3 and settings == expected and len(set(settings)) == 3
                            and not ladder["duplicate_settings"] and error <= threshold(rows[i]["task"]["backend"], value(i, o))):
                        evidence[i, axis, o] = error

    groups = defaultdict(lambda: {"tasks": [], "records": []})
    for t in tasks:
        if t["backend"] in ("qp", "dmrg"):
            axis = "cutoff" if t["backend"] == "qp" else "chi"
            groups[_key(_without(t, axis))]["tasks"].append(t)
    for i, r in rows.items():
        t = r["task"]
        if t["backend"] in ("qp", "dmrg"):
            groups[_key(_without(t, "cutoff" if t["backend"] == "qp" else "chi"))]["records"].append(i)

    def exact(i):
        t = rows[i]["task"]
        if t["backend"] in ("ed", "qp"):
            return t["backend"] == "ed" or t["cutoff"] is None
        if t["backend"] == "dmrg":
            sites = t["bath"]["levels"] + 1
            rank = 4 ** (sites // 2) if t.get("group_size", 2) == 2 else 2 ** sites
            return t["chi"] >= rank
        return t.get("numerical", {}).get("untruncated", False)

    independent, qualified = {}, {}
    for i, r in rows.items():
        backend = r["task"]["backend"]
        if backend not in ("ed", "qp", "dmrg"):
            continue
        axis = "cutoff" if backend == "qp" else "chi"
        for o in OBSERVABLES:
            if valid(i, o):
                error = 0.0 if exact(i) else evidence.get((i, axis, o))
                if error is not None:
                    independent[i, o] = (error, "unrestricted finite space" if exact(i) else "two solver refinements", None)
    qualified.update(independent)

    def terminal(i):
        t = rows[i]["task"]
        if t["backend"] == "ed":
            return True
        axis = "cutoff" if t["backend"] == "qp" else "chi"
        return _order(t[axis], axis) >= max(planned_axes[axis, _key(_without(t, axis))], default=-1)

    for i, r in rows.items():
        backend = r["task"]["backend"]
        if backend not in ("qp", "dmrg"):
            continue
        for o in OBSERVABLES:
            if not valid(i, o) or (i, o) in qualified:
                continue
            counterparts = []
            for j, s in rows.items():
                other = s["task"]["backend"]
                # DMRG may never bootstrap bond convergence from truncated QP.
                if (other == backend or other not in ("qp", "dmrg", "ed") or (j, o) not in independent
                        or not terminal(j) or s["task"].get("bath") != r["task"].get("bath")
                        or (backend == "dmrg" and not exact(j))):
                    continue
                error = abs(value(i, o) - value(j, o)) + independent[j, o][0]
                if error <= threshold(backend, value(i, o)):
                    counterparts.append((error, "independent same-bath counterpart", j))
            if counterparts:
                qualified[i, o] = min(counterparts)

    endpoints = {}
    for gkey, group in groups.items():
        t = json.loads(gkey)
        axis = "cutoff" if t["backend"] == "qp" else "chi"
        rank = lambda i, axis=axis: _order(rows[i]["task"][axis], axis)
        planned = max((_order(s[axis], axis) for s in group["tasks"]), default=-1)
        for o in OBSERVABLES:
            available = [i for i in group["records"] if value(i, o) is not None]
            good = [i for i in available if (i, o) in qualified]
            i = max(good or available, key=rank) if available else None
            error, method, counterpart = qualified.get((i, o), (measured.get((i, axis, o)), "missing solver refinement evidence", None))
            endpoint = dict(task=t, observable=o, record=i, value=value(i, o) if i is not None else None,
                            empirical_uncertainty=error, qualified=(i, o) in qualified,
                            planned_endpoint_reached=i is not None and rank(i) >= planned,
                            reason=method, counterpart=counterpart)
            endpoints[gkey, o] = endpoint
            report["finite_bath_endpoints"].append(endpoint)

    inconsistent = set()
    def consistency(kind, ids, o, delta, limit):
        passed = math.isfinite(delta) and delta <= limit
        report["consistency_checks"].append(dict(kind=kind, records=ids, observable=o,
                                                absolute_delta=delta, threshold=limit, passed=passed))
        if not passed:
            inconsistent.add(o)

    selected_ids = {e["record"] for e in endpoints.values() if e["qualified"] and e["planned_endpoint_reached"]}
    selected_ids.update(i for i, r in rows.items() if r["task"]["backend"] == "ed")
    failed_prerequisite = False
    for check in report["same_bath_checks"]:
        a, b = check["records"]
        strict = all(exact(i) or i in selected_ids for i in (a, b)) and all(control(rows[i]["task"]) for i in (a, b))
        check["prerequisite"] = strict
        for o, comparison in check["observables"].items():
            if strict:
                passed = valid(a, o) and valid(b, o) and comparison["absolute_delta"] is not None and comparison["absolute_delta"] <= 1e-8
                comparison.update(threshold=1e-8, passed=passed)
                failed_prerequisite |= not passed
            elif all(i in selected_ids and (i, o) in qualified for i in (a, b)):
                limit = sum(max(qualified[i, o][0], threshold(rows[i]["task"]["backend"], value(i, o))) for i in (a, b))
                consistency("same_bath", [a, b], o, abs(value(a, o) - value(b, o)), limit)

    missing_controls = []
    for gkey in groups:
        if control(json.loads(gkey)) and any(not endpoints[gkey, o]["qualified"] or not endpoints[gkey, o]["planned_endpoint_reached"] for o in OBSERVABLES):
            missing_controls.append(gkey)
    for entry in inventory:
        t, i = entry["task"], entry["record"]
        if control(t) and t["backend"] in ("ed", "nrg"):
            if (entry["status"] != "completed" or records[i].get("finite_converged") is not True
                    or records[i].get("error") or records[i].get("success") is False
                    or records[i].get("status") in ("failed", "error", "timeout")
                    or (not t.get("clean") and (i not in rows or not all(valid(i, o) for o in OBSERVABLES)))):
                missing_controls.append(_key(t))
            elif t.get("clean"):
                error = records[i].get("values", {}).get("clean.vacuum_error")
                failed_prerequisite |= isinstance(error, bool) or not isinstance(error, (int, float)) or not math.isfinite(error) or abs(error) > 1e-8
    if missing_controls:
        report["missing_gates"].append("planned_controls")

    for backend in ("qp", "dmrg"):
        series = defaultdict(list)
        for gkey in groups:
            t = json.loads(gkey)
            if t["backend"] == backend and not control(t):
                series[_key(_without(t, "levels"))].append(gkey)
        families = defaultdict(list)
        for skey, keys in series.items():
            keys.sort(key=lambda k: json.loads(k)["bath"]["levels"])
            families[json.loads(skey)["bath"]["kind"]].append(skey)
        # Select by the frozen plan's coverage, never by agreement with a value.
        primary = {kind: max(keys, key=lambda s: (len(series[s]), json.loads(series[s][-1])["bath"]["levels"], s))
                   for kind, keys in families.items()}
        report["final_reference"][backend] = {}
        for o in OBSERVABLES:
            missing, selected, bath_errors, changes = [], [], [], []
            for kind, skey in primary.items():
                points = [endpoints[k, o] for k in series[skey]]
                available = [p for p in points if p["value"] is not None]
                if available:
                    selected.append(available[-1])
                last = points[-3:]
                stable = len(last) == 3 and all(p["qualified"] and p["planned_endpoint_reached"] for p in last)
                if not stable:
                    missing.append(f"{kind}: three qualified planned bath endpoints")
                differences = []
                for a, b in zip(last, last[1:], strict=False):
                    if a["empirical_uncertainty"] is not None and b["empirical_uncertainty"] is not None:
                        delta = abs(a["value"] - b["value"])
                        differences.append(delta + a["empirical_uncertainty"] + b["empirical_uncertainty"])
                        changes.append(dict(records=[a["record"], b["record"]], absolute_delta=delta,
                                            solver_uncertainty=a["empirical_uncertainty"] + b["empirical_uncertainty"]))
                if differences:
                    bath_errors.append(max(differences))
                report["bath_series"].append(dict(backend=backend, observable=o, task=json.loads(skey),
                                                  records=[p["record"] for p in last], qualified=stable))
            reference = max(selected, key=lambda p: (p["task"]["bath"]["levels"], _key(p["task"]))) if selected else None
            solver = max(p["empirical_uncertainty"] for p in selected) if selected and all(p["empirical_uncertainty"] is not None for p in selected) else None
            independent_error = None
            if set(primary) >= {"cosh", "surrogate"} and len(selected) == len(primary) and solver is not None:
                independent_error = max(p["value"] + p["empirical_uncertainty"] for p in selected) - min(p["value"] - p["empirical_uncertainty"] for p in selected)
            else:
                missing.append("independent cosh and surrogate families")
            window_error = 0.0
            for kind, skeys in families.items():
                if len(skeys) < 2:
                    continue
                # Different construction settings are sensitivity variants, not
                # independent families and never consecutive bath sizes.
                baseline = {endpoints[k, o]["task"]["bath"]["levels"]: endpoints[k, o] for k in series[primary[kind]]}
                variants = defaultdict(list)
                for skey in skeys:
                    if skey == primary[kind]:
                        continue
                    for k in series[skey]:
                        p = endpoints[k, o]
                        q = baseline.get(p["task"]["bath"]["levels"])
                        if q is None or not all(e["qualified"] and e["planned_endpoint_reached"] for e in (p, q)):
                            missing.append("planned construction/window sensitivity")
                            window_error = None
                        elif window_error is not None:
                            variants[p["task"]["bath"]["levels"]].extend((p, q))
                if window_error is not None:
                    for points in variants.values():
                        window_error = max(window_error, max(p["value"] + p["empirical_uncertainty"] for p in points)
                                           - min(p["value"] - p["empirical_uncertainty"] for p in points))
            estimates = dict(solver=solver, bath=max(bath_errors, default=None), independent_bath=independent_error, window=window_error)
            error = solver + max(estimates[k] for k in ("bath", "independent_bath", "window")) if all(v is not None for v in estimates.values()) else None
            x = reference["value"] if reference else None
            limit = threshold(backend, x) if x is not None else None
            if error is None or not math.isfinite(error) or error > limit:
                missing.append("combined empirical error exceeds target or is undetermined")
            if o == "singlet.moment":
                missing.append("symmetry zero excluded")
            missing += report["missing_gates"]
            report["final_reference"][backend][o] = dict(value=x, reference_record=reference["record"] if reference else None,
                records=[p["record"] for p in selected], estimates=estimates, bath_changes=changes,
                empirical_uncertainty=error, threshold=limit, missing_gates=sorted(set(missing)),
                status="unresolved" if missing else "empirically_converged",
                reason="; ".join(sorted(set(missing))) if missing else "qualified solver/bath endpoints and independent constructions agree")

    for o in OBSERVABLES:
        refs = [report["final_reference"][b][o] for b in ("qp", "dmrg")]
        if all(r["status"] == "empirically_converged" for r in refs):
            consistency("final_reference", [r["reference_record"] for r in refs], o, abs(refs[0]["value"] - refs[1]["value"]),
                        sum(max(r["empirical_uncertainty"], r["threshold"]) for r in refs))

    def grid_key(t):
        return _key(_without(t, "z"))

    grids = {grid_key(g["task"]): g for g in report["z_grids"]}
    planned_grids = {_key(_without(t, "z")): _without(t, "z") for t in tasks
                     if t["backend"] == "nrg" and t["units"] == "gap" and not control(t)}
    anchors, points, reached = {}, {}, {}
    for lam in LAMBDAS:
        candidates = [t for t in planned_grids.values() if t["numerical"]["lambda"] == lam]
        if not candidates:
            continue
        t = max(candidates, key=lambda t: tuple(_order(t["numerical"][a], a) for a in NRG_AXES))
        anchors[lam] = t
        g = grids.get(grid_key(t))
        reached[lam] = g is not None and g["grids"]["8"]["complete"] and all(
            _order(t["numerical"][a], a) == max(_order(c["numerical"][a], a) for c in candidates) for a in NRG_AXES)
        if g is None or not g["grids"]["8"]["complete"]:
            available = [grids[grid_key(c)] for c in candidates if grid_key(c) in grids and grids[grid_key(c)]["grids"]["8"]["complete"]]
            g = max(available, key=lambda g: tuple(_order(g["task"]["numerical"][a], a) for a in NRG_AXES)) if available else None
        points[lam] = g
    signatures = []
    for t in anchors.values():
        for axis in ("lambda", "nmax"):
            t = _without(t, axis)
        signatures.append(_key(t))
    for g in points.values():
        if g is not None:
            signatures.append(_key(_without(_without(g["task"], "lambda"), "nmax")))
    series_signatures = set()
    for t in planned_grids.values():
        for axis in (*NRG_AXES, "lambda"):
            t = _without(t, axis)
        series_signatures.add(_key(t))
    same_anchors = len(set(signatures)) == 1 and len(series_signatures) == 1
    required_units, passed_units = set(), set()
    for t in anchors.values():
        t = deepcopy(t)
        t["numerical"]["z"] = 1.0
        required_units.add(_key(t))
    for check in report["unit_checks"]:
        pair = [rows[i] for i in check["records"]]
        gap = next(r for r in pair if r["task"]["units"] == "gap")
        relevant = _key(gap["task"]) in required_units or control(gap["task"])
        check["prerequisite"] = relevant
        passed = True
        for o, comparison in check["observables"].items():
            ok = comparison["absolute_delta"] is not None and comparison["absolute_delta"] <= 1e-8 and all(valid(i, o) for i in check["records"])
            comparison.update(threshold=1e-8, passed=ok)
            passed &= ok
        if relevant:
            failed_prerequisite |= not passed
            if passed:
                passed_units.add(_key(gap["task"]))
    if required_units - passed_units:
        report["missing_gates"].append("selected NRG bandwidth-unit checks")

    for g in report["z_grids"]:
        ids = g["grids"]["8"]["records"]
        g["all_axes_passed"] = {o: bool(g["grids"]["8"]["complete"] and g["comparison"][o]["passed"]
            and all((i, a, o) in evidence for i in ids for a in NRG_AXES)) for o in OBSERVABLES}
    report["final_reference"]["nrg"] = {}
    for o in OBSERVABLES:
        finite, missing, missing_axes = [], list(report["missing_gates"]), set()
        for lam, g in sorted(points.items()):
            if not reached[lam]:
                missing.append(f"Lambda={lam}: planned solver anchor not reached")
            if g is None or g["grids"]["8"]["values"][o] is None:
                missing.append(f"Lambda={lam}: complete planned z8 anchor")
                continue
            ids = g["grids"]["8"]["records"]
            estimates = {a: sum(measured[i, a, o] for i in ids) / 8 if all((i, a, o) in measured for i in ids) else None for a in NRG_AXES}
            estimates["z4_vs_z8"] = g["comparison"][o]["absolute_delta"]
            error = sum(estimates.values()) if all(v is not None for v in estimates.values()) else None
            absent = [a for a in NRG_AXES if any((i, a, o) not in evidence for i in ids)]
            if not g["comparison"][o]["passed"]:
                absent.append("z4_vs_z8")
            missing_axes.update(absent)
            missing.extend(f"Lambda={lam}: {a}" for a in absent)
            finite.append(dict(lambda_value=lam, value=g["grids"]["8"]["values"][o], records=ids,
                               task=g["task"], planned_endpoint_reached=reached[lam], estimates=estimates,
                               empirical_uncertainty=error, all_axes_passed=not absent))
        report["finite_lambda"][o] = finite
        if len(finite) != 5:
            missing.append("complete five-Lambda grid: 4, 3, 2.5, 2, 1.8")
        if not same_anchors:
            missing.append("identical numerical anchors across Lambda (except nmax)")
        fits, trials, leaveouts, high_omit = {}, [], [], []
        if len(finite) >= 3 and same_anchors and o != "singlet.moment":
            lam = np.array([p["lambda_value"] for p in finite])
            y = np.array([p["value"] for p in finite])
            errors = [p["empirical_uncertainty"] for p in finite]
            for form, x in (("log_lambda", np.log(lam)), ("lambda_minus_one", lam - 1)):
                for degree in (1, 2):
                    name = f"{form}_degree_{degree}"
                    masks = [("full", np.ones(len(lam), dtype=bool))]
                    if len(finite) == 5:
                        masks += [("leave_one_out", np.arange(5) != i) for i in range(5)]
                        masks += [("high_lambda_omit", lam < 4)]
                    for kind, mask in masks:
                        weights = np.linalg.pinv(np.vander(x[mask], degree + 1, increasing=True))[0]
                        prediction = float(weights @ y[mask])
                        uncertainty = float(np.abs(weights) @ np.array(errors)[mask]) if all(e is not None for e in errors) else None
                        trial = dict(form=name, value=prediction, lambdas=lam[mask].tolist(), weights=weights.tolist(),
                                     propagated_uncertainty=uncertainty)
                        trials.append(trial)
                        if kind == "full":
                            fits[name] = prediction
                        elif kind == "leave_one_out":
                            leaveouts.append(trial | {"omitted_lambda": float(lam[~mask][0])})
                        else:
                            high_omit.append(trial)
        x = fits.get("log_lambda_degree_1")
        spread = max(t["value"] for t in trials) - min(t["value"] for t in trials) if trials else None
        propagated = max(t["propagated_uncertainty"] for t in trials) if trials and all(t["propagated_uncertainty"] is not None for t in trials) else None
        error = propagated + spread if propagated is not None else None
        limit = threshold("nrg", x) if x is not None else None
        if error is None or not math.isfinite(error) or error > limit:
            missing.append("propagated error plus model/range envelope exceeds target or is undetermined")
        if o == "singlet.moment":
            missing.append("symmetry zero excluded")
        accepted = not missing
        reason = "; ".join(sorted(set(missing))) if missing else "empirical fit stability and propagated error meet target; not a rigorous bound"
        report["nrg_extrapolation"][o] = dict(value=x, status="empirically_converged" if accepted else "unresolved",
            accepted=accepted, reason=reason, lambdas=[p["lambda_value"] for p in finite], fits=fits,
            fit_trials=trials, leave_one_out=leaveouts, high_lambda_omit=high_omit, sensitivity=spread,
            propagated_uncertainty=propagated, empirical_uncertainty=error, threshold=limit,
            all_axes_passed=bool(finite) and all(p["all_axes_passed"] for p in finite), missing_gates=sorted(set(missing)),
            missing_axes=[a for a in (*NRG_AXES, "z4_vs_z8") if not finite or a in missing_axes],
            sensitivity_passed=spread is not None and math.isfinite(spread) and spread <= limit)
        report["final_reference"]["nrg"][o] = dict(value=x if x is not None else finite[0]["value"] if finite else None,
            candidate_kind="lambda_extrapolation" if x is not None else "finite_lambda" if finite else None, reference_record=None,
            records=[i for p in finite for i in p["records"]], finite_lambda_candidate=finite[0] if finite else None,
            status="empirically_converged" if accepted else "unresolved", reason=reason,
            estimates=dict(propagated=propagated, model_range_envelope=spread), empirical_uncertainty=error,
            threshold=limit, missing_gates=sorted(set(missing)))

    report["prerequisites_passed"] = not failed_prerequisite and not report["missing_gates"]
    for o in OBSERVABLES:
        if failed_prerequisite or o in inconsistent:
            reason = "failed control/unit prerequisite" if failed_prerequisite else "inconsistent qualified QP/DMRG evidence"
            for backend in ("qp", "dmrg", "nrg"):
                ref = report["final_reference"][backend][o]
                report["status_changed_count"] += ref["status"] == "empirically_converged"
                ref.update(status="unresolved", reason=reason)
                ref["missing_gates"].append(reason)
            report["nrg_extrapolation"][o].update(status="unresolved", accepted=False, reason=reason)
            report["nrg_extrapolation"][o]["missing_gates"].append(reason)
        for a, b in (("qp", "dmrg"), ("qp", "nrg"), ("dmrg", "nrg")):
            x, y = (report["final_reference"][backend][o]["value"] for backend in (a, b))
            report["cross_solver_discrepancy"][f"{a}_vs_{b}"][o] = abs(x - y) if x is not None and y is not None else None
    return _safe(report)
