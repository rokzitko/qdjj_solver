"""Offline targeted-policy-v1: explicit NRG anchors, empirical errors, no solvers.

All record indices address the original, unmodified payloads. Auxiliary probes
never replace the reference, including a bracketing keep check above it. Values
are already in gap units, even for bandwidth-unit input decks. No QP comparison
enters fit selection or acceptance; singlet.moment is only a symmetry check.
"""

from copy import deepcopy
import math

import numpy as np

from NRG_comparisons.test1.analysis import OBSERVABLES


LAMBDAS = (1.8, 2.0, 2.5, 3.0, 4.0)
NRG_AXES = ("keep", "keepenergy", "keepmin", "nmax")


def _key(value):
    if isinstance(value, dict):
        return tuple(sorted((k, _key(v)) for k, v in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_key(v) for v in value)
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("physical and task/plan parameters must be finite")
    return value


def _number(value):
    return float(value) if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) else None


def _safe(value):
    if isinstance(value, dict):
        return {k: _safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


def assess(records, plan, targets, inventory=None):
    """Assess complete k/Nz grids (Nz=8,16,32,64) at the supplied anchors.

    ``targets['nrg']`` contains nonnegative absolute/relative tolerances. Optional
    inventory entries have task/status/record, with original indices for completed
    tasks and null otherwise. Only completed inventory entries supply evidence.
    Identical duplicate tasks invalidate membership; inconsistent duplicates are
    rejected. Runtime/provenance metadata is neither rewritten nor used to rank.

    Each point error is the mean of four per-twist ladder errors plus the LATEST
    nested-grid change. A requested minimum-twist tail check replaces the nmax
    term by max(grid nmax error, tail error), not a diluted single-twist correction.
    These are conservative empirical estimates, not rigorous continuum bounds.
    """
    _key(plan)
    target = targets["nrg"]
    if any(_number(target[k]) is None or target[k] < 0 for k in ("absolute", "relative")):
        raise ValueError("targets must be finite and nonnegative")
    anchors = sorted(deepcopy(plan["anchors"]), key=lambda a: a["lambda"])
    if len({a["lambda"] for a in anchors}) != len(anchors):
        raise ValueError("duplicate Lambda anchors")
    for a in anchors:
        if a["lambda"] not in LAMBDAS or type(a["nz"]) is not int or a["nz"] not in (8, 16, 32, 64):
            raise ValueError("anchors require a supported Lambda and Nz=8,16,32,64")
        if "keep_checks" in a:
            a["checks"]["keep"] = a["keep_checks"]
        for axis in NRG_AXES:
            settings = a["checks"][axis]
            if (any(_number(v) is None or v <= 0 for v in [a[axis], *settings])
                    or len(settings) < 3 or sorted(set(settings)) != settings
                    or a[axis] not in settings[-3:] or (axis != "keep" and a[axis] != settings[-1])
                    or (axis == "keep" and a[axis] == settings[-3])):
                raise ValueError("checks need three ordered settings ending at or bracketing the reference")
        if a.get("tail_checks") and a["tail_checks"] != [a["nmax"] + 4, a["nmax"] + 8]:
            raise ValueError("tail_checks must be N+4,N+8")

    def validate_task(task):
        _key(task)
        n = task["numerical"]
        if (task["backend"] != "nrg" or task["units"] not in ("gap", "bandwidth")
                or type(task["clean"]) is not bool or type(n["untruncated"]) is not bool
                or any(_number(n[k]) is None for k in (*NRG_AXES, "lambda", "z"))
                or n["lambda"] < 1.8 or not 0 < n["z"] <= 1):
            raise ValueError("invalid NRG task parameters")

    by_task = {}
    for i, row in enumerate(records):
        if not isinstance(row.get("physical"), dict) or row["physical"] != records[0].get("physical"):
            raise ValueError("records must have one identical physical dict")
        _key(row["physical"])
        validate_task(row["task"])
        ids = by_task.setdefault(_key(row["task"]), [])
        if ids and any(_safe(row.get(k)) != _safe(records[ids[0]].get(k))
                       for k in ("values", "finite_converged", "error", "success", "status")):
            raise ValueError("inconsistent duplicate task records")
        ids.append(i)
    active, seen = set(range(len(records))) if inventory is None else set(), set()
    counts = dict(total=len(records) if inventory is None else len(inventory), completed=0, failed=0, unstarted=0)
    for entry in inventory or []:
        task, status, i = entry["task"], entry["status"], entry["record"]
        validate_task(task)
        key = _key(task)
        if status not in ("completed", "failed", "unstarted") or key in seen:
            raise ValueError("invalid or duplicate inventory entry")
        seen.add(key)
        counts[status] += 1
        if status == "completed":
            if type(i) is not int or not 0 <= i < len(records) or records[i]["task"] != task:
                raise ValueError("completed inventory record must match the exact task")
            active.add(i)
        elif i is not None:
            raise ValueError("unfinished inventory entry must have a null record")
    if inventory is None:
        counts["completed"] = len(records)
    counts.update(records=len(records), active_records=len(active))

    def lookup(task):
        ids = by_task.get(_key(task), [])
        return ids[0] if len(ids) == 1 and ids[0] in active else None

    def value(i, observable):
        if i is None:
            return None
        r = records[i]
        if (r.get("finite_converged") is not True or r.get("error") or r.get("success") is False
                or r.get("status") in ("failed", "error", "timeout")):
            return None
        return _number((r.get("values") or {}).get(observable))

    def task(a, z, units="gap", **changes):
        return dict(backend="nrg", clean=False, units=units, numerical={
            "lambda": a["lambda"], "z": z, "untruncated": False,
            **{k: a[k] for k in NRG_AXES}, **changes})

    def threshold(x):
        return _number(target["absolute"] + target["relative"] * abs(x)) if x is not None else None

    report = dict(schema_version=1, policy="targeted-policy-v1", uncertainty_kind="empirical_not_rigorous",
                  energy_units="gap", symmetry_only=["singlet.moment"], inventory_counts=counts,
                  unit_checks=[], clean_checks=[], missing_gates=[], finite_lambda={}, nrg_extrapolation={})
    clean_tasks = {_key(r["task"]): r["task"] for r in records if r["task"]["clean"]}
    clean_tasks.update({_key(e["task"]): e["task"] for e in inventory or [] if e["task"]["clean"]})
    coverage = {(t["numerical"]["nmax"], t["units"]) for t in clean_tasks.values()}
    if not {(n, u) for n in (2, 4) for u in ("gap", "bandwidth")} <= coverage:
        report["missing_gates"].append("CLEAN controls at nmax=2,4 in both units")
    for t in clean_tasks.values():
        i = lookup(t)
        error = value(i, "clean.vacuum_error")
        passed = error is not None and abs(error) <= 1e-8
        report["clean_checks"].append(dict(task=t, record=i, vacuum_error=error, threshold=1e-8, passed=passed))
        if not passed:
            report["missing_gates"].append("failed or missing CLEAN vacuum control")
    for a in anchors:
        ids = [lookup(task(a, 1., u)) for u in ("gap", "bandwidth")]
        checks = {}
        for o in OBSERVABLES:
            x, y = (value(i, o) for i in ids)
            delta = _number(abs(x - y)) if x is not None and y is not None else None
            checks[o] = dict(absolute_delta=delta, threshold=1e-8, passed=delta is not None and delta <= 1e-8)
        passed = all(c["passed"] for c in checks.values())
        report["unit_checks"].append(dict(lambda_value=a["lambda"], records=ids, observables=checks, passed=passed))
        if not passed:
            report["missing_gates"].append(f"Lambda={a['lambda']}: exact reference z=1 unit check")
    report["prerequisites_passed"] = bool(anchors) and not report["missing_gates"]

    def ladder(ids, o, reference):
        values = [value(i, o) for i in ids]
        changes = [abs(x - y) if x is not None and y is not None else None
                   for x, y in zip(values, values[1:], strict=False)]
        error = _number(max(changes)) if all(c is not None for c in changes) else None
        limit = threshold(reference)
        return dict(records=ids, changes=changes, empirical_uncertainty=error,
                    passed=error is not None and limit is not None and error <= limit)

    for o in OBSERVABLES:
        finite, missing = [], list(report["missing_gates"])
        for a in anchors:
            nz, grids, axis_checks = a["nz"], {}, {}
            for size in (4, 8, 16, 32, 64):
                if size > nz:
                    break
                ids = [lookup(task(a, k / size)) for k in range(1, size + 1)]
                vals = [value(i, o) for i in ids]
                mean = _number(sum(v / size for v in vals)) if all(v is not None for v in vals) else None
                grids[str(size)] = dict(records=ids, complete=mean is not None, value=mean)
            grid, coarse = grids[str(nz)], grids[str(nz // 2)]
            x = grid["value"]
            changes = [dict(coarse_nz=int(s) // 2, fine_nz=int(s), absolute_delta=(
                _number(abs(g["value"] - grids[str(int(s) // 2)]["value"]))
                if g["value"] is not None and grids[str(int(s) // 2)]["value"] is not None else None))
                for s, g in grids.items() if int(s) > 4]
            estimates, absent = {}, []
            for axis in NRG_AXES:
                checks = []
                for k in range(1, nz + 1):
                    z = k / nz
                    ids = [lookup(task(a, z, **{axis: s})) for s in a["checks"][axis][-3:]]
                    checks.append(dict(z=z, settings=a["checks"][axis][-3:],
                                       **ladder(ids, o, value(grid["records"][k - 1], o))))
                axis_checks[axis] = checks
                errors = [c["empirical_uncertainty"] for c in checks]
                estimates[axis] = _number(sum(e / nz for e in errors)) if all(e is not None for e in errors) else None
                if not all(c["passed"] for c in checks):
                    absent.append(axis)
            grid_nmax = estimates["nmax"]
            tail = None
            if a.get("tail_checks"):
                ids = [lookup(task(a, 1 / nz, nmax=n)) for n in [a["nmax"], *a["tail_checks"]]]
                tail = dict(z=1 / nz, settings=[a["nmax"], *a["tail_checks"]],
                            **ladder(ids, o, value(ids[0], o)))
                tail_error = tail["empirical_uncertainty"]
                estimates["nmax"] = max(grid_nmax, tail_error) if grid_nmax is not None and tail_error is not None else None
                if not tail["passed"]:
                    absent.append("tail_checks")
            estimates["z"] = changes[-1]["absolute_delta"]
            limit = threshold(x)
            if estimates["z"] is None or limit is None or estimates["z"] > limit:
                absent.append("z")
            if not grid["complete"] or not coarse["complete"]:
                absent.append("complete selected grid")
            error = _number(sum(estimates.values())) if all(e is not None for e in estimates.values()) else None
            finite.append(dict(lambda_value=a["lambda"], nz=nz, task=task(a, 1.), value=x, records=grid["records"],
                grids=grids, z_changes=changes, axis_checks=axis_checks, tail_check=tail, grid_nmax_uncertainty=grid_nmax,
                estimates=estimates, empirical_uncertainty=error, threshold=limit,
                all_axes_passed=not absent, missing_gates=absent, candidate_kind="finite_lambda"))
            missing.extend(f"Lambda={a['lambda']}: {gate}" for gate in absent)
        report["finite_lambda"][o] = finite
        complete = tuple(p["lambda_value"] for p in finite) == LAMBDAS and all(p["value"] is not None for p in finite)
        if not complete:
            missing.append("complete five-Lambda grid: 4, 3, 2.5, 2, 1.8")
        trials, fits = [], {}
        if complete and o != "singlet.moment":
            lam = np.array(LAMBDAS)
            for form, coordinates in (("log_lambda", np.log(lam)), ("lambda_minus_one", lam - 1)):
                for degree in (1, 2):
                    masks = [("full", np.ones(5, dtype=bool))]
                    masks += [("leave_one_out", np.arange(5) != i) for i in range(5)]
                    masks += [("high_lambda_omit", lam < 4)]
                    for kind, mask in masks:
                        points = [p for p, keep in zip(finite, mask, strict=True) if keep]
                        weights = np.linalg.pinv(np.vander(coordinates[mask], degree + 1, increasing=True))[0]
                        prediction = _number(float(weights @ [p["value"] for p in points]))
                        contributions = [{"lambda_value": p["lambda_value"], **{
                            axis: _number(float(abs(w) * e)) if e is not None else None
                            for axis, e in p["estimates"].items()},
                            "total": _number(float(abs(w) * p["empirical_uncertainty"]))
                            if p["empirical_uncertainty"] is not None else None} for p, w in zip(points, weights, strict=True)]
                        totals = {axis: _number(sum(c[axis] for c in contributions))
                                  if all(c[axis] is not None for c in contributions) else None for axis in (*NRG_AXES, "z", "total")}
                        trial = dict(form=f"{form}_degree_{degree}", kind=kind, value=prediction,
                                     lambdas=lam[mask].tolist(), weights=weights.tolist(),
                                     propagated_uncertainty=totals.pop("total"), uncertainty_contributions=totals,
                                     lambda_contributions=contributions)
                        if kind != "full":
                            trial["omitted_lambda"] = float(lam[~mask][0])
                        else:
                            fits[trial["form"]] = prediction
                        trials.append(trial)
        x = fits.get("log_lambda_degree_1")
        envelope = _number(max(t["value"] for t in trials) - min(t["value"] for t in trials)) if trials and all(t["value"] is not None for t in trials) else None
        qualified = [t for t in trials if t["propagated_uncertainty"] is not None]
        worst = max(qualified, key=lambda t: t["propagated_uncertainty"]) if qualified else None
        propagated = worst["propagated_uncertainty"] if worst and len(qualified) == len(trials) else None
        error = _number(propagated + envelope) if propagated is not None and envelope is not None else None
        limit = threshold(x)
        if error is None or limit is None or error > limit:
            missing.append("propagated error plus model/range envelope exceeds target or is undetermined")
        if o == "singlet.moment":
            missing.append("symmetry zero excluded")
        worst_axis, worst_lambda = {}, {}
        for i, t in enumerate(trials):
            for key, contribution in t["uncertainty_contributions"].items():
                if contribution is not None and contribution > worst_axis.get(key, {}).get("contribution", -1):
                    worst_axis[key] = dict(trial=i, contribution=contribution)
            for c in t["lambda_contributions"]:
                key = str(c["lambda_value"])
                if c["total"] is not None and c["total"] > worst_lambda.get(key, {}).get("contribution", -1):
                    worst_lambda[key] = dict(trial=i, contribution=c["total"])
        missing = sorted(set(missing))
        report["nrg_extrapolation"][o] = dict(value=x, empirical_uncertainty=error, threshold=limit,
            status="unresolved" if missing else "empirically_converged", accepted=not missing,
            reason="; ".join(missing) if missing else "conditional empirical stability; not a rigorous continuum bound",
            missing_gates=missing, candidate_kind="lambda_extrapolation", lambdas=[p["lambda_value"] for p in finite],
            fits=fits, fit_trials=trials,
            worst_trial=worst, worst_trial_by_axis=worst_axis, worst_trial_by_lambda=worst_lambda,
            propagated_uncertainty=propagated, model_range_envelope=envelope, sensitivity=envelope,
            estimates=dict(propagated=propagated, model_range_envelope=envelope),
            all_axes_passed=complete and all(p["all_axes_passed"] for p in finite))
    report["final_reference"] = {"nrg": report["nrg_extrapolation"]}
    return _safe(report)
