#!/usr/bin/env python3
"""Read-only analysis of test3 controls and cost pilots (standard library only)."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from NRG_comparisons.test3.manual_run.physics import OBSERVABLES, PRIMARY
from NRG_comparisons.test3 import generate


OUTPUT_NAMES = generate.OUTPUT_NAMES
COMPARISON_OBSERVABLES = (*OBSERVABLES, "excitation")
OBSERVABLE_COLUMNS = generate.OBSERVABLE_COLUMNS
COMPARISON_COLUMNS = generate.COMPARISON_COLUMNS


def digest(value):
    """The runner's canonical request digest, distinct from the file-byte hash."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def read_json(path):
    """Hash exactly the bytes parsed, even when a live runner replaces a receipt."""
    data = path.read_bytes()
    return json.loads(data), hashlib.sha256(data).hexdigest()


def compact_result(result, backend):
    """Keep every saved observable row, but summarize full final-shell spectra."""
    result = dict(result)
    diagnostics = dict(result.get("diagnostics", {}))
    spectra = diagnostics.pop("all_sector_spectra", {})
    states = []
    for original in result.get("states", []):
        row = dict(original)
        primary = row["root"] == 0 and row["twice_sz"] in PRIMARY
        finite = None if backend == "nrg" else (
            row.get("finite_converged") is True and row.get("residual") is not None
            and not row["degenerate"]
        )
        row["analysis"] = dict(
            primary=primary, role="primary" if primary else "exclusion_diagnostic",
            finite_problem_qualified=finite,
            qualification="unavailable_nrg" if finite is None else "qualified" if finite else "not_qualified",
            gap_edge_distance=abs(row["excitation"] - 1.),
            gap_edge_margin=1. - row["excitation"],
        )
        states.append(row)
    result["states"] = states
    summaries = {}
    for sector, spectrum in spectra.items():
        energies = spectrum["ground_relative_energy_over_gap"]
        spin = spectrum["twice_sz"]
        saved = {row["root"] for row in states if row["twice_sz"] == spin}
        omitted = [(root, energy) for root, energy in enumerate(energies) if root not in saved]
        above = [(root, energy) for root, energy in enumerate(energies) if energy >= 1.]
        nearest = min(enumerate(energies), key=lambda pair: abs(pair[1]-1.), default=None)
        summaries[sector] = dict(
            twice_sz=spin, root_count=len(energies), saved_roots=sorted(saved),
            subgap_root_count=sum(energy < 1. for energy in energies),
            minimum_excitation=min(energies, default=None),
            first_two_excitations=energies[:2],
            first_omitted_root=None if not omitted else dict(root=omitted[0][0], excitation=omitted[0][1]),
            nearest_excluded_root=None if not above else dict(root=above[0][0], excitation=above[0][1]),
            nearest_gap_edge_root=None if nearest is None else dict(
                root=nearest[0], excitation=nearest[1], distance=abs(nearest[1]-1.)),
        )
    if spectra:
        diagnostics["sector_spectrum_summary"] = summaries
        diagnostics["full_spectra_location"] = "original result.json: diagnostics.all_sector_spectra"
    result["diagnostics"] = diagnostics
    return result


def comparison_rows(kind, old, new, targets, matched_bath):
    """Target scales always use the new backend and abs(new), not the baseline."""
    pair_issues = []
    if all(task["status"] == "completed" and task["request"] is not None for task in (old, new)):
        if old["request"]["physical"] != new["request"]["physical"]:
            pair_issues.append("physical inputs differ")
        if matched_bath:
            if "nrg" not in (old["backend"], new["backend"]):
                baths = [task["request"].get("bath", {}) for task in (old, new)]
                matched_bath = all(key in baths[0] and baths[0][key] == baths[1].get(key)
                                   for key in ("xi", "weights", "delta", "bandwidth"))
            else:
                chains = [task["result"]["diagnostics"].get("wilson_chain", {}) if task["backend"] == "nrg"
                          else task["request"].get("bath", {}).get("metadata", {}).get("chain", {})
                          for task in (old, new)]
                matched_bath = all(
                    key in chains[0] and key in chains[1] and len(chains[0][key]) == len(chains[1][key])
                    and all(math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
                            for a, b in zip(chains[0][key], chains[1][key], strict=True))
                    for key in ("xi_over_gap", "zeta_over_gap", "delta_over_gap", "kappa_over_gap")
                )
            if not matched_bath:
                pair_issues.append("matched bath not verified from explicit request/chain coefficients")
    old_states = {(r["twice_sz"], r["root"]): r for r in old["result"]["states"]}
    new_states = {(r["twice_sz"], r["root"]): r for r in new["result"]["states"]}
    target_backend = "qp" if new["backend"] == "ed" else new["backend"]
    target = targets[target_backend]
    rows = []
    for spin, state in PRIMARY.items():
        a, b = old_states.get((spin, 0)), new_states.get((spin, 0))
        for name in COMPARISON_OBSERVABLES:
            x = None if a is None else (a.get(name) if name == "excitation" else a["values"].get(name))
            y = None if b is None else (b.get(name) if name == "excitation" else b["values"].get(name))
            reasons = list(pair_issues)
            for task in (old, new):
                if task["status"] != "completed":
                    reasons.append(f"{task['id']}: {task['status']}")
                if task["provenance"]["issues"]:
                    reasons.append(f"{task['id']}: receipt integrity not established")
            if x is None or y is None:
                reasons.append("missing state or observable")
            unavailable = bool(reasons)
            ambiguous = name not in ("energy", "excitation") and any(
                row is not None and row["degenerate"] for row in (a, b))
            if ambiguous:
                reasons.append("same-sector degeneracy: diagonal expectation is basis dependent")
            qualified = not unavailable and not ambiguous and all(
                task["backend"] == "nrg" or (row.get("finite_converged") is True
                                              and row.get("residual") is not None)
                for task, row in ((old, a), (new, b))
            )
            if not unavailable and not qualified and not ambiguous:
                reasons.append("finite-problem state qualification failed")
            if kind == "unmatched_bath_single_z":
                reasons.append("unmatched bath, single-z NRG; not a reference or uncertainty")
            difference = None if x is None or y is None else y-x
            tolerance = None if y is None else target["absolute"] + target["relative"]*abs(y)
            normalized = None if difference is None else abs(difference)/tolerance
            rows.append(dict(
                kind=kind, old_task=old["id"], new_task=new["id"], state=state, twice_sz=spin,
                observable=name, status="unavailable" if unavailable else "basis_dependent" if ambiguous
                else "compared" if qualified else "unqualified", reason="; ".join(reasons),
                old=x, new=y, difference=difference,
                absolute_difference=None if difference is None else abs(difference),
                target_backend=target_backend, target_absolute=target["absolute"],
                target_relative=target["relative"], target_tolerance=tolerance,
                normalized_change=normalized,
                within_target=None if unavailable or ambiguous else normalized <= 1.,
                finite_control_tolerance=1e-8 if kind == "finite_control" else None,
                finite_control_passed=(abs(difference) < 1e-8 if qualified else None)
                if kind == "finite_control" else None,
                matched_bath=matched_bath, qualified_pair=qualified,
            ))
    return rows


def analyze(campaign):
    """Build a snapshot without running solvers, changing receipts, or using other cases."""
    campaign = Path(campaign).resolve()
    manifest, manifest_hash = read_json(campaign / "campaign.json")
    identity = manifest["identity"]
    inputs = identity["inputs"]
    inventory = []
    for stage in ("controls", "pilot"):
        for task in identity["tasks"][stage]:
            taskid = task["id"]
            if Path(taskid).name != taskid or taskid in (".", ".."):
                raise ValueError("task id must be a directory basename")
            directory = campaign / taskid
            request, result, request_hash, result_hash = None, {}, None, None
            issues = []
            for name in ("request", "result"):
                path = directory / f"{name}.json"
                if path.exists():
                    try:
                        value, checksum = read_json(path)
                    except (OSError, ValueError) as error:
                        issues.append(f"{name}.json unreadable: {error}")
                        continue
                    if name == "request":
                        request, request_hash = value, checksum
                    else:
                        result, result_hash = value, checksum
            status = result.get("status", "incomplete" if directory.exists() else "unstarted")
            expected_physical = inputs["physical"] | task.get("physical_changes", {})
            request_digest = None if request is None else digest(request)
            if request is not None and (request.get("task") != task or request.get("physical") != expected_physical):
                issues.append("request differs from campaign task/physical inputs")
            if status == "completed" and (request is None or result.get("request_sha256") != request_digest):
                issues.append("completed result lacks a matching request digest")
            result = compact_result(result, task["backend"])
            primary = [row for row in result["states"] if row["analysis"]["primary"]]
            finite_qualified = None if task["backend"] == "nrg" else (
                status == "completed" and not issues and len(primary) == len(PRIMARY)
                and all(row["analysis"]["finite_problem_qualified"] for row in primary)
                and result.get("qualified") is True
            )
            inventory.append(dict(
                id=taskid, stage=stage, backend=task["backend"], status=status, task=task,
                request=request, result=result,
                qualification=dict(finite_problem_qualified=finite_qualified,
                                   reported_qualified=result.get("qualified"),
                                   scope="three primary roots; exclusion roots qualified separately",
                                   continuum_target_established=False),
                provenance=dict(request_path=f"{taskid}/request.json", result_path=f"{taskid}/result.json",
                                request_file_sha256=request_hash, request_sha256=request_digest,
                                result_sha256=result_hash, recorded_request_sha256=result.get("request_sha256"),
                                source_hashes=identity.get("source_hashes", {}), issues=issues),
            ))
    by_id = {task["id"]: task for task in inventory}
    if len(by_id) != len(inventory):
        raise ValueError("duplicate task id in campaign inventory")
    control_path = campaign / "controls.json"
    controls, controls_hash = None, None
    control_issues = []
    if control_path.exists():
        try:
            controls, controls_hash = read_json(control_path)
        except (OSError, ValueError) as error:
            control_issues.append(f"controls.json unreadable: {error}")
    if controls is None or controls.get("passed") is not True or not controls.get("checks") or not all(
        check.get("passed") is True for check in controls["checks"]
    ):
        control_issues.append("controls.json does not record all checks passed")
    for task in inventory:
        if task["stage"] == "controls" and (task["status"] != "completed" or task["provenance"]["issues"]):
            control_issues.append(f"{task['id']}: {task['status']}; control receipt incomplete or unverified")

    pairs = []
    for label in ("atomic", "chain", "zero", "reverse"):
        for backend in ("nrg", "qp", "dmrg"):
            if f"{label}-{backend}" in by_id and f"{label}-ed" in by_id:
                pairs.append(("finite_control", f"{label}-ed", f"{label}-{backend}", True))
    if {"chain-nrg", "chain-bandwidth"} <= by_id.keys():
        pairs.append(("finite_control", "chain-nrg", "chain-bandwidth", True))
    qp = sorted((task for task in inventory if task["stage"] == "pilot" and task["backend"] == "qp"
                 and task["task"].get("bath") == "cosh8"),
                key=lambda task: float("inf") if task["task"].get("cutoff") is None else task["task"]["cutoff"])
    dmrg = sorted((task for task in inventory if task["stage"] == "pilot" and task["backend"] == "dmrg"
                   and task["task"].get("bath") == "cosh8"), key=lambda task: task["task"]["chi"])
    nrg = sorted((task for task in inventory if task["stage"] == "pilot" and task["backend"] == "nrg"),
                 key=lambda task: task["task"]["numerical"]["keep"])
    for kind, series in (("qp_cutoff", qp), ("dmrg_bond", dmrg), ("nrg_retention", nrg)):
        pairs.extend((kind, old["id"], new["id"], True) for old, new in zip(series, series[1:], strict=False))
    full = next((task for task in qp if task["task"].get("cutoff") is None), None)
    if full is not None:
        pairs.extend(("qp_full_finite_bath", task["id"], full["id"], True) for task in qp[:-2])
        pairs.extend(("dmrg_vs_full_qp", full["id"], task["id"], True) for task in dmrg)
        available = [task for task in nrg if task["status"] == "completed" and not task["provenance"]["issues"]]
        if available and full["status"] == "completed":
            pairs.append(("unmatched_bath_single_z", full["id"], available[-1]["id"], False))
    comparisons = []
    for kind, old, new, matched in pairs:
        rows = comparison_rows(kind, by_id[old], by_id[new], inputs["targets"], matched)
        comparisons.extend(rows)
        if kind == "finite_control" and all(by_id[key]["status"] == "completed" for key in (old, new)):
            if any(row["status"] == "unavailable" for row in rows):
                control_issues.append(f"{new} vs {old}: primary finite-control comparison unavailable")
    observed = [task["result"]["monitored_seconds"] for task in inventory if "monitored_seconds" in task["result"]]
    memory = [task["result"]["peak_rss_gib"] for task in inventory if "peak_rss_gib" in task["result"]]
    return dict(
        format_version=1, scope="controls and cost pilots only; no production results",
        campaign=manifest,
        provenance=dict(campaign_path=str(campaign), campaign_sha256=manifest_hash,
                        controls_sha256=controls_hash, analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        raw_artifacts="hashes copied from receipts, not revalidated by this analyzer"),
        gates=dict(controls="passed" if not control_issues else "blocked_controls", controls_issues=control_issues,
                   pilot_escalation_allowed=not control_issues, production="deferred",
                   production_approved=False, continuum_target_established=False,
                   spectral_completeness_established=False, uncertainty_kind="not_established"),
        controls=controls, tasks=inventory, comparisons=comparisons,
        costs=dict(monitored_task_count=len(observed), total_monitored_seconds=sum(observed) if observed else None,
                   maximum_peak_rss_gib=max(memory, default=None), limits=inputs.get("limits", {}),
                   scope="sum of saved task monitoring times, not total campaign elapsed time"),
        limitations=[
            "Three primary branches are selected, not proof of spectral completeness. expected_three_subgap_states describes returned roots only.",
            "Excitations use each task's lowest saved energy; untested sectors and roots can change that minimum.",
            "Same-sector degeneracies can make diagonal expectations basis dependent; cross-sector energy ties do not.",
            "Finite-problem residual/uniqueness qualification is separate from QP cutoff, MPS bond and bath convergence.",
            "NRG residuals are unavailable. Retention changes at fixed Lambda, z and chain length are not continuum uncertainty.",
            "Targets normalize differences by absolute + relative*abs(new), using the new backend (ED uses QP targets).",
            "Finite controls use a separate absolute 1e-8 comparison; recorded analytic and derivative checks have their own tolerances.",
            "Single-z NRG versus cosh8 full QP is an unmatched-bath diagnostic, not a reference, z average or uncertainty estimate.",
            "No uncertainty or convergence claim is inherited from test1, test2, or another campaign. Production requires separate approval.",
        ],
    )


def report_markdown(report):
    return generate.report_markdown(report)


def write_report(report, output_dir):
    """Write only the four stable presentation assets; never modify the campaign."""
    output_dir = Path(output_dir).resolve()
    if "campaign_path" in report["provenance"]:
        campaign = Path(report["provenance"]["campaign_path"])
        if output_dir == campaign or output_dir.is_relative_to(campaign):
            raise ValueError("output directory must be outside the input campaign")
    generate.write_report(report, output_dir)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output-dir", type=Path, help="explicitly allow replacement of presentation assets here")
    args = parser.parse_args(argv)
    output = args.output_dir or HERE / "runs/presentation"
    if args.output_dir is None and any((output / name).exists() for name in OUTPUT_NAMES):
        parser.error("default presentation already exists; choose an explicit --output-dir to replace it")
    write_report(analyze(args.campaign), output)
    print(f"Wrote {output.resolve() / 'report.md'}; production deferred.")


if __name__ == "__main__":
    main()
