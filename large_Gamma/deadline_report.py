"""Generate the final, deadline-qualified scientific summary from measured data."""

from datetime import datetime
from pathlib import Path

from qdjj_solver.common.io import write_json
from large_Gamma import analyze, study

BEGIN = "<!-- deadline-results:start -->"
END = "<!-- deadline-results:end -->"


def write_report(directory, output, summary, *, final=False):
    output = Path(output)
    manifest, _, rows, failures, _ = analyze.collect(directory)
    config, profile = manifest["input"], manifest["deadline_profile"]
    comparisons, _ = analyze.compare(rows, config)
    now = datetime.now().astimezone()
    all_finite = all(v["finite_accepted_grid_points"] == summary["expected_base_grid_points"]
                     for v in summary["coverage"].values())
    summary.update(deadline_final=final, four_method_finite_grid_complete=all_finite,
                   report_updated=now.isoformat(), report_deadline=profile["report_deadline"],
                   report_deadline_met=final and now <= datetime.fromisoformat(profile["report_deadline"]))
    summary["status"] = ("final deadline-limited report" if final else "deadline-limited calculations in progress")
    if final and not summary["complete"]:
        summary["status"] += "; remaining bath-resolution limits are explicitly reported"
    write_json(output/"summary.json", summary)
    total = summary["expected_base_grid_points"]
    lines = ["# Deadline-limited large-hybridization results", "",
             f"**Status:** {summary['status']}.", "",
             f"Updated: {now.strftime('%Y-%m-%d %H:%M %Z')}. Requested completion: "
             f"{datetime.fromisoformat(profile['report_deadline']).strftime('%Y-%m-%d %H:%M %z')}.", "",
             "## Revised accuracy and coverage", "",
             "The deadline profile targets **$10^{-5}$ finite-bath observable accuracy** and "
             "**$10^{-3}$ empirical bath resolution**. The original tighter measurements are retained. "
             "Exact ED/BdG references retain their $10^{-8}$ reporting floor.", "",
             "| Method | Measured base-grid points | Accepted finite-bath observables | All observables bath-resolved |",
             "|---|---:|---:|---:|"]
    for method, counts in summary["coverage"].items():
        lines.append(f"| {method} | {counts['measured_grid_points']}/{total} | "
                     f"{counts['finite_accepted_grid_points']}/{total} | "
                     f"{counts['all_observables_bath_converged_grid_points']}/{total} |")
    lines += ["", "A completed deadline-limited report and a bath-converged calculation are distinct statuses. "
              "The last column and `convergence.csv` identify the achieved physical resolution, "
              "including cases whose target remains unresolved.", "",
              "## DMRG accuracy evidence", "",
              "Where exact same-bath ED/BdG is available, DMRG stops when every reported observable "
              "agrees within $10^{-5}$. ED also identifies the lowest eta sector of each parity, "
              "so unnecessary excited eta sectors can be omitted from that DMRG run. "
              "Degenerate eta minima are retained. This is an observable-specific exact-reference "
              "check; the full-state residual and its acceptance flag remain separate diagnostics.", "",
              "Without an exact reference, both eta signs are solved. Acceptance requires two stable "
              "bond refinements, the finest-state sweep check, and a full residual below "
              "$10^{-3}\\Delta$. Only accepted measurements can define reported breakdown boundaries.", ""]
    accepted = [r for r in comparisons if r["method"] == "DMRG" and r["scope"] == "matched_bath" and r["finite_accepted"]]
    if accepted:
        lines.append(f"Largest accepted DMRG discrepancy from a same-bath reference in this snapshot: "
                     f"**{max(r['absolute_error'] for r in accepted):.4g}**, in the normalization of the corresponding observable.")
        lines.append("")
    lines += ["## Signed-gap breakdown brackets", "",
              "The tables give the **first sampled upward crossing of an absolute signed-gap error "
              "of $10^{-3}\\Delta$**. Each interaction/phase row uses the largest cosh bath with "
              "a complete matched-reference Gamma grid for all three QP cutoffs. These are "
              "fixed-bath cutoff-error brackets. Nonmonotonic re-entry and other observables/thresholds "
              "are retained in `breakdown.csv`; continuum claims additionally require the stated bath checks.", ""]
    brackets = analyze.threshold_brackets(comparisons, [1e-3])
    expected_gammas = set(config["gamma_values"])
    for phi in config["phases"]:
        lines += [f"### Phase difference {phi:.12g} radians", "",
                  "| U/Delta | Signed levels per lead | 2QP Gamma/Delta | 4QP Gamma/Delta | 6QP Gamma/Delta |",
                  "|---:|---:|---|---|---|"]
        for u in config["u_values"]:
            selected = [r for r in comparisons if r["u"] == u and r["phi"] == phi
                        and r["scope"] == "matched_bath" and r["family"] == "cosh-grid"
                        and r["bandwidth"] == 100 and r["observable"] == "signed_gap"]
            levels = sorted({r["levels"] for r in selected})
            full = [n for n in levels if all(expected_gammas <= {r["gamma"] for r in selected
                    if r["levels"] == n and r["method"] == method and r["finite_accepted"]}
                    for method in ("2QP", "4QP", "6QP"))]
            if not full:
                lines.append(f"| {u:g} | — | unresolved | unresolved | unresolved |")
                continue
            n = full[-1]
            entries = []
            for method in ("2QP", "4QP", "6QP"):
                candidates = [r for r in brackets if r["u"] == u and r["phi"] == phi
                              and r["method"] == method and r["levels"] == n
                              and r["scope"] == "matched_bath" and r["family"] == "cosh-grid"
                              and r["bandwidth"] == 100 and r["observable"] == "signed_gap"
                              and r["status"] == "below_to_above"]
                if candidates:
                    row = min(candidates, key=lambda r: r["gamma_upper"])
                    entries.append(f"{row['gamma_lower']:g}–{row['gamma_upper']:g}")
                else:
                    entries.append("no resolved upward crossing")
            lines.append(f"| {u:g} | {n} | " + " | ".join(entries) + " |")
        lines.append("")
    lines += ["## Limits and retained evidence", "",
              f"Recorded numerical/resource interruptions: **{len(failures)}**. "
              "A deadline or resource stop is recorded explicitly and is never treated as convergence.", "",
              "The complete input and original-run provenance are in `manifest.json`. "
              "Scalar checkpoints are retained in the ignored run directories. "
              "`observables.csv`, `comparisons.csv`, `convergence.csv`, `bond_convergence.csv` "
              "and `failures.json` preserve the measured values and qualification flags. "
              "The input files preserve all eight interactions, both phases, and the full initial "
              "Gamma grid through 10. No large state files are retained.", ""]
    (output/"report.md").write_text("\n".join(lines), encoding="utf-8")
    readme = study.ROOT/"large_Gamma/README.md"
    text = readme.read_text(encoding="utf-8")
    if BEGIN not in text or END not in text:
        raise ValueError("README is missing the managed deadline-results block")
    before, remainder = text.split(BEGIN, 1)
    _, after = remainder.split(END, 1)
    body = ["", f"**Current status:** {summary['status']}.", "",
            f"Snapshot: {now.strftime('%Y-%m-%d %H:%M %Z')}. "
            "[Current numerical report and breakdown tables](output/report.md).", "",
            "| Method | Accepted finite-bath points | Fully bath-resolved points |",
            "|---|---:|---:|"]
    for method, counts in summary["coverage"].items():
        body.append(f"| {method} | {counts['finite_accepted_grid_points']}/{total} | "
                    f"{counts['all_observables_bath_converged_grid_points']}/{total} |")
    body.append("")
    replacement = before+BEGIN+"\n"+"\n".join(body)+END+after
    temporary = readme.with_name(".README.deadline.tmp")
    temporary.write_text(replacement, encoding="utf-8")
    temporary.replace(readme)
    archive = study.read_json(output/"manifest.json")
    archive["deadline_report_sha256"] = study.file_hash(__file__)
    for path in (output/"report.md", output/"summary.json"):
        archive["artifacts"][path.name] = dict(bytes=path.stat().st_size, sha256=study.file_hash(path))
    write_json(output/"manifest.json", archive)
