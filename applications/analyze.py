"""Summarize the archived convergence ladders without running any solver.

Run from the source root after the selected cases' three profiles have completed:
python -B -m applications.analyze
"""

import argparse
from pathlib import Path

from applications._common import read_csv, read_json, write_json

ROOT = Path(__file__).resolve().parent
CASES = ("zitko_2023_knight_shift", "zonda_2023_double_dot", "zalom_2024_multiterminal",
         "paaske_2023_surrogate_spectrum", "bargerbos_2022_parity_diagram", "choi_2004_kondo_josephson",
         "zonda_2016_unequal_gaps", "kadlecova_2017_asymmetry", "hecht_2008_gap_edge",
         "bobok_2025_chain_expansion")


def differences(rows, left, right, keys, values):
    def select(setting):
        return {tuple(row[k] for k in keys): row for row in rows if row["setting"] == setting}
    a, b = select(left), select(right)
    if not a or not a.keys() <= b.keys():
        raise ValueError(f"incomplete/mismatched convergence checkpoints: {left}, {right}")
    return {value: max(abs(float(a[key][value])-float(b[key][value])) for key in a)
            for value in values}


def analyze(case):
    if case == "large_Gamma":
        # This source-root benchmark publishes its own convergence assessment.
        return read_json(ROOT.parent/"large_Gamma/output/summary.json")
    if case == "bobok_2025_chain_expansion":
        from applications.bobok_2025_chain_expansion.analyze import analyze as chain_analysis
        return chain_analysis()
    output = ROOT/case/"output"
    rows = read_csv(output/"convergence"/"convergence.csv")
    results = {"paper": read_json(output/"paper"/"summary.json")}
    if case == CASES[0]:
        for a, b in (("L4", "L8"), ("L6", "L8"), ("L8", "L10"),
                     ("L8_w200", "L8"), ("L8_q4", "L8"), ("cosh6_q6", "L10")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, ("gamma_over_u", "phi_over_pi"), ("kappa",))
        phase = read_csv(output/"paper"/"phase.csv")
        kappa = [float(r["kappa"]) for r in phase]
        results["phase_modulation"] = max(kappa)-min(kappa)
    elif case == CASES[1]:
        for a, b in (("L2", "L4"), ("L3", "L4"), ("L5_q6", "L4"), ("L5_full", "L5_q6"), ("L5_full", "L4"),
                     ("L4_q4", "L4"), ("L4_w200", "L4"), ("L3_dmrg", "L3")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, ("curve", "phi_over_pi"),
                                                      ("signed_gap", "current", "even_current", "odd_current"))
        comparison = read_csv(output/"paper"/"comparison.csv")
        for curve in ("weak", "strong"):
            selected = [r for r in comparison if r["curve"] == curve]
            # A different transition position gives a finite jump error, not a
            # smooth-current amplitude error. Report both explicitly.
            same_phase = [r for r in selected if float(r["current"])*float(r["nrg_current"]) > 0]
            opposite = [r for r in selected if float(r["current"])*float(r["nrg_current"]) < 0
                        and abs(float(r["nrg_current"])) > .00015]
            results[f"{curve}_same_sign_max_current_error"] = max(abs(float(r["difference"])) for r in same_phase)
            results[f"{curve}_opposite_sign_reference_phases"] = [float(r["phi_over_pi"]) for r in opposite]
    elif case == CASES[2]:
        roots = {row["setting"]: float(row["chi_critical"]) for row in rows}
        results["critical_chi_by_setting"] = roots
        results["L8_versus_L10"] = abs(roots["L8"]-roots["L10"])
        results["fit_window_change"] = abs(roots["L8_w2000"]-roots["L8"])
        results["cutoff4_change"] = abs(roots["L8_q4"]-roots["L8"])
    elif case == CASES[3]:
        for a, b in (("L4", "L8"), ("L6", "L8"), ("L7", "L8"), ("L8", "L10"),
                     ("L8_w20", "L8"), ("L8_q4", "L8")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, ("gamma",), ("signed_gap", "doublet_dot_spin"))
        published = read_csv(output/"paper"/"published_surrogates.csv")
        results["published_surrogate_checkpoints"] = len(published)
        results["maximum_published_surrogate_difference"] = max(abs(float(r["difference"])) for r in published)
        comparison = read_csv(output/"paper"/"comparison.csv")
        results["largest_nrg_error_checkpoint"] = max(comparison, key=lambda r: abs(float(r["difference"])))
    elif case == CASES[4]:
        for a, b in (("L2", "L4"), ("L3", "L4"), ("L4", "L5"), ("L4_w20", "L4"), ("L4_q4", "L4")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, ("phi_over_pi",), ("detuning_critical_over_u",))
        boundary = read_csv(output/"paper"/"boundary.csv")
        results["boundary_at_phase_endpoints"] = [r for r in boundary if float(r["phi_over_pi"]) in (0., 1.)]
    elif case == CASES[5]:
        for a, b in (("L6", "L8"), ("L8", "L10"), ("L8_w2D", "L8"), ("L8_q4", "L8"),
                     ("L12_q6", "L12_q8"), ("L12_q8", "L10")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, ("delta_over_tk", "phi_over_pi"), ("signed_gap", "current"))
        comparison = read_csv(output/"paper"/"comparison.csv")
        results["reference_comparison_by_ratio"] = {}
        for ratio in sorted({r["delta_over_tk"] for r in comparison}, key=float):
            data = [r for r in comparison if r["delta_over_tk"] == ratio]
            away = [r for r in data if .02 < float(r["phi_over_pi"]) < .95]
            results["reference_comparison_by_ratio"][ratio] = dict(
                max_current_difference=max(abs(float(r["difference"])) for r in data),
                max_difference_away_from_endpoints=max(abs(float(r["difference"])) for r in away),
                opposite_sign_phases=[float(r["phi_over_pi"]) for r in away
                                      if float(r["current"])*float(r["nrg_current"]) < 0])
        if (output/"dmrg"/"convergence.csv").exists():
            backend_rows = read_csv(output/"dmrg"/"convergence.csv")
            results["L8_dmrg_versus_L8"] = differences(rows+backend_rows, "L8_dmrg", "L8",
                ("delta_over_tk", "phi_over_pi"), ("even_energy", "odd_energy", "signed_gap", "current"))
            results["dmrg_maximum_residual"] = read_json(output/"dmrg"/"summary.json")["maximum_residual"]
            results["dmrg_runtime_seconds"] = read_json(output/"dmrg"/"manifest.json")["elapsed_seconds"]
    elif case in CASES[6:8]:
        keys, values = (("gap_ratio", "detuning"), ("current", "signed_gap", "charge")) if case == CASES[6] else (
            ("u_meV", "asymmetry", "phi_over_pi"), ("tilde_epsilon",))
        for a, b in (("L2", "L4"), ("L3", "L4"), ("L4", "L5"), ("L4_q4", "L4"),
                     ("L4_w200", "L4"), ("L4_D200", "L4")):
            results[f"{a}_versus_{b}"] = differences(rows, a, b, keys, values)
    else:
        results["quadratic_error_budget"] = [r for r in rows if r["kind"] == "quadratic"]
        krylov = read_csv(output/"convergence"/"krylov.csv")
        results["interacting_resolution"] = [r for r in krylov if r["spin"] == "sum"]
        data = read_csv(output/"convergence"/"comparison.csv")
        for width in (.2, .4):
            for left, right in (("L4_k240", "L6_k240"), ("L6_k240", "L8_k240"), ("L8_k240", "L8_k480")):
                selected = []
                for row in data:
                    for label in (left, right):
                        if row["kind"] == f"interacting-{label}-b{width}":
                            selected.append(dict(row, setting=label))
                results[f"{left}_versus_{right}_b{width}"] = differences(
                    selected, left, right, ("delta_over_D", "omega_prime_over_D"), ("pi_gamma_A",))
        results["interacting_continuum_converged"] = False
    results["runtimes_seconds"] = {profile: read_json(output/profile/"manifest.json")["elapsed_seconds"]
                                  for profile in ("quick", "paper", "convergence")}
    write_json(output/"validation.json", results)
    return results


def main():
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=(*CASES, "large_Gamma"))
    args = parser.parse_args()
    for case in (args.case,) if args.case else CASES:
        print(case)
        print(json.dumps(analyze(case), indent=2))


if __name__ == "__main__":
    main()
