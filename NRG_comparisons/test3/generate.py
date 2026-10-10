#!/usr/bin/env python3
"""Regenerate Test3's public report and CSVs from shipped JSON, using only stdlib.

This module also defines the allowlisted publication boundary for raw analysis.
It never imports a solver, the historical analyzer, or a campaign runner.
"""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
OUTPUT_NAMES = ("results.json", "observables.csv", "comparisons.csv", "report.md")
OBSERVABLES = ("energy", "P0", "P1", "P2", "P_up", "P_down", "m_dot", "m_bath",
               "dot_z2", "bath_z2", "Czz", "Czz_connected")
OBSERVABLE_COLUMNS = (
    "task", "stage", "backend", "status", "state", "parity", "twice_sz", "root", "primary",
    "subgap", "degenerate", "residual", "finite_converged", "finite_problem_qualified",
    "qualification", "excitation", "gap_edge_distance", *OBSERVABLES,
)
COMPARISON_COLUMNS = (
    "kind", "old_task", "new_task", "state", "twice_sz", "observable", "status", "reason",
    "old", "new", "difference", "absolute_difference", "target_backend", "target_absolute",
    "target_relative", "target_tolerance", "normalized_change", "within_target",
    "finite_control_tolerance", "finite_control_passed", "matched_bath", "qualified_pair",
)


def fields(names="", **nested):
    return dict.fromkeys(names.split()) | nested


# None means scalar or scalar array, never an unrestricted object. Every object
# has an explicit field allowlist; tuples allow only specified dynamic map keys.
PHYSICAL = fields("geometry gap bandwidth u gamma detuning field")
LIMITS = fields("campaign_seconds cpus mathematica_instances memory_gib task_seconds")
NUMERICAL = fields("keep keepenergy keepmin lambda nmax untruncated z")
TASK = fields("backend bath bath_source chi cutoff id roots spins units",
              numerical=NUMERICAL, physical_changes=PHYSICAL)
CHAIN = fields("delta_over_gap impurity_hopping_over_gap kappa_over_gap n_bath_sites n_bonds "
               "source unused_last_bond_coefficient xi_over_gap zeta_over_gap")
BATH = fields("bandwidth delta weights xi", metadata=fields(
    "kind measure_error positive_pairs weights_renormalized", chain=CHAIN))
SEED = fields("product random_seed randomization_steps")
ROOT_INFO = fields("chi_actual lower_root_overlap norm_error sweep_converged sweep_energy_converged sweeps",
                   seed=SEED, seed_trials=[fields("energy lower_root_overlap sweep_converged", seed=SEED)],
                   sweep_statistics=fields("Delta_E Delta_S E N_updates S max_E_trunc max_S max_chi "
                                           "max_trunc_err norm_err sweep time"))
SOLVER = fields(
    "backend basis_bytes bath_convergence_checked bath_modes bath_reference bath_vacuum_shift "
    "bcs_vacuum_energy coefficient_cleanup_bound coefficient_tolerance coordinate_id decoupled_qp_energies "
    "detuning dimension direct energy_reference estimated_solver_bytes field finite_problem_converged "
    "gamma geometry gram_error groups impurity_modes method mode_order mpo_dimensions mpo_original_dimensions "
    "nnz operator_norm_bound paired_mode_compression phases phi qp_cutoff qp_vacuum real_arithmetic "
    "representation residual_method restriction rho_wn rho_ws ritz_refinement roots_energy_ordered "
    "search_energy_shift search_shift_bound tenpy tunneling u",
    bath=[BATH], gram=fields("imag real"), impurity=fields("detuning field kind u"), roots=[ROOT_INFO],
    coordinates=fields("kind restriction version", electron_annihilators=[[fields("imag operators real")]]),
    options=fields("calculate_residuals chi_max eigenpairs energy_tolerance entropy_tolerance excitation_operator "
                   "group_size lanczos_maxiter lanczos_probability_tolerance max_dimension max_memory_gib max_nnz "
                   "max_sweeps maxiter method min_sweeps mixer mixer_amplitude mixer_sweeps mode_order ncv "
                   "observables orthogonality_tolerance require_convergence residual_tolerance seed "
                   "seed_randomization_steps seed_trials svd_min threads tolerance"),
    sector=fields("eta parity particle_number twice_sz"),
)
SPECTRUM_ROOT = fields("root excitation distance")
DIAGNOSTICS = fields(
    "all_returned_states_converged analytic_bath_vacuum_over_gap backend bath_convergence_checked bath_observables "
    "deck_gap_scale energy_reference energy_units final_energy_offset_over_gap final_h5_iteration "
    "full_spectra_location qualification_scope raw_total_energy_over_gap residual_scope "
    "spectral_completeness_established spin_convention",
    primary_qualified=fields("even odd_down odd_up"), wilson_chain=CHAIN,
    param=fields(extra=fields("B"), param=fields(
        "Gamma Lambda Ninit Nmax T U band bandrescale bcsgap broaden calc0 delta diag_mode diagratio diagth "
        "discretization dm fdm finite h5U h5all h5last h5ops h5raw h5vectors keep keepall keepenergy keepmin "
        "lastall model ops options removefiles savebins silent strategy symtype tri variant z")),
    sector_spectrum_summary=(r"-?\d+", fields(
        "twice_sz root_count saved_roots subgap_root_count minimum_excitation first_two_excitations",
        first_omitted_root=SPECTRUM_ROOT, nearest_excluded_root=SPECTRUM_ROOT,
        nearest_gap_edge_root=SPECTRUM_ROOT)),
    sectors=[fields("complete_sector_diagonalized dimension energy_reference method parity requested_roots "
                    "returned_roots solved_roots twice_sz", solver_metadata=SOLVER,
                    timings=fields("basis diagonalization hamiltonian observables_and_residuals preparation total"))],
)
STATE = fields("degenerate excitation finite_converged id parity residual root subgap twice_sz",
               values=fields(" ".join(OBSERVABLES)), analysis=fields(
                   "finite_problem_qualified gap_edge_distance gap_edge_margin primary qualification role"))
HASHES = (r"[\w./+-]+", None)
PUBLIC = fields(
    "format_version scope limitations",
    inputs=fields("format_version scope", physical=PHYSICAL, limits=LIMITS,
                  targets={key: fields("absolute relative") for key in ("qp", "dmrg", "nrg")}),
    controls=fields("passed", checks=[fields("maximum_error missing name passed tolerance")]),
    costs=fields("maximum_peak_rss_gib monitored_task_count scope total_monitored_seconds", limits=LIMITS),
    gates=fields("continuum_target_established controls controls_issues pilot_escalation_allowed production "
                 "production_approved spectral_completeness_established uncertainty_kind"),
    comparisons=[fields(" ".join(COMPARISON_COLUMNS))],
    provenance=fields("campaign_sha256 controls_sha256 analyzer_sha256 raw_artifacts hash_scope",
                      recorded_source_hashes=HASHES,
                      software=fields("machine python system numpy scipy qdjj_solver qpsolver",
                                      packages=HASHES),
                      solver=fields("sha256 version", files=HASHES),
                      external_executables={key: fields("sha256 version") for key in ("math", "nrg", "nrginit")},
                      external_source_hashes=HASHES),
    tasks=[fields("backend id stage status", task=TASK,
                  request=fields(physical=PHYSICAL, bath=BATH, task=TASK),
                  qualification=fields("continuum_target_established finite_problem_qualified reported_qualified scope"),
                  provenance=fields("issues recorded_request_sha256 request_file_sha256 request_path request_sha256 "
                                    "result_path result_sha256 public_request_sha256 public_result_sha256"),
                  result=fields("status qualified monitored_seconds peak_rss_gib error",
                                physical=PHYSICAL, numerical=NUMERICAL, states=[STATE], diagnostics=DIAGNOSTICS,
                                artifacts=fields("DONE data nrg.log nrginit.log param raw.h5"),
                                coverage=fields("continuum_completeness_established degenerate_subgap_ids "
                                                "diagonal_observables_unique expected_three_subgap_states ground_energy "
                                                "ground_ids same_sector_degeneracy_tolerance_over_gap scope sectors subgap_ids")))],
)


def privacy_guard(value):
    """Reject private payloads even in otherwise allowed scientific text fields."""
    text = json.dumps(value, allow_nan=False)
    if re.search(r"(?<![\w.])/[A-Za-z_][\w.-]*/[\w./-]+|[A-Za-z]:\\\\|"
                 r"\b\d{4}-\d{2}-\d{2}\b|\b[\w.-]+@[\w.-]+\b|"
                 r'"(?:hostname|host|uid|username|job_id|node|session|execution|identity|started_unix)"\s*:', text):
        raise ValueError("private execution metadata is not publishable")


def project(value, schema):
    if value is None:
        return None
    if schema is SEED and isinstance(value, list):
        return project(value, None)  # A product seed or a randomized-seed record.
    if isinstance(schema, dict):
        if not isinstance(value, dict):
            raise ValueError("expected an object in a publication record")
        return {key: project(value[key], child) for key, child in schema.items() if key in value}
    if isinstance(schema, tuple):
        return {key: project(child, schema[1]) for key, child in value.items() if re.fullmatch(schema[0], key)}
    if isinstance(schema, list):
        return [project(child, schema[0]) for child in value]
    if isinstance(value, dict):
        raise ValueError("unexpected object in a scalar publication field")
    if isinstance(value, list):
        return [project(child, None) for child in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def public_projection(report):
    """Keep science and original hashes, never recursive execution identities.

    Original request/result hashes refer to archived bytes. Public hashes use
    sorted-key JSON of the projected nested objects, not the original receipts.
    """
    source = dict(report)
    if "campaign" in source:
        campaign = source["campaign"]
        identity = campaign["identity"]
        source["inputs"] = identity["inputs"]
        external = {}
        for path, checksum in identity.get("external_sources", {}).items():
            name = Path(path).name
            if name in external:
                raise ValueError("external source basename collision")
            external[name] = checksum
        source["provenance"] = report["provenance"] | dict(
            recorded_source_hashes=identity.get("source_hashes", {}),
            software=campaign.get("environment", {}), solver=campaign.get("solver", {}),
            external_executables=identity.get("tools", {}), external_source_hashes=external,
            hash_scope="Original hashes identify archived inputs/receipts, not current public files. "
                       "public_request_sha256 and public_result_sha256 hash sorted-key JSON of the projected objects.",
        )
    public = project(source, PUBLIC)
    for task in public["tasks"]:
        for key in ("request", "result"):
            task["provenance"][f"public_{key}_sha256"] = digest(task[key])
    privacy_guard(public)
    return public


def report_markdown(report):
    report = public_projection(report)

    def display(value):
        if value is None:
            return "n/a"
        if isinstance(value, float):
            return f"{value:.10g}"
        return str(value).replace("|", "\\|").replace("\n", " ")

    def table(columns, rows):
        return ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |",
                *("| " + " | ".join(display(value) for value in row) + " |" for row in rows), ""]

    by_id = {task["id"]: task for task in report["tasks"]}
    selected = [(by_id[key], state) for key in ("cosh8-qp-qfull", "cosh8-dmrg-chi128", "nrg-k2000")
                if key in by_id for name in ("odd_down", "odd_up", "even", "sz+0-r1")
                for state in by_id[key]["result"]["states"] if state["id"] == name]
    text = ["# Test3: Finite-Field Single-Reservoir Results", "", "## TL;DR", "",
            "**PILOT ONLY. Production deferred; no continuum extrapolation or uncertainty established.** "
            "The full-QP, highest-bond DMRG and highest-retention single-z endpoints are shown first when saved, "
            "including the second even root. "
            "No empirical NRG uncertainty from Test1 is transferred.", "",
            f"Controls: **{report['gates']['controls']}**. This regenerates saved data, not a new calculation.", "",
            "## Model", "",
            "Single superconducting reservoir: gap Delta=1, normal half-bandwidth D=100, U=2, total Gamma=0.4, "
            "detuning=0.5, local field B=0.03. Positive B raises spin up. Energies subtract the corresponding "
            "disconnected BCS vacuum and retain the centered impurity constant. "
            "`even` is parity-even, M=0, not an imposed singlet; the odd states have M=-1/2 and +1/2.", "",
            "## Energies and Gap", "",
            "Excitations are relative to each task's lowest saved energy, in gap units. "
            "QP and DMRG use the same explicit cosh8 bath; NRG uses Lambda=2, z=1, keep=2000, Nmax=42. "
            "The NRG comparison is an unmatched-bath diagnostic, not a reference or z average. "
            "Displayed digits are rounded stored values, not established accuracy.", ""]
    text += table(["Endpoint", "State", r"$E/\Delta$", r"Excitation / $\Delta$", "Subgap", "Degenerate", "Residual", "Qualification"], (
        [task["id"], row["id"], row["values"]["energy"], row["excitation"], row["subgap"], row["degenerate"], row["residual"],
         row["analysis"]["qualification"]] for task, row in selected))
    extra = [(task, row) for task in report["tasks"] if task["stage"] == "pilot"
             for row in task["result"]["states"] if row["subgap"] and not row["analysis"]["primary"]]
    if extra:
        text += ["**Additional Subgap States:** the saved pilot exclusion roots include subgap states. "
                 "`sz+0-r1` labels the second even root, not a discarded observation. This does not establish "
                 "an additional continuum bound state or spectral completeness.", ""]
    text += [
             "## Probabilities and Moment", "",
             "P0/P1/P2 are empty/singly/doubly occupied impurity probabilities. P_up/P_down are singly occupied "
             "spin probabilities, not densities. Moments are in hbar units. Each row is state-specific: "
             "same-sector degeneracy makes diagonal expectations basis dependent. Qualified full-QP rows concern "
             "only this finite bath; DMRG values remain unqualified and NRG residuals are unavailable, not zero.", ""]
    for names in (("P0", "P1", "P2", "P_up", "P_down", "m_dot", "m_bath"),
                  ("dot_z2", "bath_z2", "Czz", "Czz_connected")):
        text += table(["Endpoint", "State", *names], (
            [task["id"], row["id"], *(row["values"][name] for name in names)] for task, row in selected))
    text += ["## Achieved Accuracy", "",
             "Residual qualification, QP cutoff/MPS bond stability, and bath convergence are separate. "
             "Close DMRG energies do not override failed residuals. NRG retention changes at one z and Lambda "
             "are not a continuum uncertainty. No backend has established the requested continuum accuracy.", "",
             "Targets remain absolute 1e-9 plus relative 1e-6 for QP/DMRG; 1e-6 plus 1e-3 for NRG. "
             "Normalized changes use the new backend and abs(new). Finite controls instead use their recorded "
             "absolute bounds; these targets are not error bars.", "",
             "## Data and Reproduction", "",
             "[Results](results.json), [all state observables](observables.csv), and "
             "[comparison operands](comparisons.csv) are self-contained. Run "
             "`python -B -S NRG_comparisons/test3/generate.py --output-dir PATH` from the repository root. "
             "No ignored runs, installed backend, or raw HDF5 file is needed. "
             "See the [case README](../README.md) for fresh numerical reproduction.", "",
             "## Supplemental material", "", "### Refinement Checks", "",
             "All recorded primary comparisons are retained, including unavailable, unqualified and "
             "basis-dependent rows. A small normalized change is not an error bound.", ""]
    groups = {}
    for row in report["comparisons"]:
        groups.setdefault((row["kind"], row["old_task"], row["new_task"]), []).append(row)
    text += table(["Kind", "Old", "New", "Compared / Total", "Max Normalized Change", "Other Statuses"], (
        [*key, f"{sum(row['status'] == 'compared' for row in rows)}/{len(rows)}",
         max((row["normalized_change"] for row in rows if row["status"] == "compared"), default=None),
         ", ".join(sorted({row["status"] for row in rows if row["status"] != "compared"}))]
        for key, rows in groups.items()))
    text += ["### Residuals", "", "Every saved root, including exclusions, is retained below and in the CSV. "
             "QP residuals refer to the projected cutoff Hamiltonian. DMRG uses the untruncated finite Hamiltonian. "
             "An unqualified exclusion root does not automatically invalidate a primary root.", ""]
    text += table(["Task", "State", "Excitation", "Subgap", "Degenerate", "Residual", "Qualification"], (
        [task["id"], row["id"], row["excitation"], row["subgap"], row["degenerate"], row["residual"],
         row["analysis"]["qualification"]] for task in report["tasks"] for row in task["result"]["states"]))
    text += ["### Controls", "", *(f"- {issue}" for issue in report["gates"]["controls_issues"]), ""]
    if report["controls"] is not None:
        text += table(["Recorded Control", "Maximum Error", "Tolerance", "Passed"], (
            [check["name"], check.get("maximum_error"), check.get("tolerance"), check.get("passed")]
            for check in report["controls"].get("checks", [])))
    text += ["### Task Status and Costs", "", "Elapsed task times and memory are measurements, not scaling estimates.", ""]
    text += table(["Task", "Stage", "Status", "Primary Qualified", "Seconds", "Peak GiB", "Issue"], (
        [task["id"], task["stage"], task["status"], task["qualification"]["finite_problem_qualified"],
         task["result"].get("monitored_seconds"), task["result"].get("peak_rss_gib"),
         "; ".join([task["result"].get("error", ""), *task["provenance"]["issues"]]).strip("; ")]
        for task in report["tasks"]))
    text += ["### Limitations", "", *(f"- {note}" for note in report["limitations"]), "",
             "### Provenance", "", "The public projection keeps explicit baths, chain coefficients, numerical inputs, "
             "all saved state values, control checks, status and qualifications. Original receipt/source hashes "
             "identify archived originals, not the projected objects or current source files. Public request/result "
             "hashes separately identify sorted-key JSON of the projected nested objects. Raw artifact hashes are "
             "copied, not revalidated here. Full NRG spectra remain in the original receipts; compact sector summaries "
             "retain omitted roots and gap-edge evidence. Private execution identities are not publication data.", ""]
    return "\n".join(text)


def write_report(report, output_dir):
    """Validate the complete public payload before creating or replacing files."""
    report = public_projection(report)
    observations = []
    for task in report["tasks"]:
        for row in task["result"]["states"]:
            observations.append(dict(
                task=task["id"], stage=task["stage"], backend=task["backend"], status=task["status"], state=row["id"],
                **{key: row.get(key) for key in ("parity", "twice_sz", "root", "subgap", "degenerate", "residual",
                                               "finite_converged", "excitation")},
                **{key: row["analysis"][key] for key in ("primary", "finite_problem_qualified", "qualification", "gap_edge_distance")},
                **row["values"],
            ))
    assets = {"results.json": json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
              "report.md": report_markdown(report)}
    for name, columns, rows in (("observables.csv", OBSERVABLE_COLUMNS, observations),
                                ("comparisons.csv", COMPARISON_COLUMNS, report["comparisons"])):
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
        assets[name] = stream.getvalue()
    privacy_guard(assets)
    output_dir = Path(output_dir).resolve()
    if any((output_dir / name).is_symlink() for name in OUTPUT_NAMES):
        raise ValueError("publication assets must not be symlinks")
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, text in assets.items():
        (output_dir / name).write_bytes(text.encode())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "runs/presentation")
    args = parser.parse_args(argv)
    published = (HERE / "output").resolve()
    output = args.output_dir.resolve()
    if any(path == published or published in path.parents for path in
           [output, *((output / name).resolve() for name in OUTPUT_NAMES)]):
        parser.error("use a separate directory; published outputs are read-only")
    source = (published / "results.json").read_bytes()
    report = json.loads(source)
    if public_projection(report) != report:
        raise ValueError("shipped dataset is not an allowlisted public projection")
    write_report(report, output)
    if (output / "results.json").read_bytes() != source:
        raise ValueError("shipped dataset is not in canonical publication format")
    print(f"Regenerated {len(report['tasks'])} tasks; production deferred.")


if __name__ == "__main__":
    main()
