"""Validate the shipped scalar snapshot and regenerate tables, using only stdlib.

No numerical calculation or historical analyzer is imported. Checksums of private
originals are references, not hashes of this public payload. Recorded solver
checks are not a public revalidation of MPS, HDF5, or executable artifacts.
"""

import argparse
from collections import Counter
from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import re

if __package__:
    from .physics import OBSERVABLES, PRIMARY, finish_states, paper_offset, validate_physical
else:
    from physics import OBSERVABLES, PRIMARY, finish_states, paper_offset, validate_physical


HERE = Path(__file__).resolve().parent
CORE = ((1, -1, 0), (1, 1, 0), (0, 0, 0), (0, 0, 1))
QUANTITIES = (*OBSERVABLES, "excitation")
TARGETS = {"nrg": (1e-6, .001), "qp": (1e-9, 1e-6), "dmrg": (1e-9, 1e-6)}
TIERS = ("requested", "one_percent", "five_percent")
CHAIN_KEYS = ("xi_over_gap", "zeta_over_gap", "delta_over_gap", "kappa_over_gap")
ARTIFACTS = ("results.json", "report.md", "observables.csv", "comparisons.csv", "z_averages.csv")
CLAIMS = dict(production_approved=False, continuum_target_established=False,
              spectral_completeness_established=False, uncertainty_established=False)

# Explicit recursive allowlist for scalar numerical metadata. Coordinates and
# repeated bath/model declarations are replaced by content-addressed references.
METADATA_FIELDS = set("""
backend bath_convergence_checked bath_leads bath_levels_per_lead bath_modes_total
bath_weight_sum_per_lead chi chi_convergence_checked continuum_finite_version
coordinate_id coverage_scope current_is_real current_units degeneracy_scope
detuning dimension energy_reference energy_units full_qp_comparison_checked
group_size hamiltonian_sha256 identical_bath_record_for_both_leads layout
layout_description lowest_root_completeness_established mode_order norm_error
norm_error_tolerance options paper_energy_offset primary_qualified
primary_residual_qualified primary_unique qualification_scope requested_roots
requested_spins residual_scope residual_tolerance_over_gap returned_roots
root_selection_seed_independence sector spectral_completeness_established
trial_data_scope weights_renormalized all_returned_states_converged
calculate_residuals chi_max eigenpairs energy_tolerance entropy_tolerance
excitation_operator lanczos_maxiter lanczos_probability_tolerance max_sweeps
min_sweeps mixer mixer_amplitude mixer_sweeps observables orthogonality_tolerance
require_convergence residual_tolerance seed seed_randomization_steps seed_trials
svd_min threads eligible_energy_spread_over_gap eligible_trials
independent_trial_agreement root selected_seed trial_count product random_seed
randomization_steps eta parity particle_number twice_sz complete_sector_diagonalized
solved_roots bath_modes finite_problem_converged gram gram_error groups
impurity_modes mpo_dimensions mpo_original_dimensions operator_norm_bound qp_cutoff
residual_method restriction roots roots_energy_ordered search_energy_shift
search_shift_bound tenpy imag real chi_actual lower_root_overlap sweep_converged
sweep_energy_converged sweep_statistics sweeps energy Delta_E Delta_S E N_updates
S max_E_trunc max_S max_chi max_trunc_err norm_err sweep time even odd_down odd_up
analytic_bath_vacuum_over_gap bath_observables current_normalization deck_gap_scale
final_energy_offset_over_gap final_h5_iteration paper_offset_over_gap
raw_total_energy_over_gap representation spin_convention root_count saved_roots
subgap_root_count unsaved_subgap_roots minimum_excitation first_three_excitations
first_omitted_root nearest_excluded_root nearest_gap_edge_root excitation distance
channel delta_over_gap kappa_over_gap n_bath_sites n_bonds
unused_last_bond_coefficient xi_over_gap zeta_over_gap kind levels_per_lead
identical_physical_leads measure_error positive_pairs qp_weights_scope qp_weights
basis_bytes bath_reference bath_vacuum_shift bcs_vacuum_energy
coefficient_cleanup_bound coefficient_tolerance decoupled_qp_energies
estimated_solver_bytes method nnz paired_mode_compression phases phi qp_vacuum
real_arithmetic rho_wn rho_ws ritz_refinement tunneling u field gamma geometry
impurity direct leads matrix max_dimension max_memory_gib max_nnz maxiter ncv
tolerance version root_prefix format_version requested_roots coordinate_sha256
complete elapsed_seconds largest_block_dimension largest_block_dimension_reason
max_pretruncation_dimension mechanism_counts shells solver_vmpeak_gib timings_seconds
available_space_exhausted energy_window keep_cap keepmin_floor_or_boundary
lastall_override exact_window_count exact_window_count_reason iteration
kept_physical_states kept_reduced_states mechanism pretruncation_physical_dimension
pretruncation_reduced_dimension safeguard_added_reduced_states scale_over_gap
selection_boundary_over_gap selection_boundary_shell window_over_gap All Other diag
recalc checkpoint_seconds helper_seconds initial_load_seconds timings basis
diagonalization observables_and_residuals preparation total hamiltonian
""".split())
SPEC_FIELDS = set("""backend bath bath_source chi cutoff layout max_sweeps min_sweeps
roots roots_by_spin seed_trials spins twice_sz max_dimension max_memory_gib max_nnz
method numerical representation units role control""".split())
NUMERICAL_FIELDS = {"keep", "keepenergy", "keepmin", "lambda", "nmax", "untruncated", "z"}


def digest(value):
    """Public canonical JSON SHA-256: sorted keys, compact separators, ASCII."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def metadata(value):
    """Reject unrecognized metadata rather than treating arbitrary dicts as safe."""
    if isinstance(value, dict):
        require(all(k in METADATA_FIELDS or re.fullmatch(r"-?\d+", k) for k in value),
                "metadata key outside scientific allowlist")
        require(value.get("spectral_completeness_established", False) is False,
                "recorded spectral completeness is not established")
        for v in value.values():
            metadata(v)
    elif isinstance(value, list):
        for v in value:
            metadata(v)
    elif isinstance(value, float):
        require(math.isfinite(value), "nonfinite metadata")
    elif isinstance(value, str):
        require(not re.search(r"(?:/project/|/home/|/tmp/|/share/|SLURM_|SBATCH_|Traceback)", value),
                "private metadata string")
    else:
        require(value is None or type(value) in (bool, int), "invalid metadata scalar")


def fields(value, allowed, label):
    require(isinstance(value, dict) and set(value) <= set(allowed.split()), "unexpected "+label+" field")


def family(data, record):
    spec = data["specifications"][record["specification"]]
    runtime = data["runtimes"][record["runtime"]]
    return [record["model"], spec.get("representation", "paper"), spec.get("units", "gap"),
            {k: v for k, v in spec.get("numerical", {}).items() if k not in NUMERICAL_FIELDS},
            [r for r in runtime["sources"] if r["identity"].startswith("external/")
             or r["identity"].endswith(("/model.m", "/modeloperators.m"))]]


def state_key(row):
    return row["parity"], row["twice_sz"], row["root"]


def state_map(record):
    return {state_key(r): r for r in record["states"]}


def scalar(record, key, observable):
    rows = state_map(record)
    if observable == "odd_splitting":
        if not {(1, 1, 0), (1, -1, 0)} <= rows.keys():
            return None
        return rows[1, 1, 0]["values"]["energy"]-rows[1, -1, 0]["values"]["energy"]
    row = rows.get(key)
    return None if row is None else row.get("excitation") if observable == "excitation" else row["values"][observable]


def field_bound(data, record):
    physical = data["models"][record["model"]]
    value = scalar(record, (1, 1, 0), "odd_splitting")
    applicable = physical["ws"] == 0
    limit = abs(physical["field"])/physical["gap"]
    return dict(applicable=applicable, splitting=value, limit=limit,
                passed=None if not applicable or value is None else abs(value) <= limit+1e-8)


def finite_qualified(data, record, row):
    """Versioned recorded finite-problem checks; never cutoff/bath certification."""
    if record["backend"] == "nrg":
        return None
    residual = row["residual"]
    if (row.get("finite_converged") is not True or type(residual) not in (float, int)
            or not 0 <= residual <= (1e-7 if record["backend"] == "dmrg" else 1e-8)):
        return False
    if record["backend"] != "dmrg":
        return True
    calls = [c for c in record["calls"] if c["twice_sz"] == row["twice_sz"]]
    if len(calls) != 1:
        return False
    call = calls[0]
    m = data["solvers"][call["solver"]]
    roots = m.get("roots", [])
    r = row["root"]
    if r >= len(roots) or call["root_prefix"] != list(range(len(roots))):
        return False
    root = roots[r]
    gram = m.get("gram", {})
    count = len(roots)
    if any(len(gram.get(k, [])) != count or any(len(row) != count for row in gram[k]) for k in ("real", "imag")):
        return False
    gram_error = max(math.hypot(gram["real"][i][j]-int(i == j), gram["imag"][i][j])
                     for i in range(count) for j in range(count))
    return (m.get("finite_problem_converged") is True and m.get("roots_energy_ordered") is True
            and m.get("gram_error", math.inf) <= 1e-8 and gram_error <= 1e-8
            and root.get("norm_error", math.inf) <= 1e-8
            and root.get("lower_root_overlap", math.inf) <= 1e-8
            and root.get("sweep_converged") is True and root.get("sweep_energy_converged") is True)


def issues(data, record, *, four=True):
    """Whole-inventory gate for means, not an ordinary state comparison gate."""
    rows, errors = state_map(record), []
    required = CORE if four else ((1, -1, 0), (1, 1, 0), (0, 0, 0))
    if record["status"] != "completed" or not record["recorded_integrity_ok"]:
        errors.append("recorded completion/integrity unavailable")
    if not set(required) <= rows.keys():
        errors.append("incomplete primary/root inventory")
    if any(k not in rows or rows[k]["degenerate"] or not rows[k].get("subgap", False) for k in required):
        errors.append("missing, degenerate or above-gap corresponding state")
    if field_bound(data, record)["passed"] is False:
        errors.append("exact WS=0 field bound violated")
    if record["states"]:
        ground = min(r["values"]["energy"] for r in record["states"])
        if min((rows[k]["values"]["energy"] for k in required if k in rows and k[2] == 0), default=math.inf) > ground+1e-8:
            errors.append("available ground lies outside primary reference")
        if any(r["degenerate"] or finite_qualified(data, record, r) is False
               for r in record["states"] if r["values"]["energy"] <= ground+1e-8):
            errors.append("unqualified available ground reference")
    for key in required:
        if key in rows and finite_qualified(data, record, rows[key]) is False:
            errors.append("finite-problem checks failed")
    if any(s["unsaved_subgap_roots"] for s in record["spectrum_summary"]):
        errors.append("unsaved subgap roots")
    return sorted(set(errors))


def reference_issues(data, record):
    rows = state_map(record)
    errors = []
    if not set(CORE[:3]) <= rows.keys():
        errors.append("incomplete primary excitation reference")
    if record["states"]:
        ground = min(r["values"]["energy"] for r in record["states"])
        if any(r["degenerate"] or finite_qualified(data, record, r) is False
               for r in record["states"] if r["values"]["energy"] <= ground+1e-8):
            errors.append("unqualified available ground reference")
    else:
        errors.append("missing available ground reference")
    return errors


def state_issues(data, record, key, observable):
    """Absolute scalars need their own state, not four qualified global branches."""
    errors = []
    if record["status"] != "completed" or not record["recorded_integrity_ok"]:
        errors.append("recorded completion/integrity unavailable")
    if field_bound(data, record)["passed"] is False:
        errors.append("exact WS=0 field bound violated")
    keys = ((1, -1, 0), (1, 1, 0)) if observable == "odd_splitting" else (key,)
    rows = state_map(record)
    for selected in keys:
        row = rows.get(selected)
        if row is None:
            errors.append("missing corresponding state")
            continue
        if row["degenerate"]:
            errors.append("basis-dependent degenerate state")
        if row.get("subgap") is False:
            errors.append("corresponding state above gap")
        if finite_qualified(data, record, row) is False:
            errors.append("finite-problem state checks failed")
        if observable == "excitation" and "excitation" not in row:
            errors.append("single-sector excitation unavailable")
    if observable == "excitation":
        errors.extend(reference_issues(data, record))
    return sorted(set(errors))


def verify_call(data, record, call):
    """A valid registry hash does not establish a valid call-to-solver binding."""
    m, d = data["solvers"][call["solver"]], record["diagnostics"]
    spec = data["specifications"][record["specification"]]
    spin, meta = call["twice_sz"], call["metadata"]
    require(spec["backend"] == record["backend"] == d.get("backend", record["backend"]), "call backend mismatch")
    sector = {"parity": spin % 2, "twice_sz": spin}
    require(all(meta.get(k, v) == v for k, v in sector.items()), "call sector metadata mismatch")
    if not m:  # Electron ED has no QP/MPS solver object.
        require(record["backend"] == "ed", "missing finite solver metadata")
        return
    require(m.get("backend", record["backend"]) == record["backend"], "solver backend mismatch")
    require(all(m["sector"].get(k) == v for k, v in sector.items()), "solver sector mismatch")
    if "twice_sz" in spec:
        require(spin == spec["twice_sz"] and d["sector"] == m["sector"], "diagnostic sector mismatch")
    options = m["options"]
    for k in ("max_sweeps", "min_sweeps", "seed_trials", "max_dimension", "max_memory_gib", "max_nnz", "method"):
        if k in spec and k in options:
            require(spec[k] == options[k], "actual solver option differs: "+k)
    if "options" in d:
        require(d["options"] == options, "diagnostic solver options mismatch")
    for k in ("coordinate_id", "hamiltonian_sha256", "mode_order", "layout", "group_size"):
        if k in d and k in m:
            require(d[k] == m[k], "diagnostic solver "+k+" mismatch")
    if "dimension" in m and "dimension" in meta:
        require(m["dimension"] == meta["dimension"], "solver dimension mismatch")
    require(m["bath_modes"] == 4*len(data["baths"][record["bath"]]["xi"]), "solver bath dimension mismatch")
    if record["backend"] == "dmrg":
        order = m["mode_order"]
        require(sorted(order) == list(range(m["bath_modes"]+m["impurity_modes"])), "invalid mode ordering")
        declared = options["mode_order"]
        require(order == (list(range(len(order))) if declared is None else declared), "option mode ordering mismatch")
        group_size = options["group_size"]
        require(m["groups"] == [order[i:i+group_size] for i in range(0, len(order), group_size)], "solver groups/order mismatch")
        require(d.get("group_size", group_size) == group_size, "diagnostic group size mismatch")
        require(len(m["roots"]) == len(call["root_prefix"]) == options["eigenpairs"], "solver root-count mismatch")
        require(d.get("coordinate_id", m["coordinate_id"]) == m["coordinate_id"], "diagnostic coordinate mismatch")
        if "norm_error" in d:
            require(d["norm_error"] == [r["norm_error"] for r in m["roots"]], "diagnostic root norms mismatch")


def matched_baths(data, left, right):
    a, b = (data["baths"][r["bath"]] for r in (left, right))
    if left["backend"] != "nrg" and right["backend"] != "nrg":
        return all(a.get(k) == b.get(k) for k in ("xi", "weights", "delta", "bandwidth"))
    ac, bc = a.get("chains", []), b.get("chains", [])
    if len(ac) != 2 or len(bc) != 2:
        return False
    sds = any(data["specifications"][r["specification"]].get("representation") == "sds" for r in (left, right))
    sds = sds and all(data["models"][r["model"]][k] == 0 for r in (left, right) for k in ("wn", "ws"))
    return all(len(x[k]) == len(y[k]) and all(math.isclose(abs(v) if sds and k == "delta_over_gap" else v,
               abs(w) if sds and k == "delta_over_gap" else w, rel_tol=1e-12, abs_tol=1e-12)
               for v, w in zip(x[k], y[k], strict=True))
               for x, y in zip(ac, bc, strict=True) for k in CHAIN_KEYS)


def verify_data(data):
    """Validate public hashes and scalar algebra without opening private originals."""
    require(set(data) == {"schema", "checksum", "models", "baths", "specifications", "solvers",
                         "runtimes", "records", "outcomes", "control", "profiles", "claims",
                         "provenance", "recorded_control_checks"}, "unexpected dataset fields")
    require(data["schema"] == "test4-public-scalars-v1", "unsupported schema")
    require(data["claims"] == CLAIMS and all(v is False for v in data["claims"].values()), "unsupported coverage/production claim")
    require(data["checksum"] == digest({k: v for k, v in data.items() if k != "checksum"}), "public checksum mismatch")
    for table in ("models", "baths", "specifications", "solvers", "runtimes"):
        for key, value in data[table].items():
            require(key == digest(value), table+" public checksum mismatch")
    for physical in data["models"].values():
        validate_physical(physical)
    for spec in data["specifications"].values():
        require(set(spec) <= SPEC_FIELDS, "unrecognized scientific specification")
        require(set(spec.get("numerical", {})) <= NUMERICAL_FIELDS, "unrecognized NRG option")
        if "numerical" in spec:
            n = spec["numerical"]
            require(set(NUMERICAL_FIELDS)-{"untruncated"} <= n.keys() and n["lambda"] >= 4 and 0 < n["z"] <= 1,
                    "invalid NRG declaration")
    for runtime in data["runtimes"].values():
        fields(runtime, "software sources verification", "runtime")
        require(runtime["verification"] == "recorded", "runtime verification overclaim")
        fields(runtime["software"], "python numpy scipy qdjj_solver qpsolver machine system packages", "software")
        fields(runtime["software"].get("packages", {}), "h5py numpy physics-tenpy qdjj-solver scipy threadpoolctl", "packages")
        for source in runtime["sources"]:
            require(set(source) == {"identity", "sha256"} and not Path(source["identity"]).is_absolute()
                    and ".." not in Path(source["identity"]).parts and re.fullmatch(r"[0-9a-f]{64}", source["sha256"]),
                    "invalid canonical source identity/hash")
    for solver in data["solvers"].values():
        metadata(solver)
    for bath in data["baths"].values():
        require(set(bath) <= {"xi", "weights", "delta", "bandwidth", "metadata", "chains"}, "unrecognized bath field")
        metadata(bath.get("metadata", {}))
        if "xi" in bath:
            require(len(bath["xi"]) == len(bath["weights"]) and all(w > 0 for w in bath["weights"]), "invalid bath measure")
        if "chains" in bath:
            require(len(bath["chains"]) == 2, "both Wilson leads required")
            for i, chain in enumerate(bath["chains"], 1):
                metadata(chain)
                require(chain["channel"] == i and all(len(chain[k]) == chain["n_bath_sites"] for k in CHAIN_KEYS), "Wilson channel/array mismatch")
    ids, identities = set(), set()
    for record in data["records"]:
        require(set(record) == {"id", "backend", "role", "status", "model", "bath", "specification", "runtime",
                               "states", "calls", "diagnostics", "spectrum_summary", "retention", "checkpoint",
                               "cost", "private_original", "recorded_integrity_ok", "reported_qualified"}, "unexpected record field")
        for table, ref in (("models", "model"), ("baths", "bath"), ("specifications", "specification"), ("runtimes", "runtime")):
            require(record[ref] in data[table], "dangling "+ref)
        require(record["id"] not in ids, "ID collision")
        ids.add(record["id"])
        identity = (record["private_original"]["result_sha256"], record["model"], record["bath"], record["specification"])
        require(identity not in identities, "duplicate receipt/scientific-input identity")
        identities.add(identity)
        require(record["private_original"]["label"] == "PRIVATE_ORIGINAL_REFERENCE", "original checksum mislabeled")
        metadata(record["diagnostics"])
        metadata(record["retention"])
        fields(record["checkpoint"], "sha256 chi twice_sz roots coordinate_id hamiltonian_sha256 companion_sha256 parent_sha256", "checkpoint")
        fields(record["private_original"], "label result_sha256 request_file_sha256 request_sha256 retention_sha256", "original reference")
        fields(record["cost"], "monitored_seconds peak_rss_gib peak_vm_gib", "measured cost")
        for summary in record["spectrum_summary"]:
            metadata(summary)
        rows = record["states"]
        require(bool(rows), "completed record without scientific rows")
        checked = deepcopy(rows)
        coverage = finish_states(checked)
        physical = data["models"][record["model"]]
        require(abs(coverage["paper_energy_offset"]-paper_offset(physical)/physical["gap"]) <= 1e-8, "paper offset mismatch")
        for row, check in zip(rows, checked, strict=True):
            require(set(row) <= {"id", "parity", "twice_sz", "root", "degenerate", "finite_converged", "residual", "excitation", "subgap", "values"}, "unexpected state field")
            require(row["degenerate"] == check["degenerate"], "unrecorded same-sector degeneracy")
            if "excitation" in row:
                require(row["excitation"] == check["excitation"] and row["subgap"] == check["subgap"], "excitation reference mismatch")
            if record["backend"] == "nrg":
                require(row["residual"] is None and "finite_converged" not in row, "invented NRG residual")
            else:
                require(type(row.get("finite_converged")) is bool, "finite-converged must be recorded boolean")
        for spin in {r["twice_sz"] for r in rows}:
            roots = sorted(r["root"] for r in rows if r["twice_sz"] == spin)
            require(roots == list(range(len(roots))), "non-prefix roots")
            energies = [r["values"]["energy"] for r in sorted(rows, key=lambda r: r["root"]) if r["twice_sz"] == spin]
            require(all(a <= b+1e-8 for a, b in zip(energies, energies[1:], strict=False)), "roots are not energy ordered")
        spec = data["specifications"][record["specification"]]
        if record["backend"] == "dmrg" and "twice_sz" in spec:
            require(all("excitation" not in r and "subgap" not in r for r in rows), "single-sector global excitation")
            require({(r["twice_sz"], r["root"]) for r in rows} == {(spec["twice_sz"], r) for r in range(spec["roots"])}, "DMRG same-call prefix changed")
        for call in record["calls"]:
            fields(call, "solver twice_sz root_prefix metadata private_original_sha256", "solver call")
            metadata(call["metadata"])
            require(call["solver"] in data["solvers"], "missing solver record")
            prefix = [r["root"] for r in rows if r["twice_sz"] == call["twice_sz"]]
            require(call["root_prefix"] == prefix, "call root prefix mismatch")
            m = data["solvers"][call["solver"]]
            if record["backend"] == "qp" and "qp_cutoff" in m:
                require(m["qp_cutoff"] == spec["cutoff"] or spec["cutoff"] is None and m["qp_cutoff"] == m["bath_modes"],
                        "actual QP cutoff differs from declaration")
            if record["backend"] == "dmrg":
                require(m["options"]["chi_max"] == spec["chi"], "actual MPS chi differs from declaration")
            if "root_prefix" in m:
                solved = m["root_prefix"]
                require(solved == list(range(len(solved))) and solved[:len(prefix)] == prefix
                        and (record["backend"] != "dmrg" or solved == prefix), "recorded solved-root prefix differs")
            verify_call(data, record, call)
        require(record["backend"] == "nrg" or sorted(c["twice_sz"] for c in record["calls"]) ==
                sorted({r["twice_sz"] for r in rows}), "missing/duplicate sector calls")
        if record["checkpoint"]:
            checkpoint = record["checkpoint"]
            require(all(checkpoint[k] == record["diagnostics"][k] for k in ("hamiltonian_sha256", "coordinate_id")),
                    "checkpoint Hamiltonian/coordinate mismatch")
            require(checkpoint["roots"] == len(rows), "checkpoint prefix mismatch")
            require(checkpoint["twice_sz"] == spec["twice_sz"] and checkpoint["chi"] == spec["chi"], "checkpoint call mismatch")
        for value in record["cost"].values():
            require(type(value) in (int, float) and math.isfinite(value) and value >= 0, "invalid measured cost")
    require(len({p["id"] for p in data["profiles"]}) == len(data["profiles"]), "duplicate profile identity")
    for profile in data["profiles"]:
        verify_profile(data, profile)
    for outcome in data["outcomes"]:
        fields(outcome, "id specification model status stage recorded_stages controlled_stages reason_code private_original_result_sha256", "outcome")
    fields(data["provenance"], "original_reference_label source_reports verification finite_check_policy", "provenance")
    require(data["provenance"]["finite_check_policy"] == "recorded-finite-v1", "unknown finite-check policy")
    for source in data["provenance"]["source_reports"]:
        require(set(source) == {"identity", "sha256"}, "unexpected original source reference")
    for check in data["recorded_control_checks"]:
        fields(check, "name maximum_error measured tolerance passed missing", "recorded control")
    fields(data["control"], "bath spectra ed qp private_original_sha256 verification finite_diagnostics", "full control")
    for sector in data["control"]["spectra"]:
        require(set(sector) == {"twice_sz", "ground_relative_energy_over_gap", "indices_in_energy_order"}, "unexpected full spectrum field")
    for backend in ("ed", "qp"):
        for row in data["control"][backend]:
            fields(row, "id parity twice_sz root degenerate finite_converged residual excitation subgap values", "full control row")
        if "finite_diagnostics" in data["control"]:
            info = data["control"]["finite_diagnostics"][backend]
            fields(info, "diagnostics calls reported_qualified", "full finite control")
            metadata(info["diagnostics"])
            for call in info["calls"]:
                fields(call, "twice_sz root_prefix solver metadata", "full finite control call")
                metadata(call["metadata"])
                require(call["solver"] in data["solvers"], "missing full-control solver")
    # A checksum is integrity, not a privacy policy: deny operational keys even
    # when an edited payload has been resealed by its author.
    forbidden = {"path", "hostname", "host", "user", "uid", "account", "job_id", "command", "environment",
                 "authorization", "parent_campaign", "source_path", "timestamp", "release", "traceback", "slurm"}

    def privacy(value):
        if isinstance(value, dict):
            require(not forbidden.intersection(value), "private operational metadata")
            for v in value.values():
                privacy(v)
        elif isinstance(value, list):
            for v in value:
                privacy(v)
        elif isinstance(value, str):
            require(not re.search(r"/project/|/home/|/tmp/|/share/|SLURM_|SBATCH_|Traceback", value), "private operational value")
    privacy(data)
    return data


def average(data, profile):
    """All quantities, including connected Czz and excitations, average per z."""
    verify_profile(data, profile)
    index = {r["id"]: r for r in data["records"]}
    nz = profile["nz"]
    expected = [j/nz for j in range(1, nz+1)]
    members = [index.get(name) for name in profile["members"]]
    errors = []
    if len(members) != nz or any(t is None for t in members):
        errors.append("missing exact uniform z grid")
    present = [t for t in members if t is not None]
    settings = [data["specifications"][t["specification"]]["numerical"] for t in present]
    if ([n["z"] for n in settings] != expected or len({t["model"] for t in present}) != 1
            or any(any(n[k] != profile[k] for k in ("lambda", "nmax", "keepenergy", "keepmin")) for n in settings)
            or [n["keep"] for n in settings] != profile["keep_by_z"]):
        errors.append("different numerical/physical grid family")
    if present and any(family(data, t) != family(data, present[0]) for t in present):
        errors.append("different physical/model/tool/representation family")
    if "reference_profile" in profile:
        reference = next((p for p in data["profiles"] if p["id"] == profile["reference_profile"]), None)
        if reference is None or len(reference["members"]) != len(members) or any(
                t is None or name not in index or not matched_baths(data, t, index[name])
                for t, name in zip(members, reference["members"], strict=True)):
            errors.append("endpoint replacement must match both reference Wilson leads")
    for task in present:
        errors.extend(task["id"]+": "+e for e in issues(data, task))
    rows = []
    for key in CORE:
        for observable in (*QUANTITIES, *(('odd_splitting',) if key == (1, 1, 0) else ())):
            values = [scalar(t, key, observable) for t in present]
            complete = len(values) == nz and all(v is not None for v in values)
            rows.append(dict(profile=profile["id"], parity=key[0], twice_sz=key[1], root=key[2], observable=observable,
                mean=math.fsum(values)/nz if complete else None, minimum=min(values) if complete else None,
                maximum=max(values) if complete else None, spread=max(values)-min(values) if complete else None,
                per_z_values=values if complete else None, grid_complete=complete and not any("grid" in e for e in errors),
                qualified=complete and not errors, reason="; ".join(sorted(set(errors)))))
    return rows


def verify_profile(data, profile):
    require(set(profile)-{"reference_profile"} == {"id", "kind", "nz", "lambda", "nmax", "keepenergy", "keepmin", "keep_by_z", "members"}, "unexpected profile field")
    require(type(profile["nz"]) is int and profile["nz"] in (2, 4, 8), "unsupported exact z grid")
    require(len(profile["members"]) == len(profile["keep_by_z"]) == profile["nz"], "profile grid length mismatch")
    if profile["kind"] == "uniform":
        require("reference_profile" not in profile and len(set(profile["keep_by_z"])) == 1, "uniform profile cannot vary retention")
        return
    require(profile["kind"] == "mixed_endpoint", "unknown averaging profile kind")
    reference = next((p for p in data["profiles"] if p["id"] == profile.get("reference_profile")), None)
    require(reference is not None and reference["kind"] == "uniform", "mixed profile requires uniform reference")
    verify_profile(data, reference)
    require(all(profile[k] == reference[k] for k in ("nz", "lambda", "nmax", "keepenergy", "keepmin")), "mixed profile numerical family mismatch")
    require(profile["nz"] == 8 and profile["lambda"] == 8 and reference["keep_by_z"] == [4000]*8
            and profile["keep_by_z"] == [5000]+[4000]*7 and profile["members"][0] == "nrg-l8-z1-k5000"
            and reference["members"][0] == "nrg-l8-z1-k4000" and profile["members"][1:] == reference["members"][1:],
            "only the declared z=.125 endpoint replacement is permitted")


def compare(kind, left, right, key, observable, x, y, errors, *, backend):
    available = x is not None and y is not None
    difference = y-x if available else None
    absolute, relative = TARGETS[backend]
    row = dict(kind=kind, old_id=left, new_id=right, parity=key[0], twice_sz=key[1], root=key[2],
               observable=observable, old=x, new=y, difference=difference,
               relative_percent=100*abs(difference)/abs(y) if available and abs(y) > 1e-6 else None,
               qualified=available and not errors, reason="; ".join(sorted(set(errors))),
               target_backend=backend, target_absolute=absolute, target_relative=relative,
               target_tolerance=None if y is None else absolute+relative*abs(y))
    for tier, a, r in (("requested", absolute, relative), ("one_percent", 1e-6, .01), ("five_percent", 1e-6, .05)):
        observed = abs(difference) <= a+r*abs(y) if available else None
        row["observed_"+tier] = observed
        row[tier] = observed if row["qualified"] else None
    row["normalized_change"] = abs(difference)/row["target_tolerance"] if available else None
    return row


def assemble(data, names, *, require_qualified=True):
    """Assemble compatible entire sector calls; never splice even roots."""
    index = {r["id"]: r for r in data["records"]}
    tasks = [index[name] for name in names]
    require(tasks and all(t["backend"] == "dmrg" for t in tasks), "assembly requires DMRG calls")
    fingerprints = []
    spins = set()
    rows = []
    for t in tasks:
        d = t["diagnostics"]
        fingerprint = [t["model"], t["bath"], *(d.get(k) for k in
                       ("hamiltonian_sha256", "coordinate_id", "layout", "mode_order", "group_size"))]
        require(all(v is not None for v in fingerprint), "missing assembly fingerprint")
        fingerprints.append(fingerprint)
        sector = {r["twice_sz"] for r in t["states"]}
        require(len(sector) == 1 and not sector & spins, "cannot splice same-sector calls")
        spins |= sector
        require(sorted(r["root"] for r in t["states"]) == list(range(len(t["states"]))), "non-prefix assembly call")
        if require_qualified:
            for r in t["states"]:
                require(finite_qualified(data, t, r) is True and not r["degenerate"], "unqualified DMRG root")
        rows.extend(deepcopy(t["states"]))
    require(all(f == fingerprints[0] for f in fingerprints), "different assembly Hamiltonian/bath/coordinates")
    require({-1, 0, 1} <= spins, "incomplete assembly ground reference")
    finish_states(rows)
    return rows


def replay_controls(data):
    """Transplanted scalar control algebra; no old controller/analyzer dependency."""
    records = {r["id"]: r for r in data["records"]}
    checks = []

    def rows(name):
        return {(r["twice_sz"], r["root"]): r for r in records[name]["states"]}

    def check(name, error, tolerance=1e-8):
        checks.append(dict(name=name, maximum_error=error, tolerance=tolerance, passed=error < tolerance,
                           verification="scalar-replayed"))

    controls = [r for r in data["records"] if r["role"] == "control" and r["id"] != "nrg-l8-control"]
    require(len(controls) == 37, "37 original controls required")
    pairs = [("junction-bandwidth", "junction-nrg"), ("junction-sds", "junction-nrg")]
    for t in controls:
        spec = data["specifications"][t["specification"]]
        if t["backend"] == "nrg":
            continue
        expected = {(s, r) for s in spec["spins"] for r in range(min(
            spec.get("roots_by_spin", {}).get(str(s), spec["roots"]), math.comb(10, 5+s)))}
        check(t["id"]+" inventory", len(expected-rows(t["id"]).keys()), .5)
        check(t["id"]+" recorded finite checks", int(not all(finite_qualified(data, t, r) for r in t["states"])), .5)
        checks[-1]["verification"] = "recorded checks and scalar residuals"
        if not t["id"].startswith("hf-"):
            pairs.append((spec["bath_source"], t["id"]) if t["backend"] == "ed" else
                         (t["id"], t["id"].rsplit("-", 1)[0]+"-ed"))
    for left, right in pairs:
        a, b = rows(left), rows(right)
        check(left+" covers "+right, len(b.keys()-a.keys()), .5)
        require(records[left]["model"] == records[right]["model"], "control physical parameters differ")
        require(matched_baths(data, records[left], records[right]), "control Wilson leads differ")
        check(left+" vs "+right, max(abs(a[k]["values"][o]-b[k]["values"][o])
              for k in a.keys() & b.keys() for o in (("energy", "energy_paper")
              if a[k]["degenerate"] or b[k]["degenerate"] else OBSERVABLES)))
    atomic = rows("atomic-nrg")
    check("atomic energies, paper offset, Zeeman sign, current", max(
        abs(atomic[key]["values"][name]-value) for key, energy in (((0, 0), .5), ((0, 1), 1.5), ((-1, 0), -.015), ((1, 0), .015))
        for name, value in (("energy", energy), ("energy_paper", energy+.0625), ("current", 0.))), 1e-9)
    for label in ("junction", "bridge", "direct"):
        for step in (2e-4, 1e-4):
            a, b, c = rows(f"hf-{label}-{step:g}--1"), rows(f"hf-{label}-{step:g}-+1"), rows(label+"-nrg")
            check(f"{label} phase HF h={step:g}", max(abs((b[s, 0]["values"]["energy"]-a[s, 0]["values"]["energy"])/(2*step)
                  -c[s, 0]["values"]["current"]) for s in PRIMARY), 1e-7)
    check("Gamma=0 direct bridge current nonzero", int(max(abs(r["values"]["current"]) for r in rows("direct-nrg").values()) <= 1e-7), .5)
    for target, flip in (("reverse", True), ("phase-reverse", False)):
        errors = []
        for spin in PRIMARY:
            a, b = rows("bridge-nrg")[spin, 0]["values"], rows(target+"-nrg")[-spin if flip else spin, 0]["values"]
            errors.extend(abs(a[o]-b[o]) for o in ("energy", "energy_paper", "P0", "P1", "P2", "Czz", "Czz_connected", "dot_z2", "bath_z2"))
            errors.extend(abs(a[o]+b[o] if flip else a[o]-b[o]) for o in ("m_dot", "m_bath"))
            errors.append(abs(a["current"]-b["current"] if flip else a["current"]+b["current"]))
        check(target+" symmetry", max(errors))
    return checks


def replay_lambda_control(data):
    control = data["control"]
    nrg = next(r for r in data["records"] if r["id"] == "nrg-l8-control")
    bath = data["baths"][control["bath"]]
    require(bath["chains"] == data["baths"][nrg["bath"]]["chains"], "Lambda8 explicit bath/chain mismatch")
    chains = bath["chains"]
    require(all(c[k] == chains[0][k] for c in chains for k in CHAIN_KEYS), "Lambda8 two leads differ")
    require(all(c["delta_over_gap"] == [-1., -1.] and c["kappa_over_gap"] == [0., 0.]
                and c["n_bonds"] == 1 and c["unused_last_bond_coefficient"] is True for c in chains),
            "Lambda8 control chain convention changed")
    a, d = chains[0]["zeta_over_gap"]
    distance = math.hypot(a-d, 2*chains[0]["xi_over_gap"][0])
    expected = [(a+d-distance)/2, (a+d+distance)/2, (1-(a-d)/distance)/2, (1+(a-d)/distance)/2]
    require(all(math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-12) for x, y in zip(
        bath["xi"]+bath["weights"], expected, strict=True)), "Lambda8 chain-to-star measure changed")
    tables = {backend: {(r["twice_sz"], r["root"]): r for r in control[backend]} for backend in ("ed", "qp")}
    tables["nrg"] = {(r["twice_sz"], r["root"]): r for r in nrg["states"]}
    ground = min(r["values"]["energy"] for r in nrg["states"])
    energies = {}
    for sector in control["spectra"]:
        spin, values = sector["twice_sz"], sector["ground_relative_energy_over_gap"]
        require(len(values) == math.comb(10, 5+spin) and values == sorted(values), "Lambda8 sector multiplicity")
        require(sorted(sector["indices_in_energy_order"]) == list(range(len(values))), "Lambda8 spectrum indices")
        energies.update({(spin, root): ground+value for root, value in enumerate(values)})
    require(len(energies) == 1024 and set(tables["ed"]) == set(tables["qp"]) == set(energies), "full 1024-state inventory required")
    require(set(tables["nrg"]) == {k for k, e in energies.items() if k[1] < 3 or e-ground < 1}, "Lambda8 saved roots incomplete")
    numerical = data["specifications"][nrg["specification"]]["numerical"]
    require(numerical == dict(nmax=1, keep=1024, keepmin=0, keepenergy=-1., untruncated=True, z=.5, **{"lambda": 8.}),
            "Lambda8 control is not the declared full-space calculation")
    require(len(nrg["retention"]["shells"]) == 1 and all(nrg["retention"]["shells"][0][k] == 1024 for k in
        ("kept_reduced_states", "kept_physical_states", "pretruncation_reduced_dimension", "pretruncation_physical_dimension")), "Lambda8 truncated shell")
    groups = []
    for spin in range(-5, 6):
        block = []
        for root in range(math.comb(10, 5+spin)):
            key = spin, root
            if block and energies[key]-energies[block[-1]] > 1e-8:
                groups.append(block)
                block = []
            block.append(key)
        groups.append(block)
    unique = {g[0] for g in groups if len(g) == 1}
    maximum = 0.
    for backend in ("ed", "qp"):
        rows = control[backend]
        checked = deepcopy(rows)
        finish_states(checked)
        info = control["finite_diagnostics"][backend]
        diagnostics = info["diagnostics"]
        require(info["reported_qualified"] is True and diagnostics["all_returned_states_converged"] is True
                and diagnostics["bath_modes_total"] == 8 and diagnostics["bath_leads"] == 2
                and diagnostics["identical_bath_record_for_both_leads"] is True, "full finite control not qualified")
        require(math.isclose(diagnostics["bath_weight_sum_per_lead"], math.fsum(bath["weights"]), rel_tol=1e-12), "control measure changed")
        for call in info["calls"]:
            dim = math.comb(10, 5+call["twice_sz"])
            require(call["root_prefix"] == list(range(dim)) and call["metadata"]["dimension"] == dim
                    and call["metadata"]["complete_sector_diagonalized"] is True, "full finite sector incomplete")
        for r, c in zip(rows, checked, strict=True):
            require(r["degenerate"] == c["degenerate"] and r["excitation"] == c["excitation"], "Lambda8 finite state algebra")
            require(r["finite_converged"] is True and r["residual"] <= 1e-8, "Lambda8 residual check")
        maximum = max(maximum, *(abs(tables[backend][k]["values"]["energy"]-v) for k, v in energies.items()))
        for block in groups:
            require(all(tables[backend][k]["degenerate"] == (len(block) > 1) for k in block), "control multiplet inventories differ")
    for left, right in (("nrg", "ed"), ("nrg", "qp"), ("ed", "qp")):
        for key in tables[left].keys() & tables[right].keys():
            a, b = tables[left][key], tables[right][key]
            maximum = max(maximum, abs(a["excitation"]-b["excitation"]), *(abs(a["values"][o]-b["values"][o])
                for o in (OBSERVABLES if key in unique else ("energy", "energy_paper"))))
    for block in groups:
        if len(block) == 1:
            continue
        centroids = []
        for backend in ("ed", "qp"):
            c = {o: math.fsum(tables[backend][k]["values"][o] for k in block)/len(block) for o in OBSERVABLES}
            c["Czz_connected"] = c["Czz"]-c["m_dot"]*c["m_bath"]
            centroids.append(c)
        maximum = max(maximum, *(abs(centroids[0][o]-centroids[1][o]) for o in OBSERVABLES))
    return dict(states_per_backend=1024, unique_multiplets=len(unique), multiplets=len(groups),
                maximum_error=maximum, passed=maximum <= 1e-8, verification="scalar-replayed")


def scientific_groups(data):
    """Select whole sector calls, retaining both available and qualified evidence."""
    assemblies, groups = {}, []
    for levels in (4, 6, 8, 10, 12):
        bath = f"cosh{levels}"
        candidates = [t for t in data["records"] if t["backend"] == "dmrg" and t["role"] == "study"
                      and data["specifications"][t["specification"]].get("bath") == bath]
        for selection in ("best_available", "highest_qualified"):
            selected = []
            for spin, roots in ((0, 2), (-1, 1), (1, 1)):
                pool = [t for t in candidates if data["specifications"][t["specification"]].get("twice_sz") == spin
                        and len(t["states"]) >= roots]
                if selection == "highest_qualified":
                    pool = [t for t in pool if all(finite_qualified(data, t, r) is True and not r["degenerate"] for r in t["states"])]
                if pool:
                    selected.append(max(pool, key=lambda t: (data["specifications"][t["specification"]]["chi"], len(t["states"]), t["id"])))
            name = f"{bath}-dmrg-{selection}"
            group = dict(id=name, bath=bath, selection=selection, source_tasks=[t["id"] for t in selected],
                         complete=len(selected) == 3, qualified=False, reference_qualified=False,
                         reason="missing qualified sector/root prefix" if len(selected) != 3 else "",
                         scope="four selected states; not spectral completeness")
            if group["complete"]:
                rows = assemble(data, group["source_tasks"], require_qualified=False)
                first = selected[0]
                assembled = first | dict(id=name, states=rows, calls=[c for t in selected for c in t["calls"]],
                    recorded_integrity_ok=all(t["recorded_integrity_ok"] for t in selected), source_tasks=group["source_tasks"])
                assemblies[name] = assembled
                group["qualified"] = not issues(data, assembled)
                group["reference_qualified"] = not reference_issues(data, assembled)
                group["reason"] = "; ".join(issues(data, assembled))
            groups.append(group)
    return assemblies, groups


def point_comparisons(data, kind, old, new, *, matched_bath):
    errors = []
    if old["model"] != new["model"]:
        errors.append("physical inputs differ")
    if matched_bath and not matched_baths(data, old, new):
        errors.append("actual bath coefficients differ")
    if old["backend"] == new["backend"] == "nrg" and family(data, old) != family(data, new):
        errors.append("physical/model/tool/representation family differs")
    rows = []
    for key in CORE:
        if key not in state_map(old) or key not in state_map(new):
            continue
        for observable in QUANTITIES:
            x, y = scalar(old, key, observable), scalar(new, key, observable)
            if x is None or y is None:
                continue
            reasons = errors+state_issues(data, old, key, observable)+state_issues(data, new, key, observable)
            row = compare(kind, old["id"], new["id"], key, observable, x, y, reasons, backend=new["backend"])
            row.update(old_primary_complete=set(CORE[:3]) <= state_map(old).keys(),
                       new_primary_complete=set(CORE[:3]) <= state_map(new).keys(),
                       subgap_scope="available_spectrum" if all("subgap" in state_map(t)[key] for t in (old, new)) else "sector_only")
            rows.append(row)
    return rows


def analyze(data):
    verify_data(data)
    controls = replay_controls(data)
    control = replay_lambda_control(data)
    records = {r["id"]: r for r in data["records"]}
    assemblies, groups = scientific_groups(data)
    means = [row for profile in data["profiles"] for row in average(data, profile)]
    comparisons = []
    before, after = (records[k] for k in ("nrg-l8-z1-k4000", "nrg-l8-z1-k5000"))
    errors = issues(data, before)+issues(data, after)
    if not matched_baths(data, before, after):
        errors.append("both Wilson leads must match")
    if family(data, before) != family(data, after):
        errors.append("endpoint scientific family differs")
    for key in CORE:
        for observable in (*QUANTITIES, *(('odd_splitting',) if key == (1, 1, 0) else ())):
            comparisons.append(compare("endpoint_cap_4000_to_5000", before["id"], after["id"], key, observable,
                scalar(before, key, observable), scalar(after, key, observable), errors, backend="nrg"))
    uniform = {(state_key(r), r["observable"]): r for r in means if r["profile"] == "uniform-l8-k4000-nz8"}
    mixed = [r for r in means if r["profile"] == "mixed-l8-endpoint-nz8"]
    finite = records["cosh12-qp-q6"]
    for r in mixed:
        key, observable = state_key(r), r["observable"]
        old = uniform[key, observable]
        errors = [] if old["qualified"] and r["qualified"] else ["unqualified z mean"]
        comparisons.append(compare("endpoint_mean_impact", old["profile"], r["profile"], key, observable, old["mean"], r["mean"], errors, backend="nrg"))
        comparisons.append(compare("finite_bath_QP_diagnostic", finite["id"], r["profile"], key, observable,
            scalar(finite, key, observable), r["mean"], errors+state_issues(data, finite, key, observable), backend="nrg"))
        comparisons.append(compare("finite_bath_QP_uniform_diagnostic", finite["id"], old["profile"], key, observable,
            scalar(finite, key, observable), old["mean"], ([] if old["qualified"] else ["unqualified z mean"])+state_issues(data, finite, key, observable), backend="nrg"))
    # Explicit essential finite-bath roster, not an all-pairs Cartesian product.
    best_qp = {4: "cosh4-qp-qfull", 6: "cosh6-qp-q10", 8: "cosh8-qp-q8", 10: "cosh10-qp-q6", 12: "cosh12-qp-q6"}
    pairs = [("finite_bath_size", "cosh2-qp-full", best_qp[4], False)]
    for levels, cutoffs in ((4, (2, 4, 6, "full")), (6, (4, 6, 8, 10)), (8, (4, 6, 8)), (10, (4, 6)), (12, (4, 6))):
        pairs.extend(("qp_cutoff", f"cosh{levels}-qp-q{a}", f"cosh{levels}-qp-q{b}", True)
                     for a, b in zip(cutoffs, cutoffs[1:], strict=False))
    pairs.extend(("finite_bath_best", best_qp[a], best_qp[b], False) for a, b in ((4, 6), (6, 8), (8, 10), (10, 12)))
    pairs.append(("finite_bath_historical", best_qp[4], "cosh6-qp-q8", False))
    for levels in best_qp:
        for selection in ("best_available", "highest_qualified"):
            name = f"cosh{levels}-dmrg-{selection}"
            if name in assemblies:
                pairs.append(("dmrg_vs_full_qp" if levels == 4 else "dmrg_vs_qp", best_qp[levels], name, True))
        for spin in (0, -1, 1):
            calls = [t for t in data["records"] if t["backend"] == "dmrg"
                     and data["specifications"][t["specification"]].get("bath") == f"cosh{levels}"
                     and data["specifications"][t["specification"]].get("twice_sz") == spin]
            # Every recorded call is compared with its same-bath QP states;
            # same-sector bond steps retain only corresponding returned roots.
            calls.sort(key=lambda t: (data["specifications"][t["specification"]]["chi"], len(t["states"])))
            pairs.extend(("dmrg_state_vs_qp", best_qp[levels], t["id"], True) for t in calls)
            pairs.extend(("dmrg_bond", a["id"], b["id"], True) for a, b in zip(calls, calls[1:], strict=False)
                         if data["specifications"][a["specification"]]["chi"] != data["specifications"][b["specification"]]["chi"])
    # Adjacent NRG one-axis probes are discovered from exact settings, including
    # the early cap/floor pilots; all other numerical axes remain fixed.
    nrg = [t for t in data["records"] if t["backend"] == "nrg" and t["role"] == "study"]
    for axis in ("keep", "keepmin", "nmax", "keepenergy"):
        series = {}
        for t in nrg:
            n = data["specifications"][t["specification"]]["numerical"]
            key = digest([family(data, t), {k: v for k, v in n.items() if k != axis}])
            series.setdefault(key, []).append((n[axis], t["id"]))
        for members in series.values():
            ordered = sorted(members)
            pairs.extend(("nrg_"+axis, a[1], b[1], axis != "nmax") for a, b in zip(ordered, ordered[1:], strict=False))
    all_records = records | assemblies
    for kind, left, right, matched in pairs:
        comparisons.extend(point_comparisons(data, kind, all_records[left], all_records[right], matched_bath=matched))
    # Nested complete z grids and cross-Lambda means remain separate from
    # retention comparisons. Invalid/missing grids never gain qualified status.
    profiles = {p["id"]: p for p in data["profiles"] if p["id"] != "mixed-l8-endpoint-nz8"}
    mean_index = {(r["profile"], state_key(r), r["observable"]): r for r in means}
    for name, p in profiles.items():
        for other, q in profiles.items():
            fixed = ("lambda", "nmax", "keepenergy", "keepmin")
            kind = None
            if all(p[k] == q[k] for k in fixed) and p["keep_by_z"][0] == q["keep_by_z"][0] and p["nz"]*2 == q["nz"]:
                kind = "nrg_Nz"
            elif all(p[k] == q[k] for k in fixed) and p["nz"] == q["nz"] and p["keep_by_z"][0] == 3000 and q["keep_by_z"][0] == 4000:
                kind = "nrg_mean_cap"
            elif p["lambda"] == 4 and q["lambda"] == 8 and p["nz"] == q["nz"] == 8 and p["keep_by_z"] == q["keep_by_z"]:
                kind = "cross_lambda_discretization_and_finite_bath"
            if kind is None:
                continue
            for key in CORE:
                for observable in QUANTITIES:
                    x, y = (mean_index[n, key, observable] for n in (name, other))
                    if x["mean"] is None or y["mean"] is None:
                        continue
                    errors = [] if x["qualified"] and y["qualified"] else ["unqualified z mean"]
                    if kind == "nrg_mean_cap" and any(not matched_baths(data, records[a], records[b]) for a, b in zip(p["members"], q["members"], strict=True)):
                        errors.append("both Wilson leads must match at every z")
                    comparisons.append(compare(kind, name, other, key, observable, x["mean"], y["mean"], errors, backend="nrg"))
    observations = [dict(task=t["id"], backend=t["backend"], **{k: v for k, v in r.items() if k != "values"},
                         **r["values"], finite_problem_qualified=finite_qualified(data, t, r)) for t in data["records"] for r in t["states"]]
    return dict(observables=observations, comparisons=comparisons, z_averages=means, controls=controls,
                lambda_control=control, assemblies=assemblies, groups=groups)


def report_markdown(data, tables):
    def table(columns, rows):
        def cell(value):
            return "n/a" if value is None else f"{value:.10g}" if isinstance(value, float) else str(value).replace("|", "\\|")
        return ["| "+" | ".join(columns)+" |", "| "+" | ".join("---" for _ in columns)+" |",
                *("| "+" | ".join(cell(v) for v in r)+" |" for r in rows), ""]

    records = {r["id"]: r for r in data["records"]}
    comparisons = tables["comparisons"]
    text = ["# Test4 Working Results", "", "Coarse working evidence, not a continuum reference or error bar. "
        "Production approval, continuum accuracy and spectral completeness remain **false**.", "",
        "The separate mixed-retention Lambda8 Nz8 mean uses seven keep4000 points and keep5000 only at z=0.125. "
        "Both profiles use Nmax=16, keepenergy=16 and keepmin=400, with identical twists in the two physical leads. "
        "The original uniform keep4000 mean is unchanged. Finite cosh12 QP q6 is a bath- and cutoff-unconverged comparator, not truth.", "",
        "Energies/excitations are in gap units; current is $I/(2e\\Delta/\\hbar)$. Gamma=0.4 is total hybridization, "
        "D=100 is half-bandwidth, U=2, nu=0.75, B=0.03, phi=3*pi/5, WN=WS=0. The paper energy offset is 0.0625.", "",
        "Requested comparison tolerances remain 1e-9 + 1e-6*abs(new) for QP/DMRG and "
        "1e-6 + 0.001*abs(new) for NRG. The reported discrepancies and 1%/5% diagnostic tiers "
        "are not independently established continuum error bars.", ""]
    uniform = {(state_key(r), r["observable"]): r for r in tables["z_averages"] if r["profile"] == "uniform-l8-k4000-nz8"}
    text += table(["State (parity, twice Sz, root)", "Quantity", "Uniform NRG", "Mixed NRG", "QP cosh12 q6", "Discrepancy %"], (
        [state_key(r), r["observable"], uniform[state_key(r), r["observable"]]["mean"], r["new"], r["old"], r["relative_percent"]]
        for r in comparisons if r["kind"] == "finite_bath_QP_diagnostic"))
    text += ["Discrepancy = 100*abs(NRG-QP)/abs(NRG); near zero (<=1e-6), only absolute differences are meaningful. "
             "These are empirical deltas, not statistical errors. Absolute-energy percentages depend on the stated energy zero. "
             "Each excitation and connected Czz is formed at each z before averaging.", "",
             "## Supplemental material", "", "### Controls", "",
             "The 37 original controls retain all returned roots, including above-gap and degenerate rows. "
             "Atomic, unit, representation, phase-Hellmann-Feynman and reversal checks below are scalar-replayed. "
             "Solver norm/Gram/sweep/root checks are recorded evidence, not newly checked raw MPS/HDF5.", ""]
    text += table(["Check", "Maximum Error", "Tolerance", "Passed"], (
        [r["name"], r["maximum_error"], r["tolerance"], r["passed"]] for r in tables["controls"]))
    c = tables["lambda_control"]
    text += [f"Lambda8 untruncated control: {c['states_per_backend']} energies per backend (NRG/ED/full QP); "
             f"{c['multiplets']} multiplets, {c['unique_multiplets']} unique. Maximum scalar-replayed error: {c['maximum_error']:.10g}. "
             "Degenerate individual expectations are excluded; ED/QP multiplet means are compared. Both actual Wilson leads and unnormalized star weights are retained.", "",
             "### Convergence", "", "Requested targets use the new backend: QP/DMRG 1e-9 + 1e-6*abs(new), "
             "NRG 1e-6 + 0.001*abs(new). ED controls retain their separate absolute tolerance. "
             "The additional 1%/5% diagnostics use 1e-6 + relative*abs(new); they do not replace a finite-backend target. "
             "Passing a step is not a convergence certificate.", ""]
    for kind in ("endpoint_cap_4000_to_5000", "endpoint_mean_impact", "finite_bath_QP_diagnostic"):
        rows = [r for r in comparisons if r["kind"] == kind and r["observable"] != "odd_splitting"]
        text += [f"- {kind}: "+" / ".join(str(sum(r[t] is True for r in rows)) for t in TIERS)+f" of {len(rows)} scalars (requested / 1% / 5%)."]
    text += ["", "QP residuals qualify only the projected finite problem, not cutoff or bath convergence. "
             "DMRG qualification requires recorded finite-converged, norm, Gram, sweep and same-call root checks. "
             "Single-sector calls have no global excitations; compatible assemblies must share physics, bath, full Hamiltonian fingerprint and coordinates, with no even-root splicing.", "",
             ""]
    text += table(["Convergence Diagnostic", "Scalar Pairs", "Qualified", "Requested", "1%", "5%"], (
        [kind, len(rows), sum(r["qualified"] for r in rows), *(sum(r[tier] is True for r in rows) for tier in TIERS)]
        for kind in sorted({r["kind"] for r in comparisons if r["kind"].startswith(("nrg_", "qp_", "dmrg_", "finite_bath_size", "cross_lambda"))})
        if (rows := [r for r in comparisons if r["kind"] == kind])))
    text += ["Finite evidence below separates raw observed target agreement from qualified agreement. "
             "Ordinary sector observables use that state's checks, not another sector's residual. "
             "Excitations additionally require a qualified lowest AVAILABLE ground reference. "
             "Above-gap or degenerate comparison states are excluded individually; partial coverage remains explicit.", ""]
    text += table(["Finite Comparison", "Old", "New", "Scalars", "Qualified", "Observed Requested", "Qualified Requested"], (
        [kind, old, new, len(rows), sum(r["qualified"] for r in rows), sum(r["observed_requested"] is True for r in rows),
         sum(r["requested"] is True for r in rows)]
        for kind, old, new in sorted({(r["kind"], r["old_id"], r["new_id"]) for r in comparisons
            if r["kind"] in ("qp_cutoff", "finite_bath_size", "finite_bath_best", "finite_bath_historical", "dmrg_vs_full_qp", "dmrg_vs_qp")})
        if (rows := [r for r in comparisons if (r["kind"], r["old_id"], r["new_id"]) == (kind, old, new)])))
    text += table(["Bath", "DMRG Selection", "Complete Four States", "Qualified", "Ground Qualified", "Whole Source Calls"], (
        [g["bath"], g["selection"], g["complete"], g["qualified"], g["reference_qualified"], ", ".join(g["source_tasks"])]
        for g in tables["groups"]))
    text += ["### Negative Findings", ""]
    text += table(["Task", "Signed Odd Splitting", "Exact WS=0 Limit", "Qualified"], (
        [t["id"], b["splitting"], b["limit"], False] for t in data["records"] if (b := field_bound(data, t))["passed"] is False))
    text += ["Historical truncation-invalid points remain invalid even when their disposition allowed later work. "
             "Incomplete z grids are not renormalized; complete but invalid grids retain descriptive means marked unqualified. "
             "No NRG residual is invented.", "",
             f"Preserved: {len(records)} completed scientific records, {len(tables['observables'])} original root rows, "
             f"plus 2048 Lambda8 ED/QP control roots; {len(data['outcomes'])} compact historical noncompletion outcomes.", ""]
    text += table(["Historical Stage", "Outcome", "Count"], ([stage, status, count] for (stage, status), count in sorted(
        Counter((r["stage"], r["status"]) for r in data["outcomes"]).items())))
    text += ["Private-original checksums are reference labels only. Public content-addressed tables and the dataset checksum "
             "describe the new payload. Runtime versions and canonical relative source identities are recorded provenance, not executable rehashes. "
             "Measured work-duration sums are not parallel wall time or CPU time; RSS/VM are measured peaks.", ""]
    text += table(["Backend", "Completed Records", "Measured Work Seconds", "Peak RSS GiB", "Peak VM GiB"], (
        [backend, len(group), math.fsum(t["cost"].get("monitored_seconds", 0.) for t in group),
         max((t["cost"].get("peak_rss_gib", 0.) for t in group), default=None),
         max((t["cost"].get("peak_vm_gib", 0.) for t in group), default=None)]
        for backend in ("nrg", "qp", "dmrg", "ed") if (group := [t for t in data["records"] if t["backend"] == backend])))
    return "\n".join(text)


def write_presentation(data, output_dir, *, overwrite=False, publish=False):
    """Explicit internal publication is the only way to replace shipped tables."""
    output = Path(output_dir)
    require(not any(p.is_symlink() for p in (output, *output.parents)), "output symlink")
    output = output.resolve()
    require(publish or not output.is_relative_to(HERE / "output"), "shipped output requires internal publish")
    require(not any((output / name).is_symlink() for name in ARTIFACTS), "artifact symlink")
    if not overwrite and any((output / name).exists() for name in ARTIFACTS):
        raise FileExistsError("presentation exists; use --overwrite")
    tables = analyze(data)
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    (output / "report.md").write_text(report_markdown(data, tables), encoding="utf-8", newline="\n")
    for name in ("observables", "comparisons", "z_averages"):
        rows = tables[name]
        with (output / (name+".csv")).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=sorted({k for r in rows for k in r}))
            writer.writeheader()
            writer.writerows({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (list, dict)) else v for k, v in r.items()} for r in rows)
    return tables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables-only", action="store_true", help="stdlib-only tables (also the default)")
    parser.add_argument("--output-dir", type=Path, default=HERE / "runs/presentation")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    try:
        data = json.loads((HERE / "output/results.json").read_bytes())
        write_presentation(data, args.output_dir, overwrite=args.overwrite)
    except (ValueError, KeyError, TypeError, OSError) as error:
        parser.error(str(error))
    print("Scalar snapshot validated; tables generated. No numerical calculations performed.")


if __name__ == "__main__":
    main()
