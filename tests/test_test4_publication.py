"""Offline public-snapshot regressions, independent of frozen controllers/tests."""

from copy import deepcopy
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


CASE = Path(__file__).resolve().parents[1] / "NRG_comparisons/test4"
SPEC = importlib.util.spec_from_file_location("NRG_comparisons.test4.generate", CASE / "generate.py")
publication = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publication)


@pytest.fixture(scope="module")
def snapshot():
    return json.loads((CASE / "output/results.json").read_bytes())


@pytest.fixture(scope="module")
def tables(snapshot):
    return publication.analyze(snapshot)


def record(data, name):
    return next(r for r in data["records"] if r["id"] == name)


def profile(data, name="mixed-l8-endpoint-nz8"):
    return next(r for r in data["profiles"] if r["id"] == name)


def seal(data):
    data["checksum"] = publication.digest({k: v for k, v in data.items() if k != "checksum"})


def replace_reference(data, task, table, field, value):
    key = publication.digest(value)
    data[table][key] = value
    task[field] = key
    seal(data)


def test_preserved_union_and_false_claims(snapshot, tables):
    assert len(snapshot["records"]) == 157
    assert len(tables["observables"]) == 2616
    assert sum(len(snapshot["control"][b]) for b in ("ed", "qp")) == 2048
    assert sum(r["role"] == "control" for r in snapshot["records"]) == 38
    assert snapshot["claims"] == publication.CLAIMS
    assert all(v is False for v in snapshot["claims"].values())
    assert (CASE / "output/results.json").stat().st_size < 8*1024**2
    earlier = {"nrg-k200", "nrg-k500", "nrg-k1000", "nrg-length22", "cosh2-qp-full",
               "cosh4-qp-q2", "cosh4-qp-q4", "cosh4-qp-q6", "nrg-k1500-m300", "nrg-k2000-m300",
               "nrg-k2000-m400", "nrg-k3000-m400", "nrg-k4000-m400", "nrg-length22-k4000"}
    earlier |= {f"cosh4-dmrg-{s}-chi{chi}" for s in ("even", "odd_down", "odd_up") for chi in (64, 128)}
    assert len(earlier) == 20 and earlier <= {r["id"] for r in snapshot["records"]}
    assert {r["status"] for r in snapshot["outcomes"]} == {"failed", "unstarted", "deferred", "not_needed", "incomplete"}
    assert all(set(row["values"]) == set(publication.OBSERVABLES) for t in snapshot["records"] for row in t["states"])
    assert any(not r.get("subgap", True) for t in snapshot["records"] for r in t["states"])


def test_private_baseline_final_golden(tables):
    # Small independent fixtures captured from the private final assessment.
    expected = {"endpoint_cap_4000_to_5000": [1, 18, 38], "endpoint_mean_impact": [15, 53, 60],
                "finite_bath_QP_diagnostic": [18, 53, 60]}
    for kind, counts in expected.items():
        rows = [r for r in tables["comparisons"] if r["kind"] == kind and r["observable"] != "odd_splitting"]
        assert len(rows) == 60 and all(r["qualified"] for r in rows)
        assert [sum(r[tier] is True for r in rows) for tier in publication.TIERS] == counts
    for name in ("uniform-l8-k4000-nz8", "mixed-l8-endpoint-nz8"):
        assert len([r for r in tables["z_averages"] if r["profile"] == name]) == 61
    mixed = {(publication.state_key(r), r["observable"]): r["mean"] for r in tables["z_averages"]
             if r["profile"] == "mixed-l8-endpoint-nz8"}
    assert [mixed[k, "excitation"] for k in publication.CORE] == pytest.approx(
        [0., .0250869197, .158071090, .765661770], abs=5e-10)
    assert [mixed[k, "current"] for k in publication.CORE] == pytest.approx(
        [-.009417139, -.009716108, .059654293, -.119541100], abs=5e-10)
    assert publication.digest(sorted((r["profile"], publication.state_key(r), r["observable"], r["mean"], r["per_z_values"])
        for r in tables["z_averages"] if r["profile"] in ("uniform-l8-k4000-nz8", "mixed-l8-endpoint-nz8"))) == (
        "d3a6e6ce7248d6fc0807f769a6efaaea3d0c43c6cb53c882ac95ec5ccfc6fb6a")


def test_all_original_scientific_rows_golden(snapshot):
    # Exact floats, residuals, degeneracy, exclusions and absence of single-call
    # global excitations, independently compared by the private exporter.
    assert publication.digest(sorted((t["id"], t["states"]) for t in snapshot["records"])) == (
        "877c7fc87d2ebcbb3b373f5ee410f4073a852f525ab5772d66889b6eb0544084")
    assert publication.digest(sorted((t["id"], snapshot["models"][t["model"]], snapshot["baths"][t["bath"]],
        snapshot["specifications"][t["specification"]], [snapshot["solvers"][c["solver"]] for c in t["calls"]])
        for t in snapshot["records"])) == "bcbd2ed057e7f7c6e000903dd636904ad08f39b9f99cefa4857c2ce7fa96336e"


def test_scalar_control_replay(tables):
    assert len(tables["controls"]) == 94
    assert all(r["passed"] for r in tables["controls"])
    assert tables["lambda_control"] == dict(states_per_backend=1024, unique_multiplets=256,
        multiplets=576, maximum_error=4.206579529153487e-12, passed=True, verification="scalar-replayed")


def test_negative_evidence_and_incomplete_grids(snapshot, tables):
    for name in ("nrg-k200", "nrg-k500", "nrg-z1-k3000"):
        t = record(snapshot, name)
        assert publication.field_bound(snapshot, t)["passed"] is False
        assert "exact WS=0 field bound violated" in publication.issues(snapshot, t)
    assert abs(publication.field_bound(snapshot, record(snapshot, "nrg-z1-k3000"))["splitting"]) > .031688
    mean = [r for r in tables["z_averages"] if r["profile"] == "lambda4-n22-ke16-floor400-k3000-nz8"]
    assert len(mean) == 61 and all(r["mean"] is not None and not r["qualified"] for r in mean)
    assert any(not r["grid_complete"] and r["mean"] is None for r in tables["z_averages"])
    assert publication.field_bound(snapshot, record(snapshot, "spin-nrg"))["applicable"] is False


@pytest.mark.parametrize("lead", [0, 1])
@pytest.mark.parametrize("coefficient", publication.CHAIN_KEYS)
def test_either_wilson_lead_invalidates_matched_comparison(snapshot, lead, coefficient):
    data = deepcopy(snapshot)
    t = record(data, "nrg-l8-z1-k5000")
    bath = deepcopy(data["baths"][t["bath"]])
    bath["chains"][lead][coefficient][0] += .001
    replace_reference(data, t, "baths", "bath", bath)
    assert not publication.matched_baths(data, record(data, "nrg-l8-z1-k4000"), t)
    result = publication.analyze(data)
    rows = [r for r in result["comparisons"] if r["kind"] == "endpoint_cap_4000_to_5000"]
    assert rows and all(not r["qualified"] for r in rows)
    assert all(not r["qualified"] for r in result["z_averages"] if r["profile"] == "mixed-l8-endpoint-nz8")


def test_checksum_not_private_original_rehash(snapshot):
    data = deepcopy(snapshot)
    data["records"][0]["states"][0]["residual"] += 1e-4
    with pytest.raises(ValueError, match="public checksum"):
        publication.verify_data(data)
    assert all(r["private_original"]["label"] == "PRIVATE_ORIGINAL_REFERENCE" for r in snapshot["records"])
    assert snapshot["checksum"] not in {r["private_original"]["result_sha256"] for r in snapshot["records"]}


@pytest.mark.parametrize("key", ["user", "uid", "hostname", "account", "job_id", "command", "release", "parent_campaign", "authorization", "arbitrary_new_key"])
def test_metadata_allowlist_even_after_reseal(snapshot, key):
    data = deepcopy(snapshot)
    data["records"][0]["diagnostics"][key] = "private"
    seal(data)
    with pytest.raises(ValueError, match="metadata"):
        publication.verify_data(data)


@pytest.mark.parametrize("claim", publication.CLAIMS)
def test_formal_claims_cannot_be_promoted(snapshot, claim):
    data = deepcopy(snapshot)
    data["claims"][claim] = True
    seal(data)
    with pytest.raises(ValueError, match="claim"):
        publication.verify_data(data)


def test_root_prefix_and_single_sector_scope(snapshot):
    data = deepcopy(snapshot)
    t = record(data, "cosh4-dmrg-even-pair-chi256")
    t["states"][1]["root"] = 2
    seal(data)
    with pytest.raises(ValueError, match="prefix"):
        publication.verify_data(data)
    data = deepcopy(snapshot)
    t = record(data, "cosh4-dmrg-even-pair-chi256")
    t["states"][0].update(excitation=0., subgap=True)
    seal(data)
    with pytest.raises(ValueError, match="single-sector"):
        publication.verify_data(data)


def test_degenerate_and_missing_z_means_unqualified(snapshot):
    data = deepcopy(snapshot)
    t = record(data, "nrg-l8-z1-k5000")
    next(r for r in t["states"] if publication.state_key(r) == (1, -1, 0))["degenerate"] = True
    mean = publication.average(data, profile(data))
    assert all(not r["qualified"] and r["mean"] is not None for r in mean)
    assert any("ground" in r["reason"] for r in mean)
    data = deepcopy(snapshot)
    data["records"] = [r for r in data["records"] if r["id"] != profile(data)["members"][1]]
    assert all(not r["qualified"] and r["mean"] is None for r in publication.average(data, profile(data)))
    data = deepcopy(snapshot)
    profile(data)["members"][1] = "nrg-z2-k4000"
    with pytest.raises(ValueError, match="endpoint replacement"):
        publication.average(data, profile(data))


@pytest.mark.parametrize("levels", [6, 10, 12])
def test_bad_qp_cutoff_step_keeps_original_finite_target(tables, levels):
    rows = [r for r in tables["comparisons"] if r["kind"] == "qp_cutoff"
            and r["old_id"] == f"cosh{levels}-qp-q4" and r["new_id"] == f"cosh{levels}-qp-q6"]
    assert len(rows) == 60 and all(r["qualified"] for r in rows)
    assert sum(r["requested"] is True for r in rows) == 1
    assert sum(r["observed_requested"] is True for r in rows) == 1
    assert all((r["target_backend"], r["target_absolute"], r["target_relative"]) == ("qp", 1e-9, 1e-6) for r in rows)
    assert sum(r["one_percent"] is True for r in rows) == 60
    assert all(r["target_tolerance"] == 1e-9+1e-6*abs(r["new"]) for r in rows)


def test_essential_bath_and_cross_backend_roster(tables):
    expected = [("cosh4-qp-qfull", "cosh6-qp-q10"), ("cosh6-qp-q10", "cosh8-qp-q8"),
                ("cosh8-qp-q8", "cosh10-qp-q6"), ("cosh10-qp-q6", "cosh12-qp-q6")]
    for before, after in expected:
        rows = [r for r in tables["comparisons"] if r["kind"] == "finite_bath_best" and (r["old_id"], r["new_id"]) == (before, after)]
        assert len(rows) == 60 and all(r["qualified"] for r in rows)
        assert sum(r["requested"] is True for r in rows) == 1
    for selection in ("best_available", "highest_qualified"):
        rows = [r for r in tables["comparisons"] if r["kind"] == "dmrg_vs_full_qp" and r["new_id"] == f"cosh4-dmrg-{selection}"]
        assert len(rows) == 60 and all(r["qualified"] and r["requested"] for r in rows)
        group = next(g for g in tables["groups"] if g["id"] == f"cosh4-dmrg-{selection}")
        assert group["source_tasks"] == ["cosh4-dmrg-even-pair-chi256", "cosh4-dmrg-odd_down-chi256", "cosh4-dmrg-odd_up-chi256"]
        assert group["complete"] and group["qualified"] and group["reference_qualified"]
    assert len(tables["groups"]) == 10
    for levels in (6, 8, 10, 12):
        group = next(g for g in tables["groups"] if g["id"] == f"cosh{levels}-dmrg-highest_qualified")
        assert not group["complete"] and not group["qualified"]
    for levels in (10, 12):
        rows = [r for r in tables["comparisons"] if r["kind"] == "dmrg_vs_qp" and r["new_id"] == f"cosh{levels}-dmrg-best_available"]
        assert len(rows) == 60 and all(r["observed_requested"] and not r["qualified"] for r in rows)


def test_comparison_roster_is_bounded(tables):
    assert len(tables["comparisons"]) == 6162
    assert len({(r["kind"], r["old_id"], r["new_id"]) for r in tables["comparisons"]}) == 156
    assert len(tables["assemblies"]) == 6
    assert all(publication.state_key(r) in publication.CORE for r in tables["comparisons"])
    assert (CASE / "output/comparisons.csv").stat().st_size < 2*1024**2


def test_above_gap_extra_does_not_poison_primary_comparisons(tables):
    rows = [r for r in tables["comparisons"] if r["kind"] == "finite_bath_size" and r["old_id"] == "cosh2-qp-full"]
    assert len(rows) == 60
    assert sum(r["qualified"] for r in rows) == 45
    assert all(r["qualified"] == (publication.state_key(r) != (0, 0, 1)) for r in rows)
    extra = next(r for r in rows if publication.state_key(r) == (0, 0, 1) and r["observable"] == "excitation")
    assert extra["old"] == pytest.approx(1.13077, abs=1e-5)


def test_sector_and_excitation_qualification_are_distinct(snapshot, tables):
    rows = [r for r in tables["comparisons"] if r["kind"] == "dmrg_state_vs_qp" and r["new_id"] == "cosh6-dmrg-odd_down-chi1024"]
    assert len(rows) == 14 and all(r["qualified"] and r["requested"] for r in rows)
    assert not any(r["observable"] == "excitation" for r in rows)
    assert all(r["old_primary_complete"] and not r["new_primary_complete"] and r["subgap_scope"] == "sector_only" for r in rows)
    assert any(r["kind"] == "dmrg_bond" and r["observed_requested"] is True and not r["qualified"] for r in tables["comparisons"])
    data = deepcopy(snapshot)
    t = record(data, "cosh12-qp-q6")
    ground = min(t["states"], key=lambda r: r["values"]["energy"])
    ground["finite_converged"] = False
    assert publication.state_issues(data, t, (0, 0, 0), "current") == []
    assert "unqualified available ground reference" in publication.state_issues(data, t, (0, 0, 0), "excitation")
    t["states"] = [r for r in t["states"] if r["twice_sz"] == 0]
    assert publication.state_issues(data, t, (0, 0, 0), "energy") == []
    assert "incomplete primary excitation reference" in publication.state_issues(data, t, (0, 0, 0), "excitation")
    data = deepcopy(snapshot)
    t = record(data, "nrg-l8-z1-k5000")
    next(r for r in t["states"] if publication.state_key(r) == (0, 0, 1))["degenerate"] = True
    assert publication.state_issues(data, t, (0, 0, 0), "current") == []
    assert all(not r["qualified"] for r in publication.average(data, profile(data)))


@pytest.mark.parametrize("backend", ["qp", "dmrg"])
def test_resealed_valid_solver_registry_swap_rejected(snapshot, backend):
    data = deepcopy(snapshot)
    if backend == "qp":
        t = record(data, "cosh12-qp-q6")
        left, right = (next(c for c in t["calls"] if c["twice_sz"] == s) for s in (-1, 1))
    else:
        left = record(data, "cosh4-dmrg-odd_down-chi256")["calls"][0]
        right = record(data, "cosh4-dmrg-odd_up-chi256")["calls"][0]
    left["solver"], right["solver"] = right["solver"], left["solver"]
    seal(data)
    with pytest.raises(ValueError, match="solver sector"):
        publication.verify_data(data)


@pytest.mark.parametrize("field", ["coordinate_id", "mode_order", "group_size", "hamiltonian_sha256", "sector", "options"])
def test_resealed_diagnostic_solver_incoherence_rejected(snapshot, field):
    data = deepcopy(snapshot)
    t = record(data, "cosh4-dmrg-even-pair-chi256")
    if field == "mode_order":
        t["diagnostics"][field] = list(reversed(t["diagnostics"][field]))
    elif field == "group_size":
        t["diagnostics"][field] = 1
    elif field == "sector":
        t["diagnostics"][field]["twice_sz"] = -1
    elif field == "options":
        t["diagnostics"][field]["seed"] += 1
    else:
        t["diagnostics"][field] = "f"*64
    seal(data)
    with pytest.raises(ValueError, match="mismatch"):
        publication.verify_data(data)


@pytest.mark.parametrize("field", ["groups", "mode_order", "coordinate_id"])
def test_rehashed_solver_incoherence_rejected(snapshot, field):
    data = deepcopy(snapshot)
    t = record(data, "cosh4-dmrg-even-pair-chi256")
    call = t["calls"][0]
    m = deepcopy(data["solvers"][call["solver"]])
    m[field] = "f"*64 if field == "coordinate_id" else list(reversed(m[field]))
    replace_reference(data, call, "solvers", "solver", m)
    with pytest.raises(ValueError, match="mismatch"):
        publication.verify_data(data)


def test_uniform_profile_cannot_hide_a_mixed_cap(snapshot):
    data = deepcopy(snapshot)
    p = profile(data, "uniform-l8-k4000-nz8")
    p["members"][0] = "nrg-l8-z1-k5000"
    p["keep_by_z"][0] = 5000
    seal(data)
    with pytest.raises(ValueError, match="uniform profile"):
        publication.verify_data(data)
    with pytest.raises(ValueError, match="uniform profile"):
        publication.average(data, p)
    p["id"] = "renamed-profile"
    with pytest.raises(ValueError, match="uniform profile"):
        publication.average(data, p)
    p["kind"] = "mixed_endpoint"
    with pytest.raises(ValueError, match="uniform reference"):
        publication.average(data, p)
    assert all(r["qualified"] for r in publication.average(snapshot, profile(snapshot)))


def test_unqualified_finite_ground_and_no_invented_nrg_residual(snapshot):
    data = deepcopy(snapshot)
    t = record(data, "cosh12-qp-q6")
    ground = min(t["states"], key=lambda r: r["values"]["energy"])
    ground["finite_converged"] = False
    assert "unqualified available ground reference" in publication.issues(data, t)
    data = deepcopy(snapshot)
    record(data, "nrg-l8-z1-k5000")["states"][0]["residual"] = 0.
    seal(data)
    with pytest.raises(ValueError, match="invented NRG residual"):
        publication.verify_data(data)


def test_connected_observable_averaged_before_mean(snapshot):
    rows = publication.average(snapshot, profile(snapshot))
    key = (0, 0, 0)
    means = {r["observable"]: r["mean"] for r in rows if publication.state_key(r) == key}
    assert means["Czz_connected"] != means["Czz"]-means["m_dot"]*means["m_bath"]
    values = [publication.scalar(record(snapshot, name), key, "Czz_connected") for name in profile(snapshot)["members"]]
    assert means["Czz_connected"] == math.fsum(values)/8


def test_dmrg_recorded_checks_and_assembly(snapshot):
    names = ["cosh4-dmrg-even-pair-chi256", "cosh4-dmrg-odd_down-chi256", "cosh4-dmrg-odd_up-chi256"]
    rows = publication.assemble(snapshot, names)
    assert len(rows) == 4 and all("excitation" in r for r in rows)
    with pytest.raises(ValueError, match="splice"):
        publication.assemble(snapshot, names+["cosh4-dmrg-even-chi256"])
    for field in ("hamiltonian_sha256", "coordinate_id", "mode_order", "group_size"):
        data = deepcopy(snapshot)
        record(data, names[0])["diagnostics"][field] = "different"
        with pytest.raises(ValueError, match="different assembly"):
            publication.assemble(data, names)
    for field in ("finite_problem_converged", "gram_error"):
        data = deepcopy(snapshot)
        t = record(data, names[0])
        m = data["solvers"][t["calls"][0]["solver"]]
        m[field] = False if field == "finite_problem_converged" else .1
        assert not publication.finite_qualified(data, t, t["states"][0])
    for field, value in (("norm_error", .1), ("lower_root_overlap", .1), ("sweep_converged", False), ("sweep_energy_converged", False)):
        data = deepcopy(snapshot)
        t = record(data, names[0])
        data["solvers"][t["calls"][0]["solver"]]["roots"][0][field] = value
        assert not publication.finite_qualified(data, t, t["states"][0])
    data = deepcopy(snapshot)
    t = record(data, names[0])
    data["solvers"][t["calls"][0]["solver"]]["gram"]["real"][0][1] = .1
    assert not publication.finite_qualified(data, t, t["states"][0])


def test_full_control_scalar_mutation_not_saved_passed_boolean(snapshot):
    data = deepcopy(snapshot)
    data["control"]["ed"][0]["values"]["current"] += .01
    assert not publication.replay_lambda_control(data)["passed"]
    data = deepcopy(snapshot)
    data["control"]["spectra"][0]["ground_relative_energy_over_gap"].pop()
    with pytest.raises(ValueError, match="multiplicity"):
        publication.replay_lambda_control(data)
    data = deepcopy(snapshot)
    for r in record(data, "hf-direct-0.0001-+1")["states"]:
        r["values"]["energy"] += .001
    assert any(not r["passed"] for r in publication.replay_controls(data))


def test_artifacts_deterministic_and_relocatable(snapshot, tmp_path):
    relocated = tmp_path / "case"
    (relocated / "output").mkdir(parents=True)
    for filename in ("generate.py", "physics.py", "output/results.json"):
        shutil.copyfile(CASE / filename, relocated / filename)
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    output = tmp_path / "presentation"
    command = [sys.executable, "-B", "-S", str(relocated / "generate.py"), "--tables-only", "--output-dir", str(output)]
    result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    for filename in publication.ARTIFACTS:
        assert (output / filename).read_bytes() == (CASE / "output" / filename).read_bytes()
    assert json.loads((output / "results.json").read_bytes()) == snapshot
    assert subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, timeout=60).returncode != 0
    assert subprocess.run(command+["--overwrite"], cwd=tmp_path, env=env, capture_output=True, timeout=60).returncode == 0
    protected = command[:-1]+[str(relocated / "output"), "--overwrite"]
    assert subprocess.run(protected, cwd=tmp_path, env=env, capture_output=True, timeout=60).returncode != 0
    default = subprocess.run([sys.executable, "-B", "-S", str(relocated / "generate.py"), "--tables-only"],
                             cwd=tmp_path, env=env, capture_output=True, timeout=60)
    assert default.returncode == 0
    assert (relocated / "runs/presentation/report.md").is_file()
    main, supplemental = (output / "report.md").read_text().split("## Supplemental material\n", 1)
    assert all(text in main for text in ("Nmax=16", "keepenergy=16", "keepmin=400", "1e-9 + 1e-6*abs(new)",
                                        "1e-6 + 0.001*abs(new)", "Uniform NRG", "Mixed NRG"))
    assert "### Controls" in supplemental and "### Negative Findings" in supplemental


def test_publication_bytes_with_windows_default_newlines(snapshot, tmp_path, monkeypatch):
    original_open = Path.open

    def windows_open(path, mode="r", *args, **kwargs):
        if mode == "w" and kwargs.get("newline") is None:
            kwargs["newline"] = "\r\n"
        return original_open(path, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", windows_open)
    publication.write_presentation(snapshot, tmp_path)
    for filename in publication.ARTIFACTS:
        generated = (tmp_path / filename).read_bytes()
        assert generated == (CASE / "output" / filename).read_bytes()
        if filename.endswith((".json", ".md")):
            assert b"\r\n" not in generated and generated.endswith(b"\n")
        else:
            assert generated.count(b"\n") == generated.count(b"\r\n")


def test_output_symlink_protection(snapshot, tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    link = tmp_path / "link"
    link.symlink_to(target, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        publication.write_presentation(snapshot, link)
