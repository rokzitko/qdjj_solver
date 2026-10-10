"""Offline regression checks for the detuned test2 presentation, without solvers."""

import csv
import hashlib
import json
import math
from pathlib import Path
import re
import runpy
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest


CASE = Path(__file__).resolve().parents[1] / "NRG_comparisons/test2"
OUTPUT = CASE / "output"
TABLES = ("results.json", "observables.csv", "comparisons.csv")
FIGURES = ("comparison.svg", "convergence.svg")


@pytest.fixture(scope="module")
def results():
    data = json.loads((OUTPUT / "results.json").read_text())
    json.dumps(data, allow_nan=False)
    return data


@pytest.fixture(scope="module")
def generator():
    return runpy.run_path(str(CASE / "generate.py"))


def test_detuned_scientific_contract_and_identities(results, generator):
    generator["verify_data"](results)
    assert len(results["finite_results"]) == 14 and len(results["nrg_results"]) == 7
    assert len(results["observables"]) == 10 and "singlet.moment" not in results["observables"]
    physical = json.loads((CASE / "input/physical.json").read_text())
    assert physical["physical"] == results["model"]
    assert physical["targets"] == results["targets"]
    assert results["model"]["detuning"] == .5 and results["conventions"]["compress"] is False
    for name, row in (results["finite_results"] | results["nrg_results"]).items():
        values = row["values"]
        # Independent identities allow only roundoff, not requested-target errors.
        assert values["signed_gap"] == pytest.approx(values["doublet.energy"] - values["singlet.energy"], abs=2e-15, rel=0)
        for branch in ("singlet", "doublet"):
            p0, p1, p2 = (values[f"{branch}.P{i}"] for i in range(3))
            assert 0 <= p2 < p0 <= 1 and 0 <= p1 <= 1
            assert p0 + p1 + p2 == pytest.approx(1, abs=1e-11, rel=0)
        assert 0 < values["doublet.moment"] <= values["doublet.P1"] / 2
        assert abs(results["symmetry_diagnostics"][name]["singlet.moment"]) < 1e-10
        assert row["continuum_target_established"] is False
    assert results["controls"]["passed"] is True
    assert max(results["controls"]["discrepancies"].values()) < results["controls"]["tolerance"] == 1e-8


def test_residual_scopes_checked_roots_and_mixed_final_bonds(results):
    finite = results["finite_results"]
    for name, row in finite.items():
        qp = row["settings"]["backend"] == "qp"
        assert row["residual_scope"] == ("projected_cutoff_hamiltonian" if qp else "untruncated_finite_hamiltonian")
        assert row["residual_tolerance"] == (1e-8 if qp else 1e-7)
        count = 2 if name in {"cosh12-q6", "cosh16-q6"} else 1
        for branch in ("singlet", "doublet"):
            residuals = row["residuals"][branch]
            assert len(residuals) == row["checked_root_counts"][branch] == count
            assert row["root_separation_checked"][branch] is (count > 1)
            trials = 3 if name == "cosh16-dmrg-chi256" or (name == "cosh16-dmrg-chi512" and branch == "singlet") else 1
            assert row["seed_trial_counts"][branch] == trials
            assert row["branch_finite_problem_qualified"][branch] is all(r <= row["residual_tolerance"] for r in residuals)
        assert row["finite_problem_qualified"] is (qp or name == "cosh16-dmrg-final")
    final, previous = finite["cosh16-dmrg-final"], finite["cosh16-dmrg-chi768"]
    assert final["settings"]["bond_dimensions"] == dict(singlet=1024, doublet=768)
    assert final["branch_finite_problem_qualified"] == dict(singlet=True, doublet=True)
    assert previous["branch_finite_problem_qualified"] == dict(singlet=False, doublet=True)
    assert final["residuals"]["doublet"] == previous["residuals"]["doublet"]
    for name in results["observables"]:
        if name.startswith("doublet."):
            assert final["values"][name] == previous["values"][name]


def test_nrg_settings_are_finite_flow_not_empirical_uncertainty(results):
    assert results["nrg_comparison_id"] == "NRG-L2.0-K8000-Nz8"
    for name, row in results["nrg_results"].items():
        lam, keep, nz = re.fullmatch(r"NRG-L([\d.]+)-K(\d+)-Nz(\d+)", name).groups()
        assert row["settings"] == {"backend": "nrg", "lambda": float(lam), "keep": int(keep), "nz": int(nz),
                                   "keepenergy": 14, "keepmin": 1000, "nmax": {"1.8": 48, "2.0": 42, "2.5": 34}[lam]}
        assert row["finite_flow_completed"] is True and row["continuum_target_established"] is False
        assert row["uncertainty_kind"] == "not_established" and row["empirical_uncertainties"] is None


def test_csv_tables_retain_settings_and_all_target_normalized_operands(results):
    records = results["finite_results"] | results["nrg_results"]
    with (OUTPUT / "observables.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({(r["result_id"], r["observable"]) for r in rows}) == 210
    for row in rows:
        name, observable = row["result_id"], row["observable"]
        record = records[name]
        settings = record["settings"]
        target = results["targets"][settings["backend"]]
        assert float(row["value"]) == record["values"][observable]
        assert float(row["requested_observable_target"]) == target["absolute"] + target["relative"] * abs(float(row["value"]))
        for key, value in settings.items():
            if key == "bond_dimensions":
                for branch, bond in value.items():
                    assert row[f"{branch}_bond_dimension"] == str(bond)
            else:
                assert row[key] == ("" if value is None else str(value))
        nrg = settings["backend"] == "nrg"
        assert row["empirical_uncertainty"] == "" and row["uncertainty_kind"] == ("not_established" if nrg else "")
        assert row["finite_problem_qualified"] == ("" if nrg else str(record["finite_problem_qualified"]).lower())
        assert row["finite_flow_completed"] == ("true" if nrg else "") and row["continuum_target_established"] == "false"
        assert row["role"] == ("comparison_point" if name == results["nrg_comparison_id"] else "final" if name in results["final_result_ids"] else "refinement")
    pairs = {pair["id"]: pair for pair in results["comparisons"]}
    assert len(pairs) == len(results["comparisons"]) == 18
    with (OUTPUT / "comparisons.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({(r["comparison_id"], r["observable"]) for r in rows}) == 180
    for row in rows:
        pair, name = pairs[row["comparison_id"]], row["observable"]
        for key in ("kind", "new_id", "old_id", "denominator_type"):
            assert row[key] == pair[key]
        new, old = (records[pair[f"{side}_id"]] for side in ("new", "old"))
        assert float(row["new_value"]) == new["values"][name] and float(row["old_value"]) == old["values"][name]
        target = results["targets"][new["settings"]["backend"]]
        difference = new["values"][name] - old["values"][name]
        scale = target["absolute"] + target["relative"] * abs(new["values"][name])
        expected = dict(difference=difference, scale=scale, ratio=abs(difference) / scale)
        assert pair["observables"][name] == {key: float(row[key]) for key in expected} == expected
        assert row["denominator_type"] == "requested_observable_target"
        assert "new backend requested target, not uncertainty or error bound" in row["normalization_meaning"]


def test_local_provenance_baths_and_branch_recovery_disclosures(results):
    provenance = json.loads((CASE / "reference/provenance.json").read_text())
    hashes = provenance["published_file_sha256"]
    assert set(hashes) == {"input/physical.json", "input/baths.json", "manual_run/reference.json"} | {f"output/{name}" for name in TABLES + FIGURES}
    for name, digest in hashes.items():
        assert hashlib.sha256((CASE / name).read_bytes()).hexdigest() == digest
    assert hashes["input/baths.json"] == provenance["source_file_sha256"]["runs/main-v2/baths.json"]
    baths = json.loads((CASE / "input/baths.json").read_text())
    sources, branches = provenance["finite_endpoint_sources"], provenance["dmrg_branch_sources"]
    assert set(sources) | set(branches) == set(results["finite_results"])
    for name, row in results["finite_results"].items():
        bath = baths[sources[name]["bath_key"] if name in sources else "cosh16"]
        assert len(bath["xi"]) == len(bath["weights"]) == row["settings"]["bath_levels"]
        assert bath["bandwidth"] == 100 and bath["delta"] == 1
        assert all(math.isfinite(x) and abs(x) <= 100 for x in bath["xi"])
        assert all(math.isfinite(w) and w > 0 for w in bath["weights"])
        assert bath["metadata"].get("frequency_cutoff") == row["settings"]["frequency_window"]
    receipts = provenance["original_accepted_receipts"] | provenance["completion_accepted_receipts"]
    for receipt in receipts.values():
        assert receipt["verified"] is True
        assert receipt["result_sha256"] == provenance["source_file_sha256"][receipt["result"]]
    assert set(provenance["nrg_mean_sources"]) == set(results["nrg_results"])
    for name, source in provenance["nrg_mean_sources"].items():
        assert len(source["twist_records"]) == len(set(source["twist_records"])) == results["nrg_results"][name]["settings"]["nz"]
        assert source["contains_recovered_flow"] is (name == "NRG-L1.8-K8000-Nz4")
        assert set(source["twist_records"]) <= set(receipts) | {"nrg-l18-z4of4-k8000"}
    assert "timed out" in branches["cosh16-dmrg-chi512"]["singlet"]["evidence"]
    for branch, chi in (("singlet", "1024"), ("doublet", "0768")):
        source = branches["cosh16-dmrg-final"][branch]
        assert source["result"] == receipts[f"cosh16-{branch}-chi{chi}"]["result"]
        assert "one-root, one-warm-trial" in source["evidence"] and "no new independent-seed/root-separation claim" in source["evidence"]
    recovery = provenance["recovery"]
    assert recovery["receipt_verified"] is True and recovery["monitored_seconds"] is None
    assert recovery["worker_exit"] == "failed after NRG completed" and "did not rerun" in recovery["disclosure"]
    for name, digest in recovery["artifact_hashes_verified"].items():
        assert digest == provenance["source_file_sha256"][f"{recovery['source']}/{name}"]


@pytest.mark.parametrize("name,title", [("comparison", "Detuned impurity: finite-bath comparison"),
                                        ("convergence", "Refinement and residual diagnostics")])
def test_svg_figures_are_readable_standalone_and_match_data(results, generator, name, title):
    text = (OUTPUT / f"{name}.svg").read_text()
    root = ET.fromstring(text)
    assert root.tag == "{http://www.w3.org/2000/svg}svg" and root.get("viewBox")
    labels = " ".join("".join(node.itertext()) for node in root.findall(".//{http://www.w3.org/2000/svg}text"))
    assert title in labels and "residual" in labels and "continuum" in labels
    payload = root.find(".//{http://purl.org/dc/elements/1.1/}description")
    assert json.loads(payload.text) == generator["figure_data"](results)[name]
    assert "@import" not in text and "<?xml-stylesheet" not in text
    assert all(value.strip(" '\"").startswith("#") for value in re.findall(r"url\(([^)]+)\)", text))
    for node in root.iter():
        assert node.tag.rsplit("}", 1)[-1] not in {"script", "image", "foreignObject"}
        for key, value in node.attrib.items():
            assert not key.lower().startswith("on")
            if key.rsplit("}", 1)[-1] == "href":
                assert value.startswith("#")


def test_stdlib_tables_regenerate_identically_and_refuse_published_aliases(tmp_path):
    published = {name: (OUTPUT / name).read_bytes() for name in TABLES + FIGURES}
    command = [sys.executable, "-B", "-S", str(CASE / "generate.py"), "--tables-only", "--output-dir"]
    destination = tmp_path / "tables"
    run = subprocess.run([*command, str(destination)], cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stdout + run.stderr
    assert {p.name for p in destination.iterdir()} == set(TABLES)
    for name in TABLES:
        assert (destination / name).read_bytes() == published[name]
    alias = tmp_path / "published-alias"
    alias.symlink_to(OUTPUT, target_is_directory=True)
    file_alias = tmp_path / "file-alias"
    file_alias.mkdir()
    (file_alias / "results.json").symlink_to(OUTPUT / "results.json")
    for forbidden in (OUTPUT, alias, file_alias):
        run = subprocess.run([*command, str(forbidden)], cwd=tmp_path, capture_output=True, text=True, timeout=30)
        assert run.returncode == 2 and "published outputs are read-only" in run.stderr
    assert {name: (OUTPUT / name).read_bytes() for name in published} == published


def test_documented_local_file_links_exist():
    text = (CASE / "README.md").read_text()
    links = re.findall(r"\[[^\]]*\]\(([^\s)]+)\)", text)
    assert {"output/results.json", "reference/provenance.json"} <= set(links)
    for link in links:
        if not re.match(r"[a-zA-Z][\w+.-]*:", link) and not link.startswith("#"):
            assert (CASE / link.split("#", 1)[0]).exists(), link


def test_publication_privacy_and_results_first_structure():
    provenance = json.loads((CASE / "reference/provenance.json").read_text())
    assert "release" not in provenance["numerical_environment"]
    assert "relocation" not in provenance and "external_version_evidence" not in provenance
    assert provenance["recorded_original_file_sha256"]["manual_run/reference.json"] != provenance[
        "published_file_sha256"]["manual_run/reference.json"]
    for path in [*(OUTPUT / name for name in TABLES + FIGURES), CASE / "reference/provenance.json",
                 CASE / "manual_run/reference.json", CASE / "README.md", CASE / "manual_run/README.md"]:
        text = path.read_text()
        assert not re.search(r"/(?:project|home|Users)/|qdjj_solver_old|\b\d{4}-\d{2}-\d{2}\b|"
                             r'"(?:hostname|uid|username|session|execution)"\s*:', text), path
    text = (CASE / "README.md").read_text()
    headings = ["## TL;DR", "## Model", "## Energies and Gap", "## Probabilities and Moment",
                "## Achieved Accuracy", "## Data and Provenance", "## Reproduction", "## Supplemental material"]
    assert [text.index(heading) for heading in headings] == sorted(text.index(heading) for heading in headings)
    assert text.index("| Cosh 16, QP cutoff 6 to 7") > text.index("## Supplemental material")
