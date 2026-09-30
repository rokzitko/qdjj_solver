"""Offline checks of the static QP/DMRG/NRG tables and figures, without solvers."""

import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest


OUTPUT = Path(__file__).resolve().parents[1] / "NRG_comparisons/test1/output"


@pytest.fixture(scope="module")
def results():
    data = (OUTPUT / "results.json").read_text()
    assert len(data) < 60_000
    for excluded in ("campaign", "sha256", "runs/", "monitored_seconds", "manifest"):
        assert excluded not in data
    result = json.loads(data)
    json.dumps(result, allow_nan=False)
    return result


def test_static_results_preserve_definitions_and_physical_identities(results):
    assert results["format_version"] == 1
    assert results["model"] == dict(bandwidth=100., detuning=0., field=0., gamma=.4,
                                    gap=1., geometry="single", u=2.)
    assert results["conventions"]["compress"] is False
    assert len(results["observables"]) == 10
    assert "singlet.moment" not in results["observables"]
    assert len(results["finite_results"]) == 12
    assert results["final_result_ids"] == ["cosh26-qp-q6", "surrogate24-window400-qp-q6",
                                           "cosh12-dmrg-chi512", "cosh12-qp-q8"]
    rows = [results["nrg_reference"], *results["finite_results"].values()]
    for row in rows:
        values = row["values"]
        assert set(values) == set(results["observables"])
        # Extrapolating the gap directly permits last-bit differences in subtraction.
        assert values["signed_gap"] == pytest.approx(values["doublet.energy"] - values["singlet.energy"], abs=1e-14, rel=0)
        for branch in ("singlet", "doublet"):
            probabilities = [values[f"{branch}.P{i}"] for i in range(3)]
            assert all(0 <= p <= 1 for p in probabilities)
            assert sum(probabilities) == pytest.approx(1, abs=1e-10, rel=0)
            assert probabilities[0] == pytest.approx(probabilities[2], abs=1e-10, rel=0)
        assert 0 <= values["doublet.moment"] <= .5


def test_empirical_uncertainties_and_residual_failures_are_not_promoted(results):
    nrg = results["nrg_reference"]
    assert nrg["uncertainty_kind"] == "conditional_empirical_not_statistical"
    assert nrg["settings"]["lambdas"] == [1.8, 2., 2.5, 3., 4.]
    assert nrg["settings"]["nz"] == 32 and nrg["settings"]["reference_keep"] == 4000
    assert results["targets"] == dict(nrg=dict(absolute=1e-6, relative=1e-3),
                                      qp=dict(absolute=1e-9, relative=1e-6),
                                      dmrg=dict(absolute=1e-9, relative=1e-6))
    for name, value in nrg["values"].items():
        assert 0 < nrg["empirical_uncertainties"][name] <= 1e-6 + 1e-3 * abs(value)
    for row in results["finite_results"].values():
        qp = row["settings"]["backend"] == "qp"
        assert row["residual_scope"] == ("projected_cutoff_hamiltonian" if qp else "untruncated_finite_hamiltonian")
        assert row["residual_tolerance"] == (1e-8 if qp else 1e-7)
        assert set(row["residuals"]) == {"singlet", "doublet"}
        assert all(len(values) == 2 for values in row["residuals"].values())
        residuals = [r for values in row["residuals"].values() for r in values]
        assert all(r >= 0 for r in residuals)
        assert row["finite_problem_qualified"] is all(r <= row["residual_tolerance"] for r in residuals)
        assert row["continuum_target_established"] is False
    final = results["finite_results"]["cosh12-dmrg-chi512"]
    assert final["finite_problem_qualified"] is False
    assert sum(r <= 1e-7 for values in final["residuals"].values() for r in values) == 1


def test_comparisons_reproduce_all_operands_and_use_named_denominators(results):
    rows = dict(results["finite_results"], nrg_reference=results["nrg_reference"])
    comparisons = results["comparisons"]
    assert len(comparisons) == len({row["id"] for row in comparisons}) == 14
    for comparison in comparisons:
        new, old = rows[comparison["new_id"]], rows[comparison["old_id"]]
        assert set(comparison["observables"]) == set(results["observables"])
        for name, entry in comparison["observables"].items():
            difference = new["values"][name] - old["values"][name]
            if comparison["denominator_type"] == "nrg_empirical_uncertainty":
                assert comparison["new_id"] == "nrg_reference"
                scale = new["empirical_uncertainties"][name]
            else:
                assert comparison["denominator_type"] == "requested_observable_target"
                target = results["targets"][new["settings"]["backend"]]
                scale = target["absolute"] + target["relative"] * abs(new["values"][name])
            assert entry == dict(difference=difference, scale=scale, ratio=abs(difference) / scale)
    pairs = {(c["new_id"], c["old_id"]): c["observables"] for c in comparisons}
    for endpoint in results["final_result_ids"][:2]:
        assert max(v["ratio"] for v in pairs["nrg_reference", endpoint].values()) < .074
    assert pairs["nrg_reference", "cosh12-dmrg-chi512"]["signed_gap"]["ratio"] > 2
    assert pairs["nrg_reference", "cosh12-dmrg-chi512"]["singlet.P0"]["ratio"] > 12
    assert pairs["cosh26-qp-q6", "cosh26-qp-q5"]["signed_gap"]["ratio"] > 2.19
    assert abs(pairs["cosh12-dmrg-chi512", "cosh12-qp-q8"]["signed_gap"]["difference"]) < 4e-14


def test_csv_tables_match_the_complete_json_dataset(results):
    records = dict(results["finite_results"], nrg_reference=results["nrg_reference"])
    with (OUTPUT / "observables.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({(r["result_id"], r["observable"]) for r in rows}) == 130
    for row in rows:
        record = records[row["result_id"]]
        name = row["observable"]
        assert float(row["value"]) == record["values"][name]
        nrg = row["result_id"] == "nrg_reference"
        backend = "nrg" if nrg else record["settings"]["backend"]
        target = results["targets"][backend]
        assert row["backend"] == backend
        assert float(row["requested_observable_target"]) == target["absolute"] + target["relative"] * abs(float(row["value"]))
        if nrg:
            assert float(row["empirical_uncertainty"]) == record["empirical_uncertainties"][name]
            assert row["uncertainty_kind"] == record["uncertainty_kind"]
            assert row["role"] == "reference"
            assert row["finite_problem_qualified"] == row["continuum_target_established"] == ""
        else:
            assert row["empirical_uncertainty"] == row["uncertainty_kind"] == ""
            assert row["finite_problem_qualified"] == str(record["finite_problem_qualified"]).lower()
            assert row["continuum_target_established"] == "false"
            assert row["role"] == ("final" if row["result_id"] in results["final_result_ids"] else "refinement")
            for key in ("bath_family", "bath_levels", "qp_cutoff", "bond_dimension", "frequency_window"):
                value = record["settings"][key]
                assert row[key] == ("" if value is None else str(value))
    with (OUTPUT / "comparisons.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({(r["comparison_id"], r["observable"]) for r in rows}) == 140
    comparisons = {row["id"]: row for row in results["comparisons"]}
    for row in rows:
        comparison = comparisons[row["comparison_id"]]
        for key in ("kind", "new_id", "old_id", "denominator_type"):
            assert row[key] == comparison[key]
        for key in ("difference", "scale", "ratio"):
            assert float(row[key]) == comparison["observables"][row["observable"]][key]
        for side in ("new", "old"):
            assert float(row[f"{side}_value"]) == records[row[f"{side}_id"]]["values"][row["observable"]]


@pytest.mark.parametrize("name,title", [("comparison.svg", "Finite-bath results against NRG"),
                                        ("convergence.svg", "Refinement and residual diagnostics")])
def test_figures_are_standalone_svg_with_readable_labels(name, title):
    data = (OUTPUT / name).read_bytes()
    assert len(data) < 150_000
    root = ET.fromstring(data)
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("viewBox")
    text = " ".join(root.itertext())
    assert title in text and "residual" in text
    for element in root.iter():
        assert element.tag.rsplit("}", 1)[-1] not in {"script", "image", "foreignObject"}
        for key, value in element.attrib.items():
            if key.endswith("href"):
                assert value.startswith("#"), "figures must not load external resources"
