"""Reference contracts, executable handbook assets and publication-specific behaviour."""

import ast
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

import numpy as np
import pytest

pytest.importorskip("markdown_it")
pytestmark = pytest.mark.docs
ROOT = Path(__file__).resolve().parents[1]


def tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/"tools"/f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run_script(path, *args, cwd):
    environment = dict(os.environ)
    environment.update({key: "1" for key in tool("run_doc_examples").THREADS})
    return subprocess.run([sys.executable, "-B", str(path), *map(str, args)], cwd=cwd,
                          env=environment, capture_output=True, text=True, check=True, timeout=120)


def test_shared_reference_input_preserves_reviewed_model():
    from qdjj_solver import cosh_grid, reference_model
    from qdjj_solver.common.io import model_from_config
    config = json.loads((ROOT/"examples/01_subgap_states/input/reference.json").read_text())
    expected = reference_model(cosh_grid(pairs=4, delta=1., bandwidth=10.), u=1., gamma=.5,
                               phi=.8, detuning=.2, field=0., symmetry=False, compress=False)
    assert model_from_config(config).fingerprint() == expected.fingerprint()


def test_portable_static_workflow_and_cli(tmp_path):
    from qdjj_solver import load_result
    source = tmp_path/"handbook"
    source.mkdir()
    shutil.copy(ROOT/"examples/handbook/static.py", source/"static.py")
    inputs = tmp_path/"01_subgap_states/input"
    inputs.mkdir(parents=True)
    shutil.copy(ROOT/"examples/01_subgap_states/input/reference.json", inputs/"reference.json")
    output = tmp_path/"results"
    run_script(source/"static.py", "--output", output, cwd=tmp_path)
    result = json.loads((output/"summary.json").read_text())["values"]
    expected = json.loads((ROOT/"examples/handbook/reference.json").read_text())["values"]
    for name, entry in expected.items():
        assert abs(result[name]-entry["value"]) <= entry["atol"], name
    subprocess.run([sys.executable, "-B", "-m", "qdjj_solver.cli", "solve", str(inputs/"reference.json"),
                    "--output", str(output/"cli.json"), "--save-states"], cwd=tmp_path,
                   check=True, capture_output=True, timeout=60)
    np.testing.assert_allclose(load_result(output/"cli.json").energies,
                               load_result(output/"even.json").energies, atol=1e-10)
    moved = tmp_path/"moved"
    moved.mkdir()
    scalar = json.loads((output/"even.json").read_text())
    for name in ("even.json", scalar["state_artifact"]["file"]):
        shutil.copy(output/name, moved/name)
    restored = load_result(moved/"even.json")
    assert restored.matrix_elements("charge").shape == (2, 2)


def test_scan_analysis_does_not_resolve_states(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    run_script(ROOT/"examples/handbook/scan.py", "--points", 3, "--output", tmp_path, cwd=tmp_path)
    spec = importlib.util.spec_from_file_location("handbook_scan", ROOT/"examples/handbook/scan.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    def forbidden(*args, **kwargs):
        raise AssertionError("analysis must not invoke a solver")
    monkeypatch.setattr(module, "solve", forbidden)
    module.plot(tmp_path)
    record = json.loads((tmp_path/"scan.json").read_text())
    assert len(record["points"]) == 3
    assert all(row["status"] == "ok" for row in record["points"])
    assert record["points"][0]["doublet"] > record["points"][0]["singlet"]
    assert record["points"][-1]["doublet"] < record["points"][-1]["singlet"]


def test_matched_general_model(tmp_path):
    from qdjj_solver import load_result
    for name in ("custom.py", "custom.json"):
        shutil.copy(ROOT/"examples/handbook"/name, tmp_path/name)
    run_script(tmp_path/"custom.py", "--output", tmp_path/"python.json", cwd=tmp_path)
    subprocess.run([sys.executable, "-B", "-m", "qdjj_solver.cli", "solve", str(tmp_path/"custom.json"),
                    "--output", str(tmp_path/"cli.json"), "--save-states"], check=True, cwd=tmp_path,
                   capture_output=True, timeout=60)
    np.testing.assert_allclose(load_result(tmp_path/"python.json").energies,
                               load_result(tmp_path/"cli.json").energies, atol=1e-12)


def test_parameter_inventory_rejects_omission(tmp_path, monkeypatch):
    generator = tool("gen_docs")
    (tmp_path/"docs").mkdir()
    text = (ROOT/"docs/parameters.toml").read_text()
    text = text.replace('pairs = "Positive integer number', 'missing_pairs = "Positive integer number', 1)
    (tmp_path/"docs/parameters.toml").write_text(text)
    monkeypatch.setattr(generator, "ROOT", tmp_path)
    with pytest.raises(ValueError, match="missing=.*pairs"):
        generator.parameter_entries()


def test_format_inventory_matches_envelopes(tmp_path):
    from qdjj_solver import Sector, cosh_grid, reference_model, save_result, solve
    formats = tomllib.loads((ROOT/"docs/formats.toml").read_text())
    documented = set(formats["configuration"][0]["fields"])
    source = ast.parse((ROOT/"src/qdjj_solver/common/io.py").read_text())
    validate = next(node for node in source.body if isinstance(node, ast.FunctionDef) and node.name == "validate_config")
    allowed = next(node.value for node in validate.body if isinstance(node, ast.Assign)
                   and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "allowed")
    assert documented == ast.literal_eval(allowed) | {"backend"}
    h = reference_model(cosh_grid(1, bandwidth=5.), compress=False)
    state = solve(h, cutoff=2, sector=Sector(1, 1))
    save_result(state, tmp_path/"state.json", save_states=True)
    record = json.loads((tmp_path/"state.json").read_text())
    assert set(record) <= set(formats["record"][0]["fields"])
    qp = next(section for section in formats["record"] if section["title"] == "QP calculation metadata")
    assert set(qp["fields"]) <= set(record["calculation"])


def test_handbook_selection_is_explicit():
    pages = tool("run_doc_examples").handbook_pages()
    assert len(pages) == len(set(pages))
    assert "docs/getting_started.md" in pages and "docs/reference/api.md" in pages
    assert not any(page.startswith(("applications/", "large_Gamma/", "NRG_comparisons/")) for page in pages)
    assert "docs/bath_representations.md" not in pages  # complete research report remains linked
    for page in pages:
        assert (ROOT/page).is_file()


def test_generated_listing_is_executed_in_document_order():
    listings = tool("run_manual_listings").listings()
    assert listings[0][0].name == "practical_code.tex"
    assert "reference_model" in listings[0][2]
    assert any("save_result(even" in code for _, _, code in listings)
