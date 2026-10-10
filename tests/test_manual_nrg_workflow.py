"""Bounded standalone-controller tests using fake executables, never real NRG."""

import csv
from datetime import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "NRG_comparisons/test2"
MANUAL = CASE / "manual_run"
ARTIFACTS = ("param", "data", "raw.h5", "DONE", "nrginit.log", "nrg.log", "annotated.dat")
OBSERVABLES = ("signed_gap", "singlet.energy", "doublet.energy", "singlet.P0",
               "singlet.P1", "singlet.P2", "doublet.P0", "doublet.P1",
               "doublet.P2", "doublet.moment")
REAL_LOOPER = os.environ.get("TEST_REAL_LOOPER") or shutil.which("looper")
pytestmark = pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX shells and flock")


# Each executable records its actual argv/cwd/environment outside the run folder.
# Fake raw.h5 is intentionally not HDF5: extraction has separate numerical tests.
FAKE_TOOL = r'''
import json
import os
from pathlib import Path
import re
import sys

tool = Path(__file__).name
cwd = Path.cwd()
keys = ("NRG_WORKDIR", "TMPDIR", "MKL_THREADING_LAYER", "MKL_NUM_THREADS",
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_DYNAMIC")
with open(os.environ["TEST_NRG_EVENTS"], "a") as stream:
    stream.write(json.dumps(dict(tool=tool, argv=sys.argv[1:], cwd=str(cwd),
                                env={key: os.environ.get(key) for key in keys})) + "\n")
mode = os.environ.get("TEST_NRG_FAILURE", "")
if os.environ.get("TEST_NRG_FAIL_AT", cwd.name) != cwd.name:
    mode = ""
print(tool + " stdout", flush=True)
print(tool + " stderr", file=sys.stderr, flush=True)
if tool == "looper":
    assert sys.argv[1:] == ["param.loop"]
    text = Path("param.loop").read_text()
    nz = int(re.search(r"^#PRELUDE: \$Nz = ([0-9]+);$", text, re.M)[1])
    if mode != "silent-looper":
        for i in range(1, nz + 1):
            if mode == "missing-directory" and i == nz:
                continue
            directory = Path(str(i))
            directory.mkdir()
            z = .5 if mode == "wrong-z" and i == nz else i / nz
            param = "[param]" + text.split("[param]", 1)[1]
            param = param.replace("$Nz", str(nz)).replace("$z", str(z))
            if mode == "mixed-settings" and i == nz:
                param = param.replace("keep=8000", "keep=4000")
            (directory / "param").write_text(param)
        if mode == "extra-directory":
            Path(str(nz + 1)).mkdir()
elif tool == "nrginit":
    if mode != "initializer-no-data":
        Path("data").write_text("fake Wilson chain\n")
    if mode != "initializer-no-success":
        print("Success!", flush=True)
    if mode == "initializer-exit":
        sys.exit(17)
elif tool == "nrg":
    if mode != "solver-no-raw":
        Path("raw.h5").write_bytes(b"fake raw output\n")
    if mode != "solver-no-done":
        Path("DONE").touch()
    Path("annotated.dat").write_text("fake spectrum\n")
    if mode == "solver-param-change":
        with Path("param").open("a") as stream:
            stream.write("# changed during execution\n")
    if mode == "solver-exit":
        sys.exit(23)
else:
    raise AssertionError("fake Mathematica must not be invoked")
'''


def _module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def tools(tmp_path, monkeypatch):
    directory = tmp_path / "fake tools"
    directory.mkdir()
    paths = {"events": tmp_path / "events.jsonl"}
    for name, variable in (("looper", "LOOPER"), ("nrginit", "NRGINIT"), ("nrg", "NRG"), ("math", "MATH")):
        path = directory / name
        path.write_text(f"#!{sys.executable}\n" + FAKE_TOOL)
        path.chmod(0o755)
        monkeypatch.setenv(variable, str(path))
        paths[name] = path
    monkeypatch.setenv("TEST_NRG_EVENTS", str(paths["events"]))
    monkeypatch.delenv("TEST_NRG_FAILURE", raising=False)
    monkeypatch.delenv("TEST_NRG_FAIL_AT", raising=False)
    monkeypatch.delenv("NRG_THREADS", raising=False)
    monkeypatch.setenv("PYTHON", sys.executable)
    return paths


@pytest.fixture
def manual(tmp_path, monkeypatch, tools):
    folder = tmp_path / "copied manual"
    folder.mkdir()
    for name in ("process.py", "extract.py", "param.loop", "reference.json", "1_run", "2_proc", "run_one"):
        shutil.copy2(MANUAL / name, folder / name)
    monkeypatch.setitem(sys.modules, "extract", _module(folder / "extract.py", "copied_manual_extract"))
    return _module(folder / "process.py", "copied_manual_process")


def _events(tools):
    path = tools["events"]
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def _prepare(manual, nz=8, **settings):
    path = manual.HERE / "param.loop"
    text = re.sub(r"^#PRELUDE:.*$", f"#PRELUDE: $Nz = {nz};", path.read_text(), flags=re.M)
    for key, value in settings.items():
        text, count = re.subn(rf"^{key}=.*$", f"{key}={value}", text, flags=re.M)
        assert count == 1
    path.write_text(text)
    return manual.prepare(manual.HERE)


def _completed_sweep(manual, nz=8, **settings):
    """Synthetic receipts isolate averaging from execution and HDF5 extraction."""
    manifest = _prepare(manual, nz, **settings)
    for row in manifest["runs"]:
        directory = manual.HERE / row["directory"]
        for name in ARTIFACTS:
            if name != "param":
                (directory / name).write_text(f"{name}: twist {row['directory']}\n")
        receipt = {"param_sha256": _digest(directory / "param"),
                   "artifacts": {name: _digest(directory / name) for name in ARTIFACTS}}
        (directory / "completed.json").write_text(json.dumps(receipt))
    return manifest


def _result(directory):
    # Nonlinear twist dependence distinguishes arithmetic means and subgrids.
    i = int(directory.name)
    return {"values": {key: (j - 4) / 8 + i * i / 128 for j, key in enumerate(OBSERVABLES)},
            "diagnostics": {"synthetic_twist": i}}


def test_prepare_default_grid_and_shared_renderer_contract(manual, tools, tmp_path):
    manifest = manual.prepare(manual.HERE)
    assert manifest == json.loads((manual.HERE / "sweep.json").read_text())
    assert manifest["template_sha256"] == _digest(manual.HERE / "param.loop")
    assert {p.name for p in manual.HERE.iterdir() if p.is_dir()} >= set(map(str, range(1, 9)))
    assert [row["directory"] for row in manifest["runs"]] == list(map(str, range(1, 9)))
    assert [row["z"] for row in manifest["runs"]] == [i / 8 for i in range(1, 9)]
    assert {key: manifest["parameters"][key] for key in ("Nz", "Lambda", "Nmax", "keep", "keepmin", "keepenergy")} == {
        "Nz": "8", "Lambda": "2", "Nmax": "42", "keep": "8000", "keepmin": "1000", "keepenergy": "14"}
    for row in manifest["runs"]:
        directory = manual.HERE / row["directory"]
        assert {p.name for p in directory.iterdir()} == {"param"}
        assert row["param_sha256"] == _digest(directory / "param")
        params = manual.read_param(directory / "param")
        assert float(params.pop("z")) == row["z"]
        assert params == manifest["parameters"]
    assert [(e["tool"], e["argv"], e["cwd"]) for e in _events(tools)] == [
        ("looper", ["param.loop"], str(manual.HERE))]

    shared = _module(ROOT / "NRG_comparisons/test1/nrg.py", "shared_nrg_renderer")
    physical = dict(gap=1., bandwidth=100., u=2., gamma=.4, detuning=.5, field=0., geometry="single")
    numerical = {"lambda": 2., "z": .125, "nmax": 42, "keep": 8000, "keepmin": 1000, "keepenergy": 14}
    rendered = tmp_path / "shared.param"
    rendered.write_text(shared.render_param(physical, numerical))
    expected = manual.read_param(rendered)
    actual = manual.read_param(manual.HERE / "1/param")
    assert set(actual) == set(expected) | {"Nz"}
    for key, value in expected.items():
        try:
            number = float(value)
        except ValueError:
            assert actual[key] == value
        else:
            assert float(actual[key]) == pytest.approx(number, rel=0, abs=1e-14)
    template = (manual.HERE / "param.loop").read_text()
    assert "@ $i=1; $i<=$Nz; $i++\n!($z=$i/$Nz)\n[param]" in template


@pytest.mark.skipif(not REAL_LOOPER, reason="optional real looper input-generation check")
@pytest.mark.parametrize("nz", [3, 8])
def test_real_looper_only_generates_inputs(manual, tools, monkeypatch, nz):
    monkeypatch.setenv("LOOPER", REAL_LOOPER)
    manifest = _prepare(manual, nz)
    assert len(manifest["runs"]) == nz
    for i, entry in enumerate(manifest["runs"], 1):
        directory = manual.HERE / entry["directory"]
        assert {path.name for path in directory.iterdir()} == {"param"}
        assert float(manual.read_param(directory / "param")["z"]) == pytest.approx(i / nz, abs=1e-14)
    assert _events(tools) == []  # No initializer, NRG, or Mathematica calls.


@pytest.mark.parametrize("existing", ["directory", "file", "manifest"])
def test_prepare_refuses_existing_work_without_launching_looper(manual, tools, existing):
    path = manual.HERE / ("sweep.json" if existing == "manifest" else "1")
    if existing == "directory":
        path.mkdir()
        (path / "user-data").write_text("keep me")
    else:
        path.write_text("keep me")
    with pytest.raises(ValueError, match="existing"):
        manual.prepare(manual.HERE)
    assert _events(tools) == []
    assert (path / "user-data" if path.is_dir() else path).read_text() == "keep me"


@pytest.mark.parametrize("failure,match", [
    ("silent-looper", "expected numbered"), ("missing-directory", "expected numbered"),
    ("extra-directory", "expected numbered"), ("wrong-z", "wrong Nz or z"),
    ("mixed-settings", "all parameters except z"),
])
def test_prepare_rejects_silent_or_malformed_looper_output(manual, tools, monkeypatch, failure, match):
    monkeypatch.setenv("TEST_NRG_FAILURE", failure)
    with pytest.raises(ValueError, match=match):
        manual.prepare(manual.HERE)
    assert not (manual.HERE / "sweep.json").exists()
    assert [e["tool"] for e in _events(tools)] == ["looper"]


@pytest.mark.parametrize("old,new", [
    ("$Nz = 8;", "$Nz = 4 + 4;"),
    ("#AUTOLOOP: true", "#AUTOLOOP: nrg"),
    ("#AUTOLOOP: true", "#AUTOLOOP: true\n#OVERWRITE: true"),
])
def test_prepare_rejects_unsafe_template_before_looper(manual, tools, old, new):
    path = manual.HERE / "param.loop"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match="param.loop"):
        manual.prepare(manual.HERE)
    assert _events(tools) == []
    assert not (manual.HERE / "sweep.json").exists()


@pytest.mark.parametrize("edit", ["setting", "comment"])
def test_prepared_param_changes_block_execution_and_averaging(manual, tools, edit):
    _prepare(manual)
    path = manual.HERE / "8/param"
    text = path.read_text()
    path.write_text(text.replace("Lambda=2\n", "Lambda=2.5\n") if edit == "setting" else text + "# edited\n")
    with pytest.raises(ValueError, match="parameters changed"):
        manual.run_one(manual.HERE / "1")
    with pytest.raises(ValueError, match="parameters changed"):
        manual.average(manual.HERE)
    assert [e["tool"] for e in _events(tools)] == ["looper"]
    assert not (manual.HERE / "1/started.json").exists()
    assert not (manual.HERE / "processed").exists()


@pytest.mark.parametrize("threads", [None, "3"])
def test_run_one_logs_environment_and_verified_receipt(manual, tools, monkeypatch, threads):
    _prepare(manual)
    directory = manual.HERE / "3"
    if threads is not None:
        monkeypatch.setenv("NRG_THREADS", threads)
    for key in ("NRG_WORKDIR", "TMPDIR", "MKL_THREADING_LAYER", "MKL_NUM_THREADS",
                "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_DYNAMIC"):
        monkeypatch.setenv(key, "wrong inherited value")
    manual.run_one(directory)
    events = _events(tools)[1:]
    assert [e["tool"] for e in events] == ["nrginit", "nrg"]
    assert [e["argv"] for e in events] == [[str(tools["math"])], []]
    for event in events:
        assert event["cwd"] == str(directory)
        assert event["env"] == {
            "NRG_WORKDIR": str(directory / "scratch"), "TMPDIR": str(directory / "scratch"),
            "MKL_THREADING_LAYER": "INTEL", "MKL_NUM_THREADS": threads or "4",
            "OMP_NUM_THREADS": threads or "4", "OPENBLAS_NUM_THREADS": "1", "MKL_DYNAMIC": "FALSE"}
    assert (directory / "scratch").is_dir()
    for name in ("nrginit", "nrg"):
        log = (directory / f"{name}.log").read_text()
        assert f"{name} stdout" in log and f"{name} stderr" in log
    started = json.loads((directory / "started.json").read_text())
    receipt = json.loads((directory / "completed.json").read_text())
    assert receipt.items() >= started.items()
    assert started["threads"] == int(threads or "4")
    assert started["param_sha256"] == _digest(directory / "param")
    assert started["executables"] == {key: str(tools[key]) for key in ("nrginit", "nrg", "math")}
    assert datetime.fromisoformat(receipt["finished_utc"]) >= datetime.fromisoformat(started["started_utc"])
    assert receipt["artifacts"] == {name: _digest(directory / name) for name in ARTIFACTS}
    assert not (manual.HERE / "1/started.json").exists()
    before = {name: (directory / name).read_bytes() for name in (*ARTIFACTS, "completed.json")}
    with pytest.raises(ValueError, match="prior outputs"):
        manual.run_one(directory)
    assert _events(tools)[1:] == events
    assert before == {name: (directory / name).read_bytes() for name in before}


@pytest.mark.parametrize("failure", [
    "initializer-exit", "initializer-no-success", "initializer-no-data",
    "solver-exit", "solver-no-done", "solver-no-raw", "solver-param-change",
])
def test_failed_tools_never_get_completion_receipt(manual, tools, monkeypatch, failure):
    _prepare(manual)
    directory = manual.HERE / "1"
    monkeypatch.setenv("TEST_NRG_FAILURE", failure)
    error = subprocess.CalledProcessError if failure.endswith("-exit") else ValueError
    with pytest.raises(error):
        manual.run_one(directory)
    assert (directory / "started.json").is_file()
    assert not (directory / "completed.json").exists()
    initializer_failed = failure.startswith("initializer-")
    assert [e["tool"] for e in _events(tools)] == ["looper", "nrginit"] + ([] if initializer_failed else ["nrg"])
    assert (directory / "nrg.log").exists() is not initializer_failed
    assert "nrginit stderr" in (directory / "nrginit.log").read_text()
    if not initializer_failed:
        assert "nrg stderr" in (directory / "nrg.log").read_text()
    before = _events(tools)
    with pytest.raises(ValueError, match="prior outputs|parameters changed"):
        manual.run_one(directory)
    assert _events(tools) == before


@pytest.mark.parametrize("stale", ["DONE", "raw.h5", "nrginit.log"])
def test_run_one_refuses_unreceipted_stale_outputs(manual, tools, stale):
    _prepare(manual)
    directory = manual.HERE / "1"
    path = directory / stale
    path.write_text("old attempt")
    with pytest.raises(ValueError, match="prior outputs"):
        manual.run_one(directory)
    assert path.read_text() == "old attempt"
    assert not (directory / "started.json").exists()
    assert not (directory / "completed.json").exists()
    assert [e["tool"] for e in _events(tools)] == ["looper"]


def test_main_run_is_sequential_and_fail_fast(manual, tools, monkeypatch, tmp_path, capsys):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("TEST_NRG_FAILURE", "solver-exit")
    monkeypatch.setenv("TEST_NRG_FAIL_AT", "2")
    with pytest.raises(SystemExit) as error:
        manual.main(["run"])
    assert error.value.code == 1
    assert "manual NRG:" in capsys.readouterr().err
    assert [(e["tool"], Path(e["cwd"]).name) for e in _events(tools)] == [
        ("looper", manual.HERE.name), ("nrginit", "1"), ("nrg", "1"), ("nrginit", "2"), ("nrg", "2")]
    assert (manual.HERE / "1/completed.json").is_file()
    assert (manual.HERE / "2/started.json").is_file()
    assert not (manual.HERE / "2/completed.json").exists()
    for i in range(3, 9):
        assert {p.name for p in (manual.HERE / str(i)).iterdir()} == {"param"}
    assert not (manual.HERE / "processed").exists()


@pytest.mark.parametrize("nz,mean_squares", [(4, {4: 7.5}), (8, {8: 25.5, 4: 30.})])
def test_complete_arithmetic_means_report_comparisons_and_csv(manual, monkeypatch, capsys, nz, mean_squares):
    _completed_sweep(manual, nz)
    calls = []

    def extract(directory):
        calls.append(directory.name)
        return _result(directory)

    monkeypatch.setattr(manual, "extract_result", extract)
    summary = manual.average(manual.HERE)
    assert calls == list(map(str, range(1, nz + 1)))
    assert set(manual.OBSERVABLES) == set(OBSERVABLES)
    assert [row["settings"]["nz"] for row in summary["averages"]] == list(mean_squares)
    assert summary["uncertainty_kind"] == "not_established"
    assert summary["empirical_uncertainties"] is None
    assert summary["continuum_target_established"] is False
    reference = json.loads((manual.HERE / "reference.json").read_text())
    assert summary["reference_source"] == reference["source"]
    for row in summary["per_z"]:
        assert row["values"] == _result(manual.HERE / row["directory"])["values"]
        assert row["diagnostics"] == {"synthetic_twist": int(row["directory"])}
        assert row["completion"] == json.loads((manual.HERE / row["directory"] / "completed.json").read_text())
    for row in summary["averages"]:
        count = row["settings"]["nz"]
        assert row["directories"] == list(map(str, range(nz // count, nz + 1, nz // count)))
        assert row["reference_id"] == f"NRG-L2.0-K8000-Nz{count}"
        assert row["settings"] == reference["nrg_means"][row["reference_id"]]["settings"]
        assert row["values"] == pytest.approx({key: (j - 4) / 8 + mean_squares[count] / 128
                                               for j, key in enumerate(OBSERVABLES)}, abs=1e-15, rel=0)
        comparisons = [c for c in summary["comparisons"] if c["nz"] == count]
        assert {c["observable"] for c in comparisons} == set(OBSERVABLES)
        for comparison in comparisons:
            key = comparison["observable"]
            expected = reference["nrg_means"][row["reference_id"]]["values"][key]
            difference = row["values"][key] - expected
            assert comparison["report_id"] == row["reference_id"]
            assert comparison["mean"] == row["values"][key]
            assert comparison["report"] == expected
            assert comparison["difference"] == pytest.approx(difference, abs=1e-15)
            tolerance = reference["requested_tolerance"]
            assert comparison["difference_over_requested_tolerance"] == pytest.approx(
                abs(difference) / (tolerance["absolute"] + tolerance["relative"] * abs(row["values"][key])))
    output = manual.HERE / "processed"
    assert json.loads((output / "summary.json").read_text()) == summary
    with (output / "per_z.csv").open(newline="") as stream:
        per_z = list(csv.DictReader(stream))
    assert [row["directory"] for row in per_z] == calls
    for row in per_z:
        assert float(row["z"]) == int(row["directory"]) / nz
        assert {key: float(row[key]) for key in OBSERVABLES} == _result(manual.HERE / row["directory"])["values"]
    with (output / "averages.csv").open(newline="") as stream:
        comparisons = list(csv.DictReader(stream))
    assert len(comparisons) == 10 * len(mean_squares)
    for actual, expected in zip(comparisons, summary["comparisons"], strict=True):
        assert actual == {key: str(value) for key, value in expected.items()}
    printed = capsys.readouterr().out
    assert "NOT error estimates" in printed
    assert all(f"NRG-L2.0-K8000-Nz{count}" in printed for count in mean_squares)


def test_edited_numerical_settings_are_used_without_false_report_match(manual, monkeypatch, capsys):
    settings = dict(Lambda="2.5", Nmax="34", keep="4000", keepmin="400", keepenergy="12")
    manifest = _completed_sweep(manual, 4, **settings)
    assert all(manifest["parameters"][key] == value for key, value in settings.items())
    path = manual.HERE / "reference.json"
    reference = json.loads(path.read_text())
    # A seemingly matching Lambda/keep/Nz label is insufficient: every setting matters.
    reference["nrg_means"]["NRG-L2.5-K4000-Nz4"] = {
        "settings": {"backend": "nrg", "nz": 4, "lambda": 2.5, "nmax": 34,
                     "keep": 4000, "keepmin": 400, "keepenergy": 14.},
        "values": _result(manual.HERE / "1")["values"],
    }
    path.write_text(json.dumps(reference))
    monkeypatch.setattr(manual, "extract_result", _result)
    summary = manual.average(manual.HERE)
    assert len(summary["averages"]) == 1
    row = summary["averages"][0]
    assert row["settings"] == {"backend": "nrg", "nz": 4, "lambda": 2.5, "nmax": 34,
                               "keep": 4000, "keepmin": 400, "keepenergy": 12.}
    assert row["reference_id"] is None
    assert all(c[key] is None for c in summary["comparisons"]
               for key in ("report_id", "report", "difference", "difference_over_requested_tolerance"))
    assert "no matching stored point" in capsys.readouterr().out


@pytest.mark.parametrize("damage", [
    "missing-directory", "missing-completed.json", "missing-DONE", "missing-raw.h5",
    "changed-raw.h5", "changed-nrginit.log", "changed-annotated.dat", "receipt-artifacts", "receipt-param",
])
def test_average_refuses_incomplete_or_changed_outputs_without_partial_groups(manual, monkeypatch, damage):
    _completed_sweep(manual)
    monkeypatch.setattr(manual, "extract_result", _result)
    # Damaging an odd twist leaves the entire Nz=4 even subgrid intact.
    directory = manual.HERE / "7"
    if damage == "missing-directory":
        shutil.rmtree(directory)
    elif damage.startswith("missing-"):
        (directory / damage.removeprefix("missing-")).unlink()
    elif damage.startswith("changed-"):
        (directory / damage.removeprefix("changed-")).write_text("changed after completion")
    else:
        path = directory / "completed.json"
        receipt = json.loads(path.read_text())
        if damage == "receipt-artifacts":
            del receipt["artifacts"]["raw.h5"]
        else:
            receipt["param_sha256"] = "0" * 64
        path.write_text(json.dumps(receipt))
    with pytest.raises((ValueError, FileNotFoundError)):
        manual.average(manual.HERE)
    assert not (manual.HERE / "processed").exists()


@pytest.mark.parametrize("invalid", ["missing", "nan"])
def test_average_rejects_incomplete_or_nonfinite_extraction(manual, monkeypatch, invalid):
    _completed_sweep(manual)

    def extract(directory):
        result = _result(directory)
        if directory.name == "7":
            if invalid == "missing":
                del result["values"]["doublet.moment"]
            else:
                result["values"]["signed_gap"] = float("nan")
        return result

    monkeypatch.setattr(manual, "extract_result", extract)
    with pytest.raises(ValueError, match="incomplete or nonfinite"):
        manual.average(manual.HERE)
    assert not (manual.HERE / "processed").exists()


def test_bundled_reference_is_exact_report_subset():
    reference = json.loads((MANUAL / "reference.json").read_text())
    report = json.loads((CASE / "output/results.json").read_text())
    assert reference["physical"] == report["model"]
    assert reference["requested_tolerance"] == report["targets"]["nrg"]
    assert set(reference["nrg_means"]) == {"NRG-L2.0-K8000-Nz4", "NRG-L2.0-K8000-Nz8"}
    for key, row in reference["nrg_means"].items():
        original = report["nrg_results"][key]
        assert row["settings"] == original["settings"]
        assert row["values"] == original["values"]  # Exact bundled floats, not tolerance agreement.
        assert set(row["values"]) == set(OBSERVABLES)
        for field in ("uncertainty_kind", "empirical_uncertainties", "continuum_target_established"):
            assert reference[field] == original[field]


def test_copied_cli_prepares_without_site_packages_or_repository(manual, tools, tmp_path):
    environment = dict(os.environ, PYTHONPATH="", PYTHONNOUSERSITE="1")
    # -S excludes even editable installs and optional NumPy/HDF5 packages.
    result = subprocess.run([sys.executable, "-S", "-B", str(manual.HERE / "process.py"), "run", "--prepare-only"],
                            cwd=tmp_path, env=environment, text=True, capture_output=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert "Prepared 8 inputs" in result.stdout
    assert (manual.HERE / "8/param").is_file()
    assert not (tmp_path / "1").exists()
    assert [e["tool"] for e in _events(tools)] == ["looper"]


def test_shell_wrappers_are_executable(manual):
    missing = [name for name in ("1_run", "2_proc", "run_one") if not os.access(manual.HERE / name, os.X_OK)]
    assert not missing, f"shell wrappers lack executable permission: {missing}"


def test_shell_wrappers_resolve_their_location_and_run_one_uses_cwd(manual, tools, tmp_path):
    prepare = subprocess.run(["sh", str(manual.HERE / "1_run"), "--prepare-only"], cwd=tmp_path,
                             text=True, capture_output=True, timeout=20)
    assert prepare.returncode == 0, prepare.stderr
    assert (manual.HERE / "sweep.json").is_file()
    assert not (tmp_path / "sweep.json").exists()
    directory = manual.HERE / "3"
    run = subprocess.run(["sh", str(manual.HERE / "run_one")], cwd=directory,
                         text=True, capture_output=True, timeout=20)
    assert run.returncode == 0, run.stderr
    assert (directory / "completed.json").is_file()
    assert not (manual.HERE / "1/started.json").exists()
    assert [(e["tool"], e["cwd"]) for e in _events(tools)[1:]] == [
        ("nrginit", str(directory)), ("nrg", str(directory))]
    average = subprocess.run(["sh", str(manual.HERE / "2_proc")], cwd=tmp_path,
                             text=True, capture_output=True, timeout=20)
    assert average.returncode == 1
    assert "manual NRG:" in average.stderr
    assert str(manual.HERE / "1/completed.json") in average.stderr
    assert not (manual.HERE / "processed").exists()
    assert not (tmp_path / "processed").exists()
