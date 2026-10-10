"""Bounded offline bundle tests; no real NRG or Mathematica is ever launched."""

import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

from NRG_comparisons.test4 import manual, nrg, physics


CASE = Path(manual.__file__).resolve().parent
REAL_LOOPER = os.environ.get("TEST4_REAL_LOOPER")
pytestmark = pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX shell/flock")
FAKE_TOOL = r'''
import json
import os
from pathlib import Path
import re
import sys

name = Path(__file__).name
cwd = Path.cwd()
if name in ("nrginit", "nrg"):
    import fcntl
    with (cwd.parent / ".manual.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            raise AssertionError("controller must hold the root lock throughout execution")
with open(os.environ["TEST4_EVENTS"], "a") as stream:
    stream.write(json.dumps(dict(tool=name, argv=sys.argv[1:], cwd=str(cwd),
                                env={k: os.environ.get(k) for k in
                                     ("NRG_WORKDIR", "TMPDIR", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                      "OPENBLAS_NUM_THREADS", "MKL_DYNAMIC", "OMP_DYNAMIC")})) + "\n")
failure = os.environ.get("TEST4_FAILURE", "")
if name == "looper":
    assert sys.argv[1:] == ["param.loop"]
    text = Path("param.loop").read_text()
    nz = int(re.search(r"#PRELUDE: \$Nz = ([0-9]+);", text)[1])
    for i in range(1, nz + 1):
        if failure == "missing-directory" and i == nz:
            continue
        folder = Path(str(i))
        folder.mkdir()
        body = "[param]" + text.split("[param]", 1)[1]
        body = body.replace("$Nz", str(nz)).replace("$z", str(i / nz))
        if i == nz:
            if failure == "wrong-z":
                body = body.replace("z=1.0", "z=0.123")
            if failure == "mixed-settings":
                body = body.replace("keep=200", "keep=201")
            if failure == "extra-directory":
                Path(str(nz + 1)).mkdir()
        if failure == "uniform-change":
            body = body.replace("WN=0", "WN=0.5")
        (folder / "param").write_text(body)
elif name == "nrginit":
    assert len(sys.argv) == 2 and Path(sys.argv[1]).name == "math"
    assert Path("model.m").is_file() and Path("modeloperators.m").is_file()
    if failure != "no-data":
        Path("data").write_text("synthetic chain\n")
    if failure != "no-success":
        print("Success!", flush=True)
    if failure == "init-model-change":
        with Path("model.m").open("a") as stream:
            stream.write("(* changed *)\n")
    if failure == "init-exit":
        sys.exit(17)
elif name == "nrg":
    if failure != "no-raw":
        Path("raw.h5").write_bytes(b"synthetic output, not HDF5\n")
    if failure != "no-done":
        Path("DONE").touch()
    if failure in ("param-change", "model-change", "root-model-change"):
        path = Path("param" if failure == "param-change" else "modeloperators.m")
        if failure == "root-model-change":
            path = Path("..") / path
        with path.open("a") as stream:
            stream.write("# changed\n")
    if failure == "nrg-exit":
        sys.exit(23)
else:
    raise AssertionError("fake math must never be invoked directly")
'''


@pytest.fixture
def tools(tmp_path, monkeypatch):
    folder = tmp_path / "fake tools"
    folder.mkdir()
    paths = {"events": tmp_path / "events.jsonl"}
    for name, variable in (("looper", "LOOPER"), ("nrginit", "NRGINIT"), ("math", "MATH"), ("nrg", "NRG")):
        path = folder / name
        path.write_text(f"#!{sys.executable}\n" + FAKE_TOOL)
        path.chmod(0o755)
        monkeypatch.setenv(variable, str(path))
        paths[name] = path
    monkeypatch.setenv("TEST4_EVENTS", str(paths["events"]))
    monkeypatch.setenv("PYTHON", sys.executable)
    monkeypatch.setenv("NRG_THREADS", "1")
    monkeypatch.delenv("TEST4_FAILURE", raising=False)
    return paths


@pytest.fixture
def bundle(tmp_path, tools):
    return manual.export(tmp_path / "exported bundle")


def events(tools):
    path = tools["events"]
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def edit_template(bundle, *, nz=None, **settings):
    path = bundle / "param.loop"
    text = path.read_text()
    if nz is not None:
        text = re.sub(r"^#PRELUDE:.*$", f"#PRELUDE: $Nz = {nz};", text, flags=re.M)
    for name, value in settings.items():
        text, count = re.subn(rf"^{name}=.*$", f"{name}={value}", text, flags=re.M)
        assert count == 1
    path.write_text(text)


def result(directory):
    """Physical longitudinal identities, with deliberately nonlinear twist variation."""
    deck = nrg.read_param(directory / "param")
    z = deck["numerical"]["z"]
    states = []
    for spin, energy in ((0, .2), (-1, 0.), (1, .3)):
        values = physics.state_values(.8, .1, spin * .1, .15, spin / 2)
        values.update(energy=energy + z*z, energy_paper=energy + z*z + .0625, current=z*z - .2)
        states.append(dict(twice_sz=spin, parity=spin % 2, root=0, degenerate=False, residual=None, values=values))
    return dict(physical=deck["physical"], numerical=deck["numerical"], states=states,
                coverage=physics.finish_states(states), diagnostics={"synthetic": True})


def test_export_defaults_no_calculations_and_portable_prepare(bundle, tools, tmp_path, monkeypatch):
    assert events(tools) == []
    assert {p.name for p in bundle.iterdir()} == {
        *manual.SOURCES, "bundle.json", "param.loop", "README.md", "1_run", "run_one", "2_proc"}
    assert manual.nrg is nrg and manual.physics is physics
    for name in manual.SOURCES:
        assert (bundle / name).read_bytes() == (CASE / name).read_bytes()
    for name in ("1_run", "run_one", "2_proc"):
        assert os.access(bundle / name, os.X_OK)
        assert '$(dirname "$0")/manual.py' in (bundle / name).read_text()
    text = (bundle / "param.loop").read_text()
    assert "#PRELUDE: $Nz = 8;\n#AUTOLOOP: true\n@ $i=1; $i<=$Nz; $i++\n!($z=$i/$Nz)" in text
    assert "Nz=$Nz\n" in text and "z=$z\n" in text
    readme = (bundle / "README.md").read_text()
    assert "UNCONVERGED" in readme and "energy_paper" in readme and "2e\\Delta" in readme
    assert "no wall-time or memory limit" in readme and "report" not in readme

    moved = tmp_path / "relocated bundle"
    shutil.copytree(bundle, moved)
    shutil.rmtree(bundle)
    monkeypatch.delenv("PYTHONPATH", raising=False)
    check = subprocess.run([sys.executable, "-B", "-S", str(moved / "manual.py"), "--help"],
                           cwd=tmp_path, capture_output=True, text=True, timeout=15, check=True)
    assert "export" in check.stdout  # Direct imports work without installed packages.
    subprocess.run([str(moved / "1_run")], cwd=tmp_path, timeout=15, check=True)
    manifest = manual.load_sweep(moved)
    assert [row["z"] for row in manifest["runs"]] == [i / 8 for i in range(1, 9)]
    assert manifest["parameters"]["param"] == {
        k: v for k, v in nrg.read_param(moved / "1/param")["param"].items() if k != "z"}
    for entry in manifest["runs"]:
        folder = moved / entry["directory"]
        deck = nrg.read_param(folder / "param")
        assert deck["physical"] == pytest.approx(physics.PHYSICAL)
        assert deck["numerical"] == {"lambda": 4., "z": entry["z"], "nmax": 20, "keep": 200,
                                      "keepmin": 40, "keepenergy": 14., "untruncated": False}
        assert deck["param"]["bcsgap"] == "-0.01"
        assert set(deck["extra"]) == {"nu", "Bz", "WN", "WS", "phi"}
        assert all(nrg.LITERAL.fullmatch(v) for v in deck["extra"].values())
        assert {p.name for p in folder.iterdir()} == {"param", *manual.MODELS}
        for name in manual.MODELS:
            assert manual.sha256(folder / name) == manifest["source_sha256"][name]
    assert [row["tool"] for row in events(tools)] == ["looper"]


def test_cli_export_and_default_destination(bundle, tmp_path, tools, monkeypatch):
    target = tmp_path / "cli export"
    subprocess.run([sys.executable, "-B", str(CASE / "manual.py"), "export", str(target)],
                   cwd=tmp_path, timeout=15, check=True)
    assert manual.source_hashes(target) == manual.source_hashes(bundle)
    monkeypatch.setattr(manual, "HERE", bundle)
    assert manual.export() == bundle / "manual_run"
    assert events(tools) == []


@pytest.mark.parametrize("kind", ["source", "existing-file", "existing-dir", "symlink", "missing-source"])
def test_export_refuses_source_existing_and_missing_files(bundle, tmp_path, tools, monkeypatch, kind):
    target = tmp_path / "new export"
    if kind == "source":
        target = CASE
    elif kind == "existing-file":
        target.write_text("keep me")
    elif kind == "existing-dir":
        target.mkdir()
    elif kind == "symlink":
        target.symlink_to(tmp_path / "nonexistent", target_is_directory=True)
    else:
        monkeypatch.setattr(manual, "HERE", bundle)
        (bundle / "modeloperators.m").unlink()
    with pytest.raises(ValueError):
        manual.export(target)
    if kind == "existing-file":
        assert target.read_text() == "keep me"
    if kind == "missing-source":
        assert not target.exists()
    assert events(tools) == []


@pytest.mark.parametrize("old,new", [
    ("#AUTOLOOP: true", "#AUTOLOOP: ./run_one"),
    ("#AUTOLOOP: true", "#AUTOLOOP: true\n#BEFORE: nrg"),
    ("#AUTOLOOP: true", "#AUTOLOOP: true\n#CONTINUE"),
    ("#PRELUDE: $Nz = 8;", "#PRELUDE: $Nz = 4+4;"),
    ("!($z=$i/$Nz)", "!($z=0.5)"), ("z=$z", "z=1"), ("Nz=$Nz", "Nz=8"),
    ("WN=0", "WN=!1/2"), ("phi=", "phi=!"), ("U=2", "U=!2"),
    ("WN=0", "WN=$ENV{HOME}"), ("Lambda=4", "Lambda=3.99"),
    ("WN=0", "WN=NaN"), ("WN=0", "WN=0\nunknown=1"),
    ("bcsgap=-0.01", "bcsgap1=-0.01\nbcsgap2=-0.01"),
])
def test_reject_template_before_looper(bundle, tools, old, new):
    path = bundle / "param.loop"
    assert old in path.read_text()
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError):
        manual.prepare(bundle)
    assert events(tools) == []
    assert not (bundle / "sweep.json").exists()


@pytest.mark.parametrize("failure", ["missing-directory", "extra-directory", "wrong-z", "mixed-settings", "uniform-change"])
def test_check_looper_output_against_validated_template(bundle, tools, monkeypatch, failure):
    monkeypatch.setenv("TEST4_FAILURE", failure)
    with pytest.raises(ValueError):
        manual.prepare(bundle)
    assert [row["tool"] for row in events(tools)] == ["looper"]
    assert not (bundle / "sweep.json").exists()


def test_literal_soc_inputs_and_no_overwrite_or_resume(bundle, tools):
    edit_template(bundle, nz=3, WN=.17, WS=-.09, nu=.6, phi=1.2)
    manifest = manual.prepare(bundle)
    assert manifest["parameters"]["extra"] == {"nu": "0.6", "Bz": "0.029999999999999999",
                                               "WN": "0.17", "WS": "-0.09", "phi": "1.2"}
    for entry in manifest["runs"]:
        assert entry["physical"]["wn"] == .17 and entry["physical"]["ws"] == -.09
    with pytest.raises(ValueError, match="existing"):
        manual.prepare(bundle, execute=True)
    assert [row["tool"] for row in events(tools)] == ["looper"]


@pytest.mark.parametrize("path", ["model.m", "modeloperators.m", "nrg.py", "physics.py", "manual.py"])
def test_source_changes_rejected_before_prepare(bundle, tools, path):
    with (bundle / path).open("a") as stream:
        stream.write("\n# changed\n")
    with pytest.raises(ValueError, match="source/model hashes changed"):
        manual.prepare(bundle)
    assert events(tools) == []


@pytest.mark.parametrize("path", ["param.loop", "model.m", "1/model.m", "2/modeloperators.m", "1/param"])
def test_frozen_inputs_rejected_before_execute_and_average(bundle, tools, path):
    manual.prepare(bundle)
    with (bundle / path).open("a") as stream:
        stream.write("\n# changed\n")
    with pytest.raises(ValueError, match="changed"):
        manual.run_one(bundle / "1")
    with pytest.raises(ValueError, match="changed"):
        manual.average(bundle)
    assert [row["tool"] for row in events(tools)] == ["looper"]


def test_execute_sequential_checked_receipts_and_cwd_wrapper(bundle, tools, monkeypatch, tmp_path):
    edit_template(bundle, nz=3)
    monkeypatch.setenv("NRG_THREADS", "4")
    monkeypatch.setenv("NRG_WORKDIR", "/unrelated/work")
    manual.prepare(bundle, execute=True)
    calls = events(tools)
    assert [row["tool"] for row in calls] == ["looper", "nrginit", "nrg", "nrginit", "nrg", "nrginit", "nrg"]
    for row in calls[1:]:
        directory = Path(row["cwd"])
        assert row["env"]["NRG_WORKDIR"] == row["env"]["TMPDIR"] == str(directory / "scratch")
        assert (directory / "scratch").is_dir()
        assert row["env"]["OMP_NUM_THREADS"] == row["env"]["MKL_NUM_THREADS"] == "4"
        assert row["env"]["OPENBLAS_NUM_THREADS"] == "1"
        assert row["env"]["OMP_DYNAMIC"] == row["env"]["MKL_DYNAMIC"] == "FALSE"
        receipt = json.loads((directory / "completed.json").read_text())
        assert receipt["source_sha256"] == manual.source_hashes(bundle)
        assert receipt["artifacts"] == {name: manual.sha256(directory / name) for name in manual.ARTIFACTS}
    with pytest.raises(ValueError, match="prior outputs"):
        manual.run_one(bundle / "1")

    other = manual.export(tmp_path / "other")
    edit_template(other, nz=1)
    subprocess.run([str(other / "1_run"), "--prepare-only"], cwd=tmp_path, timeout=15, check=True)
    subprocess.run([str(other / "run_one")], cwd=other / "1", timeout=15, check=True)
    assert (other / "1/completed.json").is_file()
    failed = subprocess.run([str(other / "run_one")], cwd=other, capture_output=True, text=True, timeout=15)
    assert failed.returncode != 0


@pytest.mark.parametrize("failure", ["no-data", "no-success", "no-raw", "no-done", "init-exit", "nrg-exit",
                                     "init-model-change", "param-change", "model-change", "root-model-change"])
def test_failure_never_receipted_or_continued(bundle, tools, monkeypatch, failure):
    edit_template(bundle, nz=2)
    monkeypatch.setenv("TEST4_FAILURE", failure)
    with pytest.raises((ValueError, subprocess.CalledProcessError)):
        manual.prepare(bundle, execute=True)
    assert not (bundle / "1/completed.json").exists()
    assert not (bundle / "2/started.json").exists()
    assert (bundle / "1/started.json").exists()
    calls = events(tools)
    if failure in {"no-data", "no-success", "init-exit", "init-model-change"}:
        assert [row["tool"] for row in calls] == ["looper", "nrginit"]
    with pytest.raises(ValueError):
        manual.run_one(bundle / "1")
    assert events(tools) == calls


@pytest.mark.parametrize("stale", ["data", "DONE", "raw.h5", "scratch", "started.json"])
def test_stale_outputs_refused(bundle, tools, stale):
    manual.prepare(bundle)
    (bundle / "1" / stale).write_text("keep me")
    with pytest.raises(ValueError, match="prior outputs"):
        manual.run_one(bundle / "1")
    assert [row["tool"] for row in events(tools)] == ["looper"]
    assert (bundle / "1" / stale).read_text() == "keep me"


def test_path_resolution_thread_gate_and_root_lock(bundle, tools, monkeypatch):
    monkeypatch.setenv("PATH", str(tools["nrg"].parent) + os.pathsep + os.environ["PATH"])
    monkeypatch.delenv("NRG")
    assert manual.executable("nrg", "NRG") == str(tools["nrg"])
    alias = tools["nrg"].parent / "nrg-alias"
    alias.symlink_to(tools["nrg"])
    monkeypatch.setenv("NRG", str(alias))
    assert manual.executable("nrg", "NRG") == str(alias)
    monkeypatch.setenv("NRG", "/nonexistent/nrg")
    with pytest.raises(ValueError, match="missing nrg"):
        manual.prepare(bundle, execute=True)
    monkeypatch.setenv("NRG", str(tools["nrg"]))
    for threads in ("0", "5", "1.5"):
        monkeypatch.setenv("NRG_THREADS", threads)
        with pytest.raises(ValueError):
            manual.prepare(bundle, execute=True)
    monkeypatch.setenv("NRG_THREADS", "1")
    with manual.lock(bundle):
        with pytest.raises(ValueError, match="another manual workflow"):
            manual.prepare(bundle)
        check = subprocess.run([str(bundle / "1_run")], capture_output=True, text=True, timeout=15)
        assert check.returncode and "another manual workflow" in check.stderr
    assert events(tools) == []


def test_nominal_average_all_observables_complete_grid(bundle, tools, monkeypatch):
    manual.prepare(bundle, execute=True)
    monkeypatch.setattr(nrg, "extract_result", result)
    before = events(tools)
    summary = manual.average(bundle)
    assert len(summary["means"]) == 3 * 14
    assert summary["excluded_primary_ids"] == {}
    assert summary["expected_three_subgap_states_every_twist"]
    assert not summary["continuum_completeness_established"]
    assert not summary["twist_min_max_are_error_bars"]
    for state_id in physics.PRIMARY.values():
        rows = {row["observable"]: row for row in summary["means"] if row["id"] == state_id}
        assert set(rows) == set(physics.OBSERVABLES)
        assert rows["energy_paper"]["mean"] - rows["energy"]["mean"] == pytest.approx(.0625)
        assert rows["current"]["mean"] == pytest.approx(sum((i / 8)**2 for i in range(1, 9)) / 8 - .2)
        assert {row["nz"] for row in rows.values()} == {8}
    with (bundle / "processed/per_z.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 24 and all(row["primary"] == "True" for row in rows)
    assert json.loads((bundle / "processed/summary.json").read_text()) == summary
    assert events(tools) == before
    with pytest.raises(FileExistsError):
        manual.average(bundle)


@pytest.mark.parametrize("issue", ["missing-directory", "extra-directory", "missing-receipt", "receipt-hash",
                                   "receipt-model-hash", "raw-hash", "missing-raw", "symlink-model", "mixed-model"])
def test_incomplete_or_mixed_grid_fails_before_extraction(bundle, tools, monkeypatch, issue):
    edit_template(bundle, nz=2)
    manual.prepare(bundle, execute=True)
    folder = bundle / "2"
    if issue == "missing-directory":
        shutil.rmtree(folder)
    elif issue == "extra-directory":
        (bundle / "3").mkdir()
    elif issue == "missing-receipt":
        (folder / "completed.json").unlink()
    elif issue in {"receipt-hash", "receipt-model-hash"}:
        path = folder / "completed.json"
        receipt = json.loads(path.read_text())
        if issue == "receipt-hash":
            receipt["param_sha256"] = "incorrect"
        else:
            receipt["source_sha256"]["model.m"] = "incorrect"
        path.write_text(json.dumps(receipt))
    elif issue == "missing-raw":
        (folder / "raw.h5").unlink()
    elif issue == "symlink-model":
        (folder / "model.m").unlink()
        (folder / "model.m").symlink_to(bundle / "model.m")
    else:
        (folder / ("raw.h5" if issue == "raw-hash" else "model.m")).write_text("changed")

    def forbidden(*args):
        pytest.fail("must check the complete grid before extraction")

    monkeypatch.setattr(nrg, "extract_result", forbidden)
    before = events(tools)
    with pytest.raises((OSError, ValueError)):
        manual.average(bundle)
    assert not (bundle / "processed").exists()
    assert events(tools) == before


def test_missing_degenerate_and_extra_states_not_averaged_away(bundle, monkeypatch):
    edit_template(bundle, nz=2)
    manual.prepare(bundle, execute=True)

    def varied(directory):
        extracted = result(directory)
        states = extracted["states"]
        if directory.name == "2":
            states[:] = [row for row in states if row["twice_sz"] != 1]
            states[0]["degenerate"] = True
        extra = dict(states[0], values=dict(states[0]["values"]), root=3)
        extra["values"]["energy"] += .1
        extra["values"]["energy_paper"] += .1
        states.append(extra)
        extracted["coverage"] = physics.finish_states(states)
        return extracted

    monkeypatch.setattr(nrg, "extract_result", varied)
    summary = manual.average(bundle)
    assert set(summary["excluded_primary_ids"]) == {"even", "odd_up"}
    assert {row["id"] for row in summary["means"]} == {"odd_down"}
    assert not summary["expected_three_subgap_states_every_twist"]
    assert summary["per_z"][1]["result"]["coverage"]["missing_primary_subgap_ids"] == ["odd_up"]
    assert all(row["result"]["coverage"]["extra_subgap_ids"] == ["sz+0-r3"] for row in summary["per_z"])
    with (bundle / "processed/per_z.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    extras = [row for row in rows if row["extra_subgap"] == "True"]
    assert len(extras) == 2 and extras[1]["degenerate"] == "True"


def test_actual_h5_extraction_through_processing_wrapper(bundle, tmp_path):
    h5py = pytest.importorskip("h5py")
    edit_template(bundle, nz=2, Nmax=1)
    manifest = manual.prepare(bundle, execute=True)
    for entry in manifest["runs"]:
        folder, z = bundle / entry["directory"], entry["z"]
        (folder / "data").write_text("z\n1\n0\n999\n1\n0\n777\n1\n0\n0\n1\n0\n0\n"
                                     "Z\n1\n-1\n-1\n1\n-1\n-1\n1\n0\n0\n1\n0\n0\n")
        with h5py.File(folder / "raw.h5", "w") as h5:
            h5["stats/GS_energy"] = -4. + z*z
            h5["stats/1/energy_offset"] = -4. + z*z
            for sector, excitation in (("0", 0.), ("1", .2), ("2", .3)):
                h5[f"1/eigen/{sector}/absenergy_zero"] = [excitation]
                h5[f"1/eigen/{sector}/absenergy"] = [-4. + z*z + excitation]
                for name, value in (("n_d", .8), ("n_d_ud", .1), ("SZd", 0.), ("SZ2d", .15), ("asqJ", z*z)):
                    h5[f"1/s/{name}/{sector}/{sector}/matrix"] = [[value]]
        # Synthetic fixture provenance only; real receipts are produced by run-one.
        path = folder / "completed.json"
        receipt = json.loads(path.read_text())
        receipt["artifacts"] = {name: manual.sha256(folder / name) for name in manual.ARTIFACTS}
        path.write_text(json.dumps(receipt))
    subprocess.run([str(bundle / "2_proc")], cwd=tmp_path, timeout=15, check=True)
    summary = json.loads((bundle / "processed/summary.json").read_text())
    rows = {row["observable"]: row for row in summary["means"] if row["id"] == "odd_down"}
    assert set(rows) == set(physics.OBSERVABLES)
    assert rows["energy_paper"]["mean"] == pytest.approx(.625)
    assert rows["energy"]["mean"] == pytest.approx(.625 - .0625)
    assert rows["current"]["mean"] == pytest.approx(.625)


@pytest.mark.skipif(not REAL_LOOPER, reason="opt in with TEST4_REAL_LOOPER; preparation only")
def test_real_looper_preparation_only_in_temporary_bundle(bundle, tools, monkeypatch):
    monkeypatch.setenv("LOOPER", REAL_LOOPER)
    edit_template(bundle, nz=3)
    manifest = manual.prepare(bundle)
    assert len(manifest["runs"]) == 3
    for i in range(1, 4):
        folder = bundle / str(i)
        assert {path.name for path in folder.iterdir()} == {"param", *manual.MODELS}
        assert nrg.read_param(folder / "param")["numerical"]["z"] == pytest.approx(i / 3, abs=1e-14)
    assert events(tools) == []
