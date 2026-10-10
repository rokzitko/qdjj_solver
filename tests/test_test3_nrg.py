"""Bounded SPU1 extraction/controller checks; no NRG or Mathematica execution."""

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


MANUAL = Path(__file__).resolve().parents[1] / "NRG_comparisons/test3/manual_run"
NUMERICAL = dict(lambda_=2., z=.5, nmax=1, keep=128, keepmin=0, keepenergy=-1.)
NUMERICAL["lambda"] = NUMERICAL.pop("lambda_")


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.fixture
def adapter(monkeypatch):
    physics = module(MANUAL / "physics.py", "test3_physics")
    monkeypatch.setitem(sys.modules, "physics", physics)
    return module(MANUAL / "nrg.py", "test3_nrg")


def synthetic(directory, *, scale=1., z=.5, levels=None):
    """Independently populate the documented SPU1 layout, with unsorted roots."""
    np = pytest.importorskip("numpy")
    h5py = pytest.importorskip("h5py")
    # Use the literal manual deck, not the renderer under test.
    template = (MANUAL / "param.loop").read_text()
    param = "[param]" + template.split("[param]", 1)[1]
    param = param.replace("$Nz", "2").replace("$z", str(z))
    updates = {"Nmax": 1, "keep": 128, "keepmin": 0, "keepenergy": -1,
               "bandrescale": 100 * scale, "U": 2 * scale, "Gamma": .4 * scale,
               "delta": .5 * scale, "B": .03 * scale, "T": 1e-5 * scale}
    param = "\n".join(f"{line.split('=')[0]}={updates[line.split('=')[0]]:.17g}"
                      if "=" in line and line.split("=")[0] in updates else line for line in param.splitlines()) + "\n"
    (directory / "param").write_text(param)
    levels = levels or {"1": [2.4, .3, 1., .7], "0": [1.5, 0., 1.5], "2": [.4, 1.8],
                        "-1": [1.6, 2.], "3": [1.7], "4": []}
    with h5py.File(directory / "raw.h5", "w") as h5:
        h5["stats/GS_energy"] = -12. * scale
        h5["stats/1/energy_offset"] = -12. * scale
        for sector, energies in levels.items():
            h5[f"1/eigen/{sector}/absenergy_zero"] = np.array(energies) * scale
            h5[f"1/eigen/{sector}/absenergy"] = (np.array(energies) - 12.) * scale
            if not energies:
                continue
            charge = .8 + np.arange(len(energies)) * .1
            md = {"1": -.1, "0": -.25, "2": .25}.get(sector, .125)
            for name, values in {"n_d": charge, "n_d_ud": np.full(len(energies), .1),
                                 "SZd": np.full(len(energies), md), "SZ2d": (charge - .2) / 4}.items():
                h5[f"1/s/{name}/{sector}/{sector}/matrix"] = np.diag(values)
    # Nonzero onsite trace catches an incorrect uncentered vacuum subtraction.
    lines = ["#!9", "# symtype SPU1", "1 3 1", "e", "-0.5"]
    for tag, arrays in (("z", ([2., 999.], [1., 3.])), ("Z", ([1., 1.], [0., 0.]))):
        lines.append(tag)
        for values in arrays:
            lines += ["1", *(f"{value * scale:.17e}".replace("e", "D") for value in values)]
    (directory / "data").write_text("\n".join(lines) + "\n")
    return h5py


def test_render_source_conventions_and_manual_defaults(adapter, tmp_path):
    physics = sys.modules["physics"]
    deck = adapter.render_param(physics.PHYSICAL, NUMERICAL)
    (tmp_path / "param").write_text(deck)
    parsed = adapter.read_param(tmp_path / "param")
    p = parsed["param"]
    assert (p["model"], p["variant"], p["symtype"], p["ops"]) == ("SIAM", "MAGFIELD", "SPU1", "n_d n_d_ud SZd SZ2d")
    assert float(parsed["extra"]["B"]) == .03 and "B" not in p
    assert (float(p["bandrescale"]), float(p["bcsgap"]), float(p["delta"])) == (100., .01, .5)
    assert p["keepall"] == "" and p["Ninit"] == "0" and p["tri"] == "old"
    assert not {"spect", "sigma_d", "polarized", "Tmin"} & p.keys()
    untruncated = adapter.render_param(physics.PHYSICAL, dict(NUMERICAL, nmax=3, untruncated=True))
    assert "keepall=0,1,2\n" in untruncated
    text = (MANUAL / "param.loop").read_text()
    (tmp_path / "param").write_text("[param]" + text.split("[param]", 1)[1].replace("$Nz", "8").replace("$z", "0.125"))
    manual = adapter.read_param(tmp_path / "param")
    assert {key: manual["param"][key] for key in ("Nz", "Lambda", "Nmax", "keep", "keepmin", "keepenergy")} == {
        "Nz": "8", "Lambda": "2", "Nmax": "42", "keep": "2000", "keepmin": "400", "keepenergy": "14"}
    expected = adapter.render_param(physics.PHYSICAL, dict(NUMERICAL, z=.125, nmax=42, keep=2000, keepmin=400, keepenergy=14))
    (tmp_path / "param").write_text(expected)
    rendered = adapter.read_param(tmp_path / "param")
    for key, value in rendered["param"].items():
        try:
            assert float(manual["param"][key]) == float(value)
        except ValueError:
            assert manual["param"][key] == value
    assert float(rendered["extra"]["B"]) == float(manual["extra"]["B"])


@pytest.mark.parametrize("change", [dict(field=float("nan")), dict(geometry="junction"), dict(gap=0),
                                    dict(gamma=-1), dict(rho_ws=0)])
def test_renderer_rejects_unsupported_physics(adapter, change):
    with pytest.raises(ValueError):
        adapter.render_param(dict(sys.modules["physics"].PHYSICAL, **change), NUMERICAL)


@pytest.mark.parametrize("change", [dict(nmax=1.), dict(nmax=999), dict(keepmin=129), dict(untruncated=1),
                                    {"lambda": 1.5}, dict(z=0), dict(polarized=False)])
def test_renderer_rejects_unsupported_numerics(adapter, change):
    with pytest.raises(ValueError):
        adapter.render_param(sys.modules["physics"].PHYSICAL, dict(NUMERICAL, **change))


@pytest.mark.parametrize("old,new", [("[extra]", "[other]"), ("B=", "Bz="), ("SPU1", "SPSU2"),
                                     ("tri=old", "tri=cpp"), ("Ninit=0", "Ninit=1"),
                                     ("h5ops=true", "h5ops=false"), ("keepall=", "keepall=even"),
                                     ("[extra]", "polarized=true\n[extra]"),
                                     ("[extra]", "spect=sigma_d-sigma_d\n[extra]")])
def test_read_deck_rejects_model_and_output_changes(adapter, tmp_path, old, new):
    deck = adapter.render_param(sys.modules["physics"].PHYSICAL, NUMERICAL).replace(old, new)
    (tmp_path / "param").write_text(deck)
    with pytest.raises(ValueError):
        adapter.read_param(tmp_path / "param")


@pytest.mark.parametrize("scale", [1., .01])
@pytest.mark.parametrize("scalar_vector", [False, True])
def test_extraction_signed_sectors_sorted_roots_vacuum_and_units(adapter, tmp_path, scale, scalar_vector):
    h5py = synthetic(tmp_path, scale=scale)
    if scalar_vector:
        with h5py.File(tmp_path / "raw.h5", "a") as h5:
            for path in ("stats/GS_energy", "stats/1/energy_offset"):
                value = h5[path][()]
                del h5[path]
                h5[path] = [value]
    result = adapter.extract_result(tmp_path)
    assert result["physical"] == pytest.approx(sys.modules["physics"].PHYSICAL)
    assert result["numerical"] == dict(NUMERICAL, untruncated=False)
    diagnostics = result["diagnostics"]
    assert diagnostics["deck_gap_scale"] == scale
    assert diagnostics["final_h5_iteration"] == 1
    # Analytic eigenvalues of [[1,2],[2,3]], independent of numpy.linalg.
    vacuum = -math.hypot(2 - math.sqrt(5), 1) - math.hypot(2 + math.sqrt(5), 1)
    assert diagnostics["analytic_bath_vacuum_over_gap"] == pytest.approx(vacuum, abs=1e-13)
    assert diagnostics["wilson_chain"]["xi_over_gap"] == [2., 999.]
    assert diagnostics["wilson_chain"]["impurity_hopping_over_gap"] == pytest.approx(math.sqrt(80 / math.pi))
    states = {row["id"]: row for row in result["states"]}
    assert states["even"]["values"]["energy"] == pytest.approx(-11.7 - vacuum - .5, abs=1e-13)
    assert states["odd_down"]["values"]["energy"] == pytest.approx(-12 - vacuum - .5, abs=1e-13)
    assert states["even"]["excitation"] == pytest.approx(.3, abs=1e-13)
    assert states["even"]["values"]["P1"] == pytest.approx(.7, abs=1e-14)
    assert states["even"]["values"]["m_dot"] == -.1
    assert states["odd_down"]["values"]["m_dot"] == -.25
    assert states["odd_up"]["values"]["m_dot"] == .25
    assert states["sz-2-r0"]["twice_sz"] == -2
    assert states["sz-2-r1"]["root"] == 1 and not states["sz-2-r1"]["subgap"]
    assert states["sz-1-r1"]["degenerate"]  # Its partner is an unretained high root.
    assert "sz+0-r2" not in states  # Strict threshold: E-Eg=1 is not subgap.
    assert diagnostics["all_sector_spectra"]["1"]["indices_in_energy_order"] == [1, 3, 2, 0]
    assert "sz+0-r1" in result["coverage"]["subgap_ids"]
    assert result["coverage"]["diagonal_observables_unique"]
    for row in states.values():
        values, M = row["values"], row["twice_sz"] / 2
        assert set(values) == set(sys.modules["physics"].OBSERVABLES)
        assert row["parity"] == row["twice_sz"] % 2 and row["residual"] is None
        assert values["P_up"] - values["P_down"] == pytest.approx(2 * values["m_dot"], abs=1e-14)
        assert values["m_bath"] == pytest.approx(M - values["m_dot"], abs=1e-14)
        assert values["Czz"] == pytest.approx(M * values["m_dot"] - values["dot_z2"], abs=1e-14)
        assert values["Czz_connected"] == pytest.approx(values["m_dot"] ** 2 - values["dot_z2"], abs=1e-14)
    json.dumps(result, allow_nan=False)
    assert not (tmp_path / "DONE").exists()


def test_bandwidth_renderer_and_nonunit_physical_gap(adapter, tmp_path):
    physical = {k: v if k == "geometry" else 2.5 * v for k, v in sys.modules["physics"].PHYSICAL.items()}
    deck = adapter.render_param(physical, NUMERICAL, units="bandwidth")
    synthetic(tmp_path, scale=.01)
    (tmp_path / "param").write_text(deck)
    result = adapter.extract_result(tmp_path)
    assert result["physical"] == pytest.approx(sys.modules["physics"].PHYSICAL)
    assert result["diagnostics"]["deck_gap_scale"] == .01


def test_global_ground_missing_primary_and_subgap_degeneracy(adapter, tmp_path):
    synthetic(tmp_path, levels={"-1": [.2, 0.], "1": [.3, .3, .8], "0": [.1, 1.4]})
    result = adapter.extract_result(tmp_path)
    assert result["coverage"]["ground_ids"] == ["sz-2-r0"]
    assert "odd_up" not in result["coverage"]["subgap_ids"]
    assert {"even", "sz+0-r1"} <= set(result["coverage"]["degenerate_subgap_ids"])
    assert not result["coverage"]["diagonal_observables_unique"]
    assert "sz+0-r2" in result["coverage"]["subgap_ids"]  # All subgap roots, not just two.


@pytest.mark.parametrize("path,value", [
    ("stats/GS_energy", [0., 0.]), ("stats/1/energy_offset", -11.),
    ("1/eigen/1/absenergy", [-12., -12., -12., -12.]),
    ("1/eigen/1/absenergy_zero", [float("nan")] * 4),
    ("1/s/SZd/1/1/matrix", [[.2]]),
    ("1/s/n_d/1/1/matrix", [[float("nan") if i == 0 and j == 1 else 0 for j in range(4)] for i in range(4)]),
    ("1/s/n_d_ud/1/1/matrix", [[.9 if i == j else 0 for j in range(4)] for i in range(4)]),
    ("1/s/SZ2d/1/1/matrix", [[.5 if i == j else 0 for j in range(4)] for i in range(4)]),
    ("1/s/SZd/1/1/matrix-imag", [[.01 if i == j else 0 for j in range(4)] for i in range(4)]),
])
def test_invalid_raw_data_rejected(adapter, tmp_path, path, value):
    h5py = synthetic(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        if path in h5:
            del h5[path]
        h5[path] = value
    with pytest.raises(ValueError):
        adapter.extract_result(tmp_path)


@pytest.mark.parametrize("change", ["sector", "matrix", "final", "data_count", "pairing", "hopping"])
def test_missing_and_inconsistent_inputs_rejected(adapter, tmp_path, change):
    h5py = synthetic(tmp_path)
    if change.startswith("data") or change in ("pairing", "hopping"):
        path = tmp_path / "data"
        text = path.read_text()
        if change == "data_count":
            text = text.replace("z\n1\n", "z\n2\n")
        elif change == "pairing":
            text = text.replace("Z\n1\n1.00000000000000000D+00", "Z\n1\n2.00000000000000000D+00")
        else:
            text = text.replace("0.00000000000000000D+00", "1.00000000000000000D+00")
        path.write_text(text)
    else:
        with h5py.File(tmp_path / "raw.h5", "a") as h5:
            if change == "sector":
                h5.move("1/eigen/-1", "1/eigen/0.5")
            elif change == "matrix":
                del h5["1/s/SZd/0/0/matrix"]
            else:
                h5.create_group("10")
    with pytest.raises(ValueError):
        adapter.extract_result(tmp_path)


@pytest.fixture
def controller(adapter, tmp_path, monkeypatch):
    pytest.importorskip("fcntl", reason="standalone execution controller uses POSIX file locks")
    root = tmp_path / "copied manual"
    root.mkdir()
    for name in ("nrg.py", "physics.py", "process.py", "param.loop", "1_run", "2_proc", "run_one"):
        shutil.copy2(MANUAL / name, root / name)
    monkeypatch.setitem(sys.modules, "nrg", adapter)
    result = module(root / "process.py", "test3_process")
    for variable in ("LOOPER", "NRGINIT", "NRG", "MATH"):
        monkeypatch.setenv(variable, sys.executable)
    monkeypatch.setenv("NRG_THREADS", "1")
    return result


@pytest.fixture
def mock_run(controller, monkeypatch):
    calls = []

    def run(argv, *, cwd, check, env=None, stdout=None, stderr=None):
        cwd = Path(cwd)
        calls.append(dict(argv=argv, cwd=cwd, env=env))
        assert check
        if argv[1:] == ["param.loop"]:
            text = (cwd / "param.loop").read_text()
            nz = int(text.split("$Nz = ", 1)[1].split(";", 1)[0])
            for i in range(1, nz + 1):
                folder = cwd / str(i)
                folder.mkdir()
                param = "[param]" + text.split("[param]", 1)[1].replace("$Nz", str(nz)).replace("$z", str(i / nz))
                (folder / "param").write_text(param)
        elif len(argv) == 2:
            assert argv[1] == sys.executable
            (cwd / "data").write_text("synthetic initializer output\n")
            stdout.write("Success!\n")
        else:
            (cwd / "raw.h5").write_text("synthetic solver output\n")
            (cwd / "DONE").touch()
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(controller.subprocess, "run", run)
    return calls


@pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX flock")
def test_prepare_only_and_mock_runner_receipts(controller, mock_run):
    controller.main(["run", "--prepare-only"])
    root = controller.HERE
    manifest = controller.load_sweep(root)
    assert len(manifest["runs"]) == 8 and len(mock_run) == 1
    assert [row["z"] for row in manifest["runs"]] == [i / 8 for i in range(1, 9)]
    with pytest.raises(ValueError, match="existing"):
        controller.prepare(root)
    controller.run_one(root / "1")
    assert len(mock_run) == 3
    assert mock_run[1]["argv"] == [sys.executable, sys.executable]
    assert mock_run[2]["argv"] == [sys.executable]
    assert mock_run[2]["env"]["OMP_NUM_THREADS"] == "1"
    receipt = json.loads((root / "1/completed.json").read_text())
    assert receipt["artifacts"] == {name: controller.sha256(root / "1" / name) for name in controller.ARTIFACTS}
    with pytest.raises(ValueError, match="prior outputs"):
        controller.run_one(root / "1")
    with pytest.raises(OSError):
        controller.average(root)  # No receipt for the other seven twists.
    assert not (root / "processed").exists()


@pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX flock")
@pytest.mark.parametrize("old,new", [("#AUTOLOOP: true", "#AUTOLOOP: nrg"),
                                     ("#AUTOLOOP: true", "#AUTOLOOP: true\n#AFTER: nrg"),
                                     ("!($z=$i/$Nz)", "!(system('nrg'))"),
                                     ("keep=2000", 'keep=2000"; system("nrg"); #')])
def test_template_cannot_launch_external_solver(controller, mock_run, old, new):
    path = controller.HERE / "param.loop"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match="param.loop"):
        controller.prepare(controller.HERE)
    assert not mock_run


@pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX flock")
def test_explicit_tools_failure_and_mutated_parameters(controller, mock_run, monkeypatch):
    monkeypatch.delenv("LOOPER")
    with pytest.raises(ValueError, match="explicit absolute"):
        controller.prepare(controller.HERE)
    monkeypatch.setenv("LOOPER", sys.executable)
    controller.prepare(controller.HERE)

    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(17, args[0])

    monkeypatch.setattr(controller.subprocess, "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        controller.run_one(controller.HERE / "1")
    assert (controller.HERE / "1/started.json").exists()
    assert not (controller.HERE / "1/completed.json").exists()
    path = controller.HERE / "8/param"
    path.write_text(path.read_text() + "# changed\n")
    with pytest.raises(ValueError, match="parameters changed"):
        controller.run_one(controller.HERE / "2")


@pytest.mark.skipif(os.name != "posix", reason="manual workflow uses POSIX flock")
@pytest.mark.parametrize("issue", [None, "degenerate", "missing", "above_gap"])
def test_checked_averages_flag_extra_missing_degenerate(controller, mock_run, issue):
    root = controller.HERE
    path = root / "param.loop"
    path.write_text(path.read_text().replace("$Nz = 8;", "$Nz = 2;").replace("Nmax=42", "Nmax=1")
                    .replace("keep=2000", "keep=128").replace("keepmin=400", "keepmin=0").replace("keepenergy=14", "keepenergy=-1"))
    manifest = controller.prepare(root)
    for entry in manifest["runs"]:
        folder = root / entry["directory"]
        param = (folder / "param").read_text()
        levels = None
        if entry["directory"] == "2":
            levels = {"1": [.25, .8], "0": [0., 1.4], "2": [.6, 1.7]}
            if issue == "degenerate":
                levels["2"] = [.6, .6]
            elif issue == "missing":
                del levels["2"]
            elif issue == "above_gap":
                levels["2"] = [1.1, 1.7]
        synthetic(folder, z=entry["z"], levels=levels)
        (folder / "param").write_text(param)  # Retain the prepared bytes.
        for name in ("DONE", "nrginit.log", "nrg.log"):
            (folder / name).write_text("synthetic\n")
        controller.write_json(folder / "completed.json", {"param_sha256": entry["param_sha256"],
                              "artifacts": {name: controller.sha256(folder / name) for name in controller.ARTIFACTS}})
    summary = controller.average(root)
    assert not summary["continuum_target_established"] and not summary["expected_subgap_coverage_every_twist"]
    assert summary["reference_comparison"] is None
    assert summary["per_z"][0]["coverage"]["extra_subgap_ids"] == ["sz+0-r1"]
    assert ("odd_up" in summary["excluded_primary_ids"]) == bool(issue)
    rows = {(row["id"], row["observable"]): row for row in summary["averages"]}
    assert rows["even", "energy"]["mean"] == pytest.approx(math.fsum(
        next(row for row in record["states"] if row["id"] == "even")["values"]["energy"] for record in summary["per_z"]) / 2)
    assert (("odd_up", "energy") in rows) == (issue is None)
    with (root / "processed/per_z.csv").open() as stream:
        table = list(csv.DictReader(stream))
    assert any(row["id"] == "sz-2-r1" and row["subgap"] == "False" for row in table)
    assert len(table) == sum(len(record["states"]) for record in summary["per_z"])
    (root / "1/DONE").write_text("changed")
    with pytest.raises(ValueError, match="output changed"):
        controller.average(root)


def test_copyable_without_repository_or_optional_imports(tmp_path):
    for name in ("physics.py", "nrg.py", "process.py"):
        shutil.copy2(MANUAL / name, tmp_path / name)
    environment = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}
    code = "import nrg, physics, process; print(nrg.render_param(physics.PHYSICAL, " + repr(NUMERICAL) + "))"
    result = subprocess.run([sys.executable, "-B", "-S", "-c", code], cwd=tmp_path, env=environment,
                            capture_output=True, text=True, check=True)
    assert "symtype=SPU1" in result.stdout


def test_parent_can_import_adapter_as_package(tmp_path, monkeypatch):
    from NRG_comparisons.test3.manual_run import nrg, process

    monkeypatch.delitem(sys.modules, "physics", raising=False)
    synthetic(tmp_path)
    assert nrg.extract_result(tmp_path)["physical"]["field"] == .03
    assert process.extract_result is nrg.extract_result


@pytest.mark.skipif(os.name != "posix" or not os.environ.get("TEST3_REAL_LOOPER"),
                    reason="optional real looper input generation, never NRG")
def test_real_looper_only_writes_inputs(controller, monkeypatch):
    monkeypatch.setenv("LOOPER", os.environ["TEST3_REAL_LOOPER"])
    manifest = controller.prepare(controller.HERE)
    assert len(manifest["runs"]) == 8
    for entry in manifest["runs"]:
        folder = controller.HERE / entry["directory"]
        assert {path.name for path in folder.iterdir()} == {"param"}
        assert float(controller.read_param(folder / "param")["param"]["z"]) == entry["z"]
