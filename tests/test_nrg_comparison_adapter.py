"""Offline adapter tests: literal decks and synthetic NRG Ljubljana HDF5."""

import builtins
from configparser import ConfigParser
import importlib.util
import json
import math
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "NRG_comparisons/test1/nrg.py"
SPEC = importlib.util.spec_from_file_location("test1_nrg_adapter", PATH)
nrg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(nrg)
PHYSICAL = dict(gap=1., bandwidth=100., u=2., gamma=.4, detuning=0., field=0., geometry="single")
NUMERICAL = {"lambda": 2., "z": .5, "nmax": 3, "keep": 128, "keepenergy": -1., "keepmin": 0}


def parse(text):
    parser = ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read_string(text)
    return parser["param"]


@pytest.mark.parametrize("units,unit", [("gap", 1.), ("bandwidth", 100.)])
def test_literal_deck_units_and_static_operator_selection(units, unit):
    physical = dict(PHYSICAL, detuning=.3)
    text = nrg.render_param(physical, NUMERICAL, units)
    p = parse(text)
    assert "!" not in text and "*" not in text
    assert p["model"] == "SIAM" and p["variant"] == ""
    assert p["symtype"] == "SPSU2" and p["band"] == "flat"
    for key, value in {"bandrescale": 100 / unit, "bcsgap": .01, "U": 2 / unit,
                       "Gamma": .4 / unit, "delta": .3 / unit, "dumpEscale": 1 / unit,
                       "T": 1e-5 / unit}.items():
        assert float(p[key]) == pytest.approx(value)
    assert float(p["bandrescale"]) * float(p["bcsgap"]) == pytest.approx(1 / unit)
    assert p["ops"] == "n_d n_d_ud sigma_d"
    assert p["spect"] == "sigma_d-sigma_d"
    for key in ("finite", "h5raw", "h5last", "h5ops", "lastall"):
        assert p.getboolean(key)
    for key in ("dm", "fdm", "broaden", "savebins", "h5all", "h5vectors", "h5U"):
        assert not p.getboolean(key)
    assert p["keepall"] == "" and p["Nmax"] == "3"
    assert "Tmin" not in p and "B" not in p and "specd" not in p
    assert p["keep"] == "128" and p["keepmin"] == "0" and p["keepenergy"] == "-1"
    assert physical["detuning"] == .3 and NUMERICAL["keep"] == 128


def test_clean_and_untruncated_controls_preserve_bath():
    numerical = dict(NUMERICAL, untruncated=True, keep=2, keepenergy=8.)
    impurity = parse(nrg.render_param(PHYSICAL, numerical))
    clean = parse(nrg.render_param(PHYSICAL, numerical, clean=True))
    assert clean["model"] == "CLEAN" and clean["variant"] == ""
    assert clean["ops"] == clean["spect"] == ""
    assert not clean.getboolean("finite")
    assert clean["keepall"] == "0,1,2"  # Not 1,2,3 or an unsupported 'all'.
    for key in ("bandrescale", "bcsgap", "Gamma", "Nmax", "Lambda", "z", "discretization"):
        assert impurity[key] == clean[key]
    assert parse(nrg.render_param(dict(PHYSICAL, u=0., gamma=0.), numerical))["Gamma"] == "0"


@pytest.mark.parametrize("where,key,value,match", [
    ("numerical", "lambda", 1.799999, "Lambda < 1.8"),
    ("numerical", "lambda", float("nan"), "finite"),
    ("numerical", "z", 0., "z must"),
    ("numerical", "z", 1.1, "z must"),
    ("numerical", "nmax", 2.5, "integer"),
    ("numerical", "nmax", 999, "998"),
    ("numerical", "keep", True, "integer"),
    ("numerical", "keepmin", 129, "keepmin"),
    ("numerical", "keep", 1, "keep"),
    ("numerical", "untruncated", "false", "boolean"),
    ("physical", "field", .1, "field=0"),
    ("physical", "geometry", "junction", "single"),
    ("physical", "gap", 0., "positive"),
    ("physical", "bandwidth", float("inf"), "finite"),
    ("physical", "gamma", -.4, "nonnegative"),
])
def test_reject_unsupported_inputs(where, key, value, match):
    physical, numerical = dict(PHYSICAL), dict(NUMERICAL)
    (physical if where == "physical" else numerical)[key] = value
    with pytest.raises(ValueError, match=match):
        nrg.render_param(physical, numerical)


def test_no_hidden_numerical_defaults_or_hdf5_import(monkeypatch):
    original = builtins.__import__

    def without_h5py(name, *args, **kwargs):
        if name == "h5py":
            raise ImportError("HDF5 deliberately unavailable")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_h5py)
    module = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(module)
    assert module.render_param(PHYSICAL, dict(NUMERICAL, **{"lambda": 1.8}))
    for key in NUMERICAL:
        missing = dict(NUMERICAL)
        del missing[key]
        with pytest.raises(KeyError, match=key):
            module.render_param(PHYSICAL, missing)
    with pytest.raises(ValueError, match="units"):
        module.render_param(PHYSICAL, NUMERICAL, "unknown")
    with pytest.raises(ValueError, match="clean"):
        module.render_param(PHYSICAL, NUMERICAL, clean="false")


def make_output(work, units="gap", *, clean=False, physical=PHYSICAL, higher_ground=False):
    h5py = pytest.importorskip("h5py")
    gap_in_units = physical["gap"] / physical[units]
    (work / "param").write_text(nrg.render_param(physical, NUMERICAL, units, clean=clean))
    # Nonzero higher-spin ground tests that S/D are not silently reset to zero.
    levels = {"1": [.6, 1.1], "2": [.4, .1], "3": [0.]} if higher_ground else {"1": [.3, 1.1], "2": [.4, 0.], "3": [1.2]}
    with h5py.File(work / "raw.h5", "w") as h5:
        h5["stats/GS_energy"] = -120. * gap_in_units
        h5.create_group("2")  # Numeric shell ordering, not HDF5 iteration order.
        for step in (1, 2, 3):
            for key, value in {"energyscale": 10 / step, "scale": 10 / step, "Teff": 10 / step,
                               "abs_Egs": -1., "energy_offset": -117. - step}.items():
                h5[f"stats/{step}/{key}"] = value * gap_in_units
        for sector, energies in levels.items():
            h5[f"3/eigen/{sector}/absenergy_zero"] = [e * gap_in_units for e in energies]
            h5[f"3/eigen/{sector}/absenergy"] = [(e - 120.) * gap_in_units for e in energies]
        if not clean:
            for sector in ("1", "2"):
                h5[f"3/s/n_d/{sector}/{sector}/matrix"] = [[1.0, 0.], [0., .9]]
                h5[f"3/s/n_d_ud/{sector}/{sector}/matrix"] = [[.2, 0.], [0., .1]]
                h5[f"3/s/n_d/{sector}/{sector}/matrix-imag"] = [[0., 0.], [0., 0.]]
            h5["3/t/sigma_d/2/2/matrix"] = [[.1, 0.], [0., math.sqrt(3) * .25]]
            h5["3/t/sigma_d/2/2/matrix-imag"] = [[0., 0.], [0., 0.]]
    final = "\n".join(f"{min(values)} {sector}" for sector, values in levels.items())
    (work / "annotated.dat").write_text("0 2\n0.5 1\n\n0 2\n0.4 1\n\n0 2\n0.35 1\n\n" + final + "\n\n")
    # Same z/Z table layout emitted by nrginit/initial.m, including the extra
    # unused final hopping. Test Fortran exponents and ignore unrelated blocks.
    tables = {"z": ([50., 20., 8., 3.], [0., 0., 0., 0.]), "Z": ([1.] * 4, [0.] * 4)}
    lines = ["#!9", "# symtype SPSU2", "1 3 2", "e", "-0.5", "s n_d", "0"]
    for tag, arrays in tables.items():
        lines.append(tag)
        for values in arrays:
            lines.append("3")
            lines.extend(f"{v * gap_in_units:.17e}".replace("e", "D") for v in values)
    (work / "data").write_text("\n".join(lines) + "\n")
    return h5py


@pytest.mark.parametrize("units", ["gap", "bandwidth"])
@pytest.mark.parametrize("gap", [1., 2.])
def test_extract_units_probabilities_moment_and_full_chain(tmp_path, units, gap):
    physical = dict(PHYSICAL, gap=gap)
    make_output(tmp_path, units, physical=physical)
    result = nrg.extract_result(tmp_path, physical, NUMERICAL, units)
    json.dumps(result, allow_nan=False)
    singlet, doublet = (result["branches"][key] for key in ("singlet", "doublet"))
    assert singlet["ground_relative_energy_over_gap"] == pytest.approx(.3)
    assert doublet["ground_relative_energy_over_gap"] == 0.
    assert doublet["index"] == 1  # Select the minimum, not an assumed index 0.
    assert result["raw_total_energy_over_gap"] == pytest.approx(-120.)
    assert singlet["raw_total_energy_over_gap"] == pytest.approx(-119.7)
    assert result["doublet_minus_singlet_over_gap"] == pytest.approx(-.3)
    assert result["ground_spin_multiplicities"] == [2]
    assert [singlet[key] for key in ("P0", "P1", "P2", "moment")] == pytest.approx([.2, .6, .2, 0.])
    assert [doublet[key] for key in ("P0", "P1", "P2", "moment")] == pytest.approx([.2, .7, .1, .25])
    chain = result["wilson_chain"]
    assert chain["n_bath_sites"] == 4 and chain["n_bonds"] == 3
    assert chain["xi_over_gap"] == pytest.approx([50., 20., 8., 3.])
    assert chain["delta_over_gap"] == pytest.approx([1.] * 4)
    assert chain["zeta_over_gap"] == chain["kappa_over_gap"] == [0.] * 4
    assert chain["impurity_hopping_over_gap"] == pytest.approx(math.sqrt(80 / math.pi) / gap)
    convergence = result["shell_convergence"]
    assert convergence["last_step_gap_change_over_gap"] == pytest.approx(.05)
    assert convergence["same_parity_gap_change_over_gap"] == pytest.approx(.1)
    assert convergence["statistics"][-1]["scale_over_gap"] == pytest.approx(10 / 3)


def test_higher_spin_global_ground_and_clean_extraction(tmp_path):
    make_output(tmp_path, clean=True, higher_ground=True)
    result = nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)
    assert result["model"] == "CLEAN"
    assert result["global_ground"]["spin_multiplicity"] == 3
    assert result["ground_spin_multiplicities"] == [3]
    assert result["branches"]["singlet"]["ground_relative_energy_over_gap"] == .6
    assert result["branches"]["doublet"]["ground_relative_energy_over_gap"] == .1
    assert "moment" not in result["branches"]["doublet"]
    assert "P0" not in result["branches"]["singlet"]
    assert result["wilson_chain"]["impurity_hopping_over_gap"] == 0.


def test_free_up_spin_normalization_and_degenerate_minima(tmp_path):
    h5py = make_output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        h5["3/s/n_d/2/2/matrix"][1, 1] = 1.
        h5["3/s/n_d_ud/2/2/matrix"][1, 1] = 0.
        # An isolated physical up electron has <Sz>=1/2; its reduced spin
        # is sqrt(S(S+1))=sqrt(3)/2 in NRG's CG-only convention.
        h5["3/t/sigma_d/2/2/matrix"][1, 1] = math.sqrt(3) / 2
        h5["3/eigen/1/absenergy_zero"][:] = [.3, .3]
        h5["3/eigen/1/absenergy"][:] = [-119.7, -119.7]
    doublet = nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)["branches"]["doublet"]
    assert doublet["moment"] == pytest.approx(.5)
    assert [doublet[k] for k in ("P0", "P1", "P2")] == [0., 1., 0.]
    assert nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)["branches"]["singlet"]["degenerate_minima"] == 2


def test_optional_text_outputs_and_required_spin_matrix(tmp_path):
    h5py = make_output(tmp_path)
    for name in ("param", "data", "annotated.dat"):
        (tmp_path / name).unlink()
    result = nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)
    assert result["wilson_chain"] is None
    assert result["shell_convergence"]["branch_flow"] == []
    assert result["shell_convergence"]["last_step_gap_change_over_gap"] is None
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        del h5["3/t/sigma_d"]
    with pytest.raises(KeyError):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)


@pytest.mark.parametrize("path,value,match", [
    ("stats/GS_energy", -100., "inconsistent absolute"),
    ("3/eigen/2/absenergy_zero", [.4, .2], "inconsistent absolute"),
    ("3/s/n_d_ud/2/2/matrix", [[.2, 0.], [0., .6]], "probabilities"),
    ("3/t/sigma_d/2/2/matrix", [[.1, 0.], [0., 2.]], "moment exceeds"),
    ("3/t/sigma_d/2/2/matrix-imag", [[0., 0.], [0., .1]], "finite real"),
    ("stats/3/energy_offset", -119., "offset disagrees"),
])
def test_reject_inconsistent_hdf5(tmp_path, path, value, match):
    h5py = make_output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        h5[path][...] = value
    with pytest.raises(ValueError, match=match):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)


def test_reject_wrong_units_incomplete_shell_and_chain(tmp_path):
    h5py = make_output(tmp_path)
    with pytest.raises(ValueError, match="manifest/units"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL, "bandwidth")
    (tmp_path / "data").write_text("z\n1\n0\n0\nZ\n")
    with pytest.raises(ValueError, match="coefficient count"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        del h5["3"]
    with pytest.raises(ValueError, match="final Nmax"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)


@pytest.mark.parametrize("key,value", [
    ("keep", "256"), ("keepall", "0,1,2"), ("dumpEscale", "100"),
    ("dumpscaled", "true"), ("discretization", "Y"), ("variant", "SC"),
])
def test_reject_deck_manifest_or_convention_mismatch(tmp_path, key, value):
    make_output(tmp_path)
    deck = nrg.render_param(PHYSICAL, NUMERICAL)
    old = parse(deck)[key]
    (tmp_path / "param").write_text(deck.replace(f"\n{key}={old}\n", f"\n{key}={value}\n"))
    with pytest.raises(ValueError, match="param|variant"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)


def test_numeric_shell_order_and_crossing_ground_sectors(tmp_path):
    h5py = make_output(tmp_path)
    for name in ("param", "annotated.dat", "data"):
        (tmp_path / name).unlink()
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        h5.move("3", "10")
        h5.move("stats/3", "stats/10")
        h5.create_group("9")
        h5["10/eigen/1/absenergy_zero"][0] = 0.
        h5["10/eigen/1/absenergy"][0] = -120.
    result = nrg.extract_result(tmp_path, PHYSICAL, dict(NUMERICAL, nmax=10))
    assert result["final_h5_iteration"] == 10
    assert result["ground_spin_multiplicities"] == [1, 2]
    assert result["doublet_minus_singlet_over_gap"] == 0.


def test_reject_inconsistent_annotated_flow(tmp_path):
    make_output(tmp_path)
    (tmp_path / "annotated.dat").write_text("0 2\n0.5 1\n\n" * 4)
    with pytest.raises(ValueError, match="annotated.dat energies disagree"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)
    (tmp_path / "annotated.dat").write_text("0 2\n0.5 1\n\n")
    with pytest.raises(ValueError, match="every NRG shell"):
        nrg.extract_result(tmp_path, PHYSICAL, NUMERICAL)
