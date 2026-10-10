"""Bounded synthetic-output tests; never run NRG or Mathematica."""

import builtins
import importlib.util
import json
import math
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "NRG_comparisons/test2/manual_run/extract.py"
SPEC = importlib.util.spec_from_file_location("manual_nrg_extract", PATH)
extract = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(extract)

# Literal deck independent of the production renderer, with a short chain.
PARAM = """[param]
model=SIAM
variant=
symtype=SPSU2
band=flat
bandrescale=100
bcsgap=0.01
U=2
Gamma=0.4
delta=0.5
Lambda=2
z=0.5
Nz=8
Ninit=0
Nmax=1
keep=128
keepenergy=-1
keepmin=0
keepall=
lastall=true
discretization=Z
tri=old
options=SCHUR
diagratio=1
diag_mode=serial
diagth=1
strategy=kept
ops=n_d n_d_ud sigma_d
spect=sigma_d-sigma_d
finite=true
T=1e-5
dm=false
fdm=false
broaden=false
savebins=false
bins=50
calc0=true
dumpannotated=160
dumpscaled=false
dumpabs=false
dumpEscale=1
dumpprecision=16
dumpgroups=false
grouptol=1e-10
h5raw=true
h5all=false
h5last=true
h5ops=true
h5vectors=false
h5U=false
removefiles=true
silent=true
"""


def deck(path, **changes):
    lines = PARAM.splitlines()
    for key, value in changes.items():
        lines = [line for line in lines if not line.startswith(key + "=")]
        if value is not None:
            lines.append(f"{key}={value}")
    path.write_text("\n".join(lines) + "\n")


def output(work, *, higher_ground=False, hopping=2., onsite=(1., 3.), nmax=1):
    h5py = pytest.importorskip("h5py")
    np = pytest.importorskip("numpy")
    deck(work / "param", Nmax=nmax)
    # Both branch minima are deliberately index 1, not index 0.
    levels = {"1": [.9, .3], "2": [.4, .1 if higher_ground else 0.],
              "3": [0. if higher_ground else 1.2], "4": []}
    with h5py.File(work / "raw.h5", "w") as h5:
        h5["stats/GS_energy"] = -12.
        for step in range(nmax + 1):
            h5.create_group(str(step))
            for key, value in {"energyscale": 1., "scale": 1., "Teff": 1.,
                               "abs_Egs": -1., "energy_offset": -12.}.items():
                h5[f"stats/{step}/{key}"] = value
        for sector, energies in levels.items():
            h5[f"{nmax}/eigen/{sector}/absenergy_zero"] = energies
            h5[f"{nmax}/eigen/{sector}/absenergy"] = [e - 12. for e in energies]
        for sector in ("1", "2"):
            h5[f"{nmax}/s/n_d/{sector}/{sector}/matrix"] = np.diag([1., .9])
            h5[f"{nmax}/s/n_d_ud/{sector}/{sector}/matrix"] = np.diag([.2, .1])
        h5[f"{nmax}/t/sigma_d/2/2/matrix"] = np.diag([.1, math.sqrt(3) * .25])
    # The extra hopping 999 must never enter the finite-chain vacuum.
    tables = {"z": ([hopping] * nmax + [999.], list(onsite) + [0.] * (nmax - 1)),
              "Z": ([1.] * (nmax + 1), [0.] * (nmax + 1))}
    lines = ["#!9", "# symtype SPSU2", "1 3 2", "e", "-0.5", "s n_d", "0"]
    for tag, arrays in tables.items():
        lines.append(tag)
        for values in arrays:
            lines.append(str(nmax))
            lines.extend(f"{v:.17e}".replace("e", "D" if tag == "z" else "d") for v in values)
    (work / "data").write_text("\n".join(lines) + "\n")
    return h5py


@pytest.mark.parametrize("hopping,onsite", [(2., (1., 3.)), (.75, (-2., .5))])
def test_closed_form_centered_vacuum_minima_and_moment(tmp_path, hopping, onsite):
    output(tmp_path, hopping=hopping, onsite=onsite)
    result = extract.extract_result(tmp_path)
    assert set(result) == {"values", "diagnostics"}
    values, diagnostic = result["values"], result["diagnostics"]
    assert set(values) == {"signed_gap", "singlet.energy", "doublet.energy", "singlet.P0",
                           "singlet.P1", "singlet.P2", "doublet.P0", "doublet.P1", "doublet.P2",
                           "doublet.moment"}
    # Closed-form eigenvalues of a two-site normal chain, not numpy eigvalsh.
    center = sum(onsite) / 2
    radius = math.hypot((onsite[0] - onsite[1]) / 2, hopping)
    vacuum = -math.hypot(center - radius, 1) - math.hypot(center + radius, 1)
    assert diagnostic["analytic_bath_vacuum_over_gap"] == pytest.approx(vacuum, abs=1e-13)
    assert values["singlet.energy"] == pytest.approx(-11.7 - vacuum - .5, abs=1e-13)
    assert values["doublet.energy"] == pytest.approx(-12. - vacuum - .5, abs=1e-13)
    assert values["signed_gap"] == pytest.approx(-.3, abs=1e-13)
    assert values["doublet.energy"] - values["singlet.energy"] == pytest.approx(-.3, abs=1e-13)
    for name in ("singlet", "doublet"):
        assert [values[f"{name}.{p}"] for p in ("P0", "P1", "P2")] == pytest.approx([.2, .7, .1], abs=1e-14)
        assert diagnostic["branches"][name]["index"] == 1
    assert values["doublet.moment"] == pytest.approx(.25, abs=1e-14)
    assert diagnostic["global_ground"]["spin_multiplicity"] == 2
    assert diagnostic["wilson_chain"]["xi_over_gap"] == [hopping, 999.]
    assert diagnostic["wilson_chain"]["impurity_hopping_over_gap"] == pytest.approx(math.sqrt(80 / math.pi))
    assert diagnostic["ground_spin_multiplicities"] == [2]
    json.dumps(result, allow_nan=False)
    assert not (tmp_path / "DONE").exists()  # No controller receipts required.


def test_agrees_with_shared_extractor_and_postprocessor(tmp_path, monkeypatch):
    output(tmp_path)
    monkeypatch.syspath_prepend(str(ROOT))
    from NRG_comparisons.test1.nrg import extract_result
    from NRG_comparisons.test1.run import postprocess_nrg

    physical = dict(gap=1., bandwidth=100., u=2., gamma=.4, detuning=.5, field=0., geometry="single")
    numerical = {"lambda": 2., "z": .5, "nmax": 1, "keep": 128, "keepenergy": -1., "keepmin": 0}
    task = dict(numerical=numerical, units="gap", clean=False)
    shared = extract_result(tmp_path, physical, numerical)
    processed = postprocess_nrg(task, physical, tmp_path)
    manual = extract.extract_result(tmp_path)
    for name, value in manual["values"].items():
        assert value == pytest.approx(processed["values"][name], abs=1e-13)
    assert manual["diagnostics"]["analytic_bath_vacuum_over_gap"] == processed["diagnostics"]["analytic_bath_vacuum_over_gap"]
    for name in ("singlet", "doublet"):
        assert manual["diagnostics"]["branches"][name]["index"] == shared["branches"][name]["index"]


def test_global_ground_in_higher_sector_and_numeric_shell_order(tmp_path):
    output(tmp_path, higher_ground=True, nmax=10)
    result = extract.extract_result(tmp_path)
    assert result["diagnostics"]["final_h5_iteration"] == 10
    assert result["diagnostics"]["ground_spin_multiplicities"] == [3]
    assert result["diagnostics"]["global_ground"]["spin_multiplicity"] == 3
    assert result["diagnostics"]["branches"]["doublet"]["ground_relative_energy_over_gap"] == .1
    assert result["values"]["signed_gap"] == pytest.approx(-.2, abs=1e-13)


def test_up_electron_normalization_and_cross_sector_degeneracy(tmp_path):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        h5["1/s/n_d/2/2/matrix"][1, 1] = 1.
        h5["1/s/n_d_ud/2/2/matrix"][1, 1] = 0.
        h5["1/t/sigma_d/2/2/matrix"][1, 1] = math.sqrt(3) / 2
        h5["1/eigen/1/absenergy_zero"][1] = 0.
        h5["1/eigen/1/absenergy"][1] = -12.
    result = extract.extract_result(tmp_path)
    assert result["values"]["doublet.moment"] == pytest.approx(.5, abs=1e-14)
    assert result["values"]["signed_gap"] == 0.
    assert result["diagnostics"]["ground_spin_multiplicities"] == [1, 2]


@pytest.mark.parametrize("changes", [
    {"Lambda": "1.8", "Nmax": "998", "keep": "2", "keepmin": "2", "keepenergy": "12.5", "z": "1", "Nz": "1"},
    {"Lambda": "2.5", "keepenergy": "-1", "Nz": None, "B": "0"},
])
def test_param_accepts_numerical_edits_and_preserves_full_dict(tmp_path, changes):
    path = tmp_path / "param"
    deck(path, **changes)
    params = extract.read_param(path)
    assert params == dict(line.split("=", 1) for line in path.read_text().splitlines()[1:])
    assert "Gamma" in params and "gamma" not in params


@pytest.mark.parametrize("key,value", [
    ("model", "CLEAN"), ("variant", "SC"), ("symtype", "QS"), ("band", "other"),
    ("bandrescale", "1"), ("bcsgap", "1"), ("U", "3"), ("Gamma", ".5"), ("delta", "0"), ("B", ".1"),
    ("Ninit", "1"), ("lastall", "false"), ("keepall", "0"), ("discretization", "Y"), ("tri", "new"),
    ("strategy", "all"), ("ops", "n_d n_d_ud"), ("spect", ""), ("finite", "false"),
    ("dumpscaled", "true"), ("dumpabs", "true"), ("dumpEscale", "100"), ("dumpgroups", "true"),
    ("h5raw", "false"), ("h5last", "false"), ("h5ops", None), ("diagratio", ".5"),
    ("Lambda", "1.7999"), ("Lambda", "nan"), ("Lambda", "inf"), ("keepenergy", "-inf"),
    ("Nmax", "0"), ("Nmax", "999"), ("Nmax", "1.5"), ("Nmax", None),
    ("keep", "1"), ("keep", "true"), ("keepmin", "-1"), ("keepmin", "129"),
    ("z", "0"), ("z", "1.01"), ("Nz", "0"), ("Nz", "1.5"),
    ("perturb", "0.1 number[d[]]"), ("safeguard", "0"), ("diag_mode", "parallel"),
])
def test_param_rejects_changed_physics_conventions_or_invalid_numbers(tmp_path, key, value):
    deck(tmp_path / "param", **{key: value})
    with pytest.raises(ValueError, match=key):
        extract.read_param(tmp_path / "param")


@pytest.mark.parametrize("text", ["not an ini file", PARAM + "U=2\n", PARAM.replace("[param]", "[other]"),
                                  "[DEFAULT]\nextra=bad\n" + PARAM])
def test_malformed_param(tmp_path, text):
    (tmp_path / "param").write_text(text)
    with pytest.raises(ValueError, match="param"):
        extract.read_param(tmp_path / "param")


def test_standalone_import_and_param_read_do_not_need_optional_packages(tmp_path, monkeypatch):
    deck(tmp_path / "param")
    copied = tmp_path / "extract.py"
    copied.write_bytes(PATH.read_bytes())
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"numpy", "h5py", "qdjj_solver", "qpsolver", "NRG_comparisons"}:
            raise ImportError(f"deliberately blocked {name}")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    spec = importlib.util.spec_from_file_location("copied_extract", copied)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.read_param(tmp_path / "param")["Nz"] == "8"


def test_extraction_rereads_param_and_accepts_valid_numerical_edits(tmp_path):
    output(tmp_path)
    deck(tmp_path / "param", Lambda="1.8", keep="256", keepmin="10", keepenergy="8")
    assert extract.extract_result(tmp_path)["diagnostics"]["param"]["keep"] == "256"
    deck(tmp_path / "param", delta="0")
    with pytest.raises(ValueError, match="delta"):
        extract.extract_result(tmp_path)
    deck(tmp_path / "param", Nmax="2")
    with pytest.raises(ValueError, match="final Nmax"):
        extract.extract_result(tmp_path)


def test_copied_adapter_extracts_without_repository_imports_or_nz(tmp_path, monkeypatch):
    output(tmp_path)
    deck(tmp_path / "param", Nz=None)
    copied = tmp_path / "extract.py"
    copied.write_bytes(PATH.read_bytes())
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in {"qdjj_solver", "qpsolver", "NRG_comparisons"}:
            raise ImportError(f"deliberately blocked {name}")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    monkeypatch.chdir(tmp_path)
    spec = importlib.util.spec_from_file_location("standalone_extract", copied)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.extract_result(".")["values"]["doublet.moment"] == pytest.approx(.25, abs=1e-14)


@pytest.mark.parametrize("name", ["data", "raw.h5", "param"])
def test_required_files(tmp_path, name):
    output(tmp_path)
    (tmp_path / name).unlink()
    with pytest.raises((ValueError, FileNotFoundError), match=name):
        extract.extract_result(tmp_path)


@pytest.mark.parametrize("path,value,match", [
    ("1", 0., "missing final eigen"),
    ("1/eigen/1", 0., "missing dataset"),
    ("stats/GS_energy", -10., "inconsistent absolute"),
    ("stats/1/energy_offset", -11., "offset disagrees"),
    ("1/eigen/1/absenergy_zero", [.9, .2], "inconsistent absolute"),
    ("1/eigen/1/absenergy", [-11.1], "eigenvalue counts"),
    ("1/eigen/1/absenergy", [[-11.1, -11.7]], "dimensions"),
    ("1/eigen/1/absenergy", [float("nan"), -11.7], "finite real"),
    ("1/s/n_d/1/1/matrix", [[0., 0.], [0., float("inf")]], "finite real"),
    ("1/s/n_d_ud/2/2/matrix", [[.2, 0.], [0., .6]], "probabilities"),
    ("1/t/sigma_d/2/2/matrix", [[.1, 0.], [0., 2.]], "moment exceeds"),
    ("1/s/n_d/2/2/matrix", [[1.]], "matrix shape"),
    ("1/s/n_d/1/1/matrix", [[1., 0.], [0., .9 + .01j]], "finite real"),
    ("1/t/sigma_d/2/2/matrix-imag", [[0., 0.], [0., .1]], "finite real"),
    ("1/t/sigma_d/2/2/matrix-imag", [[0.]], "matrix shape"),
])
def test_malformed_hdf5(tmp_path, path, value, match):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        if path in h5:
            del h5[path]
        h5[path] = value
    with pytest.raises(ValueError, match=match):
        extract.extract_result(tmp_path)


@pytest.mark.parametrize("path,match", [
    ("1", "final Nmax"), ("stats/1", "final-shell statistics"),
    ("stats/1/energy_offset", "missing dataset"), ("1/eigen", "missing final eigen"),
    ("1/eigen/1", "singlet and doublet"), ("1/eigen/2", "singlet and doublet"),
    ("1/s/n_d/1/1/matrix", "missing dataset"), ("1/t/sigma_d", "missing dataset"),
])
def test_incomplete_hdf5(tmp_path, path, match):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        del h5[path]
    with pytest.raises(ValueError, match=match):
        extract.extract_result(tmp_path)


@pytest.mark.parametrize("sector", ["1", "2"])
def test_rejects_repeated_lowest_multiplets(tmp_path, sector):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        for key in ("absenergy", "absenergy_zero"):
            h5[f"1/eigen/{sector}/{key}"][0] = h5[f"1/eigen/{sector}/{key}"][1]
    with pytest.raises(ValueError, match="degenerate lowest.*basis dependent"):
        extract.extract_result(tmp_path)


def test_unmeasured_sector_degeneracy_is_only_diagnostic(tmp_path):
    h5py = output(tmp_path)
    expected = extract.extract_result(tmp_path)["values"]
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        for key in ("absenergy", "absenergy_zero"):
            value = h5[f"1/eigen/3/{key}"][0]
            del h5[f"1/eigen/3/{key}"]
            h5[f"1/eigen/3/{key}"] = [value, value]
    result = extract.extract_result(tmp_path)
    assert result["values"] == expected
    assert result["diagnostics"]["all_sector_minima"]["3"]["degenerate_minima"] == 2


@pytest.mark.parametrize("shift", [-.1, .1])
def test_rejects_nonzero_all_sector_ground(tmp_path, shift):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        for group in h5["1/eigen"].values():
            for key in ("absenergy", "absenergy_zero"):
                group[key][:] += shift
    with pytest.raises(ValueError, match="all-sector.*not zero"):
        extract.extract_result(tmp_path)


def test_optional_zero_imaginary_matrices(tmp_path):
    h5py = output(tmp_path)
    expected = extract.extract_result(tmp_path)["values"]
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        for path in ("s/n_d/1/1", "s/n_d_ud/2/2", "t/sigma_d/2/2"):
            h5[f"1/{path}/matrix-imag"] = [[0., .3], [-.3, 0.]]
    assert extract.extract_result(tmp_path)["values"] == expected


@pytest.mark.parametrize("sector", ["0", "not-a-spin"])
def test_invalid_sector_labels(tmp_path, sector):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        h5.move("1/eigen/3", f"1/eigen/{sector}")
    with pytest.raises(ValueError, match="SPSU2 spin-multiplicity sector"):
        extract.extract_result(tmp_path)


def test_branch_difference_must_match_relative_signed_gap(tmp_path):
    h5py = output(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        # Each branch passes the 1e-9 absolute-energy check separately, but
        # their opposed errors exceed that tolerance in the excitation gap.
        h5["1/eigen/1/absenergy"][1] += 8e-10
        h5["1/eigen/2/absenergy"][1] -= 8e-10
    with pytest.raises(ValueError, match="branch difference disagrees"):
        extract.extract_result(tmp_path)


@pytest.mark.parametrize("old,new,match", [
    ("\nz\n", "\ny\n", "one z"), ("\nZ\n", "\nz\n", "one z"),
    ("\nz\n1\n", "\nz\n2\n", "coefficient count"),
    ("\nz\n1\n", "\nz\nbad\n", "malformed xi"),
    ("2.00000000000000000D+00", "nan", "finite real"),
    ("1.00000000000000000d+00", "0.9", "uniform unit gap"),
    ("0.00000000000000000d+00", "0.1", "anomalous hopping"),
])
def test_malformed_chain(tmp_path, old, new, match):
    output(tmp_path)
    path = tmp_path / "data"
    assert old in path.read_text()
    path.write_text(path.read_text().replace(old, new, 1))
    with pytest.raises(ValueError, match=match):
        extract.extract_result(tmp_path)


def test_incomplete_chain(tmp_path):
    output(tmp_path)
    path = tmp_path / "data"
    path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
    with pytest.raises(ValueError, match="incomplete kappa"):
        extract.extract_result(tmp_path)


def test_optional_annotated_final_validation(tmp_path):
    output(tmp_path)
    path = tmp_path / "annotated.dat"
    path.write_text("0 2\n.5 1\n\n0 2\n.3 1\n1.2 3\n\n")
    assert extract.extract_result(tmp_path)["values"]["signed_gap"] == -.3
    path.write_text("0 2\n.5 1\n\n0 2\n.4 1\n\n")
    with pytest.raises(ValueError, match="annotated.dat energies disagree"):
        extract.extract_result(tmp_path)
    path.write_text("0 2\n.3 1\n\n")
    with pytest.raises(ValueError, match="every NRG shell"):
        extract.extract_result(tmp_path)
