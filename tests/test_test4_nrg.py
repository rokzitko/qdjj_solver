"""Offline two-channel contracts and synthetic HDF5; never execute NRG/math."""

import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest


HERE = Path(__file__).resolve().parents[1] / "NRG_comparisons/test4"
PHYSICAL = dict(geometry="junction", gap=1., bandwidth=100., u=2., gamma=.4,
                nu=.75, field=.03, phi=3 * math.pi / 5, wn=0., ws=0.)
NUMERICAL = {"lambda": 4., "z": .5, "nmax": 1, "keep": 128, "keepmin": 0, "keepenergy": -1.}


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.fixture
def adapter(monkeypatch):
    physics = module(HERE / "physics.py", "test4_physics")
    monkeypatch.setitem(sys.modules, "physics", physics)
    return module(HERE / "nrg.py", "test4_nrg")


def synthetic(directory, *, scale=1., representation="paper", levels=None, vector=False):
    """Independent literal input with unequal chains and unsorted signed sectors."""
    h5py = pytest.importorskip("h5py")
    deck = """[param]
symtype=SPU1
band=flat
Ninit=0
lastall=true
discretization=Z
tri=old
diagratio=1
diag_mode=serial
diagth=1
strategy=kept
finite=true
dm=false
fdm=false
broaden=false
savebins=false
calc0=true
h5raw=true
h5all=false
h5last=true
h5ops=true
h5vectors=false
h5U=false
removefiles=true
silent=true
Lambda=4
z=0.5
Nmax=1
keep=128
keepmin=0
keepenergy=-1
keepall=
"""
    deck += f"bandrescale={100 * scale:.17g}\nU={2 * scale:.17g}\nGamma={.4 * scale:.17g}\nT={1e-5 * scale:.17g}\n"
    if representation == "paper":
        deck += "model=model.m\nreaddir=.\noptions=COMPLEX SCHUR\nops=n_d n_d_ud SZd SZ2d asqJ\nbcsgap=-0.01\n"
        deck += f"[extra]\nnu=0.75\nBz={.03 * scale:.17g}\nWN=0\nWS=0\nphi={PHYSICAL['phi']:.17g}\n"
    else:
        deck += "model=ONE\nvariant=SDSMAGFIELD\noptions=SCHUR\nops=n_d n_d_ud SZd SZ2d hop0 hop1\n"
        deck += f"bcsgap1=-0.01\nbcsgap2=0.01\ndelta={.5 * scale:.17g}\n[extra]\nB={.03 * scale:.17g}\nphi={PHYSICAL['phi']:.17g}\n"
    (directory / "param").write_text(deck)
    levels = levels if levels is not None else {
        "1": [2.4, .3, 1., .7, .8], "0": [1.5, 0., 1.5, 1.5],
        "2": [.4, 1.8], "-1": [1.6, 2., 3., 4.], "3": [1.7], "4": [],
    }
    with h5py.File(directory / "raw.h5", "w") as h5:
        for path, value in (("stats/GS_energy", -12. * scale),
                            ("stats/1/energy_offset", -12. * scale), ("params/channels", 2)):
            h5[path] = [value] if vector else value
        for sector, energies in levels.items():
            h5[f"1/eigen/{sector}/absenergy_zero"] = np.array(energies) * scale
            h5[f"1/eigen/{sector}/absenergy"] = (np.array(energies) - 12.) * scale
            if not energies:
                continue
            charge = .8 + np.arange(len(energies)) * .05
            md = {"1": -.1, "0": -.25, "2": .25}.get(sector, .125)
            values = {"n_d": charge, "n_d_ud": np.full(len(energies), .1),
                      "SZd": np.full(len(energies), md), "SZ2d": (charge - .2) / 4}
            if representation == "paper":
                values["asqJ"] = (np.arange(len(energies)) + .2) * scale
            else:
                values.update(hop0=np.arange(len(energies)) + .2, hop1=np.arange(len(energies)) - .3)
            for name, diagonal in values.items():
                path = f"1/s/{name}/{sector}/{sector}/matrix"
                h5[path] = np.diag(diagonal)
                imaginary = np.zeros((len(energies), len(energies)))
                if len(energies) > 1:
                    imaginary[0, 1], imaginary[1, 0] = .13, -.13
                h5[path + "-imag"] = imaginary
    lines = ["#!9", "# symtype SPU1", "2 3 1", "e", "-0.5"]
    for tag, arrays in (("z", ([2., 999.], [.5, 777.], [1., 3.], [-2., 0.])),
                        ("Z", ([-1., -1.], [-1., -1.] if representation == "paper" else [1., 1.],
                               [0., 0.], [0., 0.]))):
        lines.append(tag)
        for values in arrays:
            lines += ["1", *(f"{value * scale:.17e}".replace("e", "D") for value in values)]
    (directory / "data").write_text("\n".join(lines) + "\n")
    return h5py


@pytest.mark.parametrize("representation", ["paper", "sds"])
@pytest.mark.parametrize("units", ["gap", "bandwidth"])
def test_renderer_mapping_and_roundtrip(adapter, tmp_path, representation, units):
    physical = {k: v if k in ("geometry", "nu", "phi") else 2.5 * v for k, v in PHYSICAL.items()}
    if representation == "paper":
        physical.update(wn=.7, ws=-.2)
    text = adapter.render_param(physical, dict(NUMERICAL, nmax=3, untruncated=True), units, representation)
    (tmp_path / "param").write_text(text)
    parsed = adapter.read_param(tmp_path / "param")
    p, extra, unit = parsed["param"], parsed["extra"], physical[units]
    assert p["keepall"] == "0,1,2" and float(p["Lambda"]) == 4
    assert p["finite"] == p["h5last"] == p["h5ops"] == p["lastall"] == "true"
    assert not {"spect", "specd", "polarized", "rungs", "substeps"} & p.keys()
    assert float(p["Gamma"]) == pytest.approx(physical["gamma"] / unit)
    assert parsed["representation"] == representation
    assert parsed["physical"] == pytest.approx({k: v if k in ("geometry", "nu", "phi") else v / unit
                                              for k, v in physical.items()})
    assert all(adapter.LITERAL.fullmatch(value) for value in extra.values())
    if representation == "paper":
        assert p["model"] == "model.m" and p["readdir"] == "." and p["options"] == "COMPLEX SCHUR"
        assert float(extra["WN"]) == pytest.approx(.7 / unit)  # No rho or bandwidth conversion.
        assert float(extra["WS"]) == pytest.approx(-.2 / unit)
        assert set(extra) == {"nu", "Bz", "WN", "WS", "phi"} and "asqJ" in p["ops"]
        assert float(p["bcsgap"]) == -.01
    else:
        assert (p["model"], p["variant"]) == ("ONE", "SDSMAGFIELD")
        assert float(p["delta"]) == pytest.approx(physical["u"] * (1 - physical["nu"]) / unit)
        assert float(extra["B"]) == pytest.approx(physical["field"] / unit)
        assert (float(p["bcsgap1"]), float(p["bcsgap2"])) == (-.01, .01)


@pytest.mark.parametrize("change", [dict(field=float("nan")), dict(geometry="single"), dict(gap=0),
                                    dict(gamma=-1), dict(rho_ws=0), dict(phi=True), dict(nu=float("inf"))])
def test_unsupported_physics(adapter, change):
    with pytest.raises(ValueError):
        adapter.render_param(dict(PHYSICAL, **change), NUMERICAL)


@pytest.mark.parametrize("change", [dict(nmax=1.), dict(nmax=999), dict(keepmin=129), dict(untruncated=1),
                                    {"lambda": 3.999999}, dict(z=0), dict(polarized=False), dict(keep=1)])
@pytest.mark.parametrize("representation", ["paper", "sds"])
def test_unsupported_numerics(adapter, change, representation):
    with pytest.raises(ValueError):
        adapter.render_param(PHYSICAL, dict(NUMERICAL, **change), representation=representation)


def test_invalid_representations_and_sds_bridge(adapter):
    for change in (dict(wn=.1), dict(ws=.1), dict(u=0)):
        with pytest.raises(ValueError, match="SDS"):
            adapter.render_param(dict(PHYSICAL, **change), NUMERICAL, representation="sds")
    with pytest.raises(ValueError, match="representation"):
        adapter.render_param(PHYSICAL, NUMERICAL, representation="other")
    with pytest.raises(ValueError, match="units"):
        adapter.render_param(PHYSICAL, NUMERICAL, units="other")


@pytest.mark.parametrize("old,new", [
    ("[extra]", "[other]"), ("Bz=", "B="), ("WN=", "W_N="), ("SPU1", "SPSU2"),
    ("tri=old", "tri=cpp"), ("Ninit=0", "Ninit=1"), ("h5ops=true", "h5ops=false"),
    ("keepall=", "keepall=even"), ("bcsgap=-0.01", "bcsgap=0.01"), ("readdir=.", "readdir=.."),
    ("Lambda=4", "Lambda=3.99"), ("[extra]", "polarized=true\n[extra]"),
    ("[extra]", "rungs=true\n[extra]"), ("[extra]", "substeps=true\n[extra]"),
    ("[extra]", "spect=SZd-SZd\n[extra]"), ("WN=0", "WN=!1/2"), ("WN=0", "WN=0\nU=2"),
    ("T=1.0000000000000001e-05", "T=0.1"),
])
def test_read_rejects_physics_and_output_changes(adapter, tmp_path, old, new):
    text = adapter.render_param(PHYSICAL, NUMERICAL)
    assert old in text
    (tmp_path / "param").write_text(text.replace(old, new))
    with pytest.raises(ValueError):
        adapter.read_param(tmp_path / "param")


@pytest.mark.parametrize("representation", ["paper", "sds"])
def test_lambda_floor_also_on_extraction(adapter, tmp_path, representation):
    synthetic(tmp_path, representation=representation)
    path = tmp_path / "param"
    path.write_text(path.read_text().replace("Lambda=4", "Lambda=2"))
    with pytest.raises(ValueError, match="Lambda must be >= 4"):
        adapter.extract_result(tmp_path)


@pytest.mark.parametrize("old,new", [("Lambda=4", "Lambda=!4"), ("Lambda=4\n", ""),
                                    ("[param]", "[param]\nNz=0"), ("[param]", "[param]\nNz=2\nNz=2"),
                                    ("bcsgap2=0.01", "bcsgap2=-0.01")])
def test_malformed_literal_decks(adapter, tmp_path, old, new):
    synthetic(tmp_path, representation="sds")
    path = tmp_path / "param"
    assert old in path.read_text()
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError):
        adapter.read_param(path)


def test_nonfinite_unit_conversion_rejected(adapter):
    with pytest.raises(ValueError, match="unit conversion"):
        adapter.render_param(dict(PHYSICAL, gap=1e-200, gamma=1e200), NUMERICAL)


@pytest.mark.parametrize("units,scale", [("gap", 1.), ("bandwidth", .01)])
def test_extract_canonical_nonunit_gap_and_bridge(adapter, tmp_path, units, scale):
    synthetic(tmp_path, scale=scale)
    expected = dict(PHYSICAL, wn=.7, ws=-.2)
    physical = {k: v if k in ("geometry", "nu", "phi") else 2.5 * v for k, v in expected.items()}
    (tmp_path / "param").write_text(adapter.render_param(physical, NUMERICAL, units=units))
    assert adapter.extract_result(tmp_path)["physical"] == pytest.approx(expected)


@pytest.mark.parametrize("scale", [1., .01])
@pytest.mark.parametrize("vector", [False, True])
@pytest.mark.parametrize("representation", ["paper", "sds"])
def test_extract_two_vacua_complex_matrices_units_currents(adapter, tmp_path, scale, vector, representation):
    synthetic(tmp_path, scale=scale, vector=vector, representation=representation)
    result = adapter.extract_result(tmp_path)
    assert result["physical"] == pytest.approx(PHYSICAL)
    assert result["numerical"] == dict(NUMERICAL, untruncated=False)
    diag = result["diagnostics"]
    # Analytic eigenvalues of [[1,2],[2,3]] and [[-2,.5],[.5,0]].
    vacuum = -sum(math.hypot(e, 1) for e in
                  (2 - math.sqrt(5), 2 + math.sqrt(5), -1 - math.sqrt(1.25), -1 + math.sqrt(1.25)))
    assert diag["analytic_bath_vacuum_over_gap"] == pytest.approx(vacuum, abs=1e-13)
    assert diag["deck_gap_scale"] == scale
    assert len(diag["wilson_chains"]) == 2
    assert diag["wilson_chains"][0]["xi_over_gap"] == [2., 999.]
    assert diag["wilson_chains"][1]["xi_over_gap"] == [.5, 777.]
    assert diag["wilson_chains"][1]["zeta_over_gap"] == [-2., 0.]
    assert diag["wilson_chains"][0]["delta_over_gap"] == [-1., -1.]
    states = {row["id"]: row for row in result["states"]}
    offset = .0625 if representation == "paper" else .5
    assert states["even"]["values"]["energy"] == pytest.approx(-11.7 - vacuum - offset, abs=1e-13)
    assert states["odd_down"]["values"]["energy_paper"] == pytest.approx(-12 - vacuum - offset + .0625, abs=1e-13)
    assert states["even"]["excitation"] == pytest.approx(.3, abs=1e-13)
    assert states["odd_down"]["values"]["m_dot"] == -.25
    assert states["odd_up"]["values"]["m_dot"] == .25
    assert states["sz-2-r0"]["twice_sz"] == -2
    assert states["sz-2-r2"]["root"] == 2 and not states["sz-2-r2"]["subgap"]
    assert "sz-2-r3" not in states and "sz+0-r3" not in states  # Strict E-Eg=1 boundary.
    assert states["sz-1-r2"]["degenerate"]  # Includes a discarded partner.
    expected_current = 1.2 if representation == "paper" else math.sqrt(80 / math.pi) / 4 * (
        -math.sin(PHYSICAL["phi"] / 4) * 1.2 + math.cos(PHYSICAL["phi"] / 4) * .7)
    assert states["even"]["values"]["current"] == pytest.approx(expected_current, abs=1e-13)
    assert not result["coverage"]["continuum_completeness_established"]
    assert result["coverage"]["diagonal_observables_unique"]
    for row in states.values():
        v, M = row["values"], row["twice_sz"] / 2
        assert set(v) == set(sys.modules["physics"].OBSERVABLES)
        assert row["parity"] == row["twice_sz"] % 2 and row["residual"] is None
        assert v["P_up"] - v["P_down"] == pytest.approx(2 * v["m_dot"], abs=1e-14)
        assert v["m_bath"] == pytest.approx(M - v["m_dot"], abs=1e-14)
        assert v["Czz"] == pytest.approx(M * v["m_dot"] - v["dot_z2"], abs=1e-14)
    json.dumps(result, allow_nan=False)
    assert not (tmp_path / "DONE").exists()


def test_subgap_degeneracy_all_roots_missing_primary_and_optional_channels(adapter, tmp_path):
    h5py = synthetic(tmp_path, levels={"-1": [.2, 0.], "1": [.3, .3, .8, .9], "0": [.1, 1.4]})
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        del h5["params/channels"]
    result = adapter.extract_result(tmp_path)
    coverage = result["coverage"]
    assert coverage["ground_ids"] == ["sz-2-r0"]
    assert "odd_up" not in coverage["subgap_ids"]
    assert {"even", "sz+0-r1"} <= set(coverage["degenerate_subgap_ids"])
    assert not coverage["diagonal_observables_unique"]
    assert "sz+0-r3" in coverage["subgap_ids"]


def test_real_matrices_without_imaginary_datasets(adapter, tmp_path):
    h5py = synthetic(tmp_path, levels={"1": [0.]})
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        for name in ("n_d", "n_d_ud", "SZd", "SZ2d", "asqJ"):
            del h5[f"1/s/{name}/1/1/matrix-imag"]
    assert adapter.extract_result(tmp_path)["states"][0]["values"]["current"] == .2


@pytest.mark.parametrize("path,value", [
    ("params/channels", 1), ("params/channels", [2, 2]),
    ("stats/GS_energy", [0., 0.]), ("stats/1/energy_offset", -11.),
    ("1/eigen/1/absenergy", [-12.] * 5), ("1/eigen/1/absenergy_zero", [float("nan")] * 5),
    ("1/s/SZd/1/1/matrix", [[.2]]), ("1/s/SZd/1/1/matrix-imag", [[0.]]),
    ("1/s/SZd/1/1/matrix-imag", np.eye(5) * .01),
    ("1/s/SZd/1/1/matrix-imag", np.triu(np.ones((5, 5)), 1)),
    ("1/s/asqJ/1/1/matrix", np.eye(5) * float("nan")),
    ("1/s/n_d_ud/1/1/matrix", np.eye(5) * .9), ("1/s/SZ2d/1/1/matrix", np.eye(5) * .5),
])
def test_bad_hdf5_rejected(adapter, tmp_path, path, value):
    h5py = synthetic(tmp_path)
    with h5py.File(tmp_path / "raw.h5", "a") as h5:
        if path in h5:
            del h5[path]
        h5[path] = value
    with pytest.raises(ValueError):
        adapter.extract_result(tmp_path)


@pytest.mark.parametrize("issue", ["sector", "matrix", "final", "count", "pairing", "anomalous", "extra_chain", "truncated",
                                  "missing_tag", "missing_eigen", "sector_dataset", "nonzero_ground"])
def test_incomplete_and_inconsistent_inputs(adapter, tmp_path, issue):
    h5py = synthetic(tmp_path)
    if issue in ("sector", "matrix", "final", "missing_eigen", "sector_dataset", "nonzero_ground"):
        with h5py.File(tmp_path / "raw.h5", "a") as h5:
            if issue == "sector":
                h5.move("1/eigen/-1", "1/eigen/0.5")
            elif issue == "matrix":
                del h5["1/s/asqJ/0/0/matrix"]
            elif issue == "missing_eigen":
                del h5["1/eigen"]
            elif issue == "sector_dataset":
                del h5["1/eigen/0"]
                h5["1/eigen/0"] = 0
            elif issue == "nonzero_ground":
                for sector in h5["1/eigen"]:
                    h5[f"1/eigen/{sector}/absenergy"][...] += .2
                    h5[f"1/eigen/{sector}/absenergy_zero"][...] += .2
            else:
                h5.create_group("10")
    else:
        path = tmp_path / "data"
        text = path.read_text()
        if issue == "count":
            text = text.replace("z\n1\n", "z\n2\n")
        elif issue == "pairing":
            text = text.replace("Z\n1\n-1.00000000000000000D+00", "Z\n1\n1.00000000000000000D+00")
        elif issue == "anomalous":
            text = text.replace("0.00000000000000000D+00", "1.00000000000000000D+00")
        elif issue == "extra_chain":
            text += "1\n0\n0\n"
        elif issue == "missing_tag":
            text = text.replace("Z\n", "Q\n")
        else:
            text = text.rsplit("\n", 3)[0]
        path.write_text(text)
    with pytest.raises(ValueError):
        adapter.extract_result(tmp_path)


def tiny_hamiltonian(physical, representation):
    """Independent six-mode electron ED, dot and one orbital in each lead."""
    eye = np.eye(64)
    annihilation = []
    for mode in range(6):
        matrix = np.zeros((64, 64))
        for state in range(64):
            if state & (1 << mode):
                matrix[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode) - 1)).bit_count())
        annihilation.append(matrix)
    number = [a.T @ a for a in annihilation]
    charge, md = number[0] + number[1], (number[0] - number[1]) / 2
    dot = charge - (physical["nu"] if representation == "paper" else 1) * eye
    H = (physical["u"] / 2 * dot @ dot + physical["field"] * md).astype(complex)
    if representation == "sds":
        H += physical["u"] * (1 - physical["nu"]) * charge
    J = np.zeros_like(H)
    ops = {"n_d": charge, "n_d_ud": number[0] @ number[1], "SZd": md, "SZ2d": md @ md}
    amplitude = math.sqrt(2 * physical["bandwidth"] * physical["gamma"] / math.pi)
    phi = physical["phi"]
    for lead, mode in enumerate((2, 4)):
        pairing = annihilation[mode].T @ annihilation[mode + 1].T
        H += physical["gap"] * (-1 if lead == 0 or representation == "paper" else 1) * (pairing + pairing.T)
        contact = sum(annihilation[mode + spin].T @ annihilation[spin] for spin in (0, 1))
        ops[f"hop{lead}"] = contact + contact.T
        if representation == "paper":
            sign = -1 if lead == 0 else 1
            t = amplitude / math.sqrt(2) * np.exp(sign * 1j * phi / 4)
            dt = sign * 1j / 4 * t
        else:
            t = amplitude * (math.cos(phi / 4) if lead == 0 else math.sin(phi / 4))
            dt = amplitude / 4 * (-math.sin(phi / 4) if lead == 0 else math.cos(phi / 4))
        H += t * contact + np.conjugate(t) * contact.T
        J += dt * contact + np.conjugate(dt) * contact.T
    if representation == "paper":
        for spin, sign in ((0, 1), (1, -1)):
            contact = annihilation[2 + spin].T @ annihilation[4 + spin]
            t = (physical["wn"] + sign * 1j * physical["ws"]) * np.exp(-1j * phi / 2)
            dt = -1j / 2 * t
            H += t * contact + np.conjugate(t) * contact.T
            J += dt * contact + np.conjugate(dt) * contact.T
    ops["asqJ"] = J
    spins = np.diag(sum(number[mode] * (1 if mode % 2 == 0 else -1) for mode in range(6))).astype(int)
    return H, ops, spins


@pytest.mark.parametrize("bridge", [False, True])
def test_paper_phase_derivative_and_sds_control_by_electron_ed(adapter, tmp_path, bridge):
    physical = dict(PHYSICAL, bandwidth=3., gamma=.2, wn=.17 if bridge else 0., ws=-.09 if bridge else 0.)
    results = {}
    for representation in (("paper",) if bridge else ("paper", "sds")):
        folder = tmp_path / representation
        folder.mkdir()
        h5py = synthetic(folder, representation=representation)
        (folder / "param").write_text(adapter.render_param(physical, NUMERICAL, representation=representation))
        # The second site on each chain is a disconnected spectator, energy -sqrt(5).
        arrays = ([0., 999.], [0., 777.], [0., 2.], [0., 2.], [-1., -1.],
                  [-1., -1.] if representation == "paper" else [1., 1.], [0., 0.], [0., 0.])
        lines = []
        for index, array in enumerate(arrays):
            if index in (0, 4):
                lines.append("z" if index == 0 else "Z")
            lines += ["1", *(str(x) for x in array)]
        (folder / "data").write_text("\n".join(lines) + "\n")
        H, ops, spins = tiny_hamiltonian(physical, representation)
        raw_ground = float(np.linalg.eigvalsh(H)[0]) - 2 * math.sqrt(5)
        expected = {}
        with h5py.File(folder / "raw.h5", "w") as h5:
            h5["stats/GS_energy"] = raw_ground
            h5["stats/1/energy_offset"] = raw_ground
            for spin in (-1, 0, 1):
                indices = np.flatnonzero(spins == spin)
                energies, vectors = np.linalg.eigh(H[np.ix_(indices, indices)])
                energies, vectors = energies[:3], vectors[:, :3]
                sector = str(spin + 1)
                h5[f"1/eigen/{sector}/absenergy"] = energies - 2 * math.sqrt(5)
                h5[f"1/eigen/{sector}/absenergy_zero"] = energies - 2 * math.sqrt(5) - raw_ground
                names = ("n_d", "n_d_ud", "SZd", "SZ2d", "asqJ") if representation == "paper" else (
                    "n_d", "n_d_ud", "SZd", "SZ2d", "hop0", "hop1")
                for name in names:
                    matrix = vectors.conj().T @ ops[name][np.ix_(indices, indices)] @ vectors
                    path = f"1/s/{name}/{sector}/{sector}/matrix"
                    h5[path], h5[path + "-imag"] = matrix.real, matrix.imag
                shifted = []
                for sign in (-1, 1):
                    h, _, _ = tiny_hamiltonian(dict(physical, phi=physical["phi"] + sign * 1e-5), representation)
                    shifted.append(np.linalg.eigvalsh(h[np.ix_(indices, indices)])[:3])
                expected[spin] = (shifted[1] - shifted[0]) / 2e-5
        result = adapter.extract_result(folder)
        results[representation] = result
        for row in result["states"]:
            assert row["values"]["current"] == pytest.approx(expected[row["twice_sz"]][row["root"]], abs=2e-9)
    if not bridge:
        for paper, sds in zip(results["paper"]["states"], results["sds"]["states"], strict=True):
            assert paper["id"] == sds["id"]
            assert paper["values"] == pytest.approx(sds["values"], abs=2e-12)


def test_local_model_source_contract():
    model = (HERE / "model.m").read_text()
    operators = (HERE / "modeloperators.m").read_text()
    assert model.index("def2ch[1];") < model.index("H1 =")
    assert "snegrealconstants[nu, Bz, WN, WS, phi]" in model
    assert "U/2 pow[number[d[]] - nu, 2] + Bz spinz[d[]]" in model
    assert "(gammaPol/Sqrt[2])" in model
    assert "hopphi[f[0], d[], -phi/4] + hopphi[f[1], d[], phi/4]" in model
    assert "genhop[(WN + I WS) Exp[-I phi/2], f[0], f[1], UP]" in model
    assert "genhop[(WN - I WS) Exp[-I phi/2], f[0], f[1], DO]" in model
    assert "asqCurrentSymbolic = D[Hdotlead + Hbridge, phi]" in model
    assert "/. params" not in model
    assert not re.search(r"\b(?:nu|Bz|WN|WS|phi)\s*=", model)
    for guard in ('BAND =!= "flat"', 'DISCRETIZATION =!= "Z"', 'SYMTYPE =!= "SPU1"',
                  "POLARIZED", "POL2x2", "RUNGS", '"substeps"', "lambda >= 4", 'paramnum["bcsgap"] < 0'):
        assert guard in model
    assert operators.count("/. params") == 1
    assert operators.count('mtSingletOp["asqJ",') == 1
    assert "Block[{MPVCFAST = False}" in operators and "opjc" not in model + operators
    assert operators.rstrip().endswith("]")  # The file returns the table, not Null.
    assert "ExpandFileName[$InputFileName]" in model and "ExpandFileName[$InputFileName]" in operators


def test_standalone_stdlib_import_and_package_extraction(tmp_path, monkeypatch):
    for name in ("physics.py", "nrg.py"):
        shutil.copy2(HERE / name, tmp_path / name)
    environment = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}
    code = f"import nrg, physics; print(nrg.render_param(physics.PHYSICAL, {NUMERICAL!r}))"
    result = subprocess.run([sys.executable, "-B", "-S", "-c", code], cwd=tmp_path, env=environment,
                            capture_output=True, text=True, check=True)
    assert "model=model.m" in result.stdout and "Lambda=4" in result.stdout
    from NRG_comparisons.test4 import nrg

    monkeypatch.delitem(sys.modules, "physics", raising=False)
    synthetic(tmp_path)
    assert nrg.extract_result(tmp_path)["physical"] == pytest.approx(PHYSICAL)
