"""Additional literature checks: subgap spectra, gate control and Kondo crossover."""

from io import BytesIO
import os
from pathlib import Path
import pickle
import re
import subprocess
import sys
import zipfile

import numpy as np
from numpy.testing import assert_allclose
import pytest

from qdjj_solver import DiscreteBath, reference_model
from applications import _reference
from applications._common import parity_states, read_csv, read_json
from applications.paaske_2023_surrogate_spectrum import run as spectrum
from applications.bargerbos_2022_parity_diagram import run as gate
from applications.choi_2004_kondo_josephson import run as kondo

ROOT = Path(__file__).resolve().parents[1]


def bath_from_output(module, levels, bandwidth=None, profile="quick"):
    manifest = read_json(module.CASE/"output"/profile/"manifest.json")
    records = [b for b in manifest["baths"].values() if len(b["xi"]) == levels and
               (bandwidth is None or np.isclose(b["bandwidth"], bandwidth))]
    assert len(records) == 1
    return DiscreteBath.from_record(records[0])


def test_published_surrogate_spectrum_and_uncoupled_edge():
    config = read_json(spectrum.CASE/"input"/"parameters.json")
    rows = read_csv(spectrum.CASE/"output"/"quick"/"spectrum.csv")
    for levels in config["profiles"]["quick"]["levels"]:
        bath = bath_from_output(spectrum, levels)
        for row in (r for r in rows if r["setting"] == f"L{levels}"):
            result, states = spectrum.point(bath, config["model"]["u"], float(row["gamma"]), {})
            for key in ("signed_gap", "doublet_dot_spin"):
                assert_allclose(result[key], float(row[key]), rtol=0, atol=2e-9)
            assert max(s.residuals[0] for s in states) < 2e-8
        # At zero coupling the even excitation contains a screening bath QP.
        result, _ = spectrum.point(bath, 15., 0., {})
        expected = min(15./2, min(bath.energies))
        assert_allclose(result["singlet_excitation"], expected, atol=1e-10)
        assert (expected == 1.) == bool(levels % 2)


def test_published_vector_surrogate_points():
    reference = read_csv(spectrum.CASE/"reference"/"figure3.csv")
    for levels in (1, 2, 3):
        row = next(r for r in reference if r["method"] == f"L{levels}" and
                   r["branch"] == "S" and .4 < float(r["gamma"]) < 1.)
        result, _ = spectrum.point(bath_from_output(spectrum, levels), 15., float(row["gamma"]), {})
        # Published vector coordinates have finite extraction precision.
        assert_allclose(result["singlet_excitation"], float(row["excitation"]), atol=.001, rtol=0)


def test_gerade_doublet_is_not_the_free_spectator():
    bath = bath_from_output(spectrum, 1)
    result, _ = spectrum.point(bath, 15., 10., {})
    junction = reference_model(bath, u=15., gamma=10., phi=0., symmetry=False, compress=False)
    states = parity_states(junction, {})
    assert result["doublet_excitation"] > 1.
    # The complete two-lead odd sector instead admits a free odd-channel QP.
    assert_allclose(states[1].energies[0]-states[0].energies[0], 1., atol=1e-10)
    assert_allclose(states[0].energies[0], result["even_energy"], atol=1e-10)


def test_gate_boundary_and_charge_response():
    config = read_json(gate.CASE/"input"/"parameters.json")
    settings = config["profiles"]["quick"]["settings"]
    bath = bath_from_output(gate, settings["levels"])
    for row in read_csv(gate.CASE/"output"/"quick"/"boundary.csv"):
        center = float(row["detuning_critical_over_u"])
        phi = float(row["phi_over_pi"])*np.pi
        results = [gate.point(bath, config["model"], x, phi, settings)[0]
                   for x in (center-.01, center, center+.01)]
        assert results[0]["signed_gap"] < 0 < results[2]["signed_gap"]
        assert abs(results[1]["signed_gap"]) < 3e-6
    x, phi, step = .4, np.pi/2, 1e-4
    def at(value):
        return gate.point(bath, config["model"], value, phi, settings)[0]
    center, reflected, minus, plus = at(x), at(-x), at(x-step/5), at(x+step/5)
    for parity in ("even", "odd"):
        assert_allclose(center[f"{parity}_energy"], reflected[f"{parity}_energy"], atol=1e-9)
        assert_allclose(center[f"{parity}_charge"]+reflected[f"{parity}_charge"], 2., atol=1e-9)
        derivative = (plus[f"{parity}_energy"]-minus[f"{parity}_energy"])/(2*step)
        assert_allclose(derivative, center[f"{parity}_charge"]-1, atol=2e-8)


def test_gate_boundary_against_deposited_grid():
    config = read_json(gate.CASE/"input"/"parameters.json")
    bath = bath_from_output(gate, 4, profile="paper")
    reference = read_csv(gate.CASE/"reference"/"boundary.csv")[12]
    p = float(reference["phi_over_pi"])*np.pi
    x = float(reference["detuning_critical_over_u"])
    # Bracket the continuum-method discrepancy by one published gate-grid step.
    a = gate.point(bath, config["model"], x-.01, p, {})[0]["signed_gap"]
    b = gate.point(bath, config["model"], x+.01, p, {})[0]["signed_gap"]
    assert a < 0 < b


def test_kondo_scale_and_published_parameter_baselines():
    config = read_json(kondo.CASE/"input"/"parameters.json")
    settings = config["profiles"]["quick"]["settings"]
    rows = read_csv(kondo.CASE/"output"/"quick"/"current.csv")
    for ratio in config["profiles"]["quick"]["ratios"]:
        scales = kondo.parameters(config["normal_state"], ratio)
        assert_allclose(scales["tk_over_D"], .008877583677243962, atol=1e-15)
        assert_allclose(scales["u"]/scales["bandwidth"], .2, atol=1e-15)
        assert_allclose(scales["gamma"]/scales["bandwidth"], .04, atol=1e-15)
        bath = bath_from_output(kondo, settings["levels"], scales["bandwidth"])
        for row in (r for r in rows if float(r["delta_over_tk"]) == ratio):
            result, _ = kondo.point(bath, scales, float(row["phi_over_pi"])*np.pi, settings)
            assert_allclose(result["signed_gap"], float(row["signed_gap"]), atol=2e-8, rtol=0)
            assert_allclose(result["current"], float(row["current"]), atol=2e-8, rtol=0)
        zero = kondo.point(bath, scales, 0., settings)[0]
        pi = kondo.point(bath, scales, np.pi, settings)[0]
        assert abs(zero["current"]) < 1e-8 and abs(pi["current"]) < 1e-8
        assert pi["ground_parity"] == 1  # PH-symmetric doublet at phi=pi.
    assert all(float(r["current"]) < 0 for r in rows if float(r["delta_over_tk"]) == 10.)
    assert all(float(r["current"]) > 0 for r in rows if float(r["delta_over_tk"]) == .1)


@pytest.mark.parametrize("module", (spectrum, gate, kondo), ids=lambda m: m.CASE.name)
def test_second_set_offline_runner(module, tmp_path):
    output = tmp_path/"résultats with spaces"
    env = os.environ | {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1"}
    result = subprocess.run([sys.executable, "-B", "-m", module.__name__, "--output", str(output)],
                            cwd=ROOT, env=env, capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout+result.stderr
    manifest = read_json(output/"manifest.json")
    states = read_json(output/"eigenstates.json")
    assert manifest["state_records"] == len(states) > 0
    assert max(max(r["residuals"]) for r in states) < 2e-8
    assert manifest["baths"] and manifest["scripts_sha256"]


def test_restricted_numeric_deposit_decoder():
    data = (np.arange(3.), np.ones((2, 3)), np.float64(.2))
    recovered = _reference.numeric_pickle(pickle.dumps(data, protocol=4))
    for a, b in zip(data, recovered, strict=True):
        assert_allclose(a, b)
    with pytest.raises(ValueError, match="unexpected global"):
        _reference.numeric_pickle(b"cbuiltins\nprint\n.")
    with pytest.raises(ValueError, match="nonnumeric"):
        _reference.numeric_pickle(pickle.dumps((np.array([1], dtype=object),), protocol=4))


@pytest.mark.slow
@pytest.mark.convergence
@pytest.mark.dmrg
def test_strong_coupling_independent_dmrg():
    pytest.importorskip("tenpy")
    config = read_json(kondo.CASE/"input"/"dmrg_check.json")
    scales = kondo.parameters(config["normal_state"], .1)
    bath = bath_from_output(kondo, 8, profile="dmrg")
    h = kondo.model(bath, scales, .3*np.pi)
    settings = config["profiles"]["convergence"]["settings"][0]
    qp, dmrg = parity_states(h, {}), parity_states(h, settings)
    for exact, mps in zip(qp, dmrg, strict=True):
        assert_allclose(exact.energies, mps.energies, atol=1e-8, rtol=0)
        assert_allclose(exact.observables["phase_derivative"],
                        mps.observables["phase_derivative"], atol=1e-6, rtol=0)


def test_remote_zip_uses_bounded_verified_ranges(monkeypatch):
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("target.csv", "x,y\n1,2\n")
        archive.writestr("padding", b"x"*200_000)
    data = buffer.getvalue()
    requests = []

    def urlopen(request, timeout):
        byte_range = request.headers["Range"]
        requests.append(byte_range)
        if byte_range.startswith("bytes=-"):
            start, end = max(0, len(data)-int(byte_range[7:])), len(data)-1
        else:
            start, end = map(int, re.fullmatch(r"bytes=(\d+)-(\d+)", byte_range).groups())
        response = BytesIO(data[start:end+1])
        response.status = 206
        response.headers = {"Content-Range": f"bytes {start}-{end}/{len(data)}", "ETag": '"pinned"'}
        return response

    monkeypatch.setattr(_reference, "urlopen", urlopen)
    remote = _reference.RemoteZip("https://example.invalid/archive.zip", "pinned")
    with zipfile.ZipFile(remote) as archive:
        assert archive.read("target.csv") == b"x,y\n1,2\n"
    assert len(requests) >= 2 and remote.transferred < len(data)
    with pytest.raises(ValueError, match="ETag"):
        _reference.RemoteZip("https://example.invalid/archive.zip", "wrong")
