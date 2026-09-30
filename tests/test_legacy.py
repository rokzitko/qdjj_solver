"""Independent legacy chain/DMRG checks and convention-matched synthetic NRG data."""

from copy import deepcopy
import json
import os
import subprocess
import sys

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.sparse.linalg import eigsh
from threadpoolctl import threadpool_limits

from qdjj_solver.common.baths import DiscreteBath, cosh_grid
from qdjj_solver.common.io import write_json
from qdjj_solver.dmrg_solver.legacy import chain_coefficients
from qdjj_solver.qp_solver import nrg


@pytest.mark.parametrize("xi, weights", [
    ([-7., -.3, 0., .8, 9.], [.2, .7, 1.1, .4, .3]),
    ([-.2, -.2, 1.3, 1.3], [.1, .4, .2, .6]),
], ids=["irregular-spectrum", "repeated-nodes"])
def test_chain_spectral_measure_and_independent_contact_moments(xi, weights):
    bath = DiscreteBath(xi, weights, delta=.7, bandwidth=10.)
    diagonal, hopping, norm, vacuum = chain_coefficients(bath)
    chain = np.diag(diagonal)+np.diag(hopping, 1)+np.diag(hopping, -1)
    eigenvalues, eigenvectors = np.linalg.eigh(chain)
    distinct = np.unique(bath.xi)
    assert len(diagonal) == len(distinct)
    assert len(hopping) == len(diagonal)-1 and np.all(hopping > 0)
    assert_allclose(eigenvalues, distinct, rtol=0, atol=3e-12)
    expected_weights = [bath.weights[bath.xi == energy].sum()/bath.weights.sum() for energy in distinct]
    assert_allclose(abs(eigenvectors[0])**2, expected_weights, rtol=0, atol=2e-12)
    assert_allclose(norm**2, bath.weights.sum(), rtol=0, atol=2e-15)
    # The contact spectral measure fixes every moment, independently of Lanczos.
    scale = max(1., np.max(abs(bath.xi)))
    for power in range(2*len(diagonal)+2):
        moment = norm**2*np.linalg.matrix_power(chain/scale, power)[0, 0]
        assert_allclose(moment, np.dot(bath.weights, (bath.xi/scale)**power), rtol=0, atol=4e-12)
    assert_allclose(vacuum, np.sum(distinct-np.hypot(distinct, bath.delta)), rtol=0, atol=1e-11)
    assert_allclose(bath.xi, xi, rtol=0, atol=0)


def _physical_reference(bath, geometry, twice_sz, rho_ws, u=1.3, gamma=.27, phi=.8):
    """ED in physical electrons, with no QP transform or MPO/Lanczos machinery."""
    from electron_oracle import electron_hamiltonian, electron_operators, sector_indices
    from qdjj_solver.common.algebra import number
    from qdjj_solver.common.models import Impurity

    leads = 1 if geometry == "single" else 2
    impurity = Impurity(1, u/2-u/2*(number(0)+number(1))+u*number(0)*number(1))
    tunneling = [np.sqrt(gamma/(np.pi*bath.rho*leads))*np.eye(2)]*leads
    direct = {(0, 1): np.diag([1j*rho_ws, -1j*rho_ws])/bath.rho} if leads == 2 else {}

    def matrix_at(phase):
        phases = [-phase/2, phase/2] if leads == 2 else [0.]
        return electron_hamiltonian(impurity, [bath]*leads, tunneling, phases, direct)

    matrix = matrix_at(phi)
    modes = 2+2*leads*bath.levels
    indices = sector_indices(modes, twice_sz % 2, twice_sz)
    restricted = matrix[indices][:, indices]
    with threadpool_limits(limits=1):
        energies, vectors = eigsh(restricted, k=1, which="SA", tol=1e-12,
                                  v0=np.linspace(.1, 1., len(indices)))
    vector = vectors[:, 0]
    assert np.linalg.norm(restricted @ vector-energies[0]*vector) < 2e-11
    cs = electron_operators(modes)
    up, down = cs[0].getH() @ cs[0], cs[1].getH() @ cs[1]

    def expectation(operator):
        return float(np.vdot(vector, operator[indices][:, indices] @ vector).real)

    derivative = expectation((matrix_at(phi+1e-5)-matrix_at(phi-1e-5))/(2e-5)) if leads == 2 else 0.
    return dict(energy=energies[0], impurity_spin_z=expectation((up-down)/2),
                impurity_charge=expectation(up+down), double_occupancy_0=expectation(up @ down),
                phase_derivative=derivative)


@pytest.mark.dmrg
@pytest.mark.parametrize("geometry, bath, twice_sz, rho_ws, measure", [
    pytest.param("junction", cosh_grid(1, bandwidth=3.), 1, 0., True, id="compressed-odd-junction"),
    pytest.param("junction", DiscreteBath([0.], [.8], bandwidth=3.), -1, .07, True, id="nonzero-WS-junction"),
    pytest.param("single", cosh_grid(1, bandwidth=3.), 0, 0., False, id="single-even-no-current-or-variance"),
])
def test_both_legacy_dmrg_coordinates_match_electron_ed(geometry, bath, twice_sz, rho_ws, measure):
    pytest.importorskip("tenpy")
    from qdjj_solver.dmrg_solver.legacy import dmrg_chain_reference, dmrg_reference

    exact = _physical_reference(bath, geometry, twice_sz, rho_ws)
    parameters = dict(u=1.3, gamma=.27, phi=.8, rho_ws=rho_ws, geometry=geometry,
                      twice_sz=twice_sz, chi=64, max_sweeps=22, energy_tolerance=1e-10,
                      calculate_variance=measure, threads=1)
    chain = dmrg_chain_reference(bath, **parameters)
    star = dmrg_reference(bath, **parameters, calculate_current=measure)
    for result in (chain, star):
        for name in ("energy", "impurity_spin_z", "impurity_charge", "double_occupancy_0"):
            assert_allclose(result[name], exact[name], rtol=0, atol=2e-9, err_msg=f"{result['method']}: {name}")
        assert 1 <= result["chi_actual"] <= result["chi_limit"] == 64
        assert result["norm_error"] < 2e-9
        assert abs(result["last_discarded_weight"]) < 1e-9
        assert abs(result["last_energy_change"]) < 1e-8
        assert_allclose(result["sweep_energies"][-1], result["energy"], rtol=0, atol=2e-10)
        assert result["bath"] == bath.record()
        assert result["seconds"] >= 0 and result["sweeps"] >= 10
        assert isinstance(result["sweep_converged"], bool)
        if measure:
            assert abs(result["variance"]) < 2e-9
        else:
            assert result["variance"] is None
    assert_allclose(chain["phase_derivative"], exact["phase_derivative"], rtol=0, atol=2e-9)
    if measure:
        assert_allclose(star["phase_derivative"], exact["phase_derivative"], rtol=0, atol=2e-9)
    else:
        assert star["phase_derivative"] is None
    leads = 1 if geometry == "single" else 2
    assert_allclose(chain["bcs_vacuum_energy"], leads*np.sum(bath.xi-bath.energies), rtol=0, atol=1e-13)
    assert_allclose(chain["energy_before_vacuum_subtraction"]-chain["bcs_vacuum_energy"], chain["energy"],
                    rtol=0, atol=1e-14)
    assert star["energy_reference"] == "isolated BCS reservoirs subtracted"
    assert star["paired_mode_compression"] is (geometry == "junction" and rho_ws == 0 and twice_sz != 0)


@pytest.mark.dmrg
@pytest.mark.parametrize("entry, options, message", [
    ("chain", {"geometry": "unknown"}, "geometry or spin sector"),
    ("chain", {"geometry": "single", "rho_ws": .1}, "requires two reservoirs"),
    ("chain", {"chi": 3}, "convergence parameters"),
    ("star", {"twice_sz": 2}, "unsupported reference sector"),
])
def test_legacy_dmrg_validation(entry, options, message):
    pytest.importorskip("tenpy")
    from qdjj_solver.dmrg_solver.legacy import dmrg_chain_reference, dmrg_reference

    function = dmrg_chain_reference if entry == "chain" else dmrg_reference
    with pytest.raises(ValueError, match=message):
        function(cosh_grid(1), **options)


@pytest.fixture
def nrg_records():
    parameters = dict(u=1.5, gamma=.2, phi=.4, rho_ws=.1, rho_wn=0., detuning=0., field=0.,
                      delta=1., bandwidth=10., geometry="2-reservoir", parity=1, twice_sz=1)
    calculation = {name: parameters[name] for name in ("u", "gamma", "phi", "rho_ws", "rho_wn", "detuning", "field")}
    calculation.update(bath=[dict(delta=1., bandwidth=10.)], geometry="2-reservoir",
                       sector=dict(parity=1, twice_sz=1), energy_reference="isolated BCS reservoirs subtracted")
    result = dict(format="bcs-qp-eigenstates", format_version=1, calculation=calculation,
                  energies=[-.4, .6], observables={"impurity_spin_z": [.25, -.2]})
    reference = dict(format="bcs-qp-nrg-reference", format_version=1, status="converged-reference",
                     conventions=dict(hybridization="total", phase_bias="phi_R-phi_L",
                                      impurity_convention="charge-square",
                                      energy_reference=calculation["energy_reference"]),
                     parameters=parameters,
                     nrg=dict(version="synthetic-test-v1", kept_states=1000, Lambda=2.,
                              temperature_over_delta=0., symmetry="parity-Sz", z_values=[.5, 1.],
                              convergence_notes="Synthetic interface fixture; no external NRG calculation."),
                     observables={"energy": dict(value=-.41, uncertainty=.01),
                                  "impurity_spin_z": dict(value=.24, uncertainty=.01)})
    return reference, result


def test_nrg_comparison_combines_errors_and_selects_root_without_mutation(nrg_records):
    reference, result = nrg_records
    before = deepcopy(nrg_records)
    comparison = nrg.compare(reference, result, .02)
    assert comparison["consistent"]
    assert comparison["parameters"] == reference["parameters"]
    assert comparison["nrg"] == reference["nrg"]
    for row in comparison["comparisons"]:
        assert_allclose(row["difference"], .01, rtol=0, atol=1e-16)
        assert row["combined_uncertainty"] == .03
        assert row["consistent"]
    second = nrg.compare(reference, result, .02, state=1)
    assert not second["consistent"]
    assert_allclose([row["difference"] for row in second["comparisons"]], [1.01, -.44], rtol=0, atol=3e-16)
    assert all(not row["consistent"] for row in second["comparisons"])
    assert nrg_records == before
    # The uncertainty interval is closed, including exact zero-error agreement.
    reference["observables"] = {"energy": {"value": result["energies"][0], "uncertainty": 0.}}
    assert nrg.compare(reference, result, 0.)["consistent"]


@pytest.mark.parametrize("path, value, message", [
    ("reference.format_version", 2, "unsupported NRG reference"),
    ("reference.status", "template", "not been marked converged"),
    ("result.format_version", 2, "unsupported QP result"),
    ("reference.conventions.hybridization", "per-lead", "convention is not matched"),
    ("reference.conventions.impurity_convention", "Anderson", "charge-square convention"),
    ("reference.nrg.version", "", "software version"),
    ("reference.nrg.kept_states", 0, "retained-state count"),
    ("reference.nrg.Lambda", 1., "Lambda must exceed one"),
    ("reference.nrg.temperature_over_delta", None, "temperature must be specified"),
    ("reference.nrg.symmetry", "", "symmetry sector"),
    ("reference.nrg.z_values", [0.], "z values must lie"),
    ("reference.nrg.convergence_notes", "", "describe the convergence"),
    ("result.calculation", {}, "built-in single-dot"),
    ("reference.parameters.u", 1.6, "parameter mismatch for u"),
    ("reference.observables.energy.uncertainty", -.01, "missing finite value"),
    ("reference.conventions.energy_reference", "absolute", "absolute-energy references"),
    ("reference.observables", {}, "no NRG observables"),
])
def test_nrg_rejects_unmatched_or_unvalidated_references(nrg_records, path, value, message):
    reference, result = nrg_records
    record = {"reference": reference, "result": result}
    *keys, field = path.split(".")
    for key in keys:
        record = record[key]
    record[field] = value
    with pytest.raises(ValueError, match=message):
        nrg.compare(reference, result, .02)


@pytest.mark.parametrize("uncertainty", [-1., np.nan])
def test_nrg_invalid_qp_error_estimate(nrg_records, uncertainty):
    with pytest.raises(ValueError, match="QP uncertainty must be nonnegative and finite"):
        nrg.compare(*nrg_records, uncertainty)


def test_nrg_invalid_state_index(nrg_records):
    with pytest.raises(ValueError, match="state index outside"):
        nrg.compare(*nrg_records, .02, state=2)


def test_nrg_missing_observable_does_not_silently_pass(nrg_records):
    reference, result = nrg_records
    del result["observables"]["impurity_spin_z"]
    with pytest.raises(KeyError, match="impurity_spin_z"):
        nrg.compare(reference, result, .02)


@pytest.mark.parametrize("module, consistent", [("qdjj_solver.qp_solver.nrg", True), ("qpsolver.nrg", False)])
def test_nrg_cli_output_and_exit_status_without_external_nrg(tmp_path, nrg_records, module, consistent):
    reference, result = nrg_records
    refpath, resultpath, output = tmp_path/"nrg.json", tmp_path/"qp.json", tmp_path/"comparison.json"
    write_json(refpath, reference)
    write_json(resultpath, result)
    args = [sys.executable, "-B", "-m", module, str(refpath), str(resultpath),
            "--qp-uncertainty", ".02", "--state", "0" if consistent else "1"]
    if consistent:
        args.extend(["--output", str(output)])
    completed = subprocess.run(args, cwd=tmp_path, capture_output=True, text=True, encoding="utf-8",
                               env=dict(os.environ, PYTHONIOENCODING="utf-8"), timeout=60)
    assert completed.returncode == (0 if consistent else 1), completed.stderr
    assert "energy: QP=" in completed.stdout and "impurity_spin_z: QP=" in completed.stdout
    assert "difference=" in completed.stdout and "combined uncertainty=3.000e-02" in completed.stdout
    if consistent:
        comparison = json.loads(output.read_text())
        assert comparison["consistent"] and len(comparison["comparisons"]) == 2
    else:
        assert not output.exists()
