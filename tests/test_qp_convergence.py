"""Independent physics and portable-evidence checks for the cutoff guide."""

from copy import deepcopy
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.sparse import eye

from qdjj_solver import DiscreteBath, Impurity, Sector, reference_model, solve
from qdjj_solver.common.io import fingerprint
from electron_oracle import electron_hamiltonian, electron_operators, sector_indices
from examples.qp_convergence import add_spin_observables, values_from_state
from tools import benchmark_qp_convergence as bench
from tools.plot_qp_convergence import analyze, write_tables

ARCHIVE = Path(__file__).resolve().parents[1]/"docs/benchmarks/qp_convergence"


@pytest.mark.parametrize("detuning", [0., .23])
def test_dot_probabilities_and_complete_spin_correlation_against_electron_ed(detuning):
    """Use independent electron tensors, including transverse spin exchange.

    Nonunit bath weights test that total bath spin has no contact weights.
    Detuning tests the physical occupation probabilities when P0 != P2.
    """
    bath = DiscreteBath([-.8, .8], [.27, .27], delta=.9, bandwidth=4.3)
    parameters = dict(u=1.7, gamma=.37, phi=1.1, detuning=detuning, compress=False)
    h = reference_model(bath, **parameters)
    add_spin_observables(h)
    sector = Sector(1, 1)
    result = solve(h, cutoff=None, sector=sector, options={"method": "dense"})

    impurity = Impurity.anderson(parameters["u"], detuning)
    t = np.eye(2)*np.sqrt(parameters["gamma"]/(2*np.pi*bath.rho))
    matrix = electron_hamiltonian(impurity, [bath, bath], [t, t], [-.55, .55])
    selected = sector_indices(10, 1, 1)
    energies, vectors = np.linalg.eigh(matrix[selected][:, selected].toarray())
    vector = vectors[:, 0]
    cs = electron_operators(10)
    ns = [c.getH() @ c for c in cs]
    identity = eye(2**10, format="csr")
    p0 = (identity-ns[0]) @ (identity-ns[1])
    p2 = ns[0] @ ns[1]
    p1 = ns[0]+ns[1]-2*p2
    zd = (ns[0]-ns[1])/2
    pd = cs[0].getH() @ cs[1]
    zb = sum((ns[i]-ns[i+1])/2 for i in range(2, 10, 2))
    pb = sum(cs[i].getH() @ cs[i+1] for i in range(2, 10, 2))
    correlation = zd @ zb+(pd @ pb.getH()+pd.getH() @ pb)/2

    def expectation(op):
        return float(np.vdot(vector, op[selected][:, selected] @ vector).real)

    expected = dict(energy=energies[0], q_d=2*expectation(zd), P_0=expectation(p0),
                    P_1=expectation(p1), P_2=expectation(p2), C_d_bath=expectation(correlation))
    measured = values_from_state(result)
    for key in expected:
        assert_allclose(measured[key], expected[key], rtol=0, atol=3e-11)
    assert_allclose(result.observables["dot_bath_spin"][0], expected["C_d_bath"], rtol=0, atol=3e-11)
    assert_allclose(result.observables["total_spin_squared"][0], .75, rtol=0, atol=3e-11)
    if detuning:
        assert abs(measured["P_0"]-measured["P_2"]) > .005
    else:
        assert_allclose(measured["P_0"], measured["P_2"], rtol=0, atol=3e-11)


def test_nested_cutoffs_spin_identity_dark_mode_control_and_zeeman_derivative():
    bath = DiscreteBath([-1., 1.], [.4, .4], bandwidth=5.)
    parameters = dict(u=2., gamma=.4, phi=1.3)
    full = reference_model(bath, compress=False, **parameters)
    add_spin_observables(full)
    compressed = reference_model(bath, compress=True, **parameters)
    sector = Sector(1, 1)
    energies = []
    for cutoff in (0, 1, 2, 3, 4, None):
        state = solve(full, cutoff=cutoff, sector=sector, options={"method": "dense"})
        reduced = solve(compressed, cutoff=cutoff, sector=sector, options={"method": "dense"})
        values = values_from_state(state)
        energies.append(values["energy"])
        bench.validate_values(values, state.qp_weights[0])
        assert_allclose(state.observables["dot_bath_spin"][0], values["C_d_bath"], atol=2e-11, rtol=0)
        assert_allclose(state.observables["total_spin_squared"][0], .75, atol=2e-11, rtol=0)
        for key, value in values_from_state(reduced).items():
            assert_allclose(value, values[key], atol=2e-11, rtol=0)
    assert np.all(np.diff(energies) <= 2e-11)
    # The spin is also a physical energy response in a fixed S_z branch.
    step = 2e-5
    shifted = [solve(reference_model(bath, field=sign*step, compress=False, **parameters),
                     cutoff=2, sector=sector).energies[0] for sign in (-1, 1)]
    state = solve(full, cutoff=2, sector=sector)
    assert_allclose((shifted[1]-shifted[0])/step, values_from_state(state)["q_d"], atol=2e-8, rtol=0)
    with pytest.raises(ValueError, match="uncompressed"):
        add_spin_observables(compressed)


def test_decoupled_dot_analytic_limit():
    bath = DiscreteBath([-1., 1.], [.3, .3], bandwidth=5.)
    h = reference_model(bath, u=2., gamma=0., compress=False)
    add_spin_observables(h)
    for cutoff in (0, 2, None):
        state = solve(h, cutoff=cutoff, sector=Sector(1, 1), options={"method": "dense"})
        expected = dict(energy=0., q_d=1., P_0=0., P_1=1., P_2=0., C_d_bath=0.)
        for key, value in values_from_state(state).items():
            assert_allclose(value, expected[key], atol=2e-12, rtol=0)


@pytest.mark.parametrize("checkout_newline", ["\n", "\r\n"], ids=["lf", "crlf"])
def test_reviewed_archive_observable_resolution_and_current_solver_checks(tmp_path, checkout_newline):
    manifest, baths, records = bench.load_archive(ARCHIVE)
    assert len(records) == 15
    assert all(row["bath_key"] in baths and len(row["source"]["sha256"]) == 64 for row in records)
    summary = analyze(ARCHIVE, ARCHIVE/"current")
    assert summary == bench.read_json(ARCHIVE/"summary.json")
    assert all(s["qp_bath_change_below_resolution"] for s in summary["stability"].values())
    assert not any(summary["points"][-1]["error_resolved"].values())
    assert all(summary["points"][-2]["error_resolved"].values())
    # A small energy error does not certify the relative accuracy of screening.
    q2 = summary["points"][1]
    assert q2["absolute_errors"]["q_d"]/summary["reference_values"]["q_d"] < .02
    assert q2["absolute_errors"]["C_d_bath"]/abs(summary["reference_values"]["C_d_bath"]) > .3
    assert len(summary["current_solver_verification"]["points"]) == 8
    write_tables(tmp_path, summary)
    # Exercise both Git checkout conventions, independently of the host OS.
    for name in ("convergence.csv", "table.md"):
        checked_out = tmp_path/f"checked-out-{name}"
        checked_out.write_text((ARCHIVE/name).read_text(encoding="utf-8"),
                               encoding="utf-8", newline=checkout_newline)
        # Universal-newline decoding preserves the exact table content while
        # ignoring LF/CRLF conversion performed by Git on Windows.
        assert (tmp_path/name).read_text(encoding="utf-8") == checked_out.read_text(encoding="utf-8")
    assert manifest["input"] == bench.read_json(ARCHIVE/"input.json")


@pytest.mark.parametrize("problem", ["reference", "cutoff", "missing", "bath"])
def test_report_rejects_invalid_evidence(tmp_path, problem):
    manifest, baths, records = bench.load_archive(ARCHIVE)
    if problem == "reference":
        records[-1]["diagnostics"]["sweep_converged"] = False
    elif problem == "cutoff":
        records[5]["values"]["energy"] = 0.
    elif problem == "missing":
        records.pop(5)
    else:
        row = records[5]
        bath = deepcopy(baths[row["bath_key"]])
        bath["weights"] = [2*w for w in bath["weights"]]
        row["bath_key"] = fingerprint(bath)
        baths[row["bath_key"]] = bath
    bench.save_archive(tmp_path, manifest, baths, records)
    with pytest.raises(ValueError):
        analyze(tmp_path)


def small_input():
    config = bench.read_json(ARCHIVE/"input.json")
    config.update(levels=[2, 4], cutoffs=[1, 2], bath=dict(delta=1., bandwidth=5.),
                  references=[dict(levels=2, chi_max=16), dict(levels=4, chi_max=16), dict(levels=2, chi_max=32)])
    config["dmrg_options"].update(seed_trials=2, min_sweeps=8, max_sweeps=24, residual_tolerance=1e-7)
    return config


def test_bounded_runner_resume_export_and_input_identity(tmp_path):
    output, exported = tmp_path/"run", tmp_path/"export"
    config = small_input()
    bench.run(config, "full", output)
    original = (output/"measurements.json").read_bytes()
    bench.run(config, "full", output)
    assert (output/"measurements.json").read_bytes() == original
    bench.export_archive(output, exported)
    assert bench.load_archive(output) == bench.load_archive(exported)
    changed = deepcopy(config)
    changed["model"]["u"] = 2.1
    with pytest.raises(ValueError, match="resume requires identical"):
        bench.run(changed, "full", output)
    # The check profile contains no unrestricted references and is not a full report.
    with pytest.raises(ValueError, match="incomplete"):
        analyze(output)


@pytest.mark.dmrg
def test_reference_runner_measures_full_residuals(tmp_path):
    pytest.importorskip("tenpy")
    bench.run(small_input(), "full", tmp_path, include_dmrg=True)
    _, _, records = bench.load_archive(tmp_path)
    references = [r for r in records if r["backend"] == "dmrg"]
    assert len(references) == 3
    assert all(r["diagnostics"]["finite_problem_converged"] for r in references)
    assert max(r["residual"] for r in references) < 1e-7
    assert all(r["cutoff"] is None for r in references)
    assert analyze(tmp_path)["reference_origin"] == "current"
