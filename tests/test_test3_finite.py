"""Bounded longitudinal test3 controls; never launch a campaign or NRG."""

import ast
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pytest

from NRG_comparisons.test3.finite import run_finite
from NRG_comparisons.test3.manual_run import physics
from NRG_comparisons.test3.manual_run.physics import OBSERVABLES, PHYSICAL, PRIMARY, finish_states, state_values


BATH = dict(xi=[.23], weights=[.019], delta=1., bandwidth=100.)
TWO_LEVELS = dict(xi=[-.7, 1.1], weights=[.013, .027], delta=1., bandwidth=100.)


def indexed(result):
    return {(row["twice_sz"], row["root"]): row for row in result["states"]}


def test_shared_physics_contract_is_standalone():
    assert PHYSICAL == dict(geometry="single", gap=1., bandwidth=100., u=2., gamma=.4,
                            detuning=.5, field=.03)
    assert PRIMARY == {0: "even", -1: "odd_down", 1: "odd_up"}
    assert OBSERVABLES == ("energy", "P0", "P1", "P2", "P_up", "P_down", "m_dot", "m_bath",
                           "dot_z2", "bath_z2", "Czz", "Czz_connected")
    tree = ast.parse(Path(physics.__file__).read_text())
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module)
    assert imports <= {"math", "numbers"}


def test_longitudinal_values_preserve_direct_measurements():
    values = state_values(.8, .1, .12, .15, .5,
                          m_bath=.38+1e-10, Czz=-.09+2e-10, bath_z2=.28-1e-10)
    assert "energy" not in values
    assert_allclose([values[name] for name in ("P0", "P1", "P2", "P_up", "P_down")],
                    [.3, .6, .1, .42, .18], atol=1e-15)
    assert values["m_bath"] == .38+1e-10
    assert values["Czz"] == -.09+2e-10
    assert values["bath_z2"] == .28-1e-10
    assert values["Czz_connected"] == values["Czz"]-values["m_dot"]*values["m_bath"]
    derived = state_values(.8, .1, .12, .15, .5)
    assert_allclose(derived["Czz_connected"], .12**2-.15, atol=1e-15)
    assert derived["Czz_connected"] < 0


@pytest.mark.parametrize("arguments,keywords,match", [
    ((float("nan"), .1, .1, .1, 0), {}, "finite"),
    ((.8+.01j, .1, .1, .15, 0), {}, "real"),
    ((.8, .1, .4, .15, 0), {}, "probabilities"),
    ((.8, -.1, .1, .25, 0), {}, "probabilities"),
    ((.8, .1, .1, .2, 0), {}, "dot_z2"),
    ((.8, .1, .1, .15, 0), {"m_bath": .1}, "fixed-M"),
    ((.8, .1, .1, .15, 0), {"Czz": 0.}, "fixed-M"),
    ((.8, .1, .1, .15, 0), {"bath_z2": 0.}, "fixed-M"),
    ((.8, .1, .1, .15, 0), {"Czz": float("inf")}, "finite"),
])
def test_invalid_longitudinal_values(arguments, keywords, match):
    with pytest.raises(ValueError, match=match):
        state_values(*arguments, **keywords)


@pytest.mark.parametrize("bath", [BATH, TWO_LEVELS,
                                  dict(xi=[-.8, .2, 1.3], weights=[.01, .02, .03], delta=1., bandwidth=100.)])
def test_full_electron_ed_and_qp_all_roots_and_sectors(bath):
    orbitals = len(bath["xi"])+1
    task = dict(spins=list(range(-orbitals, orbitals+1)), roots=300)
    ed = run_finite(dict(task, backend="ed"), PHYSICAL, bath)
    qp = run_finite(dict(task, backend="qp", cutoff=None), PHYSICAL, bath)
    assert len(ed["states"]) == 4**orbitals
    assert len(qp["states"]) == len(ed["states"])
    assert ed["qualified"] and qp["qualified"]
    assert ed["coverage"]["sectors"] == task["spins"]
    assert qp["diagnostics"]["all_returned_states_converged"]
    for key, expected in indexed(ed).items():
        actual = indexed(qp)[key]
        assert actual["values"]["energy"] == pytest.approx(expected["values"]["energy"], abs=3e-12)
        assert actual["degenerate"] == expected["degenerate"]
        if not actual["degenerate"]:
            # Small, complete dense sectors: roundoff, not discretization tolerance.
            assert_allclose([actual["values"][name] for name in OBSERVABLES],
                            [expected["values"][name] for name in OBSERVABLES], rtol=0, atol=3e-11)
        assert actual["residual"] < 3e-12
        assert actual["values"]["Czz_connected"] == pytest.approx(
            actual["values"]["m_dot"]**2-actual["values"]["dot_z2"], abs=2e-13)
    json.dumps(qp, allow_nan=False)


@pytest.mark.parametrize("backend", ["ed", "qp"])
def test_atomic_field_sign_energies_and_three_branch_coverage(backend):
    result = run_finite(dict(backend=backend), dict(PHYSICAL, gamma=0), BATH)
    rows = {row["id"]: row for row in result["states"]}
    assert result["coverage"]["expected_three_subgap_states"]
    assert result["coverage"]["ground_ids"] == ["odd_down"]
    assert set(result["coverage"]["subgap_ids"]) == set(PRIMARY.values())
    for name, energy, moment in (("even", .5, 0.), ("odd_down", -.015, -.5), ("odd_up", .015, .5)):
        assert rows[name]["values"]["energy"] == pytest.approx(energy, abs=1e-13)
        assert rows[name]["values"]["m_dot"] == pytest.approx(moment, abs=1e-13)
        assert rows[name]["values"]["m_bath"] == pytest.approx(0, abs=1e-13)
        assert rows[name]["values"]["Czz"] == pytest.approx(0, abs=1e-13)
    assert rows["odd_up"]["excitation"] == pytest.approx(.03, abs=1e-13)
    assert rows["odd_down"]["values"]["P_down"] == pytest.approx(1)
    assert rows["even"]["values"]["P0"] == pytest.approx(1)


@pytest.mark.parametrize("backend", ["ed", "qp"])
def test_spin_reversal_at_finite_local_field(backend):
    positive = indexed(run_finite(dict(backend=backend, roots=3), PHYSICAL, TWO_LEVELS))
    negative = indexed(run_finite(dict(backend=backend, roots=3), dict(PHYSICAL, field=-.03), TWO_LEVELS))
    for (spin, root), row in positive.items():
        a, b = row["values"], negative[-spin, root]["values"]
        for name in ("energy", "P0", "P1", "P2", "dot_z2", "bath_z2", "Czz", "Czz_connected"):
            assert a[name] == pytest.approx(b[name], abs=2e-11)
        for name in ("m_dot", "m_bath"):
            assert a[name] == pytest.approx(-b[name], abs=2e-11)
        assert a["P_up"] == pytest.approx(b["P_down"], abs=2e-11)
    # Local field polarizes even M=0; imposing a singlet/SU(2) value would be wrong.
    assert abs(positive[0, 0]["values"]["m_dot"]) > 1e-4


@pytest.mark.parametrize("step", [2e-4, 1e-4])
@pytest.mark.parametrize("parameter", ["field", "detuning", "u"])
def test_hellmann_feynman_two_step_sizes(step, parameter):
    task = dict(backend="ed", roots=2, spins=[0, -1, 1])
    measured = indexed(run_finite(dict(task, backend="qp"), PHYSICAL, TWO_LEVELS))
    minus = indexed(run_finite(task, dict(PHYSICAL, **{parameter: PHYSICAL[parameter]-step}), TWO_LEVELS))
    plus = indexed(run_finite(task, dict(PHYSICAL, **{parameter: PHYSICAL[parameter]+step}), TWO_LEVELS))
    for key, row in measured.items():
        values = row["values"]
        expected = dict(field=values["m_dot"], detuning=values["P2"]-values["P0"],
                        u=(values["P0"]+values["P2"])/2)[parameter]
        derivative = (plus[key]["values"]["energy"]-minus[key]["values"]["energy"])/(2*step)
        # O(step^2) centered differences, tested on every requested root.
        assert derivative == pytest.approx(expected, abs=3e-8)


@pytest.mark.parametrize("backend", ["ed", "qp"])
def test_degeneracy_outside_returned_slice_and_selected_sector_controls(backend):
    result = run_finite(dict(backend=backend, spins=[0], roots=1),
                        dict(PHYSICAL, gamma=0, u=.4, detuning=0), BATH)
    assert len(result["states"]) == 1
    assert result["states"][0]["degenerate"]
    assert result["coverage"]["degenerate_subgap_ids"] == ["even"]
    assert result["coverage"]["sectors"] == [0]
    assert not result["coverage"]["expected_three_subgap_states"]
    assert not result["qualified"]  # Not all three target branches were selected.
    assert result["diagnostics"]["all_returned_states_converged"]
    assert result["diagnostics"]["sectors"][0]["complete_sector_diagonalized"]


def test_finish_states_strict_threshold_cross_sector_degeneracy_and_validation():
    rows = []
    for spin, energy in ((0, 0.), (-1, 0.), (1, 1.)):
        values = state_values(0, 0, 0, 0, spin/2)
        rows.append(dict(twice_sz=spin, parity=spin % 2, root=0,
                         values=dict(values, energy=energy), degenerate=False, residual=None))
    coverage = finish_states(rows)
    assert coverage["ground_ids"] == ["even", "odd_down"]
    assert coverage["subgap_ids"] == ["even", "odd_down"]
    assert coverage["degenerate_subgap_ids"] == []  # Symmetry resolves cross-sector ties.
    assert not rows[2]["subgap"]
    assert "complete" not in coverage
    with pytest.raises(ValueError, match="at least one"):
        finish_states([])
    with pytest.raises(ValueError, match="duplicate"):
        finish_states(rows+[deepcopy(rows[0])])
    broken = deepcopy(rows)
    broken[0]["values"]["P_up"] = .01
    with pytest.raises(ValueError, match="inconsistent"):
        finish_states(broken)


def test_sparse_roots_and_cutoff_composite_observables():
    exact = run_finite(dict(backend="ed", roots=2), PHYSICAL, TWO_LEVELS)
    sparse = run_finite(dict(backend="qp", roots=2, method="sparse"), PHYSICAL, TWO_LEVELS)
    for a, b in zip(exact["states"], sparse["states"], strict=True):
        assert_allclose([a["values"][name] for name in OBSERVABLES],
                        [b["values"][name] for name in OBSERVABLES], atol=2e-11, rtol=0)
    assert sparse["diagnostics"]["sectors"][0]["solver_metadata"]["method"] == "sparse"
    for cutoff in (0, 1, 2):
        truncated = run_finite(dict(backend="qp", cutoff=cutoff), PHYSICAL, TWO_LEVELS)
        for row in truncated["states"]:
            v = row["values"]
            assert v["m_dot"]+v["m_bath"] == pytest.approx(row["twice_sz"]/2, abs=2e-13)
            assert v["Czz_connected"] == pytest.approx(v["m_dot"]**2-v["dot_z2"], abs=2e-13)
        assert truncated["diagnostics"]["residual_scope"] == "projected_cutoff_hamiltonian"
    assert truncated["qualified"]


def test_energy_gap_units_and_unnormalized_bath():
    original = run_finite(dict(backend="ed"), PHYSICAL, TWO_LEVELS)
    physical = {name: (value*2 if name != "geometry" else value) for name, value in PHYSICAL.items()}
    bath = dict(TWO_LEVELS, xi=[2*x for x in TWO_LEVELS["xi"]], delta=2., bandwidth=200.)
    for backend in ("ed", "qp"):
        scaled = run_finite(dict(backend=backend), physical, bath)
        for a, b in zip(original["states"], scaled["states"], strict=True):
            assert_allclose([a["values"][name] for name in OBSERVABLES],
                            [b["values"][name] for name in OBSERVABLES], atol=2e-12, rtol=0)
    assert TWO_LEVELS["weights"] == [.013, .027]
    assert scaled["diagnostics"]["sectors"][0]["solver_metadata"]["options"]["residual_tolerance"] == 2e-8


@pytest.mark.parametrize("task,match,error", [
    ({"backend": "ed", "cutoff": 1}, "only QP", ValueError),
    ({"backend": "qp", "cutoff": 3}, "cutoff", ValueError),
    ({"backend": "qp", "spins": [0, 0]}, "spins", ValueError),
    ({"backend": "qp", "roots": 0}, "roots", ValueError),
    ({"backend": "qp", "max_memory_gib": 9}, "resource gate", ValueError),
    ({"backend": "qp", "max_dimension": 1}, "resource gate", MemoryError),
    ({"backend": "dmrg", "chi": 256}, "DMRG pilot", ValueError),
    ({"backend": "dmrg", "max_sweeps": 41}, "DMRG pilot", ValueError),
])
def test_resource_and_task_gates(task, match, error):
    with pytest.raises(error, match=match):
        run_finite(task, PHYSICAL, BATH)


def test_ed_size_gate_and_bath_mismatch():
    with pytest.raises(ValueError, match="three signed"):
        run_finite(dict(backend="ed"), PHYSICAL,
                   dict(BATH, xi=[0., 1., 2., 3.], weights=[.01]*4))
    with pytest.raises(ValueError, match="gap/bandwidth differ"):
        run_finite(dict(backend="qp"), PHYSICAL, dict(BATH, bandwidth=2.))


@pytest.mark.dmrg
def test_dmrg_exact_small_bath_every_root():
    pytest.importorskip("tenpy")
    task = dict(roots=2)
    expected = indexed(run_finite(dict(task, backend="ed"), PHYSICAL, BATH))
    result = run_finite(dict(task, backend="dmrg", chi=32, max_sweeps=10), PHYSICAL, BATH)
    assert result["qualified"]
    assert result["diagnostics"]["all_returned_states_converged"]
    for key, row in indexed(result).items():
        assert row["finite_converged"]
        assert row["residual"] < 1e-9
        assert_allclose([row["values"][name] for name in OBSERVABLES],
                        [expected[key]["values"][name] for name in OBSERVABLES], atol=2e-9, rtol=0)
    options = result["diagnostics"]["sectors"][0]["solver_metadata"]["options"]
    assert options["seed_trials"] == 2 and options["threads"] == 1
    assert options["min_sweeps"] >= 8 and options["max_sweeps"] <= 40
    assert options["lanczos_probability_tolerance"] == 1e-22
    assert options["energy_tolerance"] == 1e-12
    assert options["residual_tolerance"] == 1e-7
    assert options["calculate_residuals"] and not options["require_convergence"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("failure", ["exclusion", "primary", "root_order"])
def test_unqualified_dmrg_states_are_preserved_with_primary_qualification(monkeypatch, failure):
    import qdjj_solver

    def fake_solve(h, *, backend, cutoff, sector, options):
        count = options["eigenpairs"]
        spin = sector.twice_sz
        md = spin/2 if abs(spin) == 1 else 0.
        charge = 1. if md else 0.
        zd2 = md*md
        obs = dict(impurity_charge=charge, double_occupancy_0=0., impurity_spin_z=md,
                   dot_z2=zd2, bathSz=spin/2-md, dotSz_bathSz=spin/2*md-zd2,
                   bathSz2=(spin/2)**2-spin*md+zd2, totSz=spin/2)
        residuals = np.array([1e-10]+[2e-4]*(count-1))
        if failure == "primary":
            residuals[0] = 2e-4
        energies = np.arange(count)+.1*abs(spin)
        if failure == "root_order":
            energies = energies[::-1]
        return SimpleNamespace(energies=energies, residuals=residuals,
                               observables={name: np.full(count, value) for name, value in obs.items()},
                               metadata=dict(roots=[dict(sweep_energy_converged=i == 0) for i in range(count)],
                                             gram=np.eye(count, dtype=complex), finite_problem_converged=False),
                               timings={})

    monkeypatch.setattr(qdjj_solver, "solve", fake_solve)
    result = run_finite(dict(backend="dmrg", chi=32), PHYSICAL, BATH)
    assert result["qualified"] == (failure == "exclusion")
    assert not result["diagnostics"]["all_returned_states_converged"]
    assert any(not row["finite_converged"] for row in result["states"])
    assert all(row["finite_converged"] == (failure == "exclusion")
               for row in result["states"] if row["id"] in PRIMARY.values())
    json.dumps(result, allow_nan=False)
