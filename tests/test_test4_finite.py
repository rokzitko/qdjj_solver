"""Bounded offline physical L/R controls; never launch NRG or a campaign."""

import ast
from copy import deepcopy
import inspect
import itertools
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pytest

from NRG_comparisons.test4 import physics
from NRG_comparisons.test4.finite import run_finite
from NRG_comparisons.test4.physics import (
    OBSERVABLES, PHYSICAL, PRIMARY, detuning, finish_states, paper_offset, state_values, validate_physical,
)


BATH = dict(xi=[.23], weights=[.019], delta=1., bandwidth=100.)
TWO_LEVELS = dict(xi=[-.7, 1.1], weights=[.013, .027], delta=1., bandwidth=100.)
COMPLEX = dict(PHYSICAL, phi=.83, wn=7., ws=5.)


def indexed(result):
    return {(row["twice_sz"], row["root"]): row for row in result["states"]}


def test_standalone_physics_and_authoritative_gate():
    assert PHYSICAL == dict(geometry="junction", gap=1., bandwidth=100., u=2., gamma=.4,
                            nu=.75, field=.03, phi=3*math.pi/5, wn=0., ws=0.)
    assert PRIMARY == {0: "even", -1: "odd_down", 1: "odd_up"}
    assert OBSERVABLES == ("energy", "P0", "P1", "P2", "P_up", "P_down", "m_dot", "m_bath",
                           "dot_z2", "bath_z2", "Czz", "Czz_connected", "energy_paper", "current")
    assert str(inspect.signature(state_values)) == "(charge, double, md, zd2, M, *, m_bath=None, Czz=None, bath_z2=None)"
    imports = set()
    for node in ast.walk(ast.parse(Path(physics.__file__).read_text())):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module)
    assert imports <= {"math", "numbers"}
    checked = validate_physical(PHYSICAL)
    assert checked == PHYSICAL and checked is not PHYSICAL
    assert detuning(checked) == .5 and paper_offset(checked) == .0625
    for u, nu in ((0., -12.), (0., 8.), (2., -2.), (2., 3.)):
        p = validate_physical(dict(PHYSICAL, u=u, nu=nu))
        assert detuning(p) == u*(1-nu)
        assert paper_offset(p) == u/2*(nu-1)**2
        if u:
            assert 1-detuning(p)/u == nu
    with pytest.raises(ValueError, match="authoritative gate"):
        validate_physical(dict(PHYSICAL, detuning=.5))


@pytest.mark.parametrize("changes", [dict(nu=float("nan")), dict(phi=float("inf")), dict(wn=1j),
                                     dict(ws=True), dict(gamma=-1), dict(gap=0), dict(u=-1),
                                     dict(geometry="single")])
def test_invalid_physical(changes):
    with pytest.raises(ValueError):
        validate_physical(dict(PHYSICAL, **changes))


def test_longitudinal_contract_retains_direct_measurements():
    values = state_values(.8, .1, .12, .15, .5,
                          m_bath=.38+1e-10, Czz=-.09+2e-10, bath_z2=.28-1e-10)
    assert_allclose([values[name] for name in ("P0", "P1", "P2", "P_up", "P_down")],
                    [.3, .6, .1, .42, .18], atol=1e-15)
    assert values["m_bath"] == .38+1e-10
    assert values["Czz"] == -.09+2e-10
    assert values["bath_z2"] == .28-1e-10
    assert values["Czz_connected"] == values["Czz"]-values["m_dot"]*values["m_bath"]
    for kwargs in (dict(m_bath=0), dict(Czz=0), dict(bath_z2=0)):
        with pytest.raises(ValueError, match="fixed-M"):
            state_values(.8, .1, .12, .15, .5, **kwargs)
    with pytest.raises(ValueError, match="probabilities"):
        state_values(.8, .1, .8, .15, .5)
    with pytest.raises(ValueError, match="dot_z2"):
        state_values(.8, .1, .12, .2, .5)


@pytest.mark.parametrize("phi", [0., .73, 3*math.pi/5, math.pi])
@pytest.mark.parametrize("wn,ws", [(0., 0.), (7., 0.), (0., 5.), (7., 5.)])
def test_complex_ed_qp_complete_one_level_spectrum_and_spin_u1(phi, wn, ws):
    p = dict(PHYSICAL, phi=phi, wn=wn, ws=ws)
    task = dict(spins=list(range(-3, 4)), roots=100)
    ed = run_finite(dict(task, backend="ed"), p, BATH)
    qp = run_finite(dict(task, backend="qp", cutoff=4), p, BATH)
    assert len(ed["states"]) == len(qp["states"]) == 64
    assert qp["diagnostics"]["bath_modes_total"] == 4
    for key, expected in indexed(ed).items():
        actual = indexed(qp)[key]
        assert actual["values"]["energy"] == pytest.approx(expected["values"]["energy"], abs=3e-12)
        assert actual["degenerate"] == expected["degenerate"]
        if not actual["degenerate"]:
            assert_allclose([actual["values"][name] for name in OBSERVABLES],
                            [expected["values"][name] for name in OBSERVABLES], rtol=0, atol=3e-10)
        v = actual["values"]
        assert isinstance(v["current"], float)
        assert v["m_dot"]+v["m_bath"] == pytest.approx(key[0]/2, abs=2e-13)
        assert v["Czz_connected"] == pytest.approx(v["m_dot"]**2-v["dot_z2"], abs=2e-13)
        assert actual["residual"] < 3e-12
    # Sum linear observables over a whole sector, also testing degenerate subspaces.
    linear = [name for name in OBSERVABLES if name != "Czz_connected"]
    for spin in task["spins"]:
        assert_allclose([sum(r["values"][name] for r in ed["states"] if r["twice_sz"] == spin) for name in linear],
                        [sum(r["values"][name] for r in qp["states"] if r["twice_sz"] == spin) for name in linear],
                        rtol=0, atol=3e-11)
    json.dumps(qp, allow_nan=False)


def test_two_levels_per_lead_all_1024_kets():
    task = dict(spins=list(range(-5, 6)), roots=300)
    ed = indexed(run_finite(dict(task, backend="ed"), COMPLEX, TWO_LEVELS))
    qp = indexed(run_finite(dict(task, backend="qp", cutoff=8), COMPLEX, TWO_LEVELS))
    assert len(ed) == len(qp) == 1024
    for key, row in ed.items():
        assert qp[key]["values"]["energy"] == pytest.approx(row["values"]["energy"], abs=4e-12)
        assert qp[key]["degenerate"] == row["degenerate"]
        if not row["degenerate"]:
            assert_allclose([qp[key]["values"][name] for name in OBSERVABLES],
                            [row["values"][name] for name in OBSERVABLES], rtol=0, atol=2e-9)


@pytest.mark.parametrize("backend", ["ed", "qp"])
@pytest.mark.parametrize("u,nu", [(2., .75), (.4, 1.3), (0., -2.), (0., 3.), (0., 1e200)])
def test_atomic_nu_charging_constant_both_bath_vacua_and_field(backend, u, nu):
    p = dict(PHYSICAL, gamma=0., u=u, nu=nu)
    result = run_finite(dict(backend=backend, spins=list(range(-3, 4)), roots=100), p, BATH)
    # Independent product spectrum: each lead contributes two spinful BCS QPs.
    atomic = [(0, 0), (1, 1), (1, -1), (2, 0)]
    energies = {spin: [] for spin in range(-3, 4)}
    for (n, md2), (nl, ml2), (nr, mr2) in itertools.product(atomic, repeat=3):
        e = u/2*(n-nu)*(n-nu)+p["field"]*md2/2+(nl+nr)*math.hypot(BATH["xi"][0], p["gap"])
        energies[md2+ml2+mr2].append(e)
    for spin, expected in energies.items():
        rows = [r for r in result["states"] if r["twice_sz"] == spin]
        assert_allclose([r["values"]["energy_paper"] for r in rows], sorted(expected), atol=5e-12, rtol=0)
        for row in rows:
            v = row["values"]
            assert v["energy_paper"]-v["energy"] == pytest.approx(u/2*(nu-1)*(nu-1), abs=2e-13)
            assert v["current"] == pytest.approx(0., abs=1e-13)
    rows = indexed(result)
    assert rows[-1, 0]["values"]["P_down"] == pytest.approx(1.)
    assert rows[1, 0]["values"]["P_up"] == pytest.approx(1.)
    assert rows[1, 0]["values"]["energy"]-rows[-1, 0]["values"]["energy"] == pytest.approx(.03)


@pytest.mark.parametrize("backend", ["ed", "qp"])
@pytest.mark.parametrize("transform", ["spin", "conjugation", "time_reversal"])
def test_phase_and_spin_reversal_with_correct_ws_sign(backend, transform):
    changes = dict(spin=dict(field=-COMPLEX["field"], ws=-COMPLEX["ws"]),
                   conjugation=dict(phi=-COMPLEX["phi"], ws=-COMPLEX["ws"]),
                   time_reversal=dict(phi=-COMPLEX["phi"], field=-COMPLEX["field"]))[transform]
    task = dict(backend=backend, roots=2)
    positive = indexed(run_finite(task, COMPLEX, BATH))
    negative = indexed(run_finite(task, dict(COMPLEX, **changes), BATH))
    spin_sign = 1 if transform == "conjugation" else -1
    current_sign = 1 if transform == "spin" else -1
    for (spin, root), row in positive.items():
        a, b = row["values"], negative[spin_sign*spin, root]["values"]
        for name in ("energy", "energy_paper", "P0", "P1", "P2", "dot_z2", "bath_z2", "Czz", "Czz_connected"):
            assert a[name] == pytest.approx(b[name], abs=2e-11)
        for name in ("m_dot", "m_bath"):
            assert a[name] == pytest.approx(spin_sign*b[name], abs=2e-11)
        assert a["current"] == pytest.approx(current_sign*b["current"], abs=2e-11)
        assert a["P_up"] == pytest.approx(b["P_down" if spin_sign == -1 else "P_up"], abs=2e-11)
    # Local Bz and the ASQ bridge preserve M, not a singlet/SU(2) assignment.
    assert abs(positive[0, 0]["values"]["m_dot"]) > 1e-4


@pytest.mark.parametrize("step", [2e-4, 1e-4])
@pytest.mark.parametrize("p", [PHYSICAL, COMPLEX, dict(COMPLEX, gamma=0., ws=0.)])
def test_current_hellmann_feynman_including_bridge_without_dot_contacts(step, p):
    task = dict(spins=[0, -1, 1], roots=2)
    measured = indexed(run_finite(dict(task, backend="qp"), p, BATH))
    minus = indexed(run_finite(dict(task, backend="ed"), dict(p, phi=p["phi"]-step), BATH))
    plus = indexed(run_finite(dict(task, backend="ed"), dict(p, phi=p["phi"]+step), BATH))
    for key, row in measured.items():
        derivative = (plus[key]["values"]["energy"]-minus[key]["values"]["energy"])/(2*step)
        # Centered finite difference O(step^2), not a continuum tolerance.
        assert derivative == pytest.approx(row["values"]["current"], abs=3e-8)
    assert abs(measured[0, 0]["values"]["current"]) > 1e-5


@pytest.mark.parametrize("step", [2e-4, 1e-4])
def test_nu_derivative_centered_and_paper_references(step):
    task = dict(backend="ed", spins=[0, -1, 1], roots=2)
    measured = indexed(run_finite(dict(task, backend="qp"), COMPLEX, BATH))
    minus = indexed(run_finite(task, dict(COMPLEX, nu=COMPLEX["nu"]-step), BATH))
    plus = indexed(run_finite(task, dict(COMPLEX, nu=COMPLEX["nu"]+step), BATH))
    for key, row in measured.items():
        v = row["values"]
        for name, expected in (("energy", -COMPLEX["u"]*(v["P2"]-v["P0"])),
                               ("energy_paper", -COMPLEX["u"]*(v["P1"]+2*v["P2"]-COMPLEX["nu"]))):
            derivative = (plus[key]["values"][name]-minus[key]["values"][name])/(2*step)
            assert derivative == pytest.approx(expected, abs=3e-7)


@pytest.mark.parametrize("backend", ["ed", "qp"])
def test_units_bath_weights_and_input_immutability(backend):
    p, bath, task = deepcopy(COMPLEX), deepcopy(TWO_LEVELS), dict(backend=backend)
    original = run_finite(task, p, bath)
    dimensional = {"gap", "bandwidth", "u", "gamma", "field", "wn", "ws"}
    scaled_p = {name: value*2 if name in dimensional else value for name, value in p.items()}
    scaled_bath = dict(bath, xi=[2*x for x in bath["xi"]], delta=2., bandwidth=200.)
    scaled = run_finite(task, scaled_p, scaled_bath)
    assert scaled_p["phi"] == p["phi"] and scaled_p["nu"] == p["nu"]
    for a, b in zip(original["states"], scaled["states"], strict=True):
        assert_allclose([a["values"][name] for name in OBSERVABLES],
                        [b["values"][name] for name in OBSERVABLES], atol=3e-11, rtol=0)
    assert p == COMPLEX and bath == TWO_LEVELS and task == dict(backend=backend)
    assert original["diagnostics"]["bath_levels_per_lead"] == 2
    assert original["diagnostics"]["bath_weight_sum_per_lead"] == .04
    # Equivalent rescaling of the unnormalized contact measure, including bridge.
    changed = run_finite(task, dict(p, gamma=p["gamma"]/2, wn=p["wn"]/2, ws=p["ws"]/2),
                         dict(bath, weights=[2*w for w in bath["weights"]]))
    assert_allclose([r["values"]["energy"] for r in changed["states"]],
                    [r["values"]["energy"] for r in original["states"]], atol=2e-12, rtol=0)


@pytest.mark.parametrize("backend", ["ed", "qp"])
def test_per_sector_roots_extra_even_subgap_and_degeneracy_lookahead(backend):
    task = dict(backend=backend, roots_by_spin={"0": 3, "-1": 2, "1": 2, "-2": 1, "2": 1})
    result = run_finite(task, dict(PHYSICAL, gamma=0., u=.4, nu=.8), BATH)
    assert len(result["states"]) == 9
    assert set(result["coverage"]["subgap_ids"]) == {*PRIMARY.values(), "sz+0-r1"}
    assert not result["coverage"]["expected_three_subgap_states"]
    assert result["coverage"]["expected_three_subgap_states_is_hypothesis"]
    assert result["qualified"]  # Extra even state does not disqualify accurate primary roots.
    assert result["coverage"]["ground_ids"] == ["odd_down"]
    tied = run_finite(dict(backend=backend, spins=[0], roots=1),
                      dict(PHYSICAL, gamma=0., u=.4, nu=1.), BATH)
    assert tied["states"][0]["degenerate"]
    assert tied["coverage"]["degenerate_subgap_ids"] == ["even"]
    assert not tied["qualified"]
    assert tied["diagnostics"]["sectors"][0]["complete_sector_diagonalized"]


def test_finish_states_threshold_same_m_only_and_offset_validation():
    rows = []
    for spin, root, energy in ((0, 0, 0.), (-1, 0, 0.), (1, 0, 1.), (0, 1, .5), (-2, 0, -.1)):
        values = dict(state_values(0., 0., 0., 0., spin/2), energy=energy,
                      energy_paper=energy+.0625, current=.1+1e-12j)
        rows.append(dict(twice_sz=spin, parity=spin % 2, root=root, values=values,
                         degenerate=False, residual=None))
    coverage = finish_states(rows)
    assert coverage["ground_ids"] == ["sz-2-r0"]  # Global among all computed sectors, not primary only.
    assert coverage["degenerate_subgap_ids"] == []
    assert "sz+0-r1" in coverage["subgap_ids"] and "odd_up" not in coverage["subgap_ids"]
    assert all(type(r["values"]["current"]) is float for r in rows)
    rows[2]["values"].update(energy=.9, energy_paper=.9625)
    assert "odd_up" not in finish_states(rows)["subgap_ids"]  # Exactly 1, strict threshold.
    tied = deepcopy(rows[0])
    tied["root"] = 2
    assert set(finish_states(rows+[tied])["degenerate_subgap_ids"]) == {"even", "sz+0-r2"}
    with pytest.raises(ValueError, match="at least one"):
        finish_states([])
    with pytest.raises(ValueError, match="duplicate"):
        finish_states(rows+[deepcopy(rows[0])])
    for name, value, match in (("energy_paper", 2., "offset"), ("current", .1j, "real"),
                               ("current", float("nan"), "finite"), ("P_up", .01, "inconsistent")):
        broken = deepcopy(rows)
        broken[0]["values"][name] = value
        with pytest.raises(ValueError, match=match):
            finish_states(broken)


def test_sparse_lookahead_and_combined_cutoff_composite_observables():
    expected = run_finite(dict(backend="ed"), COMPLEX, TWO_LEVELS)
    sparse = run_finite(dict(backend="qp", method="sparse"), COMPLEX, TWO_LEVELS)
    for a, b in zip(expected["states"], sparse["states"], strict=True):
        assert_allclose([a["values"][name] for name in OBSERVABLES],
                        [b["values"][name] for name in OBSERVABLES], atol=3e-11, rtol=0)
    for info in sparse["diagnostics"]["sectors"]:
        assert info["solved_roots"] == info["returned_roots"]+1
        assert info["solver_metadata"]["method"] == "sparse"
        assert not info["complete_sector_diagonalized"]
    for cutoff in (0, 1, 2, 3, 4):
        result = run_finite(dict(backend="qp", cutoff=cutoff), COMPLEX, BATH)
        assert result["diagnostics"]["residual_scope"] == "projected_cutoff_hamiltonian"
        for row in result["states"]:
            v = row["values"]
            assert v["m_dot"]+v["m_bath"] == pytest.approx(row["twice_sz"]/2, abs=2e-13)
            assert v["Czz_connected"] == pytest.approx(v["m_dot"]**2-v["dot_z2"], abs=2e-13)


@pytest.mark.parametrize("task,match,error", [
    ({"backend": "other"}, "backend", ValueError),
    ({"backend": "ed", "cutoff": 1}, "only QP", ValueError),
    ({"backend": "qp", "cutoff": 5}, "combined L/R", ValueError),
    ({"backend": "qp", "spins": [0, 0]}, "spins", ValueError),
    ({"backend": "qp", "roots": 0}, "roots", ValueError),
    ({"backend": "qp", "roots_by_spin": {0: 3}}, "string keys", ValueError),
    ({"backend": "qp", "roots_by_spin": {"0": 0}}, "positive integer", ValueError),
    ({"backend": "qp", "roots_by_spin": {"3": 1}}, "selected spin", ValueError),
    ({"backend": "qp", "max_dimension": 1_000_001}, "resource gate", ValueError),
    ({"backend": "qp", "max_memory_gib": 9}, "resource gate", ValueError),
    ({"backend": "qp", "max_nnz": 50_000_001}, "resource gate", ValueError),
    ({"backend": "qp", "method": "matrix-free"}, "method", ValueError),
    ({"backend": "dmrg", "chi": 256}, "DMRG control", ValueError),
    ({"backend": "dmrg", "max_sweeps": 41}, "DMRG control", ValueError),
])
def test_task_resource_gates(task, match, error):
    with pytest.raises(error, match=match):
        run_finite(task, PHYSICAL, BATH)


def test_size_preflight_before_basis_and_bath_validation(monkeypatch):
    import qdjj_solver

    def forbidden_solve(*args, **kwargs):
        pytest.fail("resource preflight must happen before solve/basis allocation")

    monkeypatch.setattr(qdjj_solver, "solve", forbidden_solve)
    with pytest.raises(MemoryError, match="resource gate"):
        run_finite(dict(backend="qp", max_dimension=1), PHYSICAL, BATH)
    three = dict(BATH, xi=[-.8, .2, 1.3], weights=[.01, .02, .03])
    with pytest.raises(ValueError, match="two signed bath levels per lead"):
        run_finite(dict(backend="ed"), PHYSICAL, three)
    for task in (dict(method="dense"), dict(roots=4000)):
        with pytest.raises(MemoryError, match="512 states"):
            run_finite(dict(task, backend="qp"), PHYSICAL, three)
    with pytest.raises(ValueError, match="gap/bandwidth differ"):
        run_finite(dict(backend="ed"), PHYSICAL, dict(BATH, bandwidth=2.))
    for changes in (dict(weights=[0.]), dict(weights=[]), dict(xi=[float("nan")])):
        with pytest.raises(ValueError, match="explicit bath"):
            run_finite(dict(backend="ed"), PHYSICAL, dict(BATH, **changes))


@pytest.mark.dmrg
def test_dmrg_three_sites_complex_primary_and_extra_even():
    pytest.importorskip("tenpy")
    task = dict(spins=[0, -1, 1], roots=1, roots_by_spin={"0": 2})
    expected = indexed(run_finite(dict(task, backend="ed"), COMPLEX, BATH))
    result = run_finite(dict(task, backend="dmrg", chi=32, max_sweeps=10), COMPLEX, BATH)
    assert result["qualified"] and result["diagnostics"]["all_returned_states_converged"]
    for key, row in indexed(result).items():
        assert row["residual"] < 1e-8
        assert_allclose([row["values"][name] for name in OBSERVABLES],
                        [expected[key]["values"][name] for name in OBSERVABLES], atol=2e-8, rtol=0)
    options = result["diagnostics"]["sectors"][0]["solver_metadata"]["options"]
    assert options["seed_trials"] == 2 and options["threads"] == 1 and options["group_size"] == 2
    assert options["min_sweeps"] >= 8 and options["max_sweeps"] <= 40
    assert options["lanczos_probability_tolerance"] == 1e-22
    assert options["energy_tolerance"] == 1e-12 and options["residual_tolerance"] == 1e-7
    assert options["calculate_residuals"] and not options["require_convergence"]
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize("failure", ["exclusion", "primary", "root_order"])
def test_unqualified_dmrg_rows_preserved(monkeypatch, failure):
    import qdjj_solver

    def fake_solve(h, *, backend, cutoff, sector, options):
        count, spin = options["eigenpairs"], sector.twice_sz
        md = spin/2 if abs(spin) == 1 else 0.
        charge, zd2 = 1. if md else 0., md*md
        assert "rho_ws_derivative" not in h.observables
        obs = dict(impurity_charge=charge, double_occupancy_0=0., impurity_spin_z=md,
                   dot_z2=zd2, bathSz=spin/2-md, dotSz_bathSz=spin/2*md-zd2,
                   bathSz2=(spin/2)**2-spin*md+zd2, totSz=spin/2, phase_derivative=.02)
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
    result = run_finite(dict(backend="dmrg"), PHYSICAL, BATH)
    assert result["qualified"] == (failure == "exclusion")
    assert not result["diagnostics"]["all_returned_states_converged"]
    assert len(result["states"]) == 10
    assert all(row["finite_converged"] == (failure == "exclusion")
               for row in result["states"] if row["id"] in PRIMARY.values())
    json.dumps(result, allow_nan=False)
