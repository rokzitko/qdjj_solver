"""Small, offline analytic controls, independent of the native QP extension."""

import json
import math
from unittest.mock import Mock

import numpy as np
from numpy.testing import assert_allclose
import pytest

from NRG_comparisons.test1 import control


PHYSICAL = dict(gap=1., bandwidth=100., u=2., gamma=.4, detuning=0., field=0., geometry="single")


def bath(physical, xi=(0.,), weights=(1.,)):
    return dict(xi=list(xi), weights=list(weights), delta=physical["gap"],
                bandwidth=physical["bandwidth"], metadata={})


@pytest.mark.parametrize("nmax", [0, 2, 4, 6, 100])
def test_wilson_chain_lengths_normalization_and_units(nmax):
    physical = dict(PHYSICAL, gap=2., bandwidth=7.)
    record = control.wilson_bath(physical, 2., nmax)
    json.dumps(record, allow_nan=False)
    assert set(record) == {"xi", "weights", "delta", "bandwidth", "metadata"}
    assert len(record["xi"]) == len(record["weights"]) == nmax+1
    assert record["delta"] == 2. and record["bandwidth"] == 7.
    assert_allclose(sum(record["weights"]), 1., atol=1e-14, rtol=0)
    assert np.all(np.array(record["weights"]) > 0)
    metadata = record["metadata"]
    assert metadata["discretization"] == "Z" and metadata["z"] == 1.
    assert metadata["nmax"] == nmax and metadata["lambda"] == 2.
    assert metadata["chain_onsite"] == [0.]*(nmax+1)
    assert len(metadata["chain_hoppings"]) == nmax
    assert_allclose(metadata["impurity_hopping"], math.sqrt(2*7*.4/math.pi), atol=1e-14)
    if nmax:
        # Closed first two Z coefficients for Lambda=2, in physical D units.
        assert_allclose(metadata["chain_hoppings"][:2],
                        [7/(math.log(2)*math.sqrt(7)), 6*7/(math.log(2)*math.sqrt(434))], rtol=2e-15)
    assert physical == dict(PHYSICAL, gap=2., bandwidth=7.)


def test_one_site_wilson_bath_does_not_call_lapack(monkeypatch):
    diagonalize = Mock(side_effect=AssertionError("one-site chain reached LAPACK"))
    monkeypatch.setattr(control, "eigh_tridiagonal", diagonalize)
    record = control.wilson_bath(PHYSICAL, 2., 0)
    diagonalize.assert_not_called()
    assert record["xi"] == [0.]
    assert record["weights"] == [1.]
    assert record["metadata"]["chain_hoppings"] == []


@pytest.mark.parametrize("Lambda", [1.8, 2., 3.])
@pytest.mark.parametrize("nmax", [2, 4])
def test_wilson_measure_matches_infinite_logarithmic_star_moments(Lambda, nmax):
    record = control.wilson_bath(PHYSICAL, Lambda, nmax)
    xi, weights = np.array(record["xi"])/PHYSICAL["bandwidth"], np.array(record["weights"])
    # Z,z=1 has nodes +/-a*q^m and weights (1-q)*q^m/2. Summing this
    # geometric series is independent of the chain-hopping implementation.
    q = 1/Lambda
    a = (1-q)/math.log(Lambda)
    for power in range(2*nmax+2):
        expected = 0. if power % 2 else a**power*(1-q)/(1-q**(power+1))
        assert_allclose(weights @ xi**power, expected, atol=2e-15, rtol=2e-13)


@pytest.mark.parametrize("detuning", [-.3, .3])
@pytest.mark.parametrize("field", [-.2, .2])
def test_decoupled_atomic_energy_reference_probabilities_and_up_moment(detuning, field):
    physical = dict(PHYSICAL, gap=1.7, u=1.6, gamma=0., detuning=detuning, field=field)
    record = bath(physical, xi=(.7, -.2), weights=(.2, .3))
    result = control.electron_ed(physical, record)
    singlet, doublet = (result["branches"][key] for key in ("singlet", "doublet"))
    assert_allclose(singlet["energy"], (.8-abs(detuning))/1.7, atol=3e-14)
    assert_allclose(doublet["energy"], field/(2*1.7), atol=3e-14)
    assert_allclose(result["signed_gap"], (field/2-.8+abs(detuning))/1.7, atol=3e-14)
    assert_allclose([singlet[k] for k in ("P0", "P1", "P2", "moment")],
                    [float(detuning > 0), 0., float(detuning < 0), 0.], atol=3e-14)
    assert_allclose([doublet[k] for k in ("P0", "P1", "P2", "moment")], [0., 1., 0., .5], atol=3e-14)
    assert singlet["spin_squared"] == pytest.approx(0., abs=1e-14)
    assert doublet["spin_squared"] == pytest.approx(.75, abs=1e-14)
    assert max(singlet["residual"], doublet["residual"]) < 1e-12


def test_coupled_quadratic_one_level_analytic_branches_without_weight_renormalization():
    physical = dict(PHYSICAL, gap=1.3, bandwidth=3., u=0.)
    # Deliberately not unit weight: finite bath records must not be normalized.
    record = bath(physical, weights=(.37,))
    result = control.electron_ed(physical, record)
    delta = physical["gap"]
    t2 = 2*physical["bandwidth"]*physical["gamma"]*.37/math.pi
    root = math.sqrt(delta**2+4*t2)
    singlet, doublet = (result["branches"][key] for key in ("singlet", "doublet"))
    assert_allclose(singlet["energy"], (delta-root)/delta, atol=2e-14)
    assert_allclose(doublet["energy"], (delta-root)/(2*delta), atol=2e-14)
    assert_allclose(result["signed_gap"], (root-delta)/(2*delta), atol=2e-14)
    singlet_p2 = .25+(delta/(2*root))**2
    doublet_p2 = (1-delta/root)/4
    assert_allclose([singlet[k] for k in ("P0", "P1", "P2", "moment")],
                    [singlet_p2, 1-2*singlet_p2, singlet_p2, 0.], atol=2e-14)
    assert_allclose([doublet[k] for k in ("P0", "P1", "P2", "moment")],
                    [doublet_p2, 1-2*doublet_p2, doublet_p2, (1-2*doublet_p2)/2], atol=2e-14)
    assert record["weights"] == [.37]
    json.dumps(result, allow_nan=False)


def test_degenerate_atomic_singlets_do_not_claim_individual_observables():
    physical = dict(PHYSICAL, u=1., gamma=0.)
    result = control.electron_ed(physical, bath(physical))
    singlet = result["branches"]["singlet"]
    assert singlet["energy"] == pytest.approx(.5)
    assert singlet["degeneracy"] == 2
    assert all(singlet[key] is None for key in ("P0", "P1", "P2", "moment"))
    assert result["branches"]["doublet"]["moment"] == pytest.approx(.5)


@pytest.mark.parametrize("u", [2., 4.])
def test_mixed_spin_atomic_minimum_is_rejected_not_called_singlet(u):
    physical = dict(PHYSICAL, u=u, gamma=0.)
    # For u=2 the empty/double impurity states also coincide with the
    # impurity-plus-bath-QP singlet/triplet. At u=4 only the latter are lowest.
    with pytest.raises(ValueError, match="not purely spin 0.*mixed-spin degeneracy"):
        control.electron_ed(physical, bath(physical))


@pytest.mark.parametrize("nmax", [2, 4])
def test_matched_small_wilson_baths_spin_residual_and_energy_unit_invariance(nmax):
    record = control.wilson_bath(PHYSICAL, 2., nmax)
    result = control.electron_ed(PHYSICAL, record)
    factor = 2.3
    scaled = {k: factor*v if k != "geometry" else v for k, v in PHYSICAL.items()}
    other = control.electron_ed(scaled, control.wilson_bath(scaled, 2., nmax))
    for name, spin in (("singlet", 0.), ("doublet", .75)):
        row = result["branches"][name]
        assert row["degeneracy"] == 1
        assert row["spin_squared"] == pytest.approx(spin, abs=1e-12)
        assert row["residual"] < 1e-10
        assert row["method"] == ("dense" if nmax == 2 else "sparse")
        assert_allclose(sum(row[k] for k in ("P0", "P1", "P2")), 1., atol=2e-14)
        assert_allclose(row["P0"], row["P2"], atol=2e-12)
        assert abs(row["moment"]) <= row["P1"]/2+1e-13
        keys = ("energy", "P0", "P1", "P2", "moment")
        assert_allclose([row[k] for k in keys], [other["branches"][name][k] for k in keys], atol=2e-12)
    assert_allclose(result["signed_gap"], other["signed_gap"], atol=2e-12)


def test_sparse_matches_dense_sector_solution(monkeypatch):
    physical = dict(PHYSICAL, detuning=.23)
    record = control.wilson_bath(physical, 2., 2)
    dense = control.electron_ed(physical, record)
    monkeypatch.setattr(control, "_DENSE_THRESHOLD", 0)
    sparse = control.electron_ed(physical, record)
    for name in ("singlet", "doublet"):
        assert sparse["branches"][name]["method"] == "sparse"
        keys = ("energy", "P0", "P1", "P2", "moment")
        assert_allclose([sparse["branches"][name][k] for k in keys],
                        [dense["branches"][name][k] for k in keys], atol=2e-12)


def test_sparse_degenerate_minimum_checks_complete_dense_subspace(monkeypatch):
    physical = dict(PHYSICAL, u=1., gamma=1e-24)
    monkeypatch.setattr(control, "_DENSE_THRESHOLD", 0)
    result = control.electron_ed(physical, bath(physical))
    singlet = result["branches"]["singlet"]
    assert singlet["method"] == "dense"
    assert singlet["degeneracy"] == 2
    assert singlet["energy"] == pytest.approx(.5, abs=1e-13)
    assert all(singlet[key] is None for key in ("P0", "P1", "P2", "moment"))


def test_bad_eigensolver_residual_is_rejected(monkeypatch):
    original = control.eigh

    def inaccurate(matrix):
        energies, vectors = original(matrix)
        return energies+.01, vectors

    monkeypatch.setattr(control, "eigh", inaccurate)
    with pytest.raises(ValueError, match="residual is too large"):
        control.electron_ed(PHYSICAL, bath(PHYSICAL))


def test_nonzero_local_field_spin_mixing_is_not_mislabeled():
    physical = dict(PHYSICAL, field=.3)
    with pytest.raises(ValueError, match="not purely spin"):
        control.electron_ed(physical, bath(physical))


@pytest.mark.parametrize("Lambda,nmax,match", [
    (1.799, 2, "Lambda < 1.8"), (float("nan"), 2, "finite"), (True, 2, "finite"),
    (2., True, "integer"), (2., 2.5, "integer"), (2., -1, "integer"), (2., 101, "integer"),
])
def test_invalid_wilson_settings(Lambda, nmax, match):
    with pytest.raises(ValueError, match=match):
        control.wilson_bath(PHYSICAL, Lambda, nmax)


@pytest.mark.parametrize("key,value,match", [
    ("gap", 0., "positive"), ("bandwidth", -1., "positive"), ("gap", float("nan"), "finite"),
    ("field", True, "finite"), ("u", -1., "nonnegative"), ("gamma", -.1, "nonnegative"),
    ("geometry", "junction", "single"),
])
def test_invalid_physical_contract(key, value, match):
    physical = dict(PHYSICAL, **{key: value})
    with pytest.raises(ValueError, match=match):
        control.electron_ed(physical, bath(PHYSICAL))
    with pytest.raises(ValueError, match=match):
        control.wilson_bath(physical, 2., 2)


@pytest.mark.parametrize("change,match", [
    ({"xi": [], "weights": []}, "nonempty"), ({"xi": [[0.]]}, "one-dimensional"),
    ({"weights": [1., 2.]}, "equal size"), ({"weights": [0.]}, "strictly positive"),
    ({"weights": [float("inf")]}, "finite"), ({"xi": [float("nan")]}, "finite"),
    ({"delta": 2.}, "delta.*match"), ({"bandwidth": 1.}, "bandwidth.*match"),
    ({"delta": True}, "delta.*match"),
    ({"xi": [0.]*6, "weights": [1/6]*6}, "five signed bath levels"),
])
def test_invalid_bath_contract_and_ed_size_limit(change, match):
    with pytest.raises(ValueError, match=match):
        control.electron_ed(PHYSICAL, bath(PHYSICAL) | change)
