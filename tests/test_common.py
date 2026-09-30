"""Independent, bounded algebra, quadrature, and symmetry checks."""

from inspect import Parameter, signature
from itertools import product
import json
from unittest.mock import Mock

from hypothesis import given, settings, strategies as st
import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.integrate import quad
from scipy.optimize import OptimizeResult

from qdjj_solver.common import baths as bath_module
from qdjj_solver.common.algebra import (
    FermionOperator, annihilate, create, hermitian_pair, number,
)
from qdjj_solver.common.baths import (
    DiscreteBath, cosh_grid, discrete_g, fit_surrogate, hybridization_g,
)
from qdjj_solver.common.problem import Hamiltonian, Sector
from electron_oracle import electron_operators, evaluate


def raw_matrix(terms, annihilators):
    """Evaluate the *input* words, without the solver's normal ordering."""
    result = np.zeros(annihilators[0].shape, dtype=complex)
    for word, coefficient in terms:
        term = coefficient * np.eye(len(result), dtype=complex)
        for index in word:
            matrix = annihilators[abs(index) - 1]
            term = term @ (matrix.conj().T if index > 0 else matrix)
        result += term
    return result


def polynomial(terms):
    return sum((FermionOperator({word: coefficient}) for word, coefficient in terms),
               FermionOperator())


def test_all_canonical_anticommutation_relations():
    cs = electron_operators(4)
    identity = np.eye(16)
    for i, j in product(range(4), repeat=2):
        ai, aj = annihilate(i), annihilate(j)
        for left, right, expected in (
            (ai, aj, np.zeros((16, 16))),
            (ai.dagger(), aj.dagger(), np.zeros((16, 16))),
            (ai, aj.dagger(), identity if i == j else np.zeros((16, 16))),
        ):
            anticommutator = left * right + right * left
            assert_allclose(evaluate(anticommutator, cs).toarray(), expected,
                            atol=0., rtol=0.)
        assert_allclose(evaluate(ai * ai.dagger(), cs).toarray(),
                        identity - (cs[i].getH() @ cs[i]).toarray(), atol=0., rtol=0.)
    # Contractions can separate two identical factors before they vanish.
    assert (annihilate(0) * create(0) * annihilate(0)).terms == annihilate(0).terms
    assert (create(0) * annihilate(1) * create(0)).terms == {}


_coefficient = st.tuples(st.integers(-4, 4), st.integers(-4, 4)).map(
    lambda pair: complex(*pair) / 4)
_word = st.lists(st.sampled_from([-3, -2, -1, 1, 2, 3]), max_size=8).map(tuple)
_terms = st.lists(st.tuples(_word, _coefficient), max_size=7)


@settings(max_examples=60, derandomize=True, deadline=None)
@given(_terms, _terms, _coefficient)
def test_complex_polynomials_against_raw_jordan_wigner_words(left, right, scalar):
    cs = electron_operators(3)
    dense_cs = [c.toarray() for c in cs]
    a, b = polynomial(left), polynomial(right)
    am, bm = raw_matrix(left, dense_cs), raw_matrix(right, dense_cs)
    for operator, expected in (
        (a, am), (a + b, am + bm), (a - b, am - bm), (a * b, am @ bm),
        (a.dagger(), am.conj().T), (scalar * a, scalar * am),
        (scalar - a, scalar * np.eye(8) - am),
        ((a + scalar) / (2 - 3j), (am + scalar * np.eye(8)) / (2 - 3j)),
        (hermitian_pair(a), am + am.conj().T),
    ):
        assert_allclose(evaluate(operator, cs).toarray(), expected, atol=2e-13, rtol=2e-13)
        for word in operator.terms:
            creators = [i for i in word if i > 0]
            annihilators = [i for i in word if i < 0]
            assert word == tuple(sorted(creators) + sorted(annihilators))
            assert len(set(word)) == len(word)
    assert a.dagger().dagger().terms == a.terms
    assert FermionOperator.from_records(a.to_records()).terms == a.terms


def test_complex_bogoliubov_substitution_preserves_car_and_polynomials():
    cs = electron_operators(2)
    theta, phase = .37, .81
    u, v = np.cos(theta), np.exp(1j * phase) * np.sin(theta)
    maps = [u * annihilate(0) + v * create(1),
            u * annihilate(1) - v * create(0)]
    matrices = [evaluate(a, cs).toarray() for a in maps]
    for i, j in product(range(2), repeat=2):
        assert_allclose(matrices[i] @ matrices[j].conj().T
                        + matrices[j].conj().T @ matrices[i],
                        np.eye(4) if i == j else np.zeros((4, 4)), atol=3e-16, rtol=0.)
        assert_allclose(matrices[i] @ matrices[j] + matrices[j] @ matrices[i],
                        0., atol=3e-16, rtol=0.)
    terms = [((), .7j), ((-1, 1, -2, 2), .3 - .2j),
             ((1, 2), .21j), ((-2, 1, -1, 2, -2), -.8)]
    result = polynomial(terms).substitute(maps)
    assert_allclose(evaluate(result, cs).toarray(), raw_matrix(terms, matrices),
                    atol=4e-16, rtol=2e-15)
    assert FermionOperator().substitute(maps).terms == {}


@pytest.mark.parametrize("coefficient", [np.nan, complex(0., np.inf)])
def test_operator_rejects_nonfinite_coefficients(coefficient):
    with pytest.raises(ValueError, match="finite"):
        FermionOperator({(): coefficient})


@pytest.mark.parametrize("index", [0, 1.2])
def test_operator_rejects_invalid_signed_indices(index):
    with pytest.raises(ValueError, match="nonzero signed integers"):
        FermionOperator({(index,): 1})


@pytest.mark.parametrize("mode", [-1, 1.5])
def test_mode_validation(mode):
    with pytest.raises(ValueError, match="nonnegative integer"):
        annihilate(mode)


def test_operator_numeric_protocol_cleanup_and_records():
    op = number(np.int64(0))
    assert (2 - op).terms == {(): 2, (1, -1): -1}
    assert (2 + op).terms == {(): 2, (1, -1): 1}
    for expression in (lambda: op + object(), lambda: op * object(), lambda: object() * op):
        with pytest.raises(TypeError):
            expression()
    dirty = FermionOperator({(): 1e-9 + 2j, (1,): 3 + 1e-9j,
                             (-1,): 1e-9 + 1e-9j, (1, -1): 1e-8})
    assert dirty.cleaned(1e-8).terms == {(): 2j, (1,): 3, (1, -1): 1e-8}
    assert dirty.cleaned(0).terms == dirty.terms
    assert number(0).hermiticity_error() == 0.
    assert FermionOperator().hermiticity_error() == 0.
    assert FermionOperator.scalar(2j).hermiticity_error() == 4.
    records = [dict(operators=[-1, 1], real=2., imag=1.),
               dict(operators=[1, -1], real=2., imag=1.),
               dict(operators=[], real=-2., imag=-1.)]
    assert FermionOperator.from_records(records).terms == {}


@pytest.mark.parametrize("kwargs, message", [
    ({"parity": 2}, "parity must"),
    ({"twice_sz": .5}, "integer or None"),
    ({"twice_sz": 0}, "parity and twice_sz"),
    ({"eta": 0}, "eta must"),
    ({"particle_number": -1}, "particle_number must"),
    ({"particle_number": 2}, "particle_number must"),
    ({"twice_sz": -3, "particle_number": 1}, "particle_number and twice_sz"),
])
def test_invalid_sector_contract(kwargs, message):
    with pytest.raises(ValueError, match=message):
        Sector(**kwargs)


@pytest.mark.parametrize("kwargs, message", [
    ({"nimp": 3}, "impurity dimension"),
    ({"spins": (1, 0)}, "spin labels"),
    ({"eta_labels": (1,)}, "eta labels"),
    ({"eta_labels": (1, 0)}, "eta labels"),
    ({"operator": number(2)}, "undeclared fermionic mode"),
    ({"operator": create(0) + annihilate(0)}, "fermion parity"),
    ({"operator": 1j * number(0)}, "not Hermitian"),
])
def test_invalid_hamiltonian_contract(kwargs, message):
    arguments = dict(operator=number(0), nimp=2, spins=(1, -1))
    arguments.update(kwargs)
    with pytest.raises(ValueError, match=message):
        Hamiltonian(**arguments)


def test_sector_missing_eta_and_excess_particle_number():
    h = Hamiltonian(number(0), 1, (1, -1))
    with pytest.raises(ValueError, match="eta reduction is not defined"):
        h.check_sector(Sector(1, 1, 1))
    with pytest.raises(ValueError, match="exceeds the mode count"):
        h.check_sector(Sector(1, 1, particle_number=3))
    # Numpy integers are part of the documented sector contract.
    h.check_sector(Sector(np.int64(1), np.int64(1), particle_number=np.int64(1)))


@pytest.mark.parametrize("word, conserved", [
    pytest.param((1, -1), (True, True, True), id="density"),
    pytest.param((1, 2), (False, True, True), id="pair-changes-number"),
    pytest.param((1, -2), (True, False, True), id="spin-flip"),
    pytest.param((1, -3), (True, True, False), id="hop-changes-eta"),
])
def test_sector_conservation_witnesses_against_dense_charge_blocks(word, conserved):
    """Each nonconserving witness breaks just one of the three optional charges."""
    spins, labels = (1, -1, 1), (1, 1, -1)
    occupations = np.array([[(state >> i) & 1 for i in range(3)] for state in range(8)])
    charges = [occupations.sum(axis=1), occupations @ spins,
               np.prod(np.where(occupations, labels, 1), axis=1)]
    sectors = [Sector(1, None, particle_number=1), Sector(1, 1), Sector(1, None, eta=1)]
    matrix = raw_matrix([(word, .7 + .2j)], [c.toarray() for c in electron_operators(3)])
    matrix += matrix.conj().T
    h = Hamiltonian(hermitian_pair(FermionOperator({word: .7 + .2j})), 1, spins, labels)
    for charge, sector, preserved in zip(charges, sectors, conserved, strict=True):
        assert np.all(matrix[charge[:, None] != charge[None, :]] == 0) == preserved
        if preserved:
            h.check_sector(sector)
        else:
            with pytest.raises(ValueError, match="does not conserve"):
                h.check_sector(sector)
    h.check_sector(Sector(1, None))


def test_hamiltonian_records_coordinates_and_empty_system():
    h = Hamiltonian(hermitian_pair(.2j * create(0) * annihilate(1)), 1, [1, -1], [1, -1],
                    {"occupation": number(0)}, {"tag": "complex"})
    restored = Hamiltonian.from_record(h.record())
    assert restored.record() == h.record()
    assert restored.fingerprint() == h.fingerprint() == h.coordinate_fingerprint()
    assert restored.bath_modes == 1
    coordinates = dict(restriction="none", electron_annihilators=[annihilate(i).to_records()
                                                                  for i in range(2)])
    h.metadata["coordinates"] = coordinates
    changed = Hamiltonian(2 * h.operator, 1, h.spins, h.eta_labels,
                          metadata={"coordinates": coordinates, "tag": "other"})
    assert h.fingerprint() != changed.fingerprint()
    assert h.coordinate_fingerprint() == changed.coordinate_fingerprint()
    assert h.physical_operator(h.operator).terms == h.operator.terms
    with pytest.raises(ValueError, match="undeclared electron mode"):
        h.physical_operator(number(2))
    empty = Hamiltonian.from_record(dict(operator=[], nimp=0, spins=[]))
    assert empty.record()["eta_labels"] is None
    assert empty.observables == empty.metadata == {}
    assert empty.bath_modes == 0
    empty.check_sector(Sector(0, 0, particle_number=0))
    with pytest.raises(ValueError, match="uncompressed declared electron map"):
        empty.physical_operator(FermionOperator.scalar(1))
    empty.metadata["coordinates"] = {"restriction": "none"}
    with pytest.raises(ValueError, match="no declared electron map"):
        empty.physical_operator(FermionOperator.scalar(1))


@pytest.mark.parametrize("xi, weights", [([], []), ([[0.]], [1.]), ([0., 1.], [1.])])
def test_bath_rejects_incompatible_arrays(xi, weights):
    with pytest.raises(ValueError, match="nonempty one-dimensional"):
        DiscreteBath(xi, weights)


@pytest.mark.parametrize("xi, weights", [([np.nan], [1.]), ([0.], [0.])])
def test_bath_rejects_invalid_nodes_or_measure(xi, weights):
    with pytest.raises(ValueError, match="finite and contact weights strictly positive"):
        DiscreteBath(xi, weights)


@pytest.mark.parametrize("field, value", [("delta", 0.), ("bandwidth", np.inf)])
def test_bath_and_grid_reject_invalid_scales(field, value):
    with pytest.raises(ValueError, match="positive and finite"):
        DiscreteBath([0.], [1.], **{field: value})
    with pytest.raises(ValueError, match="positive and finite"):
        cosh_grid(2, **{field: value})


@pytest.mark.parametrize("count", [0, 1.5])
def test_positive_integer_grid_and_fit_sizes(count):
    with pytest.raises(ValueError, match="positive integer"):
        cosh_grid(count)
    with pytest.raises(ValueError, match="positive integer"):
        fit_surrogate(count)


def test_bath_owns_readonly_arrays_and_record_roundtrip():
    xi, weights, metadata = np.array([-.7, .7]), np.array([.3, .3]), {"kind": "manual"}
    bath = DiscreteBath(xi, weights, .4, 2.5, metadata)
    xi[0], weights[0], metadata["kind"] = 99., 8., "changed"
    assert bath.paired and bath.levels == 2
    assert bath.rho == .2
    assert_allclose(bath.energies, np.sqrt(.7**2 + .4**2), atol=0., rtol=2e-16)
    assert bath.metadata == {"kind": "manual"}
    for array in (bath.xi, bath.weights):
        with pytest.raises(ValueError, match="read-only"):
            array[0] = 0.
    assert DiscreteBath.from_record(bath.record()).record() == bath.record()
    assert not DiscreteBath([-.7, .7], [.3, .4]).paired
    assert not DiscreteBath([-.7, .8], [.3, .3]).paired
    assert DiscreteBath([0.], [.5]).paired


def test_cosh_quadrature_moments_and_independent_continuum_integrals():
    delta, bandwidth = .73, 5.4
    bath = cosh_grid(np.int64(32), delta, bandwidth)
    assert bath.levels == 64 and bath.paired
    for power in range(9):
        moment = np.dot(bath.weights, (bath.xi / bandwidth)**power)
        assert_allclose(moment, 0. if power % 2 else 1 / (power + 1), atol=2e-14, rtol=0.)
    omega = np.array([[0., .07, -.7], [1.3, 4., 70.]])
    expected = np.array([quad(lambda xi, w=w: 1 / (np.pi * (w*w + delta*delta + xi*xi)),
                              -bandwidth, bandwidth, epsabs=1e-13, epsrel=1e-13)[0]
                         for w in omega.ravel()]).reshape(omega.shape)
    assert_allclose(discrete_g(bath, omega), expected, atol=3e-15, rtol=2e-14)
    assert_allclose(hybridization_g(omega, delta, bandwidth), expected, atol=3e-15, rtol=2e-14)
    assert_allclose(discrete_g(bath, 0.), expected[0, 0], atol=3e-15, rtol=2e-14)
    coarse = cosh_grid(1, delta, bandwidth)
    upper = np.arcsinh(bandwidth / delta)
    measure = upper * delta * np.cosh(upper / 2) / bandwidth
    assert_allclose(coarse.weights.sum(), measure, atol=0., rtol=2e-16)
    assert_allclose(coarse.metadata["measure_error"], measure - 1, atol=0., rtol=2e-16)
    assert abs(measure - 1) > .05


@pytest.mark.parametrize("kwargs, message", [
    ({"starts": 0}, "mesh or number of starts"),
    ({"frequency_points": 9}, "mesh or number of starts"),
    ({"frequency_min": 2., "frequency_cutoff": 1.}, "mesh or number of starts"),
    ({"delta": 0.}, "energy scales"),
    ({"bandwidth": np.inf}, "energy scales"),
    ({"relative_weight": np.nan}, "energy scales"),
])
def test_surrogate_fitting_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        fit_surrogate(2, **kwargs)


@pytest.mark.parametrize("max_nfev", [
    0, -1, np.int64(0), 1.5, 1., np.nan, np.inf, "5000", None,
    True, False, np.bool_(True), np.bool_(False),
])
def test_surrogate_rejects_invalid_evaluation_budget(monkeypatch, max_nfev):
    optimizer = Mock(side_effect=AssertionError("invalid budget reached the optimizer"))
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    with pytest.raises(ValueError, match="max_nfev must be a positive integer"):
        fit_surrogate(2, starts=1, frequency_points=10, max_nfev=max_nfev)
    optimizer.assert_not_called()


def test_surrogate_budget_preserves_signature_and_old_coefficients(monkeypatch):
    parameters = signature(fit_surrogate).parameters
    old_names = ["levels", "delta", "bandwidth", "frequency_cutoff", "frequency_min",
                 "frequency_points", "starts", "seed", "relative_weight"]
    assert list(parameters) == [*old_names, "max_nfev"]
    assert all(parameters[name].kind is Parameter.POSITIONAL_OR_KEYWORD for name in old_names)
    assert [parameters[name].default for name in old_names[1:]] == [
        1., 100., 100., 1e-3, 1000, 4, 1729, 0.]
    assert parameters["max_nfev"].kind is Parameter.KEYWORD_ONLY
    assert parameters["max_nfev"].default == 5000
    optimizer = Mock(wraps=bath_module.least_squares)
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    arguments = (2, .8, 8., 12., .003, 70, 2, 901, .35)
    with pytest.raises(TypeError, match="positional"):
        fit_surrogate(*arguments, 5000)
    optimizer.assert_not_called()
    default = fit_surrogate(*arguments)
    explicit = fit_surrogate(*arguments, max_nfev=5000)
    assert [call.kwargs["max_nfev"] for call in optimizer.call_args_list] == [5000]*4
    # Pre-budget API coefficients; allow small SciPy/platform termination differences.
    assert_allclose(default.xi, [-1.1533760078876094, 1.1533760078876094], rtol=2e-7, atol=0.)
    assert_allclose(default.weights, [.22350570425768274]*2, rtol=2e-7, atol=0.)
    assert_allclose(default.metadata["cost"], .09864683634958052, rtol=2e-12, atol=0.)
    assert_allclose(default.xi, explicit.xi, rtol=0., atol=0.)
    assert_allclose(default.weights, explicit.weights, rtol=0., atol=0.)
    for bath in (default, explicit):
        assert bath.metadata["max_nfev"] == 5000
        assert bath.metadata["relative_weight"] == .35


@pytest.mark.parametrize("levels, relative_weight", [(1, 0.), (2, .35), (5, 1.)])
def test_odd_even_surrogate_fits_and_analytic_jacobian(monkeypatch, levels, relative_weight):
    optimizer = bath_module.least_squares
    fits, budgets = [], []

    def checked_optimizer(fun, initial, *, jac, **kwargs):
        # Finite differences independently check the logarithmic-residue and
        # logarithmic-energy derivatives, including the unpaired zero level.
        step = 2e-6
        finite_difference = np.column_stack([
            (fun(initial + step * direction) - fun(initial - step * direction)) / (2 * step)
            for direction in np.eye(len(initial))])
        assert_allclose(jac(initial), finite_difference, atol=2e-9, rtol=2e-7)
        result = optimizer(fun, initial, jac=jac, **kwargs)
        fits.append(result)
        budgets.append(kwargs["max_nfev"])
        return result

    monkeypatch.setattr(bath_module, "least_squares", checked_optimizer)
    delta, bandwidth = .8, 8.
    bath = fit_surrogate(np.int64(levels), delta=delta, bandwidth=bandwidth,
                         frequency_min=.003, frequency_cutoff=12., frequency_points=70,
                         starts=2, seed=901, relative_weight=relative_weight,
                         max_nfev=np.int64(200))
    assert budgets == [200, 200]
    assert all(type(budget) is int for budget in budgets)
    assert bath.metadata["max_nfev"] == 200
    trials = bath.metadata["optimizer_trials"]
    assert len(trials) == len(fits) == 2
    for trial, fit in zip(trials, fits, strict=True):
        assert trial == dict(success=bool(fit.success), status=int(fit.status), message=fit.message,
                             evaluations=int(fit.nfev), optimality=float(fit.optimality),
                             cost=float(fit.fun @ fit.fun))
        assert trial["success"] and 0 < trial["evaluations"] <= 200
        assert_allclose(trial["cost"], 2*fit.cost, atol=0., rtol=2e-15)
    selected = bath.metadata["selected_start"]
    assert selected == min(range(len(trials)), key=lambda i: trials[i]["cost"])
    assert bath.metadata["cost"] == trials[selected]["cost"]
    assert bath.metadata["evaluations"] == trials[selected]["evaluations"]
    assert bath.metadata["fit_optimality"] == trials[selected]["optimality"]
    assert bath.levels == levels and bath.paired
    assert np.count_nonzero(bath.xi == 0) == levels % 2
    assert np.all(bath.weights > 0) and np.all(np.diff(bath.xi) > 0)
    omega = np.geomspace(.003, 12., 70)
    target = hybridization_g(omega, delta, bandwidth)
    scaling = target**relative_weight * target.max()**(1 - relative_weight)
    residual = (discrete_g(bath, omega) - target) / scaling
    assert_allclose(bath.metadata["cost"], residual @ residual, atol=2e-15, rtol=2e-13)
    assert_allclose(bath.metadata["max_relative_fit_error"],
                    np.max(abs(discrete_g(bath, omega) / target - 1)), atol=0., rtol=2e-15)
    assert_allclose(bath.metadata["measure_error"], bath.weights.sum() - 1, atol=0., rtol=0.)
    assert bath.metadata["relative_weight"] == relative_weight
    if levels == 1:
        column = 1 / ((omega**2 + delta**2) * scaling)
        residue = np.dot(column, target / scaling) / np.dot(column, column)
        assert_allclose(bath.weights / (np.pi * bath.rho), [residue], atol=2e-12, rtol=2e-12)
    else:
        # An external, denser frequency mesh tests interpolation, not just cost.
        validation = np.geomspace(.003, 12., 137)
        error = abs(discrete_g(bath, validation) - hybridization_g(validation, delta, bandwidth))
        assert np.max(error) < (.1 if levels == 2 else .035)


@pytest.mark.parametrize("extra, selected", [(2.**-26, 1), (0., 0)])
def test_surrogate_selects_recorded_cost_when_norms_tie(monkeypatch, extra, selected):
    residuals = [np.r_[1., extra, np.zeros(8)], np.r_[1., np.zeros(9)]]
    # Distinct squared costs can round to the same norm after the square root.
    assert [np.linalg.norm(residual) for residual in residuals] == [1., 1.]
    fits = [OptimizeResult(success=True, status=1, message="converged",
                           x=np.log([.5, i+1.]), fun=residual, nfev=1, optimality=0.,
                           cost=.5*np.dot(residual, residual))
            for i, residual in enumerate(residuals)]
    monkeypatch.setattr(bath_module, "least_squares", Mock(side_effect=fits))
    bath = fit_surrogate(2, starts=2, frequency_points=10)
    trials = bath.metadata["optimizer_trials"]
    assert trials[0]["cost"] == 1.+extra**2
    assert trials[1]["cost"] == 1.
    assert bath.metadata["selected_start"] == selected
    assert bath.metadata["cost"] == 1.
    assert_allclose(bath.xi, [-selected-1., selected+1.], atol=0., rtol=0.)


@pytest.mark.parametrize("success, field, value", [
    (False, "x", 0.), (True, "x", np.nan), (True, "x", np.inf),
    (True, "fun", np.nan), (True, "fun", np.inf),
])
def test_surrogate_reports_optimizer_failure(monkeypatch, success, field, value):
    fits = [OptimizeResult(success=success, x=np.zeros(2), fun=np.zeros(10),
                           status=i+1 if success else -i, message=f"termination of trial {i}",
                           optimality=0., nfev=7-i, cost=0.) for i in range(2)]
    for fit in fits:
        fit[field][0] = value
        fit.cost = .5*(fit.fun @ fit.fun)
    optimizer = Mock(side_effect=fits)
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    with pytest.raises(RuntimeError, match="did not converge") as error:
        fit_surrogate(2, starts=2, frequency_points=10, max_nfev=7)
    assert [call.kwargs["max_nfev"] for call in optimizer.call_args_list] == [7, 7]
    assert "max_nfev=7" in str(error.value)
    for i, fit in enumerate(fits):
        assert f"start {i}: status={fit.status}, nfev={fit.nfev}, {fit.message}" in str(error.value)


def test_surrogate_one_evaluation_budget_is_enforced():
    with pytest.raises(RuntimeError, match="did not converge") as error:
        fit_surrogate(2, starts=2, frequency_points=10, max_nfev=1)
    assert "max_nfev=1" in str(error.value)
    for i in range(2):
        assert f"start {i}: status=0, nfev=1" in str(error.value)
    assert str(error.value).count("The maximum number of function evaluations is exceeded.") == 2


@pytest.mark.parametrize("field, value", [
    ("x", np.nan), ("x", np.inf), ("fun", np.nan), ("fun", np.inf),
])
def test_surrogate_skips_nonfinite_start_and_records_json_safe_trials(monkeypatch, field, value):
    invalid = OptimizeResult(success=True, status=1, message="`gtol` termination condition is satisfied.",
                             x=np.zeros(2), fun=np.zeros(10), nfev=1, optimality=value, cost=0.)
    invalid[field][0] = value
    invalid.cost = .5*(invalid.fun @ invalid.fun)
    valid = OptimizeResult(success=True, status=2, message="`ftol` termination condition is satisfied.",
                           x=np.log([.5, 1.]), fun=np.full(10, .2), nfev=2, optimality=.001,
                           cost=.5*np.dot(np.full(10, .2), np.full(10, .2)))
    optimizer = Mock(side_effect=[invalid, valid])
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    bath = fit_surrogate(2, starts=2, frequency_points=10, max_nfev=3)
    record = json.loads(json.dumps(bath.record(), allow_nan=False))
    assert DiscreteBath.from_record(record).record() == bath.record()
    assert record["metadata"]["selected_start"] == 1
    assert record["metadata"]["optimizer_trials"] == [
        dict(success=True, status=1, message=invalid.message, evaluations=1,
             optimality=None, cost=None if field == "fun" else 0.),
        dict(success=True, status=2, message=valid.message, evaluations=2,
             optimality=.001, cost=float(valid.fun @ valid.fun)),
    ]
    assert_allclose(bath.xi, [-1., 1.], atol=0., rtol=0.)
    assert bath.metadata["cost"] == 2*valid.cost


@pytest.mark.filterwarnings("error::RuntimeWarning")
@pytest.mark.parametrize("success", [False, True])
def test_surrogate_skips_overflowed_cost_and_records_json_safe_trials(monkeypatch, success):
    overflowed = OptimizeResult(success=success, status=2 if success else 0,
                                message="overflowed trial termination", x=np.zeros(2),
                                fun=np.full(10, 1e200), nfev=7, optimality=1e200, cost=np.inf)
    valid = OptimizeResult(success=True, status=1, message="`gtol` termination condition is satisfied.",
                           x=np.log([.5, 2.]), fun=np.ones(10), nfev=2, optimality=0., cost=5.)
    assert np.all(np.isfinite(overflowed.fun))
    optimizer = Mock(side_effect=[overflowed, valid])
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    with np.errstate(over="warn", invalid="warn"):
        bath = fit_surrogate(2, starts=2, frequency_points=10, max_nfev=7)
    assert optimizer.call_count == 2
    record = json.loads(json.dumps(bath.record(), allow_nan=False))
    assert DiscreteBath.from_record(record).record() == bath.record()
    metadata = record["metadata"]
    assert metadata["selected_start"] == 1
    assert metadata["optimizer_trials"] == [
        dict(success=success, status=overflowed.status, message=overflowed.message,
             evaluations=7, optimality=1e200, cost=None),
        dict(success=True, status=1, message=valid.message, evaluations=2, optimality=0., cost=10.),
    ]
    assert metadata["cost"] == 2*valid.cost
    assert metadata["evaluations"] == valid.nfev
    assert_allclose(bath.xi, [-2., 2.], atol=0., rtol=0.)


@pytest.mark.filterwarnings("error::RuntimeWarning")
@pytest.mark.parametrize("success", [False, True])
def test_surrogate_reports_all_overflowed_costs_as_failure(monkeypatch, success):
    fits = [OptimizeResult(success=success, status=i+1 if success else -i,
                           message=f"overflowed trial {i} termination", x=np.zeros(2),
                           fun=np.full(10, 1e200), nfev=7-i, optimality=1e200, cost=np.inf)
            for i in range(2)]
    assert all(np.all(np.isfinite(fit.fun)) for fit in fits)
    optimizer = Mock(side_effect=fits)
    monkeypatch.setattr(bath_module, "least_squares", optimizer)
    with np.errstate(over="warn", invalid="warn"), pytest.raises(RuntimeError, match="did not converge") as error:
        fit_surrogate(2, starts=2, frequency_points=10, max_nfev=7)
    assert optimizer.call_count == 2
    assert "max_nfev=7" in str(error.value)
    for i, fit in enumerate(fits):
        assert f"start {i}: status={fit.status}, nfev={fit.nfev}, {fit.message}" in str(error.value)


def test_surrogate_ignores_failed_start_and_keeps_converged_fit(monkeypatch):
    optimizer = bath_module.least_squares
    failed_start = False

    def first_start_fails(*args, **kwargs):
        nonlocal failed_start
        result = optimizer(*args, **kwargs)
        if not failed_start:
            result.success = False
            result.status = 0
            result.message = "The maximum number of function evaluations is exceeded."
            result.nfev = kwargs["max_nfev"]
            result.fun[:] = 0.  # A failed run must not win even with a smaller cost.
            result.cost = 0.
            failed_start = True
        return result

    monkeypatch.setattr(bath_module, "least_squares", first_start_fails)
    bath = fit_surrogate(2, starts=2, frequency_points=30, bandwidth=5., frequency_cutoff=5.,
                         max_nfev=100)
    assert bath.metadata["selected_start"] == 1
    failed, selected = bath.metadata["optimizer_trials"]
    assert failed["success"] is False and failed["status"] == 0
    assert failed["evaluations"] == 100 and failed["cost"] == 0.
    assert selected["success"] is True
    assert bath.metadata["cost"] == selected["cost"]
    assert bath.metadata["evaluations"] == selected["evaluations"]
    assert bath.metadata["fit_optimality"] == selected["optimality"]
    assert bath.metadata["cost"] > 0.
    omega = np.geomspace(1e-3, 5., 30)
    target = hybridization_g(omega, bandwidth=5.)
    residual = (discrete_g(bath, omega) - target) / target.max()
    assert_allclose(bath.metadata["cost"], residual @ residual, atol=2e-15, rtol=2e-13)
