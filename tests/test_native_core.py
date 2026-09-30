"""Direct native contracts, checked against unbounded-integer Fock arithmetic."""

import gc
import itertools
import subprocess
import sys
import textwrap

from hypothesis import example, given, settings, strategies as st
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from scipy.sparse import csc_matrix

from qdjj_solver.qp_solver import _core


pytestmark = pytest.mark.native
NO_SZ = 1_000_000_000


def native_basis(nimp=2, spins=(1, -1, 1, -1), cutoff=2, parity=0,
                 sz=NO_SZ, eta=(), target_eta=0, maximum=100_000, number=-1):
    return _core.Basis(nimp, list(spins), cutoff, parity, sz, list(eta), target_eta, maximum, number)


def integer_states(basis):
    return [sum(int(word) << (64 * w) for w, word in enumerate(row))
            for row in basis.occupations()]


def occupation_oracle(nimp, spins, cutoff, parity, sz, eta, target_eta, number):
    """Enumerate bit strings, independently of native spin-group combinations."""
    result = {}
    for bits in itertools.product((0, 1), repeat=len(spins)):
        particles = sum(bits)
        q = sum(bits[nimp:])
        if q > cutoff or particles % 2 != parity:
            continue
        if sz != NO_SZ and sum(s * n for s, n in zip(spins, bits, strict=True)) != sz:
            continue
        if target_eta and np.prod([e for e, n in zip(eta, bits, strict=True) if n]) != target_eta:
            continue
        if number >= 0 and particles != number:
            continue
        result[sum(n << i for i, n in enumerate(bits))] = q
    return result


def fermion_action(state, string):
    """Right-to-left action using Python integers, including arbitrarily high bits."""
    amplitude = 1
    for op in reversed(string):
        bit = 1 << (abs(op) - 1)
        if bool(state & bit) == (op > 0):
            return None, 0
        amplitude *= -1 if (state & (bit - 1)).bit_count() % 2 else 1
        state ^= bit
    return state, amplitude


def operator_oracle(source, destination, coefficients, strings):
    rows = {state: i for i, state in enumerate(destination)}
    matrix = np.zeros((len(destination), len(source)), dtype=complex)
    for j, state in enumerate(source):
        for coefficient, string in zip(coefficients, strings, strict=True):
            target, sign = fermion_action(state, string)
            if target in rows:
                matrix[rows[target], j] += coefficient * sign
    return matrix


@st.composite
def basis_cases(draw):
    n = draw(st.integers(1, 7))
    spins = draw(st.lists(st.sampled_from((-1, 1)), min_size=n, max_size=n))
    nimp = draw(st.integers(0, n))
    cutoff = draw(st.integers(0, n - nimp))
    parity = draw(st.integers(0, 1))
    sz = draw(st.one_of(st.just(NO_SZ), st.integers(-n, n)))
    eta = draw(st.lists(st.sampled_from((-1, 1)), min_size=n, max_size=n))
    target_eta = draw(st.sampled_from((0, -1, 1)))
    number = draw(st.integers(-1, n))
    return nimp, spins, cutoff, parity, sz, eta, target_eta, number


@given(basis_cases())
@example(case=(0, (1, -1), 0, 0, NO_SZ, (), 0, 0))  # Vacuum, no impurity.
@example(case=(4, (1, -1, 1, -1), 0, 0, 0, (-1, 1, -1, 1), 1, 4))  # Filled impurity.
@example(case=(0, (1, 1, -1, -1), 2, 0, 2, (-1, 1, 1, 1), -1, 2))  # Fixed N, Sz, eta.
@settings(max_examples=160, deadline=None, derandomize=True)
def test_native_basis_independent_occupation_enumeration(case):
    nimp, spins, cutoff, parity, sz, eta, target_eta, number = case
    expected = occupation_oracle(*case)
    if not expected:
        with pytest.raises(ValueError, match="empty symmetry sector"):
            native_basis(nimp, spins, cutoff, parity, sz, eta, target_eta, number=number)
        return
    basis = native_basis(nimp, spins, cutoff, parity, sz, eta, target_eta, number=number)
    states = integer_states(basis)
    assert basis.dimension == len(expected) == len(set(states))
    assert dict(zip(states, basis.qp_counts(), strict=True)) == expected
    assert basis.occupations().dtype == np.uint64
    assert np.all(np.diff(basis.qp_counts()) >= 0)
    assert basis.bytes >= basis.occupations().nbytes + basis.qp_counts().nbytes


@pytest.mark.parametrize("kwargs, message", [
    ({"nimp": -1}, "dimensions"), ({"nimp": 5}, "dimensions"),
    ({"nimp": 0, "spins": (), "cutoff": 0}, "dimensions"),
    ({"cutoff": -1}, "dimensions"), ({"cutoff": 3}, "dimensions"),
    ({"parity": -1}, "dimensions"), ({"parity": 2}, "dimensions"),
    ({"maximum": 0}, "dimensions"), ({"number": -2}, "dimensions"),
    ({"number": 5}, "dimensions"), ({"spins": (1, 0, 1, -1)}, "spin labels"),
    ({"target_eta": 2}, "eta must"), ({"target_eta": -2}, "eta must"),
    ({"target_eta": 1}, "missing eta"),
    ({"target_eta": -1, "eta": (1,)}, "missing eta"),
    ({"eta": (1, 0, 1, -1)}, "eta labels"),
    ({"sz": 12}, "empty symmetry"),
])
def test_native_basis_validation(kwargs, message):
    with pytest.raises(ValueError, match=message):
        native_basis(**kwargs)


def test_native_dimension_limit_and_legacy_constructor():
    with pytest.raises(ValueError, match="max_dimension"):
        native_basis(maximum=1)
    legacy = _core.Basis(2, [1, -1, 1, -1], 2, 0, NO_SZ, [], 0, 8)
    assert legacy.dimension == 8
    assert_array_equal(legacy.occupations(), native_basis(maximum=8).occupations())


STRINGS = [(), (1, -1), (2, -2), (1, 2, -2, -1), (1, 3, -3, -1),
           (1, 3), (3, 4), (-4, -3), (-3, -1), (4, -1), (1, -4),
           (1, 2, -4, -3), (1, 3, -4, -2), (1,), (-2,), (1, 3, -2), (1, -4, -2)]


@pytest.mark.parametrize("complex_coefficients", [False, True])
@pytest.mark.parametrize("cutoff", [0, 1, 2])
@pytest.mark.parametrize("parity", [0, 1])
def test_native_action_sparse_diagonal_and_cross_basis(complex_coefficients, cutoff, parity):
    source = native_basis(cutoff=cutoff, parity=parity)
    destination = native_basis(cutoff=2, parity=1 - parity)
    coefficients = np.arange(1, len(STRINGS) + 1) / 7
    if complex_coefficients:
        coefficients = coefficients + 1j * coefficients[::-1]
    op = _core.Operator(source, coefficients, STRINGS)
    states = integer_states(source)
    expected = operator_oracle(states, states, coefficients, STRINGS)
    assert op.real is (not complex_coefficients)
    assert_allclose(op.diagonal(), np.diag(expected), rtol=0, atol=1e-14)
    # A strided input forces the binding's conversion path; a zero entry tests skipping.
    vector = np.arange(2 * source.dimension, dtype=float)[::2] / 3
    assert_allclose(op.matvec_complex(vector), expected @ vector, rtol=0, atol=2e-14)
    assert_allclose(op.matvec_complex(vector + 1j * vector[::-1]),
                    expected @ (vector + 1j * vector[::-1]), rtol=0, atol=5e-14)
    if complex_coefficients:
        with pytest.raises(ValueError, match="real action"):
            op.matvec_real(vector)
    else:
        assert_allclose(op.matvec_real(vector), expected @ vector, rtol=0, atol=2e-14)
    values, rows, pointers = op.sparse(10_000)
    matrix = csc_matrix((values, rows, pointers), shape=expected.shape)
    assert matrix.has_canonical_format
    assert_allclose(matrix.toarray(), expected, rtol=0, atol=1e-14)
    assert pointers[0] == 0 and pointers[-1] == len(values) == len(rows)
    cross = operator_oracle(states, integer_states(destination), coefficients, STRINGS)
    assert_allclose(op.between(vector, destination), cross @ vector, rtol=0, atol=2e-14)
    assert_allclose(op.between(vector + 1j * vector[::-1], source),
                    expected @ (vector + 1j * vector[::-1]), rtol=0, atol=5e-14)
    # Equal parity with a smaller cutoff exercises successful and missing diagonal lookups.
    reduced = native_basis(cutoff=0, parity=parity)
    rectangular = operator_oracle(states, integer_states(reduced), coefficients, STRINGS)
    assert_allclose(op.between(vector, reduced), rectangular @ vector, rtol=0, atol=2e-14)


@pytest.mark.parametrize("nmodes", [63, 64, 65, 127, 128, 129])
@pytest.mark.parametrize("complex_coefficients", [False, True])
def test_native_multiword_fermion_signs(nmodes, complex_coefficients):
    active = sorted({i for i in (0, 1, 31, 62, 63, 64, 95, 126, 127, 128, nmodes - 1) if i < nmodes})
    spins = [1 if i in active else -1 for i in range(nmodes)]
    basis = native_basis(nimp=2, spins=spins, cutoff=2, sz=2, number=2)
    states = integer_states(basis)
    assert set(states) == {sum(1 << i for i in pair) for pair in itertools.combinations(active, 2)}
    assert basis.occupations().shape == (len(states), (nmodes + 63) // 64)
    strings = [(), (1, -1), (1, nmodes, -nmodes, -1)]
    strings += [(i + 1, -(j + 1)) for i in active for j in active]
    strings += [tuple(i + 1 for i in active[:2]) + tuple(-(j + 1) for j in active[-2:][::-1]),
                (1, nmodes), (-nmodes, -1)]
    coefficients = np.sin(np.arange(len(strings)) + 1)
    if complex_coefficients:
        coefficients = coefficients + 1j * np.cos(np.arange(len(strings)))
    op = _core.Operator(basis, coefficients, strings)
    expected = operator_oracle(states, states, coefficients, strings)
    rng = np.random.default_rng(nmodes)
    vector = rng.normal(size=len(states))
    assert_allclose(op.matvec_complex(vector + 1j * vector[::-1]),
                    expected @ (vector + 1j * vector[::-1]), rtol=0, atol=3e-14)
    if not complex_coefficients:
        assert_allclose(op.matvec_real(vector), expected @ vector, rtol=0, atol=2e-14)
    assert_allclose(csc_matrix(op.sparse(100_000), shape=expected.shape).toarray(), expected, rtol=0, atol=1e-14)
    assert_allclose(op.diagonal(), np.diag(expected), rtol=0, atol=1e-14)
    destination = native_basis(nimp=2, spins=spins, cutoff=1, parity=1, sz=1, number=1)
    odd_strings = [(-(i + 1),) for i in active]
    odd_strings += [(nmodes, -(active[1] + 1), -1), (1,)]
    odd_coefficients = [(-1j if complex_coefficients else -1) ** i for i in range(len(odd_strings))]
    odd = _core.Operator(basis, odd_coefficients, odd_strings)
    cross = operator_oracle(states, integer_states(destination), odd_coefficients, odd_strings)
    assert_allclose(odd.between(vector, destination), cross @ vector, rtol=0, atol=1e-14)


@given(st.lists(st.tuples(st.sets(st.integers(1, 4)), st.sets(st.integers(1, 4)),
                          st.integers(-3, 3), st.integers(-3, 3)), max_size=16),
       st.integers(0, 2), st.integers(0, 1))
@settings(max_examples=90, deadline=None, derandomize=True)
def test_native_random_canonical_polynomials(terms, cutoff, parity):
    basis = native_basis(cutoff=cutoff, parity=parity)
    coefficients = [complex(real, imag) for _, _, real, imag in terms]
    strings = [tuple(sorted(cre)) + tuple(-i for i in sorted(ann, reverse=True))
               for cre, ann, _, _ in terms]
    expected = operator_oracle(integer_states(basis), integer_states(basis), coefficients, strings)
    op = _core.Operator(basis, coefficients, strings)
    actual = np.column_stack([op.matvec_complex(v) for v in np.eye(basis.dimension)])
    assert_allclose(actual, expected, rtol=0, atol=1e-14)
    assert_allclose(csc_matrix(op.sparse(1000), shape=expected.shape).toarray(), expected, rtol=0, atol=1e-14)


def test_native_sparse_duplicate_coalescing_cancellation_and_limits():
    basis = native_basis()
    strings = [(), (), (3, -1), (3, -1), (4, -2), (4, -2)]
    coefficients = [2, -2, 3, -3, 2, 5]
    op = _core.Operator(basis, coefficients, strings)
    expected = operator_oracle(integer_states(basis), integer_states(basis), coefficients, strings)
    nnz = np.count_nonzero(expected)
    values, rows, pointers = op.sparse(nnz)
    assert len(values) == nnz and np.all(values != 0)
    assert_allclose(csc_matrix((values, rows, pointers), shape=expected.shape).toarray(), expected, rtol=0, atol=0)
    with pytest.raises(ValueError, match="max_nnz"):
        op.sparse(nnz - 1)
    for empty in (_core.Operator(basis, [], []), _core.Operator(basis, [1, -1], [(3, -1), (3, -1)])):
        values, rows, pointers = empty.sparse(0)
        assert len(values) == len(rows) == 0
        assert_array_equal(pointers, np.zeros(basis.dimension + 1))
        assert_array_equal(empty.matvec_real(np.ones(basis.dimension)), np.zeros(basis.dimension))


@pytest.mark.parametrize("coefficients, strings, message", [
    ([1], [], "count mismatch"), ([], [(1,)], "count mismatch"),
    ([np.nan], [()], "non-finite"), ([np.inf], [()], "non-finite"),
    ([complex(0, np.nan)], [()], "non-finite"), ([complex(0, np.inf)], [()], "non-finite"),
    ([1], [(0,)], "mode outside"), ([1], [(5,)], "mode outside"),
    ([1], [(-5,)], "mode outside"), ([1], [(2**31 - 1,)], "mode outside"),
    ([1], [(-1, 1)], "canonical"), ([1], [(2, 1)], "canonical"),
    ([1], [(-1, -2)], "canonical"), ([1], [(1, 1)], "canonical"),
    ([1], [(-1, -1)], "canonical"), ([1], [(1, -2, 3)], "canonical"),
    ([1], [(1, 2, -1, -2)], "canonical"),
])
def test_native_operator_validation(coefficients, strings, message):
    with pytest.raises(ValueError, match=message):
        _core.Operator(native_basis(), coefficients, strings)


@pytest.mark.parametrize("method", ["matvec_real", "matvec_complex", "between"])
@pytest.mark.parametrize("shape", [(), (1,), (8, 1), (1, 8), (9,)])
def test_native_vector_dimension_errors(method, shape):
    basis = native_basis()
    op = _core.Operator(basis, [], [])
    args = (np.ones(shape), basis) if method == "between" else (np.ones(shape),)
    with pytest.raises(ValueError, match="state-vector dimension"):
        getattr(op, method)(*args)


@pytest.mark.parametrize("destination", [dict(nimp=1), dict(nimp=1, spins=(1, -1, 1), cutoff=2)])
def test_native_transition_mode_mismatch(destination):
    source = native_basis()
    op = _core.Operator(source, [1], [()])
    with pytest.raises(ValueError, match="different fermionic modes"):
        op.between(np.ones(source.dimension), native_basis(**destination))


def test_native_arrays_and_operator_own_their_storage():
    basis = native_basis()
    states, counts = basis.occupations(), basis.qp_counts()
    saved_states, saved_counts = states.copy(), counts.copy()
    states[:] = 0
    counts[:] = -9
    assert_array_equal(basis.occupations(), saved_states)
    assert_array_equal(basis.qp_counts(), saved_counts)
    coefficients, strings = [2., 3.], [[], [1, -1]]
    op = _core.Operator(basis, coefficients, strings)
    expected = operator_oracle(integer_states(basis), integer_states(basis), coefficients, strings)
    coefficients[:] = [100., 200.]
    strings[1][:] = [4, -4]
    diagonal = op.diagonal()
    diagonal[:] = 99
    assert_allclose(op.diagonal(), np.diag(expected), rtol=0, atol=0)
    values, rows, pointers = op.sparse(1000)
    values[:] = 0
    rows[:] = 0
    pointers[:] = 0
    del basis
    gc.collect()
    vector = np.arange(len(expected), dtype=float)
    result = op.matvec_real(vector)
    assert_allclose(result, expected @ vector, rtol=0, atol=0)
    result[:] = -1
    assert_allclose(op.matvec_real(vector), expected @ vector, rtol=0, atol=0)
    retained = op.diagonal()
    del op
    gc.collect()
    assert_allclose(retained, np.diag(expected), rtol=0, atol=0)


@pytest.mark.parametrize("expression, message", [
    ("core.Operator(None, [], [])", "non-null basis"),
    ("op.between(np.ones(basis.dimension), None)", "non-null destination"),
    ("core.Operator(basis, [1], [[-(1 << 31)]])", "mode outside"),
    ("core.Operator(basis, [1], [[1, -(1 << 31)]])", "mode outside"),
])
def test_native_unsafe_input_regressions_in_subprocess(expression, message):
    # Load the exact parent binary: isolated builds and sanitizer jobs must not
    # accidentally exercise a different extension from the editable install.
    script = textwrap.dedent("""
        import importlib.util
        import sys
        import numpy as np
        spec = importlib.util.spec_from_file_location('_core', sys.argv[1])
        core = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(core)
        basis = core.Basis(2, [1, -1, 1, -1], 2, 0, 1000000000, [], 0, 100)
        op = core.Operator(basis, [], [])
        try:
            eval(sys.argv[2])
        except ValueError as error:
            assert sys.argv[3] in str(error), str(error)
        else:
            raise AssertionError('unsafe input was accepted')
    """)
    completed = subprocess.run([sys.executable, "-c", script, _core.__file__, expression, message],
                               text=True, capture_output=True, timeout=30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
