"""Cooperative signal handling and exception-safe native cancellation."""

import signal

import numpy as np
import pytest

from qdjj_solver.qp_solver import _core


pytestmark = [pytest.mark.native,
              pytest.mark.skipif(not hasattr(signal, "setitimer"), reason="POSIX interval timer")]
NO_SZ = 1_000_000_000


@pytest.mark.parametrize("operation", ["basis", "operator", "real", "complex", "between", "sparse"])
def test_native_signal_unwinds_and_objects_remain_usable(operation):
    spins = [1, -1] * 14
    args = (2, spins, 6, 0, NO_SZ, [], 0, 2_000_000)
    if operation == "basis":
        invoke = lambda: _core.Basis(*args)
    else:
        basis = _core.Basis(*args)
        strings = [[i + 1, -(j + 1)] for i in range(28) for j in range(28) if i != j]
        coefficients = [0.01] * len(strings)
        if operation == "operator":
            invoke = lambda: _core.Operator(basis, coefficients, strings)
        else:
            operator = _core.Operator(basis, coefficients, strings)
            vector = np.ones(basis.dimension)
            invoke = {"real": lambda: operator.matvec_real(vector),
                      "complex": lambda: operator.matvec_complex(vector),
                      "between": lambda: operator.between(vector, basis),
                      "sparse": lambda: operator.sparse(100_000_000)}[operation]

    def interrupted(signum, frame):
        raise KeyboardInterrupt("native cancellation regression")

    previous = signal.signal(signal.SIGALRM, interrupted)
    try:
        with pytest.raises(KeyboardInterrupt, match="native cancellation regression"):
            signal.setitimer(signal.ITIMER_REAL, 0.001)
            invoke()
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)

    # A thrown signal must not leave the GIL or extension in an unusable state.
    small = _core.Basis(2, [1, -1], 0, 0, NO_SZ, [], 0, 10)
    identity = _core.Operator(small, [1.0], [[]])
    np.testing.assert_array_equal(identity.matvec_real(np.ones(small.dimension)), 1.0)
