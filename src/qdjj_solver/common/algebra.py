"""Exact fermionic operator algebra shared by both numerical backends.

In a monomial, ``+(i+1)`` denotes creation and ``-(i+1)`` annihilation.
Operators are written in their mathematical order (rightmost acts first).
Canonical strings have ascending creators followed by descending annihilators.
"""

from __future__ import annotations

from functools import lru_cache
from numbers import Number

import numpy as np


@lru_cache(maxsize=65536)
def _normal_order(ops: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], int], ...]:
    def rank(op):
        return (0, op) if op > 0 else (1, op)

    for i in range(len(ops) - 1):
        a, b = ops[i:i + 2]
        if a == b:
            return ()
        if rank(a) > rank(b):
            result = {}
            swapped = ops[:i] + (b, a) + ops[i + 2:]
            for key, value in _normal_order(swapped):
                result[key] = result.get(key, 0) - value
            if a < 0 and a == -b:
                for key, value in _normal_order(ops[:i] + ops[i + 2:]):
                    result[key] = result.get(key, 0) + value
            return tuple((key, value) for key, value in result.items() if value)
    return ((ops, 1),)


class FermionOperator:
    """Finite sum of fermionic monomials; scalar multiplication uses ``*``.

    Multiplication normal-orders the *full* operator, not matrices of singly
    projected creation operators. The distinction matters at the QP boundary.
    """

    def __init__(self, terms=None):
        accumulated = {}
        for ops, coefficient in (terms or {}).items():
            coefficient = complex(coefficient)
            if not np.isfinite(coefficient):
                raise ValueError("operator coefficients must be finite")
            ops = tuple(ops)
            if any(not isinstance(i, (int, np.integer)) or i == 0 for i in ops):
                raise ValueError("operator indices must be nonzero signed integers")
            for key, sign in _normal_order(ops):
                accumulated[key] = accumulated.get(key, 0) + sign * coefficient
        self.terms = {key: value for key, value in accumulated.items() if value != 0}

    @classmethod
    def scalar(cls, value):
        return cls({(): value})

    def __add__(self, other):
        if isinstance(other, Number):
            other = self.scalar(other)
        if not isinstance(other, FermionOperator):
            return NotImplemented
        result = self.terms.copy()
        for key, value in other.terms.items():
            result[key] = result.get(key, 0) + value
        return FermionOperator(result)

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, Number):
            return FermionOperator({key: value * other for key, value in self.terms.items()})
        if not isinstance(other, FermionOperator):
            return NotImplemented
        result = {}
        for a, ca in self.terms.items():
            for b, cb in other.terms.items():
                for key, sign in _normal_order(a + b):
                    result[key] = result.get(key, 0) + sign * ca * cb
        return FermionOperator(result)

    def __rmul__(self, other):
        return self * other if isinstance(other, Number) else NotImplemented

    def __truediv__(self, other):
        return self * (1 / other)

    def dagger(self):
        result = {}
        for key, value in self.terms.items():
            adjoint = tuple(-i for i in key[::-1])
            for ordered, sign in _normal_order(adjoint):
                result[ordered] = result.get(ordered, 0) + sign*value.conjugate()
        return FermionOperator(result)

    def substitute(self, annihilators):
        """Substitute canonical annihilation operators (e.g. a Bogoliubov map)."""
        result = FermionOperator()
        creators = [a.dagger() for a in annihilators]
        for key, coefficient in self.terms.items():
            term = self.scalar(coefficient)
            for op in key:
                term = term * (creators[op - 1] if op > 0 else annihilators[-op - 1])
            result += term
        return result

    def cleaned(self, tolerance=1e-14):
        """Remove only floating-point cancellation at the specified absolute scale."""
        result = {}
        for key, value in self.terms.items():
            value = complex(0 if abs(value.real) < tolerance else value.real,
                            0 if abs(value.imag) < tolerance else value.imag)
            if value:
                result[key] = value
        return FermionOperator(result)

    def hermiticity_error(self):
        return max((abs(v) for v in (self - self.dagger()).terms.values()), default=0.0)

    def to_records(self):
        return [{"operators": list(key), "real": value.real, "imag": value.imag}
                for key, value in sorted(self.terms.items())]

    @classmethod
    def from_records(cls, records):
        result = cls()
        for record in records:
            result += cls({tuple(record["operators"]): complex(record["real"], record["imag"])})
        return result


def annihilate(mode: int) -> FermionOperator:
    if not isinstance(mode, (int, np.integer)) or mode < 0:
        raise ValueError("mode must be a nonnegative integer")
    return FermionOperator({(-(int(mode) + 1),): 1})


def create(mode: int) -> FermionOperator:
    return annihilate(mode).dagger()


def number(mode: int) -> FermionOperator:
    return create(mode) * annihilate(mode)


def hermitian_pair(operator: FermionOperator) -> FermionOperator:
    return operator + operator.dagger()
