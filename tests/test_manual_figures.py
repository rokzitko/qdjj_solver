"""Independent physical checks of the analytical introductory figure curves."""

import importlib
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.docs


@pytest.fixture
def primer(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/"docs"/"manual"/"scripts"))
    return importlib.import_module("fig_primer")


@pytest.mark.parametrize("pairing", [0., .6, .3+.4j])
def test_atomic_curves_against_even_block(primer, pairing):
    # Direct matrix diagonalization also checks the complex pair-amplitude convention.
    for detuning in (-1.2, 0., .8):
        h = np.array([[1.-detuning, -np.conj(pairing)], [-pairing, 1.+detuning]])
        values, vectors = np.linalg.eigh(h)
        lower, upper, charge = primer.atomic(detuning, 2., pairing)
        np.testing.assert_allclose([lower, upper], values, atol=2e-15)
        if abs(values[1]-values[0]) > 1e-12:
            np.testing.assert_allclose(charge, 2.*abs(vectors[1, 0])**2, atol=2e-15)
        else:
            assert np.isnan(charge)  # no unique charge in a degenerate even subspace
    lower, _, _ = primer.atomic(np.sqrt(1.-abs(pairing)**2), 2., pairing)
    assert abs(lower) < 1e-14  # crossing with the odd state's zero energy


@pytest.mark.parametrize("transmission", [0., .2, .95, 1.])
def test_contact_current_is_energy_derivative(primer, transmission):
    # Current units are e Delta/hbar, so I = -2 dE_A/dphi.
    phi = np.array([.2, .8, 2.2, 4.8])
    step = 1e-5
    energy, current = primer.contact(phi, transmission)
    plus, _ = primer.contact(phi+step, transmission)
    minus, _ = primer.contact(phi-step, transmission)
    np.testing.assert_allclose(current, -(plus-minus)/step, atol=5e-10, rtol=1e-8)
    assert np.all((energy >= 0.) & (energy <= 1.))
    edge, _ = primer.contact(0., transmission)
    assert edge == 1.
    minimum, _ = primer.contact(np.pi, transmission)
    np.testing.assert_allclose(minimum, np.sqrt(1.-transmission), atol=2e-15)


def test_thermal_degeneracy_and_classical_ysr_record(primer):
    record = primer.mc.load("primer")
    thermal = record["thermal"]
    x = np.asarray(thermal["x"])
    # Explicit Boltzmann weights for three states: one singlet, two doublet components.
    weights = np.array([np.ones_like(x), np.exp(-x), np.exp(-x)])
    expected = weights[1:].sum(axis=0)/weights.sum(axis=0)
    np.testing.assert_allclose(thermal["doublet_probability"], expected, atol=2e-15)
    ysr = record["ysr"]
    alpha = np.asarray(ysr["alpha"])
    assert np.asarray(ysr["signed_energy"])[alpha == 1.].item() == 0.
    # Independent scattering phase-shift identity E/Delta = cos(2 atan(alpha)).
    np.testing.assert_allclose(ysr["signed_energy"], np.cos(2.*np.arctan(alpha)), atol=5e-16)
