"""Independent Padé moments, published coefficients, and physical bath equivalence."""

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import eigh_tridiagonal

from qdjj_solver import chain_expansion
from qdjj_solver.common.baths import discrete_g
from qdjj_solver.common.io import bath_from_config


@pytest.mark.parametrize("levels", [1, 2, 3, 4, 7, 8, 16])
@pytest.mark.parametrize("bandwidth", [.02, 10., 100.])
def test_finite_band_pade_moments(levels, bandwidth):
    bath = chain_expansion(levels, delta=.7, bandwidth=bandwidth)
    angle = np.arctan(bandwidth/.7)
    moment = angle
    for k in range(levels):
        if k:
            moment = (np.sin(angle)*np.cos(angle)**(2*k-1)+(2*k-1)*moment)/(2*k)
        finite = np.sum(bath.weights/(np.pi*bath.rho*.7)/(1+(bath.xi/.7)**2)**(k+1))
        # L coefficients in the expansion in x=(omega/Delta)^2 match the
        # independent analytic integrals (2/pi) integral cos(theta)^(2k)dtheta.
        assert_allclose(finite, 2*moment/np.pi, rtol=8e-13, atol=2e-15)
    assert bath.paired and np.all(bath.weights > 0) and bath.levels == levels
    assert (bath.xi == 0).sum() == levels % 2
    assert max(abs(bath.xi)) <= bandwidth*(1+1e-12)


def test_published_chain_coefficients_and_narrow_band_limit():
    for d in (.01, 10., 100.):
        g0 = 2/np.pi*np.arctan(d)
        factor = 2*(1+d*d)*np.arctan(d)/(d+(1+d*d)*np.arctan(d))
        assert_allclose(chain_expansion(1, bandwidth=d).metadata["chain_h"], [g0], atol=1e-14)
        assert_allclose(chain_expansion(2, bandwidth=d).metadata["chain_h"],
                        [g0*factor, factor-1], atol=2e-14)
    # Independently printed Appendix A, Eq. (46), rounded to six digits.
    assert_allclose(chain_expansion(3, bandwidth=10.).metadata["chain_h"],
                    [2.46516, 2.17024, .329662], atol=5e-6, rtol=0)


@pytest.mark.parametrize("levels", [1, 2, 3, 4, 8])
@pytest.mark.parametrize("wide_band", [False, True])
def test_star_and_chain_resolvents(levels, wide_band):
    bath = chain_expansion(levels, delta=1.3, bandwidth=70., wide_band=wide_band)
    h = np.array(bath.metadata["chain_h"])
    if wide_band:
        ell = np.arange(1, levels)
        assert_allclose(h, np.r_[levels, (levels**2-ell**2)/(4*ell**2-1)], atol=0)
    matrix = np.diag(bath.delta*np.sqrt(h[1:]), 1)+np.diag(bath.delta*np.sqrt(h[1:]), -1)
    for omega in (0., .2, 1., 10.):
        target = bath.delta*h[0]*np.linalg.inv(matrix@matrix+(bath.delta**2+omega**2)*np.eye(levels))[0, 0]
        assert_allclose(discrete_g(bath, omega), target, rtol=5e-13)


def test_wide_band_normalization_cancels_and_config_roundtrip():
    first = chain_expansion(4, bandwidth=10., wide_band=True)
    second = chain_expansion(4, bandwidth=100., wide_band=True)
    assert_allclose(first.weights/first.rho, second.weights/second.rho, atol=1e-14)
    assert abs(first.weights.sum()-1) > .1
    omega = np.array([0., .1, .5, 1., 2.])
    assert_allclose(discrete_g(first, omega), (1+omega**2/2)/(1+omega**2+omega**4/8), atol=3e-15)
    configured = bath_from_config(dict(kind="chain-expansion", levels=4, bandwidth=10., wide_band=True))
    assert configured.record() == first.record()
    assert bath_from_config(dict(kind="discrete", **first.record())).record() == first.record()
    eigs = eigh_tridiagonal(np.zeros(4), np.sqrt([5., .8, .2]), eigvals_only=True)
    assert_allclose(first.xi, eigs, atol=5e-15)


@pytest.mark.parametrize("kwargs", [dict(levels=0), dict(levels=True), dict(levels=2.5),
                                  dict(levels=2, delta=0), dict(levels=2, bandwidth=np.inf),
                                  dict(levels=2, wide_band="yes")])
def test_chain_expansion_invalid_inputs(kwargs):
    with pytest.raises(ValueError):
        chain_expansion(**kwargs)
