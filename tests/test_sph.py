import unittest
import numpy as np
from scipy.special import sph_harm_y

from pydft import spherical_harmonic
from pydft.spherical_harmonics import (
    SphericalHarmonicsCache,
    real_sph_harm_l_legendre,
    real_sph_harm_l_scipy,
    spherical_harmonic_cart,
)

class TestSphericalHarmonics(unittest.TestCase):      

    def test_spherical_harmonics(self):
        # test Ylm(0,0) (s)
        np.testing.assert_almost_equal(spherical_harmonic(0, 0, 0, 0), 1 / np.sqrt(4 * np.pi))
        vals = np.random.rand(10) * 2.0 * np.pi
        np.testing.assert_almost_equal(spherical_harmonic(0, 0, vals, 0), np.ones_like(vals) / np.sqrt(4 * np.pi))

        # test Ylm(1,0) (p_z)
        theta = np.random.rand(10) * 2.0 * np.pi # azimuthal
        np.testing.assert_almost_equal(spherical_harmonic(1, 0, vals, 0), np.ones_like(theta) * 1/2 * np.sqrt(3 / np.pi))
        phi = np.random.rand(10) * np.pi # polar
        np.testing.assert_almost_equal(np.real(sph_harm_y(1, 0, phi, 0)), 1/2 * np.sqrt(3 / np.pi) * np.cos(phi))
        np.testing.assert_almost_equal(spherical_harmonic(1, 0, 0, phi), 1/2 * np.sqrt(3 / np.pi) * np.cos(phi))

        # test Ylm(1,0) (p_z)
        theta = np.random.rand(10) * 2.0 * np.pi # azimuthal
        np.testing.assert_almost_equal(spherical_harmonic(1, 0, vals, 0), np.ones_like(theta) * 1/2 * np.sqrt(3 / np.pi))
        phi = np.random.rand(10) * np.pi # polar
        np.testing.assert_almost_equal(np.real(sph_harm_y(1, 0, phi, 0)), 1/2 * np.sqrt(3 / np.pi) * np.cos(phi))
        np.testing.assert_almost_equal(spherical_harmonic(1, 0, 0, phi), 1/2 * np.sqrt(3 / np.pi) * np.cos(phi))

    def test_spherical_harmonic_m_branches_and_cartesian_wrapper(self):
        theta = np.array([0.2, 1.1, 2.0])
        phi = np.array([0.4, 1.0, 2.4])

        np.testing.assert_allclose(
            spherical_harmonic(2, 1, theta, phi),
            np.sqrt(2) * np.real(sph_harm_y(2, 1, phi, theta)),
        )
        np.testing.assert_allclose(
            spherical_harmonic(2, -1, theta, phi),
            np.sqrt(2) * np.imag(sph_harm_y(2, 1, phi, theta)),
        )
        np.testing.assert_allclose(
            spherical_harmonic_cart(1, 0, np.array([0.0, 0.0, 1.0])),
            spherical_harmonic(1, 0, 0.0, 0.0),
        )

    def test_scipy_and_legendre_spherical_harmonics_match(self):
        theta = np.array([0.2, 1.1, 2.0])
        phi = np.array([0.4, 1.0, 2.4])

        np.testing.assert_allclose(
            real_sph_harm_l_scipy(3, theta, phi),
            real_sph_harm_l_legendre(3, theta, phi),
            atol=1e-12,
        )

    def test_spherical_harmonic_cache_lifecycle(self):
        SphericalHarmonicsCache.clear_cache()
        try:
            SphericalHarmonicsCache.precache_bulk([(1, 50)], max_workers=1)
            cached = SphericalHarmonicsCache.get_l(1, 50)

            SphericalHarmonicsCache.precache_bulk([(1, 50)], max_workers=1)
            self.assertIs(SphericalHarmonicsCache.get_l(1, 50), cached)

            SphericalHarmonicsCache.freeze()
            self.assertIs(SphericalHarmonicsCache.get_l(1, 50), cached)
            with self.assertRaisesRegex(RuntimeError, "miss after freeze"):
                SphericalHarmonicsCache.get_l(2, 50)
        finally:
            SphericalHarmonicsCache.clear_cache()

if __name__ == '__main__':
    unittest.main()
