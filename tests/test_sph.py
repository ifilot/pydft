import unittest
import sys
import os
import numpy as np
from scipy.special import sph_harm_y

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import spherical_harmonic, spherical_harmonic_cart

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

if __name__ == '__main__':
    unittest.main()
