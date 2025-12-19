# -*- coding: utf-8 -*-

from scipy.special import sph_harm_y
import numpy as np
from scipy.special import lpmv
from math import factorial, pi

class SphericalHarmonicsCache:
    """
    Smart cache for real spherical harmonics Y_lm(theta, phi)
    """

    _cache = {}  # (l, angpts_key) -> (2*l+1, npts)

    @staticmethod
    def _angpts_key(angpts):
        return tuple(tuple(p) for p in angpts)

    @classmethod
    def get_l(cls, l, angpts):
        """
        Return Y_lm for fixed l and all m in [-l, l]
        Shape: (2*l+1, npts)
        """
        key = (l, cls._angpts_key(angpts))

        if key in cls._cache:
            return cls._cache[key]

        npts = len(angpts)
        Yl = np.zeros((2*l + 1, npts))

        for i, m in enumerate(range(-l, l + 1)):
            Yl[i, :] = [
                spherical_harmonic(l, m, theta, phi)
                for theta, phi in angpts
            ]

        cls._cache[key] = Yl
        return Yl

    @classmethod
    def get_ylm(cls, lmax, angpts):
        """
        Return stacked Y_lm for l = 0..lmax
        Shape: ((lmax+1)^2, npts)
        """
        blocks = [cls.get_l(l, angpts) for l in range(lmax + 1)]
        return np.vstack(blocks)

    @classmethod
    def clear_cache(cls):
        cls._cache.clear()

def spherical_harmonic(l, m, theta, phi):
    """
    Calculate value of spherical harmonic function
    
    l:     angular quantum number
    m:     magnetic quantum number
    theta: azimuthal angle in radians
    phi:   polar angle in radians
    """
    if m < 0:
        val = np.sqrt(2) * np.imag(sph_harm_y(l, np.abs(m), phi, theta))
    elif m > 0:
        val = np.sqrt(2) * np.real(sph_harm_y(l, m, phi, theta))
    else:
        val = np.real(sph_harm_y(l, m, phi, theta))
    
    return val

def spherical_harmonic_cart(l, m, p):
    """
    Calculate the value of the spherical harmonic depending on the position
    p in Cartesian coordinates. This function assumes that the position p
    lies on the unit sphere
    
    l:    angular quantum number
    m:    magnetic quantum number
    p:    position three-vector on the unit sphere
    """
    theta = np.arctan2(p[1], p[0])  # azimuthal
    phi = np.arccos(p[2])           # polar

    return spherical_harmonic(l, m, theta, phi)

def real_sph_harm_l_scipy(l, theta, phi):
    """
    Evaluate real Spherical Harmonics for all m-values corresponding to a single
    l-value.

    Parameters
    ----------
    l : int
        Angular momentum quantum number
    theta : array_like
        Azimuthal angles
    phi : array_like
        Polar angles
    """
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    npts = theta.size
    Y = np.zeros((2*l + 1, npts))

    # m = 0
    Y[l, :] = sph_harm_y(l, 0, phi, theta).real

    # m > 0 and m < 0
    for m in range(1, l + 1):
        Yc = sph_harm_y(l, m, phi, theta)
        Y[l + m, :] = np.sqrt(2) * Yc.real
        Y[l - m, :] = np.sqrt(2) * Yc.imag

    return Y

def real_sph_harm_l_legendre(l, theta, phi):
    """
    Real spherical harmonics matching your sph_harm_y-based definition,
    but computed via associated Legendre polynomials.
    """
    theta = np.asarray(theta)
    phi = np.asarray(phi)

    x = np.cos(phi)
    npts = theta.size
    Y = np.empty((2*l + 1, npts))

    for m in range(0, l + 1):
        P = lpmv(m, l, x)

        norm = np.sqrt((2*l + 1) / (4*pi) * factorial(l - m) / factorial(l + m))

        if m == 0:
            Y[l, :] = norm * P
        else:
            Y[l + m, :] = np.sqrt(2) * norm * P * np.cos(m * theta)
            Y[l - m, :] = np.sqrt(2) * norm * P * np.sin(m * theta)

    return Y