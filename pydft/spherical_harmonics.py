# -*- coding: utf-8 -*-
"""Real spherical harmonics used for density and Hartree-potential expansions."""

import numpy as np
from scipy.special import sph_harm_y, lpmv
from math import factorial, pi
from concurrent.futures import ThreadPoolExecutor

from .angulargrid import AngularGrid

class SphericalHarmonicsCache:
    """
    Cache for real spherical harmonics Y_lm evaluated on Lebedev angular grids.

    Cached objects:
        (l, nangpts) -> ndarray of shape (2*l + 1, nangpts)

    Assumptions:
    - For a given nangpts, the Lebedev angular grid is deterministic.
    - Angular points are generated internally and never part of the cache key.
    """

    _cache = {}          # (l, nangpts) -> Y_l
    _theta_phi = {}      # nangpts -> (theta, phi)
    _frozen = False

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @classmethod
    def _get_theta_phi(cls, nangpts):
        """
        Lazily generate (theta, phi) arrays for a Lebedev grid of size nangpts.
        """
        if nangpts not in cls._theta_phi:
            ag = AngularGrid()
            coeffs = ag.get_coefficients(nangpts)

            theta = coeffs[:,0]
            phi = coeffs[:,1]

            cls._theta_phi[nangpts] = (theta, phi)

        return cls._theta_phi[nangpts]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def get_l(cls, l, nangpts):
        """
        Return Y_lm for fixed l and all m in [-l, l].

        Shape: (2*l + 1, nangpts)
        """
        key = (l, nangpts)

        if key in cls._cache:
            return cls._cache[key]

        if cls._frozen:
            raise RuntimeError(
                f"SphericalHarmonicsCache miss after freeze: {key}"
            )

        theta, phi = cls._get_theta_phi(nangpts)
        Yl = real_sph_harm_l_legendre(l, theta, phi)

        cls._cache[key] = Yl
        return Yl

    @classmethod
    def get_ylm(cls, lmax, nangpts):
        """
        Return all real spherical harmonics up to ``lmax`` on a Lebedev grid.

        The rows are ordered by increasing ``l`` and, within each ``l``, by
        ``m=-l, ..., l``. The resulting array has shape
        ``((lmax + 1)**2, nangpts)``.
        """
        return np.vstack([
            cls.get_l(l, nangpts) for l in range(lmax + 1)
        ])

    # ------------------------------------------------------------------
    # Bulk pre-caching (parallel)
    # ------------------------------------------------------------------

    @classmethod
    def precache_bulk(cls, requests, max_workers=None):
        """
        Precompute spherical harmonics for multiple (lmax, nangpts) requests.

        Parameters
        ----------
        requests : iterable of (lmax, nangpts)
        max_workers : int or None
            Passed to ThreadPoolExecutor
        """
        jobs = []
        keys = []

        for lmax, nangpts in requests:
            theta, phi = cls._get_theta_phi(nangpts)

            for l in range(lmax + 1):
                key = (l, nangpts)
                if key not in cls._cache:
                    jobs.append((l, theta, phi))
                    keys.append(key)

        if not jobs:
            return

        # Parallel numeric work (no cache mutation)
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            results = list(ex.map(job_compute_l, jobs))

        # Serial cache insertion
        for key, Yl in zip(keys, results):
            cls._cache[key] = Yl

    # ------------------------------------------------------------------
    # Cache lifecycle control
    # ------------------------------------------------------------------

    @classmethod
    def freeze(cls):
        """
        Prevent new cache entries from being created.
        """
        cls._frozen = True

    @classmethod
    def clear_cache(cls):
        """
        Clear all cached spherical harmonics and angular grids.
        """
        cls._cache.clear()
        cls._theta_phi.clear()
        cls._frozen = False

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
    l-value. This function is faster than looping over `spherical_harmonic`, yet
    remains slower than `real_sph_harm_l_legendre`, which is the preferred
    function for this task.

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
    but computed via associated Legendre polynomials. This function seems
    to perform best.

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

def job_compute_l(args):
    """
    Helper function for the ThreadPoolExecutor
    """
    l, theta, phi = args
    return real_sph_harm_l_legendre(l, theta, phi)
