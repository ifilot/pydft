# -*- coding: utf-8 -*-

from scipy.special import sph_harm
import numpy as np

def spherical_harmonic(l, m, theta, phi):
    """
    Calculate value of spherical harmonic function
    
    l:     angular quantum number
    m:     magnetic quantum number
    theta: azimuthal angle in radians
    phi:   polar angle in radians
    """
    if m < 0:
        val = np.sqrt(2) * np.imag(sph_harm(np.abs(m), l, theta, phi))
    elif m > 0:
        val = np.sqrt(2) * np.real(sph_harm(m, l, theta, phi))
    else:
        val = np.real(sph_harm(m, l, theta, phi))
    
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
    theta = np.arctan2(p[1], p[0])
    phi = np.arccos(p[2])

    return spherical_harmonic(l, m, theta, phi)