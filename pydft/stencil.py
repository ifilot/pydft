"""Numba-accelerated cubic stencil interpolation kernels."""

import numba as nb

@nb.njit(parallel=True, fastmath=True)
def stencil_interp(vals, x, weights, out):
    """
    Evaluate four-point cubic interpolation stencils for many channels.

    Parameters
    ----------
    vals : ndarray
        Source values with shape ``(radial_points, channels)``.
    x : ndarray
        Base stencil indices for each evaluation point. Each index ``j`` uses
        source rows ``j-1`` through ``j+2``.
    weights : ndarray
        Cubic interpolation weights with shape ``(points, 4)``.
    out : ndarray
        Output array with shape ``(channels, points)``.
    """
    npts = x.size
    channels = vals.shape[1]
    for c in nb.prange(channels):
        for p in range(npts):
            j = x[p]
            out[c, p] = (weights[p,0] * vals[j-1, c] +
                         weights[p,1] * vals[j,   c] +
                         weights[p,2] * vals[j+1, c] +
                         weights[p,3] * vals[j+2, c])
