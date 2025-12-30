import numba as nb

@nb.njit(parallel=True, fastmath=True)
def stencil_interp(vals, x, weights, out):
    npts = x.size
    channels = vals.shape[1]
    for c in nb.prange(channels):
        for p in range(npts):
            j = x[p]
            out[c, p] = (weights[p,0] * vals[j-1, c] +
                         weights[p,1] * vals[j,   c] +
                         weights[p,2] * vals[j+1, c] +
                         weights[p,3] * vals[j+2, c])