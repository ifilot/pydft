import numba as nb

@nb.njit(parallel=True, fastmath=True)
def stencil_interp(ulm_rev, i, w, out):
    Neval = i.size
    Nchan = ulm_rev.shape[1]
    for c in nb.prange(Nchan):
        for p in range(Neval):
            j = i[p]
            out[c, p] = (w[p,0]*ulm_rev[j-1, c] +
                        w[p,1]*ulm_rev[j,   c] +
                        w[p,2]*ulm_rev[j+1, c] +
                        w[p,3]*ulm_rev[j+2, c])