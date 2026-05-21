import numpy as np
from pydft import MoleculeBuilder, DFT
from pytessel import PyTessel

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    dft = DFT(co, basis='sto3g')
    res = dft.scf(1e-4, verbose=False)
        
    # grab a copy of the MolecularGrid object and of the coefficient matrix
    molgrid = dft.get_molgrid_copy()
    orbc = res['orbc']

    sz = 10     # size of the domain
    npts = 50   # number of sampling points per cartesian direction

    build_isosurface('1pi_pos.ply', molgrid, orbc[:,4], 0.03, sz, npts)
    build_isosurface('1pi_neg.ply', molgrid, orbc[:,4], -0.03, sz, npts)

def build_isosurface(filename, molgrid, coeff, isovalue, sz, npts):
    # produce meshgrid for the unit cell
    x = np.linspace(-sz/2, sz/2, npts)
    zz, yy, xx = np.meshgrid(x, x, x, indexing='ij')
    N = len(x)
    gridpoints = np.zeros((N, N, N, 3))
    gridpoints[:,:,:,0] = xx
    gridpoints[:,:,:,1] = yy
    gridpoints[:,:,:,2] = zz
    gridpoints = gridpoints.reshape((-1,3))

    # build isosurface
    scalarfield = molgrid.get_amplitude_at_points(gridpoints, coeff).reshape(npts, npts, npts)
    unitcell = np.diag(np.ones(3) * sz)

    pytessel = PyTessel()
    vertices, normals, indices = pytessel.marching_cubes(scalarfield.flatten(), scalarfield.shape, unitcell.flatten(), isovalue)
    pytessel.write_ply(filename, vertices, normals, indices)

if __name__ == '__main__':
    main()
