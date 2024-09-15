import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g', verbose=False)
en = dft.scf(1e-6)

# generate grid of points and calculate the electron density for these points
sz = 4      # size of the domain
npts = 100  # number of sampling points per cartesian direction

# produce meshgrid for the xz-plane
x = np.linspace(-sz/2,sz/2,npts)
zz, xx = np.meshgrid(x, x, indexing='ij')
gridpoints = np.zeros((zz.shape[0], xx.shape[1], 3))
gridpoints[:,:,0] = xx
gridpoints[:,:,2] = zz
gridpoints[:,:,1] = np.zeros_like(xx) # set y-values to 0
gridpoints = gridpoints.reshape((-1,3))

# calculate (logarithmic) scalar field and convert if back to an 2D array
gradient = np.linalg.norm(dft.get_gradient_at_points(gridpoints), axis=1)
gradient = np.log10(gradient.reshape((npts, npts)))

# build contour plot
fig, ax = plt.subplots(1,1, dpi=144, figsize=(4,4))
im = ax.contourf(x, x, gradient, levels=np.linspace(-2,4,7, endpoint=True))
ax.contour(x, x, gradient, colors='black', levels=np.linspace(-2,4,7, endpoint=True))
ax.set_aspect('equal', 'box')
ax.set_xlabel('x-coordinate [a.u.]')
ax.set_ylabel('z-coordinate [a.u.]')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax, orientation='vertical')
ax.set_title('Electron density gradient magnitude')