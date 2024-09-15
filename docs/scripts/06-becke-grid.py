from pydft import MolecularGrid
from pyqint import MoleculeBuilder
import numpy as np
import matplotlib.pyplot as plt

# construct molecule
mol = MoleculeBuilder().from_name('benzene')
cgfs, atoms = mol.build_basis('sto3g')

# construct molecular grid
molgrid = MolecularGrid(atoms, cgfs)

# produce grid of sampling points to calculate the atomic
# weight coefficients for
N = 150
sz = 8
x = np.linspace(-sz,sz,N)
xv,yv = np.meshgrid(x,x)
points = np.array([[x,y,0] for x,y in zip(xv.flatten(),yv.flatten())])

# calculate the atomic weights
mweights = molgrid.calculate_weights_at_points(points, k=3)

# plot the atomic weights
plt.figure(dpi=144)
plt.imshow(np.max(mweights,axis=0).reshape((N,N)),
           extent=(-sz,sz,-sz,sz), interpolation='bicubic')
plt.xlabel('x [a.u.]')
plt.ylabel('y [a.u.]')
plt.colorbar()
plt.grid(linestyle='--', color='black', alpha=0.5)

# add the atoms to the plot
r = np.zeros((len(atoms), 3))
for i,at in enumerate(atoms):
    r[i] = at[0]
plt.scatter(r[0:6,0], r[0:6,1], s=50.0, color='grey', edgecolor='black')
plt.scatter(r[6:12,0], r[6:12,1], s=50.0, color='white', edgecolor='black')

plt.tight_layout()