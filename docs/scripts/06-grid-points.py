from pydft import MolecularGrid
from pyqint import MoleculeBuilder
import matplotlib.pyplot as plt

# construct molecule
mol = MoleculeBuilder().from_name('co')
cgfs, atoms = mol.build_basis('sto3g')

# construct molecular grid
molgrid = MolecularGrid(atoms, cgfs, nshells=16, nangpts=74)
molgrid.initialize()

# obtain the grid points
gridpoints = molgrid.get_grid_coordinates()

# plot the atomic weights
colors = '#DD0000', '#222222'
fig = plt.figure(dpi=300, figsize=(8,6))
ax = fig.add_subplot(projection='3d')
for i in range(0, len(gridpoints)):
    ax.scatter(gridpoints[i][:,0], gridpoints[i][:,1], gridpoints[i][:,2],
               s = 1.5, alpha=0.5, color=colors[i])

# set axes
ax.set_xlim(-5,5)
ax.set_ylim(-5,5)
ax.set_zlim(-5,5)
ax.set_xlabel('x [a.u.]')
ax.set_ylabel('y [a.u.]')
ax.set_zlabel('z [a.u.]')
ax.set_box_aspect(aspect=None, zoom=0.8)