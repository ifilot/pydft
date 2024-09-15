import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os,sys

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, DFT

mol_builder = MoleculeBuilder()
mol = mol_builder.from_name('co')

# construct dft object
dft = DFT(mol, basis='sto3g')
energy = dft.scf()
C = dft.get_data()['C']

molgrid = dft.get_molgrid_copy()
she1 = molgrid.get_spherical_harmonic_expansion_of_amplitude(C[:,4], radial_factor=True)
she2 = molgrid.get_spherical_harmonic_expansion_of_amplitude(C[:,5], radial_factor=True)

fig, ax = plt.subplots(1,4,dpi=144, figsize=(10,4))

for i in range(0,4):
    she = she1 if i < 2 else she2
    limit = max(abs(np.min(she[i % 2])), abs(np.max(she[i % 2])) )
    im = ax[i].imshow(she[i % 2,:,0:16], origin='lower', vmin=-limit, vmax=limit, cmap='BrBG')
    
    title = '$\psi_{4}$' if i < 2 else '$\psi_{5}$'
    title += '/ C' if i % 2 == 0 else '/ O'
    ax[i].set_title(title)
    
    ax[i].set_ylabel(r'$r_{i}^{2} \rho_{lm} (r_{i})$')
    ax[i].set_xlabel('$Y_{lm}$')
    
    ax[i].set_xticks([-0.5,0.5,3.5,8.5])
    ax[i].set_xticklabels(['s','p','d','f'])
    ax[i].set_xticks(np.arange(0,16), minor=True)
    ax[i].grid(linestyle='--', color='black', alpha=0.5)
    
    # add colorbars
    divider = make_axes_locatable(ax[i])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

plt.tight_layout()