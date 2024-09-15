import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    dft = DFT(co, basis='sto3g', verbose=False)
    en = dft.scf(1e-6)
    
    # get the projection onto the spherical harmonics
    lmprojection = dft.get_molgrid_copy().get_rho_lm_atoms()
    
    fig, ax = plt.subplots(1, 2, dpi=144, figsize=(6,4))
    plot_spherical_harmonic_projection(ax[1], fig, lmprojection[0], title='O')
    plot_spherical_harmonic_projection(ax[0], fig, lmprojection[1], title='C')
    plt.tight_layout()
    

def plot_spherical_harmonic_projection(ax, fig, mat, minval=5e-2, title=None):
    im = ax.imshow(np.flip(mat, axis=0), origin='lower', cmap='PiYG', vmin=-minval,vmax=minval)
    ax.grid(linestyle='--', color='black', alpha=0.8)
    ax.set_xticks([-0.5,0.5,3.5,8.5,15.5], labels=['s','p','d','f','g'])
    ax.set_yticks(np.arange(0, mat.shape[0], 5)-0.5, 
                  labels=np.arange(0, mat.shape[0], 5))
    ax.set_ylabel('Radial points $r_{i}$')
    ax.set_xlim(-.5,15)
    ax.set_xlabel('Spherical harmonics (lm)')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')
    
    if title is not None:
        ax.set_title(title)

if __name__ == '__main__':
    main()