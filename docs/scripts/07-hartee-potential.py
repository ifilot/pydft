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

if __name__ == '__main__':
    main()