import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    dft = DFT(co, basis='sto3g', verbose=False)
    en = dft.scf(1e-6)
    
    # grab molecular orbital energies and coefficients
    orbc = dft.get_data()['C']
    orbe = dft.get_data()['orbe']  
    
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
    
    # grab a copy of the MolecularGrid object
    molgrid = dft.get_molgrid_copy()
    
    fig, ax = plt.subplots(2, 5, dpi=144, figsize=(16,6))
    for i in range(len(orbe)):
        m_ax = ax[i//5, i%5]
        
        # calculate field
        field = molgrid.get_amplitude_at_points(gridpoints, orbc[:,i]).reshape((npts, npts))
        
        # plot field
        levels = np.linspace(-2,2,17, endpoint=True)
        im = m_ax.contourf(x, x, field, levels=levels, cmap='PiYG')
        m_ax.contour(x, x, field, colors='black')
        m_ax.set_aspect('equal', 'box')
        m_ax.set_xlabel('x-coordinate [a.u.]')
        m_ax.set_ylabel('z-coordinate [a.u.]')
        divider = make_axes_locatable(m_ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')
        m_ax.set_title('MO %i: %6.4f Ht' % (i+1, orbe[i]))

    plt.tight_layout()

if __name__ == '__main__':
    main()