import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    dft = DFT(co, basis='sto3g', verbose=False)
    en = dft.scf(1e-6)
    
    # build list of basis functions
    labels = []
    for a in co.atoms:
        for o in ['1s', '2s', '2px', '2py', '2pz']:
            labels.append('%s - %s' % (a[0],o))
    
    fig, ax = plt.subplots(1, 2, dpi=144, figsize=(8,4))
    plot_matrix(ax[0], dft.get_data()['S'], xlabels=labels, ylabels=labels, title='Overlap matrix')
    plot_matrix(ax[1], dft.get_data()['F'], xlabels=labels, ylabels=labels, title='Hamiltonian matrix')

def plot_matrix(ax, mat, xlabels=None, ylabels=None, title = None, xlabelrot=90):
    """
    Produce plot of matrix
    """
    ax.imshow(mat, vmin=-1, vmax=1, cmap='PiYG')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(i, j, '%.2f' % mat[j,i], ha='center', va='center',
                    fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.hlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    ax.vlines(np.arange(1, mat.shape[0])-0.5, -0.5, mat.shape[0] - 0.5,
              color='black', linestyle='--', linewidth=1)
    
    ax.set_xticks(np.arange(0, mat.shape[0]))
    if xlabels is not None:
        ax.set_xticklabels(xlabels, rotation=xlabelrot)
    ax.set_yticks(np.arange(0, mat.shape[0]))
    
    if ylabels is not None:
        ax.set_yticklabels(ylabels, rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=7)
    
    if title is not None:
        ax.set_title(title)

if __name__ == '__main__':
    main()