# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os, sys

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, MolecularGrid

def main():
    # construct molecule
    mol = MoleculeBuilder().from_name('benzene')
    cgfs, atoms = mol.build_basis('sto3g')
    
    # construct molecular grid
    molgrid = MolecularGrid(atoms, cgfs)
    
    # produce grid of sampling points to calculate the atomic
    # weight coefficients for
    N = 100
    sz = 8
    x = np.linspace(-sz,sz,N)
    xv,yv = np.meshgrid(x,x)
    points = np.array([[x,y,0] for x,y in zip(xv.flatten(),yv.flatten())])
    
    # calculate the atomic weights
    mweights = molgrid.calculate_weights_at_points(points, k=3)
    
    # plot the atomic weights
    plt.imshow(np.max(mweights,axis=0).reshape((N,N)),
               extent=(-sz,sz,-sz,sz))
    plt.xlabel('x')
    plt.xlabel('y')
    plt.colorbar()
    plt.tight_layout()

if __name__ == '__main__':
    main()