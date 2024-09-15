import numpy as np
from pydft import MoleculeBuilder, DFT

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")

# use SVWM XC functional
dft = DFT(co, basis='sto3g', verbose=False, functional='svwn5')
print('SVWN: ', dft.scf(1e-5), 'Ht')

# use PBE XC functional
dft = DFT(co, basis='sto3g', verbose=False, functional='pbe')
print('PBE: ', dft.scf(1e-5), 'Ht')