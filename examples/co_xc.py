# -*- coding: utf-8 -*-
from pydft import MoleculeBuilder,DFT

#
# Example: Calculate total electronic energy for CO using different
#          XC-functionals.
#

CO = MoleculeBuilder().from_name("CO")
dft = DFT(CO, basis='sto3g', functional='svwn5')
en = dft.scf(1e-4)
print("Total electronic energy (SVWN5): %f Ht" % en)

CO = MoleculeBuilder().from_name("CO")
dft = DFT(CO, basis='sto3g', functional='pbe')
en = dft.scf(1e-4)
print("Total electronic energy (PBE): %f Ht" % en)