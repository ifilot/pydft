# -*- coding: utf-8 -*-
import os,sys

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, DFT

#
# Example: Calculate total electronic energy for CO using standard
#          settings.
#

CO = MoleculeBuilder().from_name("H2")
dft = DFT(CO, basis='sto3g')
en = dft.scf(1e-4)
print("Total electronic energy: %f Ht" % en)