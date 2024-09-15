# -*- coding: utf-8 -*-
import os,sys

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, DFT

co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g')
en = dft.scf(1e-4)
print("Total electronic energy: %f Ht" % en)

res = dft.get_data()
print(res.keys())