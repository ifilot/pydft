# -*- coding: utf-8 -*-
from pydft import MoleculeBuilder, DFT

mol = MoleculeBuilder().from_name("H2")
dft = DFT(mol, basis='sto3g')
res = dft.scf(1e-4, verbose=True)
print("Total electronic energy: %f Ht" % res['energy'])