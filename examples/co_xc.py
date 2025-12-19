# -*- coding: utf-8 -*-
from pydft import MoleculeBuilder,DFT

mol = MoleculeBuilder().from_name("CO")

dft = DFT(mol, basis='sto3g', functional='svwn5')
res = dft.scf(1e-4)
print("Total electronic energy (SVWN5): %f Ht" % res['energy'])

dft = DFT(mol, basis='sto3g', functional='pbe')
res = dft.scf(1e-4)
print("Total electronic energy (PBE): %f Ht" % res['energy'])