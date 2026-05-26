# -*- coding: utf-8 -*-
from pydft import MoleculeBuilder, DFT

mol = MoleculeBuilder().from_name("CO")
dft = DFT(mol, basis='sto3g', nshells={'C': 32, 'O': 32}, nangpts={'C': 770, 'O': 770})
res = dft.scf(1e-4, verbose=True)
print("Total electronic energy: %f Ht" % res['energy'])
dft.print_time_statistics()