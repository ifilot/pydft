# -*- coding: utf-8 -*-
from pydft import MoleculeBuilder, DFT

mol = MoleculeBuilder().from_name("h2")
dft = DFT(mol, basis='sto3g')
res = dft.scf(1e-4, verbose=True)
print("Total electronic energy: %f Ht" % res['energy'])
dft.print_time_statistics()

print(res['ekin'])
print(res['enuc'])
print(res['erepe'])
print(res['exc'])
print(res['ex'])
print(res['ec'])