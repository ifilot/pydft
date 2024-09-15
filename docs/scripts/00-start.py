import numpy as np
from pydft import MoleculeBuilder, DFT

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g', verbose=False)
en = dft.scf(1e-6)

# print total energy
print("Total electronic energy:     %12.6f Ht" % en)