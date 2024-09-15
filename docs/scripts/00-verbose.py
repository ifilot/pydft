import numpy as np
from pydft import MoleculeBuilder, DFT

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g', verbose=True)
en = dft.scf(1e-5)