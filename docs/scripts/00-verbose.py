import logging

from pydft import MoleculeBuilder, DFT

logging.basicConfig(level=logging.INFO, format="%(message)s")

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g')
res = dft.scf(1e-5, verbose=True)
