import numpy as np
import sys,os

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, DFT
import pydft

print(pydft.__version__)

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g', parallel=True)
en = dft.scf(1e-5, verbose=True)