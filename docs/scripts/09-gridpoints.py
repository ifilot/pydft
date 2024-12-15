import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    
    for n in [32,64,128]:
        dft = DFT(co, basis='sto3g', nshells=n)
        en = dft.scf(1e-6, verbose=False)

if __name__ == '__main__':
    main()