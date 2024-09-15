# -*- coding: utf-8 -*-
import os,sys
import numpy as np
import time

# add a reference to load the module
ROOT = os.path.dirname(__file__)
sys.path.insert(1, os.path.join(ROOT, '..'))

from pydft import MoleculeBuilder, DFT

#
# Example: Calculate total electronic energy for CO using standard
#          settings.
#

CO = MoleculeBuilder().from_name("CO")

rshells = [8,16,32,64,128]
ang = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

energies = np.zeros((len(rshells), len(ang)))
timestats = np.zeros((len(rshells), len(ang)))

for i,nr in enumerate(rshells):
    for j,a in enumerate(ang):
        start_time = time.time()
        dft = DFT(CO, basis='sto3g', nshells=nr, nangpts=a)
        en = dft.scf(1e-4)
        end_time = time.time()
        elapsed_time = end_time - start_time
        timestats[i,j] = elapsed_time
        energies[i,j] = en
        
        print(nr, a, en, elapsed_time)
        
np.savetxt('energies.txt', energies)
np.savetxt('timestats.txt', timestats)