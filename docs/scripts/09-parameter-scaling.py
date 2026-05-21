import numpy as np
from pydft import MoleculeBuilder, DFT
import matplotlib.pyplot as plt
import time

def main():
    # perform DFT calculation on the CO molecule
    co = MoleculeBuilder().from_name("CO")
    
    plt.figure(dpi=144)
    
    # explore scaling with number of radial shells
    t = []
    print('Spherical shells')
    nvals = [8,16,32,64,92,128,256]
    for n in nvals:
        st = time.time()
        dft = DFT(co, basis='sto3g', nshells={'C': n, 'O': n})
        res = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (n, res['energy'], dt))
    print()
    
    plt.loglog(nvals, t, 'o--', label=r'$N_{r}$')
    
    # explore scaling with number of angular points
    t = []    
    print('Angular points')
    avals = [38,74,110,194,302,590]
    for a in avals:
        st = time.time()
        dft = DFT(co, basis='sto3g',
                  nshells={'C': 64, 'O': 64},
                  nangpts={'C': a, 'O': a})
        res = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (a, res['energy'], dt))
    print()
        
    plt.loglog(avals, t, 'o--', label=r'$N_{a}$')
        
    # explore scaling with lmax value
    t = []
    print('Lmax value')
    lvals = [2,3,4,6,8,12,16,24]
    for lmax in lvals:
        st = time.time()
        dft = DFT(co, basis='sto3g',
                  nshells={'C': 64, 'O': 64},
                  nangpts={'C': 590, 'O': 590},
                  lmax={'C': lmax, 'O': lmax})
        res = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (lmax, res['energy'], dt))
    print()
    
    plt.loglog(lvals, t, 'o--', label=r'$l_{\text{max}}$')
    
    # explore scaling with number of discretization points used in the
    # finite difference scheme to construct the Hartree potential
    t = []
    print('Number of finite difference coefficients')
    fdvals = [3,5,7,9,11,13,15,17]
    for fd in fdvals:
        st = time.time()
        dft = DFT(co, basis='sto3g',
                  nshells={'C': 64, 'O': 64},
                  nangpts={'C': 194, 'O': 194},
                  fdpts=fd)
        res = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (fd, res['energy'], dt))
    print()
        
    plt.loglog(fdvals, t, 'o--', label=r'$N_{\text{fd}}$')
    plt.grid()
    plt.xlabel('Value [-]')
    plt.ylabel('Computation time [s]')
    plt.legend()

if __name__ == '__main__':
    main()
