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
        dft = DFT(co, basis='sto3g', nshells=n)
        en = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (n, en, dt))
    print()
    
    plt.loglog(nvals, t, 'o--', label=r'$N_{r}$')
    
    # explore scaling with number of angular points
    t = []    
    print('Angular points')
    avals = [38,74,110,194,302,590]
    for a in avals:
        st = time.time()
        dft = DFT(co, basis='sto3g', nshells=64, nangpts=a)
        en = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (a, en, dt))
    print()
        
    plt.loglog(avals, t, 'o--', label=r'$N_{a}$')
        
    # explore scaling with lmax value
    t = []
    print('Lmax value')
    lvals = [2,3,4,6,8,12,16,24]
    for lmax in lvals:
        st = time.time()
        dft = DFT(co, basis='sto3g', nshells=64, nangpts=590, lmax=lmax)
        en = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (lmax, en, dt))
    print()
    
    plt.loglog(lvals, t, 'o--', label=r'$l_{\text{max}}$')
    
    # explore scaling with number of discretization points used in the
    # finite difference scheme to construct the Hartree potential
    t = []
    print('Number of finite difference coefficients')
    fdvals = [3,5,7,9,11,13,15,17]
    for fd in fdvals:
        st = time.time()
        dft = DFT(co, basis='sto3g', nshells=64, nangpts=194, fdpts=fd)
        en = dft.scf(1e-6, verbose=False)
        dt = time.time() - st
        t.append(dt)
        print("%3i %12.4f %5.4f" % (fd, en, dt))
    print()
        
    plt.loglog(lvals, t, 'o--', label=r'$N_{\text{fd}}$')
    plt.grid()
    plt.xlabel('Value [-]')
    plt.ylabel('Computation time [s]')
    plt.legend()

if __name__ == '__main__':
    main()