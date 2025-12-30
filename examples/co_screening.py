# -*- coding: utf-8 -*-
import numpy as np
import time
from pydft import MoleculeBuilder, DFT
import os
import matplotlib.pyplot as plt

ROOT = os.path.dirname(__file__)
energiesfile = os.path.join(ROOT, 'energies.txt')
timefile = os.path.join(ROOT, 'timestats.txt')

rshells = [8,16,32,64,128]
ang = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 
    590, 770, 974]

if not os.path.exists(energiesfile) or not os.path.exists(timefile):

    mol = MoleculeBuilder().from_name("CO")
    energies = np.zeros((len(rshells), len(ang)))
    timestats = np.zeros((len(rshells), len(ang)))

    for i,nr in enumerate(rshells):
        for j,a in enumerate(ang):
            start_time = time.perf_counter()
            dft = DFT(mol, basis='sto3g', nshells={'C': nr, 'O': nr}, nangpts={'C': a, 'O': a})
            res = dft.scf(1e-4)
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            timestats[i,j] = elapsed_time
            energies[i,j] = res['energy']
            
            print(nr, a, res['energy'], elapsed_time)
        
    np.savetxt(energiesfile, energies)
    np.savetxt(timefile, timestats)

energies = np.loadtxt(energiesfile).reshape((len(rshells), len(ang)))
timestats = np.loadtxt(timefile).reshape((len(rshells), len(ang)))

plt.figure(dpi=144)
nfit = 20
for i,r in enumerate(rshells):
    (line,) = plt.loglog(ang, timestats[i], 'o', alpha=0.5, label=r'$N_{\text{shells}} = %i$' % r)

    # Select last nfit points
    x_fit = ang[-nfit:]
    y_fit = timestats[i, -nfit:]

    # Fit in log-log space
    logx = np.log10(x_fit)
    logy = np.log10(y_fit)
    slope, intercept = np.polyfit(logx, logy, 1)

    # Reconstruct fitted line in linear space
    xx = np.linspace(100, 6000, 20)
    y_trend = 10**intercept * xx**slope

    color = line.get_color()

    # Plot trend line
    plt.loglog(
        xx,
        y_trend,
        linestyle='--',
        color=color,
        linewidth=2,
        alpha=0.5,
        label=r'Fit: $n \approx %.2f$' % slope
    )
plt.xlabel('Number of angular points')
plt.ylabel('Wall clock time [s]')
plt.grid(linestyle='--')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'co_trend_angplot.png'))
plt.close()

plt.figure(dpi=144)
ngridpts = np.outer(rshells, ang).flatten()
times = timestats.flatten()
mask = ngridpts > 100000
plt.loglog(ngridpts, times, 'o', alpha=0.5, label='Data')
logx = np.log10(ngridpts[mask])
logy = np.log10(times[mask])
slope, intercept = np.polyfit(logx, logy, 1)
x_trend = np.linspace(10000, 1e6, 30)
y_trend = 10**intercept * x_trend**slope
plt.loglog(x_trend, y_trend, '--', alpha=0.9, color='black',
           label=r'Fit: $n \approx %.2f$' % slope)
plt.xlabel('Number of grid points')
plt.ylabel('Wall clock time [s]')
plt.grid(linestyle='--')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'co_trend_time.png'))
plt.close()

plt.figure(dpi=144)
nfit = 20
for i,r in enumerate(rshells):
    plt.loglog(ang, np.abs((energies[i] - np.min(energies[-1,-1])) / np.min(energies[-1,-1])), 'o', alpha=0.5, label=r'$N_{\text{shells}} = %i$' % r)
plt.xlabel('Number of angular points')
plt.ylabel('Relative error [-]')
plt.grid(linestyle='--')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(ROOT, 'co_energy_angplot.png'))
plt.show()