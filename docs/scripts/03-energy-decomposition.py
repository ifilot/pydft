import numpy as np
from pydft import MoleculeBuilder, DFT

co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g', verbose=False)
en = dft.scf(1e-4)
print("Total electronic energy:     %12.6f Ht" % en)
print()

# retrieve molecular matrices
res = dft.get_data()
P = res['P']
T = res['T']
V = res['V']
J = res['J']

# calculate energy terms
Et = np.einsum('ji,ij', P, T)
Ev = np.einsum('ji,ij', P, V)
Ej = 0.5 * np.einsum('ji,ij', P, J)
Ex = res['Ex']
Ec = res['Ec']
Exc = res['Exc']
Enuc = res['enucrep']

print('Kinetic energy:              %12.6f Ht' % Et)
print('Nuclear attraction:          %12.6f Ht' % Ev)
print('Electron-electron repulsion: %12.6f Ht' % Ej)
print('Exchange energy:             %12.6f Ht' % (Ex))
print('Correlation energy:          %12.6f Ht' % (Ec))
print('Exchange-correlation energy: %12.6f Ht' % (Exc))
print('Nucleus-nucleus repulsion:   %12.6f Ht' % (Enuc))
print()
print('Sum: %12.6f Ht' % (Et + Ev + Ej + Exc + Enuc))
