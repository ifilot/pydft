from pydft import MoleculeBuilder, DFT

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")

# use SVWM XC functional
dft = DFT(co, basis='sto3g', functional='svwn5')
print('SVWN: ', dft.scf(1e-5)['energy'], 'Ht')

# use PBE XC functional
dft = DFT(co, basis='sto3g', functional='pbe')
print('PBE: ', dft.scf(1e-5)['energy'], 'Ht')