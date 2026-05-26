from pydft import MoleculeBuilder, DFT

# perform DFT calculation on the CO molecule
co = MoleculeBuilder().from_name("CO")
dft = DFT(co, basis='sto3g')
res = dft.scf(1e-6, verbose=False)

# print total energy
print("Total electronic energy:     %12.6f Ht" % res['energy'])
