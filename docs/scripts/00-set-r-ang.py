from pydft import MoleculeBuilder, DFT

mol = MoleculeBuilder.from_name("CO")
dft = DFT(mol, basis='sto3g',
    nshells={'C': 32, 'O': 32},
    nangpts={'C': 230, 'O': 230}
)
res = dft.scf(1e-6, verbose=False)

print("Total electronic energy:     %12.6f Ht" % res['energy'])