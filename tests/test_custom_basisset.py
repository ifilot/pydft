import unittest
import numpy as np
from pydft import MoleculeBuilder, DFT
from pyqint import cgf

class TestCustomBasisSet(unittest.TestCase):
        
    def test_h2(self):
        """
        Test DFT calculation of hydrogen molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('h2')

        cgfs = []
        for n in mol.get_nuclei():
            _cgf = cgf(n[0])

            _cgf.add_gto(0.154329, 3.425251, 0, 0, 0)
            _cgf.add_gto(0.535328, 0.623914, 0, 0, 0)
            _cgf.add_gto(0.444635, 0.168855, 0, 0, 0)

            cgfs.append(_cgf)

        # construct dft object
        dft = DFT(mol, basis=cgfs)
        energy = dft.scf()

        answer = -1.1570136
        np.testing.assert_almost_equal(energy, answer, 4)

if __name__ == '__main__':
    unittest.main()
