import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import MoleculeBuilder, DFT

class TestParallel(unittest.TestCase):

           
    def test_benzene(self):
        """
        Test DFT calculation of benzene molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('benzene')

        # construct dft object
        dft = DFT(mol, basis='sto3g', parallel=True)
        energy = dft.scf()

        answer = -228.06536332917378
        np.testing.assert_almost_equal(energy, answer, 4)

if __name__ == '__main__':
    unittest.main()
