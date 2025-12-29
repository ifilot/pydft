import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import MoleculeBuilder, DFT

class TestFiniteDifferenceSchemes(unittest.TestCase):

    def test_co(self):
        """
        Test DFT calculation of CO molecule using various finite difference
        schemes (always diagonal matrices, but with different number of
                 non-zero diagonals)
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('CO')

        answers = [-110.961 , -111.1438, -111.1431, -111.1467, -111.1461, -111.1445]

        res = []
        for fdpts in [3,5,7,9,11,13]:
            dft = DFT(mol, basis='sto3g', fdpts=fdpts)
            res.append(dft.scf()['energy'])

        np.testing.assert_almost_equal(res, answers, 4)

if __name__ == '__main__':
    unittest.main()
