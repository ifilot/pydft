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

        answers = [
            -110.96476654905626,
            -111.14751240002687,
            -111.14683799873168,
            -111.15048475412650,
            -111.14985895386218,
            -111.14831640666813,
        ]

        for fdpts,answer in zip([3,5,7,9,11,13], answers):
            dft = DFT(mol, basis='sto3g', fdpts=fdpts)
            res = dft.scf()

            np.testing.assert_almost_equal(res['energy'], answer, 4)

if __name__ == '__main__':
    unittest.main()
