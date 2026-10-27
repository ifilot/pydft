import unittest
import numpy as np

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

        answers = [-110.5669, -110.8619, -110.8546, -110.8678, -110.8716, -110.8589]

        res = []
        for fdpts in [3,5,7,9,11,13]:
            dft = DFT(mol, basis='sto3g', fdpts=fdpts)
            res.append(dft.scf()['energy'])

        np.testing.assert_almost_equal(res, answers, 4)

if __name__ == '__main__':
    unittest.main()
