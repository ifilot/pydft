import unittest
import numpy as np

from pydft import MoleculeBuilder, DFT

class TestXC(unittest.TestCase):

    def test_svwn5(self):
        """
        Test DFT calculation of Helium atom
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('He')

        # construct dft object
        dft = DFT(mol, basis='sto3g', functional='svwn5')
        res = dft.scf()

        answer = -2.8081735181814724
        np.testing.assert_almost_equal(res['energy'], answer, 4)

    def test_pbe(self):
        """
        Test DFT calculation of Helium atom
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('He')

        # construct dft object
        dft = DFT(mol, basis='sto3g', functional='pbe')
        res = dft.scf()

        answer = -2.8287483384734955
        np.testing.assert_almost_equal(res['energy'], answer, 4)