import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

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
        energy = dft.scf()

        answer = -2.809567
        np.testing.assert_almost_equal(energy, answer, 4)

    def test_pbe(self):
        """
        Test DFT calculation of Helium atom
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('He')

        # construct dft object
        dft = DFT(mol, basis='sto3g', functional='pbe')
        energy = dft.scf()

        answer = -2.8301366521757787
        np.testing.assert_almost_equal(energy, answer, 4)