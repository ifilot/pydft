import unittest
import numpy as np

from pydft import MoleculeBuilder, DFT
from pydft.xcfunctionals import Functionals

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

        answer = -2.7704595622068138
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

    def test_invalid_functional_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Illegal XC-functional'):
            Functionals('not-a-functional')

    def test_pbe_numerical_exchange_derivative_path(self):
        functionals = Functionals('pbe')
        rho = np.array([0.2, 0.5])
        gamma = np.array([0.01, 0.02])

        derivative = functionals._Functionals__pbe_x_deriv_numerical(rho, gamma)

        self.assertEqual(derivative.shape, rho.shape)
        self.assertTrue(np.all(np.isfinite(derivative)))
