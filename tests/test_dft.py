import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import MoleculeBuilder, DFT

class TestDFT(unittest.TestCase):

    def test_helium(self):
        """
        Test DFT calculation of Helium atom
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('He')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -2.809567
        np.testing.assert_almost_equal(res['energy'], answer, 4)

    def test_h2o(self):
        """
        Test DFT calculation of water molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('H2O')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -74.9357662036297
        np.testing.assert_almost_equal(res['energy'], answer, 4)

    def test_co(self):
        """
        Test DFT calculation of CO molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('CO')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -111.14308742936126
        np.testing.assert_almost_equal(res['energy'], answer, 4)
        
    def test_bf3(self):
        """
        Test DFT calculation of BF3 molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('bf3')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -318.24725157047226
        np.testing.assert_almost_equal(res['energy'], answer, 4)
    
    def test_ch4(self):
        """
        Test DFT calculation of methane molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('ch4')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -39.802673700611706
        np.testing.assert_almost_equal(res['energy'], answer, 4)
    
    def test_co2(self):
        """
        Test DFT calculation of co2 molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('co2')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -185.0058824009409
        np.testing.assert_almost_equal(res['energy'], answer, 4)
    
    def test_ethylene(self):
        """
        Test DFT calculation of an ethylene molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('ethylene')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -77.15675419076709
        np.testing.assert_almost_equal(res['energy'], answer, 4)
        
    def test_h2(self):
        """
        Test DFT calculation of hydrogen molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('h2')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -1.1568407835459502
        np.testing.assert_almost_equal(res['energy'], answer, 4)
        
    def test_lih(self):
        """
        Test DFT calculation of LiH molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('lih')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -7.865528839441138
        np.testing.assert_almost_equal(res['energy'], answer, 4)
           
    def test_benzene(self):
        """
        Test DFT calculation of benzene molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('benzene')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -228.03953311038978
        np.testing.assert_almost_equal(res['energy'], answer, 4)

if __name__ == '__main__':
    unittest.main()
