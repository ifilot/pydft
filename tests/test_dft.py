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
        energy = dft.scf()

        answer = -2.809567
        np.testing.assert_almost_equal(energy, answer, 4)

    def test_h2o(self):
        """
        Test DFT calculation of water molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('H2O')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -74.93750649033662
        np.testing.assert_almost_equal(energy, answer, 4)

    def test_co(self):
        """
        Test DFT calculation of CO molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('CO')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -111.14683799873168
        np.testing.assert_almost_equal(energy, answer, 4)
        
    def test_bf3(self):
        """
        Test DFT calculation of BF3 molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('bf3')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -318.2753601057624
        np.testing.assert_almost_equal(energy, answer, 4)
    
    def test_ch4(self):
        """
        Test DFT calculation of methane molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('ch4')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -39.80546943950668
        np.testing.assert_almost_equal(energy, answer, 4)
    
    def test_co2(self):
        """
        Test DFT calculation of co2 molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('co2')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -185.01476596292474
        np.testing.assert_almost_equal(energy, answer, 4)
    
    def test_ethylene(self):
        """
        Test DFT calculation of an ethylene molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('ethylene')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -77.16253628136911
        np.testing.assert_almost_equal(energy, answer, 4)
        
    def test_h2(self):
        """
        Test DFT calculation of hydrogen molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('h2')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -1.1570136
        np.testing.assert_almost_equal(energy, answer, 4)
        
    def test_lih(self):
        """
        Test DFT calculation of LiH molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('lih')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        energy = dft.scf()

        answer = -7.865828015750928
        np.testing.assert_almost_equal(energy, answer, 4)
           
    # def test_benzene(self):
    #     """
    #     Test DFT calculation of benzene molecule
    #     """
    #     mol_builder = MoleculeBuilder()
    #     mol = mol_builder.from_name('benzene')

    #     # construct dft object
    #     dft = DFT(mol, basis='sto3g')
    #     energy = dft.scf()

    #     answer = -228.0695756259556
    #     np.testing.assert_almost_equal(energy, answer, 4)

if __name__ == '__main__':
    unittest.main()
