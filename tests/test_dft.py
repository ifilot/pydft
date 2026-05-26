import unittest
import numpy as np

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

        answer = -2.8081735181814724
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

        answer = -74.92583960745394
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

        answer = -111.1305117481223
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

        answer = -318.19442532820096
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

        answer = -39.790602071993334
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

        answer = -184.97691522055553
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

        answer = -77.13648424441982
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

        answer = -1.1562627795279132
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

        answer = -7.862755794011036
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

        answer = -227.9563844056453
        np.testing.assert_almost_equal(res['energy'], answer, 4)

if __name__ == '__main__':
    unittest.main()
