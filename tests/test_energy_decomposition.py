import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import DFT, MoleculeBuilder
import pydft

class TestEnergyDecomposition(unittest.TestCase):
    """
    Test existence of dictionary keys
    """

    def test_energy_decomposition(self):
        """
        Test energy decomposition of a molecule
        """
        co = MoleculeBuilder().from_name("CO")
        dft = DFT(co, basis='sto3g')
        en = dft.scf(1e-4)
        print("Total electronic energy: %f Ht" % en)

        # retrieve molecular matrices
        res = dft.get_data()
        P = res['P']
        T = res['T']
        V = res['V']
        J = res['J']

        # calculate energy terms
        Et = np.einsum('ji,ij', P, T)
        Ev = np.einsum('ji,ij', P, V)
        Ej = 0.5 * np.einsum('ji,ij', P, J)
        Ex = res['Ex']
        Ec = res['Ec']
        Exc = res['Exc']
        Enuc = res['enucrep']

        # print('Kinetic energy:              %12.6f' % Et)
        # print('Nuclear attraction:          %12.6f' % Ev)
        # print('Electron-electron repulsion: %12.6f' % Ej)
        # print('Exchange energy:             %12.6f' % (Ex))
        # print('Correlation energy:          %12.6f' % (Ec))
        # print('Exchange-correlation energy: %12.6f' % (Exc))
        # print('Nucleus-nucleus repulsion:   %12.6f' % (Enuc))

        esum = Et + Ev + Ej + Exc + Enuc
        np.testing.assert_almost_equal(esum, en, decimal=5)

    def test_key_existence(self):
        """
        Test whether all keys are in dictionary
        """
        h2 = MoleculeBuilder().from_name("H2")
        dft = DFT(h2, basis='sto3g')
        en = dft.scf(1e-4)
        mats = dft.get_data()
        
        keys = ['S', 'T', 'V', 'C', 'J', 'P', 'XC', 'F', 'Exc', 'Ex', 'Ec', 'energies', 'energy', 'orbc', 'orbe', 'enucrep']
        
        for key in keys:
            assert key in mats
            
        # also test for false positives
        assert 'fakekey' not in mats

if __name__ == '__main__':
    unittest.main()