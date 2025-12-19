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
        res = dft.scf(1e-4)
        print("Total electronic energy: %f Ht" % res['energy'])

        # retrieve molecular matrices
        res = dft.get_data()
        P = res['density']
        T = res['kinetic']
        V = res['nuclear']
        J = res['hartree']

        # calculate energy terms
        Et = np.einsum('ji,ij', P, T)
        Ev = np.einsum('ji,ij', P, V)
        Ej = 0.5 * np.einsum('ji,ij', P, J)
        Ex = res['ex']
        Ec = res['ec']
        Exc = res['exc']
        Enuc = res['enucrep']

        # print('Kinetic energy:              %12.6f' % Et)
        # print('Nuclear attraction:          %12.6f' % Ev)
        # print('Electron-electron repulsion: %12.6f' % Ej)
        # print('Exchange energy:             %12.6f' % (Ex))
        # print('Correlation energy:          %12.6f' % (Ec))
        # print('Exchange-correlation energy: %12.6f' % (Exc))
        # print('Nucleus-nucleus repulsion:   %12.6f' % (Enuc))

        esum = Et + Ev + Ej + Exc + Enuc
        np.testing.assert_almost_equal(esum, res['energy'], decimal=5)

    def test_key_existence(self):
        """
        Test whether all expected result keys are present in the data dictionary.
        """
        h2 = MoleculeBuilder().from_name("H2")
        dft = DFT(h2, basis="sto3g")
        dft.scf(1e-4)

        data = dft.get_data()

        expected_keys = [
            # system
            "mol",
            "nuclei",
            "cgfs",

            # energies
            "energy",
            "energies",
            "ekin",
            "enuc",
            "enucrep",
            "ex",
            "ec",
            "exc",

            # orbitals
            "orbc",
            "orbe",

            # matrices
            "overlap",
            "kinetic",
            "nuclear",
            "hcore",
            "density",
            "fock",
            "hartree",
            "xc",

            # timing
            "time_stats",
        ]

        for key in expected_keys:
            assert key in data, f"Missing key: {key}"

        # also test for false positives
        assert "fakekey" not in data

if __name__ == '__main__':
    unittest.main()