import unittest
import numpy as np
from pydft import MoleculeBuilder, MolecularGrid

class TestCustomBasisSet(unittest.TestCase):
        
    def test_bragg_slater_bf3(self):
        """
        Test grabbing correct Bragg-Slater radii
        """
        # construct molecule
        mol = MoleculeBuilder().from_name('bf3')
        cgfs, atoms = mol.build_basis('sto3g')
        
        # construct molecular grid
        molgrid = MolecularGrid(atoms, cgfs)
        molgrid.initialize()

        ANGSTROM2BOHR = 1.88973

        # test Bragg-Slater radius for B
        np.testing.assert_almost_equal(
            molgrid.get_atomic_grid(0).get_bragg_slater_radius(), 
            0.85 * ANGSTROM2BOHR)
        
        # test Bragg-Slater radius for F
        np.testing.assert_almost_equal(
            molgrid.get_atomic_grid(1).get_bragg_slater_radius(), 
            0.50 * ANGSTROM2BOHR)
        
    def test_bragg_slater_ch4(self):
        """
        Test grabbing correct Bragg-Slater radii
        """
        # construct molecule
        mol = MoleculeBuilder().from_name('ch4')
        cgfs, atoms = mol.build_basis('sto3g')
        
        # construct molecular grid
        molgrid = MolecularGrid(atoms, cgfs)
        molgrid.initialize()

        ANGSTROM2BOHR = 1.88973

        # test Bragg-Slater radius for B
        np.testing.assert_almost_equal(
            molgrid.get_atomic_grid(0).get_bragg_slater_radius(), 
            0.70 * ANGSTROM2BOHR)
        
        # test Bragg-Slater radius for F
        np.testing.assert_almost_equal(
            molgrid.get_atomic_grid(1).get_bragg_slater_radius(), 
            0.35 * ANGSTROM2BOHR)

if __name__ == '__main__':
    unittest.main()