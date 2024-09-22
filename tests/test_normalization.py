import unittest
import numpy as np
from pydft import MoleculeBuilder, DFT

class TestNormalization(unittest.TestCase):
    """
    Test whether normalization yields **exactly** the number of electrons of
    the system
    """

    def test_helium(self):
        """
        Test DFT calculation of Helium atom
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('He')
        Nelec = 2.0

        # construct dft object
        dft = DFT(mol, basis='sto3g', normalize=False)
        dft.scf()
        nelec = dft.get_molgrid_copy().count_electrons()

        err_no_norm = np.abs(Nelec - nelec)

        # construct dft object
        dft = DFT(mol, basis='sto3g', normalize=True)
        dft.scf()
        nelec = dft.get_molgrid_copy().count_electrons()

        err_norm = np.abs(Nelec - nelec)
        self.assertTrue(err_norm < err_no_norm)
        np.testing.assert_almost_equal(nelec, Nelec, 8)

    def test_methane(self):
        """
        Test DFT calculation of methane molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('CH4')
        Nelec = 10.0

        # construct dft object
        dft = DFT(mol, basis='sto3g', normalize=False)
        dft.scf()
        nelec = dft.get_molgrid_copy().count_electrons()

        err_no_norm = np.abs(Nelec - nelec)

        # construct dft object
        dft = DFT(mol, basis='sto3g', normalize=True)
        dft.scf()
        nelec = dft.get_molgrid_copy().count_electrons()

        err_norm = np.abs(Nelec - nelec)
        self.assertTrue(err_norm < err_no_norm)
        np.testing.assert_almost_equal(nelec, Nelec, 8)