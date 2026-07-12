import unittest
import contextlib
import io
from unittest.mock import patch

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
        with self.assertLogs('pydft.dft', level='INFO') as logs:
            res = dft.scf(verbose=True)

        answer = -2.7704595622068138
        np.testing.assert_almost_equal(res['energy'], answer, 4)
        self.assertTrue(any('| E =' in msg for msg in logs.output))

        pts = np.array([[0.1, 0.0, 0.0], [0.0, 0.2, 0.1]])
        self.assertEqual(dft.get_density_at_points(pts).shape, (2,))
        self.assertEqual(dft.get_gradient_at_points(pts).shape, (2, 3))

        with self.assertRaisesRegex(ValueError, 'Nx3'):
            dft.get_density_at_points(np.array([0.0, 0.0, 0.0]))

        with self.assertRaisesRegex(ValueError, 'Nx3'):
            dft.get_gradient_at_points(np.array([[0.0, 0.0]]))

        with contextlib.redirect_stdout(io.StringIO()) as stream:
            dft.print_time_statistics()
        self.assertIn('Construction times', stream.getvalue())

    def test_h2o(self):
        """
        Test DFT calculation of water molecule
        """
        mol_builder = MoleculeBuilder()
        mol = mol_builder.from_name('H2O')

        # construct dft object
        dft = DFT(mol, basis='sto3g')
        res = dft.scf()

        answer = -74.72985279614994
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

        answer = -110.85461362334834
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

        answer = -317.55352572879667
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

        answer = -39.60074741004088
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

        answer = -184.5409746994328
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

        answer = -76.82946120370394
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

        answer = -1.1204764001484104
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

        answer = -7.789199497213202
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

        answer = -227.14299165089105
        np.testing.assert_almost_equal(res['energy'], answer, 4)

    def test_scalar_grid_settings(self):
        mol = MoleculeBuilder().from_name('He')

        dft = DFT(mol, basis='sto3g', nshells=8, nangpts=50, lmax=5)
        res = dft.scf(tol=1e-4)

        self.assertTrue(np.isfinite(res['energy']))

    def test_eigenvalue_failure_reports_scf_iteration(self):
        mol = MoleculeBuilder().from_name('He')
        dft = DFT(mol, basis='sto3g', nshells=8, nangpts=50, lmax=5)
        dft._DFT__setup()

        with patch('pydft.dft.np.linalg.eigh', side_effect=np.linalg.LinAlgError('boom')):
            with self.assertLogs('pydft.dft', level='ERROR'):
                with self.assertRaisesRegex(np.linalg.LinAlgError, 'SCF iteration 0'):
                    dft._DFT__iterate(0)

if __name__ == '__main__':
    unittest.main()
