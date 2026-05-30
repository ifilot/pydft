import unittest
import numpy as np

from pyqint import PyQInt, Molecule
from pydft import MoleculeBuilder, MolecularGrid, DFT

class TestMolecularGrid(unittest.TestCase):

    def test_grid_points(self):
        mol = MoleculeBuilder().from_name('h2')
        cgfs, atoms = mol.build_basis('sto3g')
        
        # construct molecular grid
        molgrid = MolecularGrid([at for at in mol], 
                                cgfs,
                                nshells = {'H' : 20},
                                nangpts = {'H' : 50},
                                lmax = {'H' : 5})
        molgrid.initialize()
        mgpts = molgrid.get_molecular_grid_coordinates()
        agpts = molgrid.get_grid_coordinates()

        np.testing.assert_equal(len(mgpts), np.sum([len(g) for g in agpts]))
        np.testing.assert_equal(mgpts[0:1000], agpts[0])
        np.testing.assert_equal(mgpts[1000:], agpts[1])

        # test radial grid points
        r = molgrid.get_rgridpoints()
        np.testing.assert_equal(r.shape, (2,2000))

        # test relative coordinates
        relgpts = molgrid.get_atomic_relative_grid_coordinates()
        np.testing.assert_equal(relgpts.shape, (2,2000,3))
        np.testing.assert_equal(relgpts[0], mgpts - mol[0][1])
        np.testing.assert_equal(relgpts[1], mgpts - mol[1][1])

    def test_amplitudes(self):
        mol = Molecule('Fe')
        mol.add_atom('Fe', 0.0, 0.0, 0.0)
        mol.add_atom('Fe', 0.0, 0.0, 5.0)
        mol.add_atom('Fe', 0.0, 0.0, 10.0)
        cgfs, nuclei = mol.build_basis('sto3g')
        integrator = PyQInt()
        ao = cgfs[-7]
                
        # construct molecular grid
        molgrid = MolecularGrid([at for at in mol],
                                cgfs,
                                nshells = {'Fe' : 32},
                                nangpts = {'Fe' : 110},
                                lmax = {'Fe' : 8})
        molgrid.initialize()

        # test caching of amplitudes
        for i in range(len(mol)):
            amps = integrator.plot_wavefunction(molgrid.get_grid_coordinates()[i], [1], [ao])
            amps_cached = molgrid.get_cached_amplitudes()[i][-7]
            np.testing.assert_almost_equal(amps, amps_cached)

    def test_weight_functions(self):
        """
        Test DFT calculation of Helium atom
        """
        # construct molecule
        mol = MoleculeBuilder().from_name('bf3')
        cgfs, atoms = mol.build_basis('sto3g')
        
        # construct molecular grid
        molgrid = MolecularGrid([at for at in mol], 
                                cgfs,
                                nshells = {'B' : 32, 'F' : 32},
                                nangpts = {'B' : 110, 'F' : 110})
        
        # produce grid of sampling points to calculate the atomic
        # weight coefficients for
        N = 10
        sz = 8
        x = np.linspace(-sz,sz,N)
        xv,yv = np.meshgrid(x,x)
        points = np.array([[x,y,0] for x,y in zip(xv.flatten(),yv.flatten())])
        
        # calculate the atomic weights
        mweights = molgrid.calculate_weights_at_points(points, k=3)
        
        # verify results
        np.testing.assert_almost_equal(mweights[0,40], 
                                       7.438664028803838e-06, 4)
        np.testing.assert_almost_equal(mweights[0,50], 
                                       0.0023096748652996083, 4)
        np.testing.assert_almost_equal(mweights[0,60], 
                                       0.03962756655918886, 4)
        np.testing.assert_almost_equal(mweights[0,70], 
                                       0.08886587061825403, 4)
        np.testing.assert_almost_equal(mweights[0,80], 
                                       0.04752230028517094, 4)
        
        np.testing.assert_almost_equal(mweights[1,70],
                                       0.41031416799944265, 4)
        np.testing.assert_almost_equal(mweights[1,80],
                                       0.7991044370004883, 4)

    def test_arbitrary_point_fields_and_grid_analysis_helpers(self):
        mol = MoleculeBuilder().from_name('He')
        dft = DFT(mol, basis='sto3g', nshells=8, nangpts=50, lmax=5)
        res = dft.scf(tol=1e-4)
        molgrid = dft.get_molgrid_copy()

        pts = np.array([
            [0.1, 0.0, 0.0],
            [0.0, 0.2, 0.1],
        ])
        P = res['P']

        self.assertEqual(molgrid.get_density_at_points(pts, P).shape, (2,))
        self.assertEqual(molgrid.get_gradient_at_points(pts, P).shape, (2, 3))
        self.assertEqual(molgrid.get_amplitude_at_points(pts, res['orbc'][:, 0]).shape, (2,))
        self.assertEqual(molgrid.calculate_coulomb_potential_at_points(pts).shape, (2,))
        self.assertEqual(molgrid.get_exchange_potential_at_points(pts, P).shape, (2,))
        self.assertEqual(molgrid.get_correlation_potential_at_points(pts, P).shape, (2,))
        self.assertEqual(len(molgrid.get_becke_weights()), 1)
        self.assertEqual(molgrid.get_becke_weights()[0].shape, (400,))
        self.assertEqual(molgrid.get_densities()[0].shape, (400,))
        self.assertEqual(molgrid.get_gradients()[0].shape, (400, 3))

        with self.assertRaisesRegex(ValueError, 'Nx3'):
            molgrid.get_density_at_points(np.array([0.0, 0.0, 0.0]), P)

        with self.assertRaisesRegex(ValueError, 'Nx3'):
            molgrid.calculate_coulomb_potential_at_points(np.array([[0.0, 0.0]]))

        with self.assertRaisesRegex(ValueError, 'singular'):
            molgrid.calculate_coulomb_potential_at_points(np.array([[0.0, 0.0, 0.0]]))

        she = molgrid.get_spherical_harmonic_expansion_of_amplitude(res['orbc'][:, 0])
        she_radial = molgrid.get_spherical_harmonic_expansion_of_amplitude(
            res['orbc'][:, 0],
            radial_factor=True,
        )
        self.assertEqual(she.shape, (1, 8, 36))
        self.assertEqual(she_radial.shape, she.shape)
        self.assertFalse(np.allclose(she, she_radial))

        self.assertEqual(molgrid.get_rho_lm_atoms().shape, (1, 8, 36))
        self.assertTrue(np.isfinite(molgrid.count_electrons_from_rho_lm()))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_nuclear_attraction_local()))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_nuclear_attraction_full()))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_coulomb()))
        self.assertEqual(molgrid.get_hartree_potential().shape, (400,))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_coulomb_no_interpolation()))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_exchange()))
        self.assertTrue(np.isfinite(molgrid.calculate_dfa_kinetic()))

    def test_gga_point_exchange_and_correlation_potentials(self):
        mol = MoleculeBuilder().from_name('He')
        dft = DFT(
            mol,
            basis='sto3g',
            functional='pbe',
            nshells=8,
            nangpts=50,
            lmax=5,
        )
        res = dft.scf(tol=1e-4)
        molgrid = dft.get_molgrid_copy()
        pts = np.array([
            [0.1, 0.0, 0.0],
            [0.0, 0.2, 0.1],
        ])

        self.assertEqual(molgrid.get_exchange_potential_at_points(pts, res['P']).shape, (2,))
        self.assertEqual(molgrid.get_correlation_potential_at_points(pts, res['P']).shape, (2,))

if __name__ == '__main__':
    unittest.main()
