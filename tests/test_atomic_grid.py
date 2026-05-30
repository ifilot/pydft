import unittest
import numpy as np
from pyqint import PyQInt, Molecule
from pydft import AtomicGrid
from pydft.atomicgrid import job_build_atomic_grid


class TestAtomicGrid(unittest.TestCase):

    def test_hartee_potential(self):
        """
        Test calculation of the Hartree potential and the electronic
        self-repulsion using a regular grid
        """
        # grab a non-trivial sampling AO
        mol = Molecule('Fe')
        mol.add_atom('Fe', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')
        integrator = PyQInt()
        ao = cgfs[-7]

        # build an atomic grid
        ag = AtomicGrid(('H', [0,0,0]))
        pts = ag.get_gridpoints()
        amps = integrator.plot_wavefunction(pts, [1], [ao])

        # count total number of electrons
        ag.set_density(np.power(amps,2))
        np.testing.assert_almost_equal(ag.count_electrons(), 1.0, 4)

        # evaluate Hatree potential and calculate E-E repulsion
        ag.build_hartree_potential()
        E = ag.calculate_coulomb_energy()
        exact = integrator.repulsion(ao,ao,ao,ao)
        np.testing.assert_almost_equal(E, exact, 4)
    
    def test_hartee_potential_coarse_grid(self):
        """
        Test calculation of the Hartree potential and the electronic
        self-repulsion using a somewhat coarser grid
        """
        # grab a non-trivial sampling AO
        mol = Molecule('Fe')
        mol.add_atom('Fe', 0.0, 0.0, 0.0)
        cgfs, nuclei = mol.build_basis('sto3g')
        integrator = PyQInt()
        ao = cgfs[-7]

        # build an atomic grid
        ag = AtomicGrid(('H', [0,0,0]), nshells=20, nangpts=50, lmax=5)
        pts = ag.get_gridpoints()
        amps = integrator.plot_wavefunction(pts, [1], [ao])

        # count total number of electrons
        ag.set_density(np.power(amps,2))
        np.testing.assert_almost_equal(ag.count_electrons(), 1.0, 3)

        # evaluate Hatree potential and calculate E-E repulsion
        ag.build_hartree_potential()
        E = ag.calculate_coulomb_energy()
        exact = integrator.repulsion(ao,ao,ao,ao)
        np.testing.assert_almost_equal(E, exact, 3)

    def test_constructor_rejects_invalid_stencil_sizes(self):
        with self.assertRaisesRegex(ValueError, "at least be 3"):
            AtomicGrid(('H', [0, 0, 0]), fdpts=1)

        with self.assertRaisesRegex(ValueError, "Only odd"):
            AtomicGrid(('H', [0, 0, 0]), fdpts=4)

        ag = AtomicGrid(('H', [0, 0, 0]), nshells=3, nangpts=50, lmax=5)
        with self.assertRaisesRegex(ValueError, "Need at least 4"):
            ag.build_cubic_stencil_evaluation(np.array([0.1]))

    def test_atomic_grid_helpers_and_local_energy_terms(self):
        ag = AtomicGrid(('H', np.array([0.0, 0.0, 0.0])), nshells=8, nangpts=50, lmax=5)
        density = np.exp(-np.linalg.norm(ag.get_gridpoints(), axis=1))
        gradient = np.zeros((ag.get_nr_pts(), 3))

        ag.set_density(density)
        ag.set_gradient(gradient)

        self.assertEqual(ag.get_charge(), 1)
        self.assertEqual(ag.get_lmax(), 5)
        self.assertEqual(ag.get_nr_pts(), 400)
        self.assertEqual(ag.get_full_grid().shape, (8, 50, 3))
        self.assertEqual(ag.get_radial_grid().shape, (8,))
        self.assertEqual(ag.get_weights().shape, (400,))
        self.assertEqual(ag.get_weights_angular_points().shape, (50,))
        self.assertEqual(ag.get_becke_weights().shape, (8, 50))
        np.testing.assert_allclose(ag.get_gradient_squared(), np.zeros(400))
        self.assertGreater(ag.count_electrons(), 0.0)

        self.assertTrue(np.isfinite(ag.get_dfa_exchange()))
        self.assertTrue(np.isfinite(ag.get_dfa_kinetic()))
        self.assertTrue(np.isfinite(ag.get_dfa_nuclear_local()))
        self.assertTrue(np.isfinite(ag.calculate_nuclear_attraction()))

        ag.build_hartree_potential()
        self.assertEqual(ag.get_rho_lm().shape, (8, 36))
        self.assertEqual(ag.get_ulm().shape, (8, 36))
        self.assertEqual(ag.get_ylm().shape, (36, 50))
        self.assertEqual(ag.get_ylm_atom().shape, (36, 50))

        interpolated = ag.calculate_interpolated_ulm(ag.get_radial_grid()[:3])
        self.assertEqual(interpolated.shape, (36, 3))
        np.testing.assert_allclose(
            ag.calculate_interpolated_ulm_at_points(ag.get_radial_grid()[:3]),
            interpolated,
        )

        ag.build_cubic_stencil_evaluation(ag.get_radial_grid()[:3])
        self.assertEqual(ag.calculate_interpolated_ulm_stencil().shape, (36, 3))

        self.assertTrue(np.isfinite(ag.calculate_coulomb_energy()))
        self.assertEqual(ag.get_local_hartree_potential().shape, (400,))
        self.assertTrue(np.isfinite(ag.calculate_coulomb_energy_interpolation()))

    def test_job_build_atomic_grid(self):
        ag = job_build_atomic_grid((('H', np.array([0.0, 0.0, 0.0])), 8, 50, 5, 7))

        self.assertIsInstance(ag, AtomicGrid)
        self.assertEqual(ag.get_nr_pts(), 400)

if __name__ == '__main__':
    unittest.main()
