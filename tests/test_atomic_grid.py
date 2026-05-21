import unittest
import numpy as np
from pyqint import PyQInt, Molecule
from pydft import AtomicGrid


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

if __name__ == '__main__':
    unittest.main()