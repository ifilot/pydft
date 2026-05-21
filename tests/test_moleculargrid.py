import unittest
import numpy as np

from pyqint import PyQInt, Molecule
from pydft import MoleculeBuilder, MolecularGrid

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

if __name__ == '__main__':
    unittest.main()