import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from pydft import MoleculeBuilder, MolecularGrid

class TestMolecularGrid(unittest.TestCase):
    """
    Perform a quick test to verify whether all DFT routines
    are properly working
    """

    def test_weight_functions(self):
        """
        Test DFT calculation of Helium atom
        """
        # construct molecule
        mol = MoleculeBuilder().from_name('bf3')
        cgfs, atoms = mol.build_basis('sto3g')
        
        # construct molecular grid
        molgrid = MolecularGrid(atoms, cgfs)
        
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