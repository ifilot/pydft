import unittest

import numpy as np

from pydft import DFT, MoleculeBuilder


class TestXCMatrixDerivative(unittest.TestCase):

    def test_pbe_xc_matrix_matches_energy_derivative(self):
        mol = MoleculeBuilder().from_name('H2')
        dft = DFT(mol, basis='sto3g', functional='pbe', normalize=False)
        res = dft.scf()
        molgrid = dft.get_molgrid_copy()

        pmat = res['density']
        direction = np.array([[0.25, -0.15],
                              [-0.15, -0.10]])

        molgrid.build_density(pmat, normalize=False)
        xmat, ex = molgrid.calculate_exchange()
        cmat, ec = molgrid.calculate_correlation()
        analytical = np.einsum('ij,ij', xmat + cmat, direction)

        eps = 1e-4
        energies = []
        for sign in [1.0, -1.0]:
            molgrid.build_density(pmat + sign * eps * direction,
                                  normalize=False)
            _xmat, ex = molgrid.calculate_exchange()
            _cmat, ec = molgrid.calculate_correlation()
            energies.append(ex + ec)

        numerical = (energies[0] - energies[1]) / (2 * eps)
        np.testing.assert_allclose(analytical, numerical,
                                   rtol=1e-6, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
