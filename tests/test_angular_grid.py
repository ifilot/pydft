import unittest

import numpy as np

from pydft.angulargrid import AngularGrid


class TestAngularGrid(unittest.TestCase):

    def test_coefficients_are_cached_and_dataset_sizes_are_reported(self):
        grid = AngularGrid()

        coeffs = grid.get_coefficients(50)
        np.testing.assert_allclose(grid.get_coefficients(50), coeffs)
        self.assertIn(50, grid.get_dataset_sizes())

    def test_invalid_lebedev_order_is_rejected(self):
        grid = AngularGrid()

        with self.assertRaisesRegex(ValueError, "There is no Lebedev order"):
            grid.get_coefficients(51)


if __name__ == '__main__':
    unittest.main()
