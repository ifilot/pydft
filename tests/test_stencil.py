import unittest

import numpy as np

from pydft.stencil import stencil_interp


class TestStencilInterpolation(unittest.TestCase):

    def test_stencil_interp_matches_manual_cubic_weights(self):
        vals = np.array([
            [1.0, 10.0],
            [2.0, 20.0],
            [3.0, 30.0],
            [4.0, 40.0],
            [5.0, 50.0],
        ])
        x = np.array([1, 2])
        weights = np.array([
            [0.1, 0.2, 0.3, 0.4],
            [0.4, 0.3, 0.2, 0.1],
        ])

        expected = np.array([
            [
                0.1 * vals[0, 0] + 0.2 * vals[1, 0] + 0.3 * vals[2, 0] + 0.4 * vals[3, 0],
                0.4 * vals[1, 0] + 0.3 * vals[2, 0] + 0.2 * vals[3, 0] + 0.1 * vals[4, 0],
            ],
            [
                0.1 * vals[0, 1] + 0.2 * vals[1, 1] + 0.3 * vals[2, 1] + 0.4 * vals[3, 1],
                0.4 * vals[1, 1] + 0.3 * vals[2, 1] + 0.2 * vals[3, 1] + 0.1 * vals[4, 1],
            ],
        ])

        out = np.zeros((2, 2))
        stencil_interp.py_func(vals, x, weights, out)
        np.testing.assert_allclose(out, expected)

        compiled_out = np.zeros((2, 2))
        stencil_interp(vals, x, weights, compiled_out)
        np.testing.assert_allclose(compiled_out, expected)


if __name__ == '__main__':
    unittest.main()
