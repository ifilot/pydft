import unittest
import sys
import os
import numpy as np

# add a reference to load the pyDFT module
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pydft

class TestVersion(unittest.TestCase):

    def test_version(self):
        strver = pydft.__version__.split('.')
        self.assertTrue(strver[0].isnumeric())
        self.assertTrue(strver[1].isnumeric())
        self.assertTrue(strver[2].isnumeric())

if __name__ == '__main__':
    unittest.main()
