import unittest
import pydft

class TestVersion(unittest.TestCase):

    def test_version(self):
        strver = pydft.__version__.split('.')
        self.assertTrue(strver[0].isnumeric())
        self.assertTrue(strver[1].isnumeric())
        self.assertTrue(strver[2].isnumeric())

if __name__ == '__main__':
    unittest.main()
