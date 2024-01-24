import numpy as np
import unittest
from nose2.tools import params

import MockSZ.Constants as test_ct

class TestConstants(unittest.TestCase):
    @params("h", "k", "c", "me", "eV", "st", "Tcmb")
    def test_Constants(self, attr):
        self.assertTrue(hasattr(test_ct, attr))

if __name__ == "__main__":
    import nose2
    nose2.main()

