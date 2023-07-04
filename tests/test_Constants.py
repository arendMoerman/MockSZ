"""
@file
Test constants in MockSZ.
"""

import numpy as np
import unittest
from MockSZ.Constants import Constants

class TestConstants(unittest.TestCase):
    def test_Constants(self):
        test_ct = Constants()
        self.assertTrue(hasattr(test_ct, "h"))
        self.assertTrue(hasattr(test_ct, "k"))
        self.assertTrue(hasattr(test_ct, "c"))
        self.assertTrue(hasattr(test_ct, "me"))

if __name__ == "__main__":
    import nose2
    nose2.main()

