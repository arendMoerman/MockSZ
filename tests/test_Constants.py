"""
@file
Test constants in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.Constants as test_ct

class TestConstants(unittest.TestCase):
    def test_Constants(self):
        self.assertTrue(hasattr(test_ct, "h"))
        self.assertTrue(hasattr(test_ct, "k"))
        self.assertTrue(hasattr(test_ct, "c"))
        self.assertTrue(hasattr(test_ct, "me"))

if __name__ == "__main__":
    import nose2
    nose2.main()

