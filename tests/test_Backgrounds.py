"""
@file
Test custom logger of PyPO.
"""

import numpy as np
import unittest
from MockSZ.Backgrounds import CMB

class TestBackgrounds(unittest.TestCase):
    def test_CMB(self):
        test_cmb = CMB()
        self.assertTrue(isinstance(test_cmb, CMB))

        freq = 100e9
        test_out = test_cmb.getSpecificIntensity(freq)

        self.assertTrue(isinstance(test_out, float))
        
        freq = np.linspace(100e9, 200e9, num=100)
        test_out = test_cmb.getSpecificIntensity(freq)

        self.assertTrue(test_out.shape[0] == 100)

if __name__ == "__main__":
    import nose2
    nose2.main()
