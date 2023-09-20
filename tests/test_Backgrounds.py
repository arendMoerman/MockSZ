"""
@file
Test custom logger of PyPO.
"""

import numpy as np
import unittest
import MockSZ.Backgrounds as MBack

class TestBackgrounds(unittest.TestCase):
    def test_CMB(self):
        freq = 100e9
        test_out = MBack.getSpecificIntensityCMB(freq)

        self.assertTrue(isinstance(test_out, float))
        
        freq = np.linspace(100e9, 200e9, num=100)
        test_out = MBack.getSpecificIntensityCMB(freq)

        self.assertTrue(test_out.shape[0] == 100)

if __name__ == "__main__":
    import nose2
    nose2.main()
