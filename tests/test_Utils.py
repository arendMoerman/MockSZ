"""
@file
Test utilities file in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.Utils as MUtils

class TestUtils(unittest.TestCase):
    def test_getBetaFromVelocity(self):
        test_vel = 100 # m / s
        test_out = MUtils.getBetaFromVelocity(test_vel)

        self.assertTrue(isinstance(test_out, float))
    
    def test_getGammaFromVelocity(self):
        test_vel = 100 # m / s
        test_out = MUtils.getGammaFromVelocity(test_vel)

        self.assertTrue(isinstance(test_out, float))

    def test_getGammaFromBeta(self):
        test_beta = 0.5 # m / s
        test_out = MUtils.getGammaFromBeta(test_beta)

        self.assertTrue(isinstance(test_out, float))
if __name__ == "__main__":
    import nose2
    nose2.main()


