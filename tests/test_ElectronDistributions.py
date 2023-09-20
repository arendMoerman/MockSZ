"""
@file
Test electron distributions.
"""

import numpy as np
import unittest
import MockSZ.ElectronDistributions as EDist

class TestBackgrounds(unittest.TestCase):
    def test_relativisticMaxwellian(self):
        beta = 0.1
        Te = 1000
        test_out = EDist.relativisticMaxwellian(beta, Te)

        self.assertTrue(isinstance(test_out, float))
    
    def test_relativisticPowerLaw(self):
        beta = 0.1
        alpha = 2
        test_out = EDist.relativisticPowerlaw(beta, alpha=alpha)

        self.assertTrue(isinstance(test_out, float))
        
        test_out = EDist.relativisticPowerlaw(beta)

        self.assertTrue(isinstance(test_out, float))

if __name__ == "__main__":
    import nose2
    nose2.main()

