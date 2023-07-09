"""
@file
Test constants in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.Conversions as MConv
from MockSZ.Backgrounds import CMB

class TestConversions(unittest.TestCase):
    def test_eV_Temp(self):
        self.assertAlmostEqual(MConv.eV_Temp(1), 1.1604518e4, delta=1e-3)
    
    def test_SI_JySr(self):
        self.assertAlmostEqual(MConv.SI_JySr(1)*1e-26, 1, delta=1e-3)

    def test_SI_Temp(self):
        nu = np.linspace(1, 100)
        cmb = CMB()
        I0 = cmb.getSpecificIntensity(nu)
        self.assertAlmostEqual(np.mean(MConv.SI_Temp(I0, nu)), cmb.T ,delta=1e-3)

    def test_pc_m(self):
        self.assertAlmostEqual(MConv.pc_m(1)*1e-16, 3.0857, delta=1e-3)

if __name__ == "__main__":
    import nose2
    nose2.main()


