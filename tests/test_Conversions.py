"""
@file
Test constants in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.Conversions as MConv
import MockSZ.Backgrounds as MBack
import MockSZ.Constants as ct

class TestConversions(unittest.TestCase):
    def test_eV_Temp(self):
        self.assertAlmostEqual(MConv.eV_Temp(1), 1.1604518e4, delta=1e-3)
    
    def test_SI_JySr(self):
        self.assertAlmostEqual(MConv.SI_JySr(1)*1e-26, 1, delta=1e-3)

    def test_SI_Temp(self):
        nu = np.linspace(1, 100)
        I0 = MBack.getSpecificIntensityCMB(nu)
        self.assertAlmostEqual(np.mean(MConv.SI_Temp(I0, nu)), ct.Tcmb ,delta=1e-3)

    def test_pc_m(self):
        self.assertAlmostEqual(MConv.pc_m(1)*1e-16, 3.0857, delta=1e-3)
    
    def test_v_beta(self):
        test_vel = 100 # m / s
        test_out = MConv.v_beta(test_vel)

        self.assertTrue(isinstance(test_out, float))
    
    def test_v_gamma(self):
        test_vel = 100 # m / s
        test_out = MConv.v_gamma(test_vel)

        self.assertTrue(isinstance(test_out, float))

    def test_beta_gamma(self):
        test_beta = 0.5 # m / s
        test_out = MConv.beta_gamma(test_beta)

        self.assertTrue(isinstance(test_out, float))

    def test_freq_x(self):
        freq = 100
        test_out = MConv.freq_x(freq)
        self.assertTrue(isinstance(test_out, float))
    
    def test_x_freq(self):
        x = 1
        test_out = MConv.x_freq(x)
        self.assertTrue(isinstance(test_out, float))
    
    def test_Te_theta(self):
        Te = 100
        test_out = MConv.Te_theta(Te)
        self.assertTrue(isinstance(test_out, float))


if __name__ == "__main__":
    import nose2
    nose2.main()


