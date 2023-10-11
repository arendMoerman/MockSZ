"""
@file
Test single pointing spectral distortions in MockSZ.
"""

import numpy as np
import unittest
import matplotlib.pyplot as pt
import MockSZ.SinglePointing as SPoint 

class TestSinglePointing(unittest.TestCase):
    def test_getSpecIntensityRM(self):
        Te = 59e6
        tau_e = 1
        num_arr = 10

        mu = np.linspace(1, 1000, num=num_arr) * 1e9
        
        test_out = SPoint.getSpecIntensityRM(mu, Te, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
        
        test_out = SPoint.getSpecIntensityRM(mu, Te, tau_e, log_nu=True)
        self.assertEqual(test_out.shape, (num_arr,))
    
    def test_getSpecIntensityPL(self):
        alpha = 2.5
        tau_e = 1
        num_arr = 10

        mu = np.linspace(1, 1000, num=num_arr) * 1e9
        
        test_out = SPoint.getSpecIntensityPL(mu, alpha, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
        
        alpha = None
        test_out = SPoint.getSpecIntensityPL(mu, alpha, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
        
        test_out = SPoint.getSpecIntensityPL(mu, alpha, tau_e, log_nu=True)
        self.assertEqual(test_out.shape, (num_arr,))

    def test_getSpecIntensityKSZ(self):
        beta_z = 0.1
        tau_e = 1
        num_arr = 10

        nu = np.linspace(1, 1000, num=num_arr) * 1e9
        
        test_out = SPoint.getSpecIntensityKSZ(nu, beta_z, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
if __name__ == "__main__":
    import nose2
    nose2.main()





