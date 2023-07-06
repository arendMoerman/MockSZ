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
        num_arr = 1000

        mu = np.linspace(1, 1000, num=num_arr) * 1e9
        
        test_out = SPoint.getSpecIntensityRM(mu, Te, tau_e)
        test_out2 = SPoint.getSpecIntensityRM(mu, 3*Te, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
    
    def test_getSpecIntensityPL(self):
        alpha = 2.5
        tau_e = 1
        num_arr = 1000

        mu = np.linspace(1, 1000, num=num_arr) * 1e9
        
        test_out = SPoint.getSpecIntensityPL(mu, alpha, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))
        
        alpha = None
        test_out = SPoint.getSpecIntensityPL(mu, alpha, tau_e)
        self.assertEqual(test_out.shape, (num_arr,))

if __name__ == "__main__":
    import nose2
    nose2.main()





