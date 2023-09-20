"""
@file
Test multi-electron statistics in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.MultiStats as MMulti

class TestMultiStats(unittest.TestCase):
    def test_P1_RM(self):
        Te = 10e6
        num_arr = 10

        s = np.linspace(-1, 1, num=num_arr)
        
        test_out = MMulti.getP1_RM(s, Te)
        self.assertEqual(test_out.shape, (num_arr,))
    
    def test_MJ_parallel(self):
        ns = 10
        nb = 5

        beta = np.arange(nb) 
        beta_lim = np.ones(ns) * 0.01
        num_mu = 10
        s = np.linspace(-1, 1, ns)
        dbeta = 0.0001
        Te = 100

        out = MMulti._MJ_parallel(beta, beta_lim, num_mu, s, dbeta, Te)

        self.assertEqual(ns, out.size)

    def test_P1_PL(self):
        alpha = 2.5
        num_arr = 10

        s = np.linspace(-1, 10, num=num_arr)
        
        test_out = MMulti.getP1_PL(s, alpha)
        self.assertEqual(test_out.shape, (num_arr,))
    
    def test_PL_parallel(self):
        ns = 10
        nb = 5

        beta = np.arange(nb) 
        beta_lim = np.ones(ns) * 0.01
        num_mu = 10
        s = np.linspace(-1, 1, ns)
        dbeta = 0.0001
        alpha = 2

        out = MMulti._PL_parallel(beta, beta_lim, num_mu, s, dbeta, alpha)

        self.assertEqual(ns, out.size)

if __name__ == "__main__":
    import nose2
    nose2.main()




