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
        num_arr = 1000

        s = np.linspace(-1, 1, num=num_arr)
        
        test_out = MMulti.getP1_RM(s, Te)
        self.assertEqual(test_out.shape, (num_arr,))
    
    def test_P1_PL(self):
        alpha = 2.5
        num_arr = 1000

        s = np.linspace(-1, 10, num=num_arr)
        
        test_out = MMulti.getP1_PL(s, alpha)
        self.assertEqual(test_out.shape, (num_arr,))

if __name__ == "__main__":
    import nose2
    nose2.main()




