"""
@file
Test single-electron statistics in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.SingleStats as MSingle

class TestSingleStats(unittest.TestCase):
    def test_getPsbThomson(self):
        num_arr = 10
        s = np.linspace(-1, 1, num=num_arr)
        beta = 0.5
        
        test_out = MSingle.getPsbThomson(s, beta, num_arr)
        self.assertEqual(test_out.shape, (num_arr, 1))

        s = np.linspace(-1, 1, num=num_arr)
        beta = np.linspace(0.1, 0.5, num=num_arr * 2)
        
        test_out = MSingle.getPsbThomson(s, beta, num_arr)
        self.assertEqual(test_out.shape, (num_arr, num_arr*2))
        
        beta = np.linspace(0.01, 0.5, num=num_arr)
        s = 0.
        
        test_out = MSingle.getPsbThomson(s, beta, num_arr)
        self.assertEqual(test_out.shape, (1, num_arr))
        
        beta = 0.5 
        s = 0.
        
        test_out = MSingle.getPsbThomson(s, beta, num_arr)
        self.assertEqual(test_out.shape, (1,))

if __name__ == "__main__":
    import nose2
    nose2.main()



