"""
@file
Test multi-electron statistics in MockSZ.
"""

import numpy as np
import unittest
import MockSZ.MultiStats as MMulti

import matplotlib.pyplot as pt

class TestMultiStats(unittest.TestCase):
    def test_P1_RM(self):
        Te = 59e6
        num_arr = 1000

        s = np.linspace(-1, 1, num=num_arr)
        
        test_out = MMulti.getP1_RM(s, Te)
        self.assertEqual(test_out.shape, (num_arr,))

        #print(test_out)
        pt.plot(s, test_out)
        pt.show()

if __name__ == "__main__":
    import nose2
    nose2.main()




