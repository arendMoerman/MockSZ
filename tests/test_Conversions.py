import numpy as np
import unittest
from nose2.tools import params

import MockSZ.Conversions as test_cv

class TestConversions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_scal = 10
        cls.n_test = 1000
        cls.test_arr = np.linspace(5, 15, cls.n_test)
    
    @params(test_cv.keV_theta, test_cv.keV_Temp, test_cv.SI_JySr,
            test_cv.freq_x, test_cv.x_freq)
    def test_Conversions(self, func):
        self.assertTrue(type(func(self.test_scal)), float)
        
        out_arr = test_cv.keV_theta(self.test_arr)
        self.assertEqual(out_arr.size, self.n_test)
    
    def test_SI_Temp(self):
        self.assertTrue(type(test_cv.SI_Temp(self.test_scal, self.test_scal)), float)
        
        out_arr = test_cv.SI_Temp(self.test_arr, self.test_arr)
        self.assertEqual(out_arr.size, self.n_test)

if __name__ == "__main__":
    import nose2
    nose2.main()


