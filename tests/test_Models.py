"""
@file
Test custom logger of PyPO.
"""

import numpy as np
import unittest
from MockSZ.Models import IsoBetaModel
import MockSZ.Conversions as MConv
import matplotlib.pyplot as pt

import MockSZ.SinglePointing as MSingle

class TestModels(unittest.TestCase):
    def test_IsoBetaModel(self):
        Te = 56e6
        ne0 = 0.1
        rc = MConv.pc_m(3e6)
        Da = MConv.pc_m(99e6)

        beta = 0.7

        isob = IsoBetaModel(Te, ne0, rc, beta, Da)

        Az, El = np.mgrid[-1:1:100*1j, -1:1:100*1j]

        theta = np.sqrt(Az**2 + El**2)

        tau = isob.opticalDepths(theta)
        self.assertEqual(tau.shape, (100, 100))
        
        nu = np.linspace(1, 1000, num=1000)
        tSZ = isob.tSZMap(theta, nu)
        self.assertEqual(tSZ.shape, (100, 100, 1000))
        
        kSZ = isob.kSZMap(theta, nu)
        self.assertEqual(kSZ.shape, (100, 100, 1000))

if __name__ == "__main__":
    import nose2
    nose2.main()

