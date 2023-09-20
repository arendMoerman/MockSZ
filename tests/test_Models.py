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

        nAz = 10
        nEl = 10

        nfreq = 10

        Az, El = np.mgrid[-1:1:nAz*1j, -1:1:nEl*1j]

        theta = np.sqrt(Az**2 + El**2)

        tau = isob.opticalDepths(theta)
        self.assertEqual(tau.shape, (nAz, nEl))
        
        nu = np.linspace(1, 1000, num=nfreq)
        tSZ = isob.tSZMap(theta, nu)
        self.assertEqual(tSZ.shape, (nAz, nEl, nfreq))
        
        kSZ = isob.kSZMap(theta, nu)
        self.assertEqual(kSZ.shape, (nAz, nEl, nfreq))

if __name__ == "__main__":
    import nose2
    nose2.main()

