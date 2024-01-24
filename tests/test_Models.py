import numpy as np
import unittest
from nose2.tools import params

import MockSZ.Models as test_md

class TestModels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nAz = 10
        cls.Az = np.linspace(-10, 10, cls.nAz)
        cls.nEl = 8
        cls.El = np.linspace(-10, 10, cls.nEl)
        cls.ibeta = 0.7
        cls.ne0 = 1.2e-2
        cls.thetac = 15
        cls.Da = 1500

        cls.Te = 15.33
        cls.alpha = 2.5
        cls.v_pec = 100

        cls.n_test = 1000
        cls.nu_arr = np.linspace(100, 500, cls.n_test)

        cls.s = 1
        cls.s_arr = np.linspace(-1, 1, cls.n_test)
        
        cls.beta = 0.1
        cls.beta_arr = np.linspace(0, 0.9, cls.n_test)

    @params(True, False)
    def test_IsoBetaModel(self, CMB):
        isobObj = test_md.IsoBetaModel()
        isob_grid = isobObj.getIsoBeta(self.Az, self.El, self.ibeta, 
                                      self.ne0, self.thetac, self.Da,
                                      grid=True)

        self.assertEqual(isob_grid.shape[0], self.nAz)
        self.assertEqual(isob_grid.shape[1], self.nEl)
        
        isob_trace = isobObj.getIsoBeta(self.Az, self.Az, self.ibeta, 
                                      self.ne0, self.thetac, self.Da,
                                      grid=False)
        
        self.assertEqual(isob_trace.size, self.nAz)

        isob_cube_grid = isobObj.getIsoBetaCube(isob_grid, self.nu_arr, self.Te, self.v_pec, CMB)

        self.assertEqual(isob_cube_grid.shape[0], self.nAz)
        self.assertEqual(isob_cube_grid.shape[1], self.nEl)
        self.assertEqual(isob_cube_grid.shape[2], self.n_test)
        
        isob_cube_trace = isobObj.getIsoBetaCube(isob_trace, self.nu_arr, self.Te, self.v_pec, CMB)

        self.assertEqual(isob_cube_trace.shape[0], self.nAz)
        self.assertEqual(isob_cube_trace.shape[1], self.n_test)
    
    @params(True, False)
    def test_SinglePointing(self, CMB):
        spObj = test_md.SinglePointing()

        tSZ = spObj.getSingleSignal_tSZ(self.nu_arr, self.Te, no_CMB=CMB)
        self.assertEqual(tSZ.shape, self.nu_arr.shape)
        
        ntSZ = spObj.getSingleSignal_ntSZ(self.nu_arr, self.alpha, no_CMB=CMB)
        self.assertEqual(ntSZ.shape, self.nu_arr.shape)
        
        kSZ = spObj.getSingleSignal_kSZ(self.nu_arr, self.v_pec)
        self.assertEqual(kSZ.shape, self.nu_arr.shape)

    def test_ScatteringKernels(self):
        skObj = test_md.ScatteringKernels()

        sscatter = skObj.getSingleScattering(self.s_arr, self.beta)
        self.assertEqual(sscatter.shape, self.s_arr.shape)
        
        mj = skObj.getMaxwellJuttner(self.beta_arr, self.Te)
        self.assertEqual(mj.shape, self.beta_arr.shape)
        
        pl = skObj.getPowerlaw(self.beta_arr, self.alpha)
        self.assertEqual(pl.shape, self.beta_arr.shape)
        
        mjscatter = skObj.getMultiScatteringMJ(self.s_arr, self.Te)
        self.assertEqual(mjscatter.shape, self.s_arr.shape)
        
        plscatter = skObj.getMultiScatteringPL(self.s_arr, self.alpha)
        self.assertEqual(plscatter.shape, self.s_arr.shape)

if __name__ == "__main__":
    import nose2
    nose2.main()



