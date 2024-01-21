"""!
@file
Cluster models for making SZ maps.
"""

import numpy as np

import MockSZ.Bindings as MBind

class IsoBetaModel(object):
    """!
    Class representing an isothermal-beta model.
    Should be instantiated and serves as an interface for MockSZ when simulating these types of clusters.

    @ingroup clustermodels
    """
    def getIsoBeta(self, Az, El, ibeta, ne0, thetac, Da, grid=False):
        res = MBind.getIsoBeta(Az, El, ibeta, ne0, thetac, Da, grid)
        return res
    
    def getIsoBetaCube(self, isobeta, nu_arr, Te, v_pec, n_s=500, n_beta=500, n_mu=500, no_CMB=False):
        res_tSZ = MBind.getSinglePointing_tSZ(nu_arr, Te, n_s, n_beta, no_CMB)
        res_kSZ = MBind.getSinglePointing_kSZ(nu_arr, v_pec, n_mu)
        
        res_SZ = res_tSZ + res_kSZ

        if len(isobeta.shape) == 1:
            shape = (isobeta.size, nu_arr.size)
            res = np.ones(shape)
            for i in range(nu_arr.size):
                res[:,i] = isobeta
        else:
            shape = (isobeta.shape[0], isobeta.shape[1], nu_arr.size)
            res = np.ones(shape)
            for i in range(nu_arr.size):
                res[:,:,i] = isobeta

        res *= res_SZ
        
        return res

class SinglePointing(object):
    """! 
    Class for generating a single pointing SZ signal.

    """
    def getSingleSignal_tSZ(self, nu_arr, Te, n_s=500, n_beta=500, no_CMB=False):
        res = MBind.getSinglePointing_tSZ(nu_arr, Te, n_s, n_beta, no_CMB)
        return res
    
    def getSingleSignal_ntSZ(self, nu_arr, alpha, n_s=500, n_beta=500, no_CMB=False):
        res = MBind.getSinglePointing_ntSZ(nu_arr, alpha, n_s, n_beta, no_CMB)

        return res
    
    def getSingleSignal_kSZ(self, nu_arr, v_pec, n_mu=500):
        res = MBind.getSinglePointing_kSZ(nu_arr, v_pec, n_mu)

        return res

class ScatteringKernels(object):
    """!
    Class for investigating the scattering kernels.
    It is not to be used per se for simulations, but it provides a nice overview of available distributions.
    In addition, it is helpful for troubleshooting and debugging to have access to the scattering kernels directly.
    So, it is mostly informative, and not meant to be used during actual simulations.
    """

    def getSingleScattering(self, s_arr, beta, num_mu=1000):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getThomsonScatter(s_arr, beta, num_mu)

        return res
    
    def getMaxwellJuttner(self, beta_arr, Te):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getMaxwellJuttner(beta_arr, Te)

        return res
    
    def getPowerlaw(self, beta_arr, alpha):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getPowerlaw(beta_arr, alpha)

        return res
    
    def getMultiScatteringMJ(self, s_arr, Te, n_beta=500):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getMultiScatteringMJ(s_arr, Te, n_beta)

        return res
    
    def getMultiScatteringPL(self, s_arr, alpha, n_beta=500):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getMultiScatteringPL(s_arr, alpha, n_beta)

        return res
