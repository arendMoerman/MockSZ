"""!
@file
Cluster/single-pointing models for making SZ maps.
Additionally, this file contains a class that is useful for investigating scattering kernels and electron distributions.
"""

import numpy as np
from time import time

import MockSZ.Bindings as MBind
import MockSZ.Constants as MConst
import MockSZ.Conversions as MConv

def timer_func(func): 
    """!
    Decorator for timing functions in MockSZ.
    Useful for estimating performance.
    """

    def wrap_func(*args, **kwargs): 
        t1 = time() 
        result = func(*args, **kwargs) 
        t2 = time()
        if kwargs.get("timer") is not None and kwargs.get("timer") == True:
            return result, t2-t1
            #print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s') 
        return result 
    return wrap_func


class SinglePointing(object):
    """! 
    Class for generating a single pointing SZ signal. Can choose between tSZ, kSZ and ntSZ (powerlaw).

    Attributes:
        clib Library containing backend functions.
    
    @ingroup singlepointing
    """

    def __init__(self, param, v_pec=None, phi_cl=0, tau_e=1, no_CMB=False):
        """!
        Initialise a single-pointing model of a galaxy cluster.

        Be careful: if param is an electron temperature, using the ntSZ option leads to nonsensical results.
        Similarly, setting param to a powerlaw alpha and running the tSZ routine leads to nonsense.

        @param param Parameter (Te or alpha) governing cluster properties.
            Te should be given in keV, alpha is dimensionless.
        @param v_pec Absolute peculiar velocity of cluster, in km / s.
            Defaults to None, which puts the cluster at rest w.r.t. the CMB frame.
        @param phi_cl Angle between line-of-sight (pointing away from observer) and cluster velocity, in degrees.
            Defaults to 0 degrees, i.e., the cluster is receding from us.
        @param tau_e Electron optical depth of cluster.
            Defaults to None, which puts the optical depth at 1.
        @param no_CMB Whether or not to add the CMB signal to the distortions.
            Defaults to False (CMB added on top).
        """

        self.param = param
        self.v_pec = v_pec
        self.tau_e = tau_e
        self.no_CMB = no_CMB

        if v_pec is not None:
            self.beta_cl = self.v_pec * 1e3 / MConst.c
            self.beta_cl_z = self.beta_cl * np.cos(np.radians(phi_cl))

        else:
            self.beta_cl = None
            self.beta_cl_z = None

        self.clib = MBind.loadMockSZlib()
    
    @timer_func
    def getSingleSignal_tSZ(self, nu_arr, timer=False, acc=1e-6):
        """!
        Generate a single pointing signal of the tSZ effect.

        @param nu_arr Numpy array of frequencies for tSZ effect, in Hz.
        @param timer Time function execution. Used in decorator.
        @param acc Required relative accuracy of integration.
        
        @returns res 1D array containing tSZ effect.
        """

        res = MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, self.no_CMB, acc, 
                                            func=self.clib.MockSZ_getSignal_tSZ)

        if self.v_pec is not None:
            res += MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, self.no_CMB, self.beta_cl, 
                                            func=self.clib.MockSZ_getSignal_tSZ_beta2)

        return res
    
    @timer_func
    def getSingleSignal_ntSZ(self, nu_arr, timer=False, acc=1e-6):
        """!
        Generate a single pointing signal of the ntSZ effect, according to a powerlaw.

        @param nu_arr Numpy array of frequencies for ntSZ effect, in Hz.
        @param timer Time function execution. Used in decorator.
        @param acc Required relative accuracy of integration.
        
        @returns res 1D array containing ntSZ effect.
        """

        res = MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, self.no_CMB, acc, 
                                            func=self.clib.MockSZ_getSignal_ntSZ)

        return res
    
    @timer_func
    def getSingleSignal_kSZ(self, nu_arr, timer=False, acc=1e-6):
        """!
        Generate a single pointing signal of the kSZ effect.

        The CMB can be added by adding the kSZ signal to a tSZ or ntSZ signal containing the CMB.

        @param nu_arr Numpy array of frequencies for kSZ effect, in Hz.
        @param timer Time function execution. Used in decorator.
        @param acc Required relative accuracy of integration.
        
        @returns res 1D array containing kSZ effect.
        """

        res = MBind.getDistributionTwoParam(nu_arr, self.beta_cl_z, self.tau_e, self.no_CMB, acc, 
                                            func=self.clib.MockSZ_getSignal_kSZ)
        if self.param is not None:
            res += MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, self.no_CMB, self.beta_cl_z, 
                                            func=self.clib.MockSZ_getSignal_kSZ_betatheta)
            
            b2t_fac = self.beta_cl_z**2 * 3/2 - self.beta_cl**2/2

            res += MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, self.no_CMB, b2t_fac, 
                                            func=self.clib.MockSZ_getSignal_kSZ_betat2heta)

        return res

class IsoBetaModel(SinglePointing):
    """!
    Class representing an isothermal-beta model.
    Serves as an interface for MockSZ when simulating these types of clusters.

    Attributes:
        clib Library containing backend functions.

    @ingroup clustermodels
    """
    
    def __init__(self, param, v_pec=None, phi_cl=0, no_CMB=False):
        """!
        Initialise a single-pointing model of a galaxy cluster.
        Under the hood, calls the constructor of a single-pointing class.

        @param param Parameter (Te or alpha) governing cluster properties.
            Te should be given in keV, alpha is dimensionless.
        @param v_pec Absolute peculiar velocity of cluster, in km / s.
            Defaults to None, which puts the cluster at rest w.r.t. the CMB frame.
        @param phi_cl Angle between line-of-sight (pointing away from observer) and cluster velocity, in degrees.
            Defaults to 0 degrees, i.e., the cluster is receding from us.
        @param no_CMB Whether or not to add the CMB signal to the distortions.
            Defaults to False (CMB added on top).
        """
        
        super().__init__(param, v_pec, phi_cl, 1, True)
        self.no_CMB_cl = no_CMB
    
    def getIsoBeta(self, Az, El, ibeta, ne0, thetac, Da, grid=False):
        """!
        Get an isothermal-beta optical depth screen. 

        @param Az Numpy array containing the range of Azimuth co-ordinates, in arcseconds.
        @param El Numpy array containing the range of Elevation co-ordinates, in arcseconds.
        @param ibeta Beta parameter for isothermal model.
        @param ne0 Central electron number density, in number / cm**3.
        @param thetac Angular cluster core radius, in arcseconds.
        @param Da Angular diameter distance to cluster, in Megaparsec.
        @param grid Whether or not to evaluate the model on a 2D grid spanned by Az and El, or on a 1D trace.
            If grid=True, the screen will be of size Az.size * El.size.
            If grid=False (default), it is required that Az.size ==  El.size, and the screen will be of size Az.size = El.size.

        @returns res The optical depth screen.
        """

        res = MBind.getIsoBeta(Az, El, ibeta, ne0, thetac, Da, grid)
        return res
    
    def getIsoBetaCube(self, isobeta, nu_arr, acc=1e-6):
        """!
        Get an isothermal-beta model from an optical depth screen.

        @param isobeta An optical depth screen generated by self.getIsoBeta.
        @param nu_arr Numpy array of frequencies for SZ effect, in Hz.
        @param acc Required relative accuracy of integration.
        
        @returns res 2D or 3D grid (depending on dimensions of isobeta) containing SZ signal attenuated by optical depth in isobeta.
        """

        res_kSZ = self.getSingleSignal_kSZ(nu_arr, acc=acc)

        res_tSZ = self.getSingleSignal_kSZ(nu_arr, acc=acc)
        
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
        
        if self.no_CMB_cl == False:
            res += MBind.getCMB(nu_arr)
        
        return res

class ScatteringKernels(object):
    """!
    Class for investigating the scattering kernels.
    It is not to be used per se for simulations, but it provides a nice overview of available distributions.
    In addition, it is helpful for troubleshooting and debugging to have access to the scattering kernels directly.
    So, it is mostly informative, and not meant to be used during actual simulations.

    Attributes:
        clib Library containing backend functions.

    @ingroup scatteringkernels
    """
    
    def __init__(self):
        self.clib = MBind.loadMockSZlib()

    def getSingleScattering(self, s_arr, beta, acc=1e-6):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        This method assumes a Thomson scattering in the electron rest frame.

        @param s_arr Numpy array of logarithmic frequency shifts s.
        @param beta Dimensionless electron velocity.
        @param acc Required relative accuracy of integration.

        @returns res 1D array containing single-electron scattering probabilities.
        """
        
        res = MBind.getDistributionSingleParam(s_arr, beta, acc, func=self.clib.MockSZ_getThomsonScatter)

        return res
    
    def getMaxwellJuttner(self, beta_arr, Te):
        """!
        Obtain a Maxwell-Juttner distribution.

        @param beta_arr Numpy array of dimensionless electron velocities.
        @param Te Electron temperature in keV.

        @returns res 1D array containing Maxwell-Juttner distribution.
        """
        
        res = MBind.getDistributionSingleParam(beta_arr, Te, acc=1e-6, 
                                               func=self.clib.MockSZ_getMaxwellJuttner)

        return res
    
    def getPowerlaw(self, beta_arr, alpha):
        """!
        Obtain a relativistic powerlaw distribution.

        @param beta_arr Numpy array of dimensionless electron velocities.
        @param alpha Slope of powerlaw.

        @returns res 1D array containing relativistic powerlaw distribution.
        """
        
        res = MBind.getDistributionSingleParam(beta_arr, alpha, acc=1e-6, 
                                               func=self.clib.MockSZ_getPowerlaw)

        return res
    
    def getMultiScatteringMJ(self, s_arr, Te, acc=1e-6):
        """!
        Obtain multi-electron scattering kernel, for a range of beta.
        This kernel is calculated using a Maxwell-Juttner distribution.

        @param s_arr Numpy array of logarithmic frequency shifts s.
        @param Te Electron temperature in keV.
        @param acc Required relative accuracy of integration.

        @returns res 1D array containing mulit-electron scattering probabilities.
        """
        
        res = MBind.getDistributionSingleParam(s_arr, Te, acc, 
                                               func=self.clib.MockSZ_getMultiScatteringMJ)

        return res
    
    def getMultiScatteringPL(self, s_arr, alpha, acc=1e-6):
        """!
        Obtain multi-electron scattering kernel, for a range of beta.
        This kernel is calculated using a relativistic powerlaw distribution.

        @param s_arr Numpy array of logarithmic frequency shifts s.
        @param alpha Slope of powerlaw.
        @param acc Required relative accuracy of integration.

        @returns res 1D array containing mulit-electron scattering probabilities.
        """
        
        res = MBind.getDistributionSingleParam(s_arr, alpha, acc, 
                                               func=self.clib.MockSZ_getMultiScatteringPL)

        return res
