"""!
@file
Cluster/single-pointing models for making SZ maps.
Additionally, this file contains a class that is useful for investigating scattering kernels and electron distributions.
"""

# STL
from time import time
from typing import Callable, Optional, Sequence

# External packages
import numpy as np
import scipy.constants as const

# MockSZ-specifics
import MockSZ.Bindings as MBind
import MockSZ.Conversions as MConv

def timer_func(func : Callable) -> Callable: 
    """!
    Decorator for timing functions in MockSZ.
    Useful for estimating performance.

    @param func Function to be timed.

    @returns Function to be timed, wrapped in timer.
    """

    def wrap_func(*args, **kwargs): 
        t1 = time() 
        result = func(*args, **kwargs) 
        t2 = time()
        if kwargs.get("timer") == True:
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

    def __init__(self, param    : Optional[float] = None, 
                       v_pec    : Optional[float] = None, 
                       phi_cl   : Optional[float] = 0, 
                       tau_e    : Optional[float] = 1, 
                       no_CMB   : Optional[bool]  = False) -> None:
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
            self.beta_cl = self.v_pec * 1e3 / const.c
            self.beta_cl_z = np.cos(np.radians(phi_cl))

        else:
            self.beta_cl = None
            self.beta_cl_z = None

        self.clib = MBind.loadMockSZlib()
    
    @timer_func
    def getSingleSignal_tkSZ(self, nu_arr   : Sequence[float], 
                                   timer    : Optional[bool]  = False, 
                                   acc      : Optional[float] = 1e-6) -> Sequence[float]:
        """!
        Generate a single pointing signal of the tSZ effect.

        @param nu_arr Array of frequencies for tSZ effect, in Hz.
        @param timer Time function execution. Used in decorator.
        @param acc Required relative accuracy of integration.
        
        @returns res 1D array containing tSZ effect.
        """

        res = np.zeros(nu_arr.size)
        
        if not self.no_CMB:
            res += self.getCMB(nu_arr)
        
        if self.param is not None:
            res += MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, acc, 
                                    func=self.clib.MockSZ_getSignal_tSZ)

        if self.v_pec is not None:
                
            res += MBind.getDistributionTwoParam(nu_arr, self.beta_cl * self.beta_cl_z, self.tau_e, acc, 
                                        func=self.clib.MockSZ_getSignal_kSZ)
            if self.param is not None:
                res += self.tau_e * MBind.getDistributionTwoParam(nu_arr, self.param, self.beta_cl, self.beta_cl_z, 
                                            func=self.clib.MockSZ_getSignal_corrections)
        
        return res
    
    @timer_func
    def getSingleSignal_ntkSZ(self, nu_arr  : Sequence[float], 
                                    timer   : Optional[bool]  = False, 
                                    acc     : Optional[float] = 1e-6) -> Sequence[float]:
        """!
        Generate a single pointing signal of the ntSZ effect, according to a powerlaw.

        @param nu_arr Numpy array of frequencies for ntSZ effect, in Hz.
        @param timer Time function execution. Used in decorator.
        @param acc Required relative accuracy of integration.
        
        @returns res 1D array containing ntSZ effect.
        """
        
        res = np.zeros(nu_arr.size)
        
        if not self.no_CMB:
            res += self.getCMB(nu_arr)
        
        if self.param is not None:
            res += MBind.getDistributionTwoParam(nu_arr, self.param, self.tau_e, acc, 
                                    func=self.clib.MockSZ_getSignal_ntSZ)

        if self.v_pec is not None:
            res += MBind.getDistributionTwoParam(nu_arr, self.beta_cl_z, self.tau_e, acc, 
                                        func=self.clib.MockSZ_getSignal_kSZ)

        return res

    def getCMB(self, nu_arr : Sequence[float]) -> Sequence[float]:
        return MBind.getCMB(nu_arr)

class IsoBetaModel(SinglePointing):
    """!
    Class representing an isothermal-beta model.
    Serves as an interface for MockSZ when simulating these types of clusters.

    Attributes:
        clib Library containing backend functions.

    @ingroup clustermodels
    """
    
    def __init__(self, param    : float, 
                       v_pec    : Optional[float] = None, 
                       phi_cl   : Optional[float] = 0, 
                       no_CMB   : Optional[bool]  = False) -> None:
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
    
    def getIsoBeta(self, Az     : Sequence[float], 
                         El     : Sequence[float], 
                         ibeta  : float, 
                         ne0    : float, 
                         thetac : float, 
                         Da     : float, 
                         grid   : Optional[bool] = False) -> Sequence[float]:
        """!
        Get an isothermal-beta optical depth screen. 

        @param Az Array containing the range of Azimuth co-ordinates, in arcseconds.
        @param El Array containing the range of Elevation co-ordinates, in arcseconds.
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
    
    def getIsoBetaCube(self, isobeta : Sequence[float], 
                             nu_arr  : Sequence[float], 
                             acc     : Optional[float] = 1e-6) -> Sequence[float]:
        """!
        Get an isothermal-beta model from an optical depth screen.

        @param isobeta An optical depth screen generated by self.getIsoBeta.
        @param nu_arr Array of frequencies for SZ effect, in Hz.
        @param acc Required relative accuracy of integration.
        
        @returns res 2D or 3D grid (depending on dimensions of isobeta) containing SZ signal attenuated by optical depth in isobeta.
        """

        res_SZ = self.getSingleSignal_tkSZ(nu_arr, acc=acc)

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
        
        if not self.no_CMB_cl:
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
    
    def __init__(self) -> None:
        self.clib = MBind.loadMockSZlib()

    def getSingleScattering(self, s_arr : Sequence[float], 
                                  beta  : float, 
                                  acc   : float = 1e-6) -> Sequence[float]:
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
    
    def getMaxwellJuttner(self, beta_arr : Sequence[float], 
                                Te       : float) -> Sequence[float]:
        """!
        Obtain a Maxwell-Juttner distribution.

        @param beta_arr Numpy array of dimensionless electron velocities.
        @param Te Electron temperature in keV.

        @returns res 1D array containing Maxwell-Juttner distribution.
        """
        
        res = MBind.getDistributionSingleParam(beta_arr, Te, acc=1e-6, 
                                               func=self.clib.MockSZ_getMaxwellJuttner)

        return res
    
    def getPowerlaw(self, beta_arr : Sequence[float], 
                          alpha    : float) -> Sequence[float]:
        """!
        Obtain a relativistic powerlaw distribution.

        @param beta_arr Numpy array of dimensionless electron velocities.
        @param alpha Slope of powerlaw.

        @returns res 1D array containing relativistic powerlaw distribution.
        """
        
        res = MBind.getDistributionSingleParam(beta_arr, alpha, acc=1e-6, 
                                               func=self.clib.MockSZ_getPowerlaw)

        return res
    
    def getMultiScatteringMJ(self, s_arr : Sequence[float], 
                                   Te    : float, 
                                   acc   : Optional[float] = 1e-6) -> Sequence[float]:
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
    
    def getMultiScatteringPL(self, s_arr : Sequence[float], 
                                   alpha : float, 
                                   acc   : Optional[float] = 1e-6) -> Sequence[float]:
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
