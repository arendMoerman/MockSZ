"""!
@file
Cluster models for making SZ maps.
"""

from scipy.special import gamma
import numpy as np

import MockSZ.Utils as MUtils
import MockSZ.SinglePointing as MSingle
import MockSZ.Constants as ct

import MockSZ.Bindings as MBind

import matplotlib.pyplot as pt

class IsoBetaModel():
    """!
    Class representing an isothermal-beta model.
    Should be instantiated and serves as an interface for MockSZ when simulating these types of clusters.

    @ingroup clustermodels
    """

    def __init__(self, Te, ne0, rc, beta, Da, v_pec=0):
        """!
        Constructor: initisalise isothermal-beta model.

        @param Te Electron temperature in Kelvin.
        @param ne0 Central electron number density per cubic meter.
        @param rc Cluster core radius in meters.
        @param beta Structural parameter for model.
        @param Da Angular diameter distance to cluster in meters.
        @param v_pec Peculiar velocity of cluster, in m/s.
        """

        self.Te = Te
        self.ne0 = ne0
        self.beta = beta
        self.rc = rc
        self.Da = Da
        self.v_pec = v_pec / ct.c

        self.thetac = rc / Da

        self.tau0 = ne0 * ct.st * rc * np.sqrt(np.pi) * gamma(3/2*beta - 0.5) / gamma(3/2 * beta)

    def tSZMap(self, theta, freqHz):
        """!
        Generate a thermal SZ map over a 3D grid.
        The first two axes represent azimuth and elevation.
        The third axis is the frequency axis.

        @param theta Angular 2D grid on which to project cluster, in degrees.
        @param freqHz Frequency on which to evaluate SZ effect in Hertz.
        
        @returns tSZ A 3D numpy.ndarray of size (theta.shape[0], theta.shape[1], freqHz.size), containing the tSZ effect over the cluster, as function of frequency.
        """

        temp_tSZ = MSingle.getSpecIntensityRM(freqHz, self.Te, tau_e=1)
        tau = self.opticalDepths(theta)
        
        tSZ = tau[..., None] * temp_tSZ
        return tSZ
    
    def kSZMap(self, theta, freqHz):
        """!
        Generate a kinematic SZ map over a 3D grid.
        The first two axes represent azimuth and elevation.
        The third axis is the frequency axis.

        @param theta Angular 2D grid on which to project cluster, in degrees.
        @param freqHz Frequency on which to evaluate SZ effect in Hertz.
        
        @returns kSZ A 3D numpy.ndarray of size (theta.shape[0], theta.shape[1], freqHz.size), containing the kSZ effect over the cluster, as function of frequency.
        """

        temp_kSZ = MSingle.getSpecIntensityKSZ(freqHz, self.v_pec, tau_e=1)
        kSZ = self.opticalDepths(theta)[..., None] * temp_kSZ

        return kSZ

    def opticalDepths(self, theta):
        """!
        Calculate optical depth along sightline in cluster.

        @param theta On-sky angle, with respect to cluster center in degrees.
        """
        theta = np.radians(theta)

        tau_e = self.tau0 * (1 + (theta/self.thetac)**2)**(0.5 - 1.5 * self.beta)
        return tau_e

class SinglePointing(object):
    """! 
    Class for generating a single pointing SZ signal.

    """

    def getSingleSignal_tSZ(self, nu_arr, Te, n_s=500, n_beta=500, no_CMB=False):
        res = MBind.getSinglePointing_tSZ(nu_arr, Te, n_s, n_beta, no_CMB)

        return res
    
    def getSingleSignal_kSZ(self, nu_arr, v_pec, n_mu=500):
        res = MBind.getSinglePointing_kSZ(nu_arr, v_pec, n_mu)

        return res

class ScatteringKernels(object):
    """!
    Class for investigating the scattering kernels.
    It is not to be used per se for simulations, but it provides a nice overview of available distributions.
    So, it is mostly informative.
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
    
    def getRelatPowerlaw(self, beta_arr, alpha):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getRelatPowerlaw(beta_arr, alpha)

        return res
    
    def getMultiScatteringMJ(self, s_arr, Te, n_beta=100):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getMultiScatteringMJ(s_arr, Te, n_beta)

        return res
    
    def getMultiScatteringRP(self, s_arr, alpha, n_beta=100):
        """!
        Obtain single-electron scattering kernel, for a range of s and a single beta.
        """
        
        res = MBind.getMultiScatteringRP(s_arr, alpha, n_beta)

        return res
