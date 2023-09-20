"""!
@file
File containing expressions for single pointing spectral distortions.
"""

import numpy as np

import MockSZ.MultiStats as MStats
import MockSZ.Utils as MUtils
import MockSZ.Constants as ct
import MockSZ.Backgrounds as MBack
import MockSZ.Conversions as MConv

import matplotlib.pyplot as pt

def getSpecIntensityRM(nu, Te, tau_e):
    lims_s = [-1.5, 2.5]
    func = MStats.getP1_RM
    Inu = getSpecIntensity(nu, Te, tau_e, func, lims_s)

    return Inu
    #return MIntegral.solveTripleIntegral(Te, nu, tau_e, dist="MJ")

def getSpecIntensityPL(nu, alpha, tau_e):
    lims_s = [-1.5, 10]
    func = MStats.getP1_PL
    Inu = getSpecIntensity(nu, alpha, tau_e, func, lims_s)

    return Inu

def getSpecIntensity(nu, param, tau_e, func, lims_s):
    """!
    Calculate the specific intensity of the CMBR distortion along a single line of sight.
    Can choose between thermalised electrons (Maxwellian) or non-thermal population (power law).
    Note that this only returns the thermal SZ effect: distortion due to random scattering.

    @param nu Range of frequencies over which to evaluate the intensity in Hertz.
    @param param Parameter for electron distribution. If relativistic Maxwellian, electron temperature of the cluster gas. If power law, spectral slope.
    @param tau_e Optical depth of cluster gas along line of sight. Note that this method assumes optically thin gases, i.e. tau_e << 1.
    @param func Electron distribution to use. Can choose between relativistic Maxwellian or power law.
    
    @returns Itot Comptonised CMBR specific intensity relative to CMBR.
    """

    s_range = np.linspace(lims_s[0], lims_s[1], num=1000)
    ds = s_range[1] - s_range[1]

    I0 = MBack.getSpecificIntensityCMB(nu)

    trans_I0 = -tau_e * I0 # CMB transmitted through cluster, attenuated by tau_e

    S, NU = MUtils.getXYGrid(s_range, nu)
    S += ds / 2
    P1 = func(s_range, param)
    P1_mat = np.vstack([P1] * S.shape[1]).T

    # Now, evaluate I0 on an mu*e^(-s) grid
    I0_mat = MBack.getSpecificIntensityCMB(NU*np.exp(-S))
    scatter_I0 = tau_e * np.sum(P1_mat * I0_mat * (s_range[1] - s_range[0]), axis=0)

    Itot = (scatter_I0 + trans_I0)

    return Itot

def getSpecIntensityKSZ(nu, beta_z, tau_e):
    """!
    Calculate the specific intensity of the CMBR distortion along a single line of sight.
    Note that this only returns the kinematic SZ effect: distortion due to gas bulk motion.

    @param nu Range of frequencies over which to evaluate the intensity in Hertz.
    @param beta_z Beta parameter of cluster receding velocity.
    @param tau_e Optical depth of cluster gas along line of sight. Note that this method assumes optically thin gases, i.e. tau_e << 1.
    
    @returns Itot Comptonised CMBR specific intensity relative to CMBR.
    """
    
    gamma_z = MConv.beta_gamma(beta_z)

    mu_range = np.linspace(-1, 1, num=1000)
    
    dmu = mu_range[1] - mu_range[0]

    I0 = MBack.getSpecificIntensityCMB(nu)

    MU, NU = MUtils.getXYGrid(mu_range, nu)

    MU += dmu / 2

    X = ct.h * NU / ct.k / ct.Tcmb
    X2 = X * gamma_z**2 * (1 + beta_z) * (1 - beta_z * MU)

    Itot = tau_e * np.sum(3/8 * (1 + MU**2) * ((np.exp(X) - 1) / (np.exp(X2) - 1) - 1) * dmu, axis=0) * I0

    return Itot
