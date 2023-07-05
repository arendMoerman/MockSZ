"""!
@file
File containing expressions for single pointing spectral distortions.
"""

import numpy as np

import MockSZ.MultiStats as MStats
import MockSZ.Utils as MUtils
from MockSZ.Backgrounds import CMB

import matplotlib.pyplot as pt

def getSpecIntensityRM(mu, Te, tau_e):
    func = MStats.getP1_RM
    Imu = getSpecIntensity(mu, Te, tau_e, func)

    return Imu

def getSpecIntensityPL(mu, alpha, tau_e):
    func = MStats.getP1_PL
    Imu = getSpecIntensity(mu, alpha, tau_e, func)

    return Imu

def getSpecIntensity(mu, param, tau_e, func):
    """!
    Calculate the specific intensity of the CMBR distortion along a single line of sight.

    @param mu Range of frequencies over which to evaluate the intensity in Hertz.
    @param param Parameter for electron distribution. If relativistic Maxwellian, electron temperature of the cluster gas. If power law, spectral slope.
    @param tau_e Optical depth of cluster gas along line of sight. Note that this method assumes optically thin gases, i.e. tau_e << 1.
    @param func Electron distribution to use. Can choose between relativistic Maxwellian or power law.
    """

    s_range = np.linspace(-20, 20, num=1000)

    cmb = CMB()

    I0 = cmb.getSpecificIntensity(mu)

    trans_I0 = (1 - tau_e) * I0 # CMB transmitted through cluster, attenuated by tau_e

    S, MU = MUtils.getXYGrid(s_range, mu)

    P1 = func(s_range, param)

    P1_mat = np.vstack([P1] * S.shape[1]).T
  
    pt.plot(s_range, P1)
    pt.show()

    # Now, evaluate I0 on an mu*e^(-s) grid
    I0_mat = cmb.getSpecificIntensity(MU*np.exp(-S))
    scatter_I0 = tau_e * np.sum(P1_mat * I0_mat * (s_range[1] - s_range[0]), axis=0)

    Itot = (scatter_I0 + trans_I0 - I0)

    return Itot
