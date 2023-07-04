import numpy as np
import scipy.special as sp

from MockSZ.Constants import Constants as ct
import MockSZ.Utils as MUtils

def getDimTemp(Te):
    """
    Get dimensionless electron temperature.

    @param Te Electron temperature in Kelvin.
    
    @returns theta Dimensionless electron temperature.
    """

    theta = ct.k * Te / (ct.me * ct.c**2)
    return theta

def relativisticMaxwellian(beta, Te):
    """
    Generate an electron population from a relativistic Maxwellian.

    @param beta Range of beta factors over which to define the distribution.
    @param Te Mean electron temperature in Kelvin.
    
    @returns pe Electron probability distribution.
    """

    theta = getDimTemp(Te)
    gamma = MUtils.getGammaFromBeta(beta)

    nomi = gamma**5 * beta**2 * np.exp(-gamma / theta)
    deno = theta * sp.kn(2, 1/theta)

    return nomi / deno

def relativisticPowerlaw(beta, alpha=None):
    """
    Generate an electron population from a relativistic power law.

    @param beta Range of beta factors over which to define the distribution.
    @param alpha Slope of power law.
    
    @returns pe Electron probability distribution.
    """

    gamma = MUtils.getGammaFromBeta(beta)

    gamma1 = np.min(gamma)
    gamma2 = np.max(gamma)

    if alpha is None:
        A = np.log10(gamma2) - np.log10(gamma1)
    else:
        A = (1 - alpha) * (gamma2**(1 - alpha) - gamma1**(1 - alpha))

    return A * gamma**(-alpha)
