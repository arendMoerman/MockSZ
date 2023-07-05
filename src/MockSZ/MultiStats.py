import numpy as np

import MockSZ.SingleStats as MSingle
import MockSZ.Utils as MUtils
import MockSZ.ElectronDistributions as EDist

import matplotlib.pyplot as pt

def getP1_RM(s, Te, num_beta=100, num_mu=100):
    """!
    Generate the one-scattering ensemble scattering kernel P1.
    This kernel corresponds to a relativistic Maxwellian.

    @param s Range of log frequency shifts over which to evaluate P1.
    @param Te Electron temperature of plasma in Kelvin.
    @param num_beta Number of beta factor points to evaluate.
    @param num_mu Number of direction cosines to evaluate for single electron scattering cross-section.

    @returns P1 The P1 scattering kernel.
    """

    beta_lim = (np.exp(np.absolute(s)) - 1) / (np.exp(np.absolute(s)) + 1)

    dbeta = (1 - beta_lim) / num_beta

    P1 = np.zeros(s.shape)
    for i in range(num_beta):
        be = beta_lim + i*dbeta
        Psb = MSingle.getPsbThomson(s, be, num_mu, grid=False)
        Pe = EDist.relativisticMaxwellian(be, Te)
        P1 += Pe * Psb * dbeta
        
    return P1

def getP1_PL(s, alpha, num_beta=100, num_mu=100):
    """!
    Generate the one-scattering ensemble scattering kernel P1.
    This kernel corresponds to a power law.

    @param s Range of log frequency shifts over which to evaluate P1.
    @param alpha Slope of power law.
    @param num_beta Number of beta factor points to evaluate.
    @param num_mu Number of direction cosines to evaluate for single electron scattering cross-section.

    @returns P1 The P1 scattering kernel.
    """
    beta_lim = (np.exp(np.absolute(s)) - 1) / (np.exp(np.absolute(s)) + 1)

    dbeta = (1 - beta_lim) / num_beta
    P1 = np.zeros(s.shape)
    for i in range(num_beta):
        be = beta_lim + i*dbeta
        Psb = MSingle.getPsbThomson(s, be, num_mu, grid=False)
        Pe = EDist.relativisticPowerlaw(be, np.min(beta_lim), alpha=alpha)
        P1 += Pe * Psb * dbeta * be * (1 - be**2)**(-3/2)
        
    return P1

