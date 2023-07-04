import numpy as np

from MockSZ.Constants import Constants as ct

def getBetaFromVelocity(velocity):
    """
    Obtain the beta factor from a velocity.

    @param velocity The electron velocity in m / s. Float or numpy array.

    @returns beta The beta factor. Float or numpy array.
    """
    
    beta = velocity / ct.c

    return beta

def getGammaFromVelocity(velocity)
    """
    Obtain the gamma factor from a velocity.

    @param velocity The electron velocity in m / s. Float or numpy array.

    @returns gamma The gamma factor. Float or numpy array.
    """

    beta = getBetaFromVelocity(velocity)
    gamma = getGammaFromBeta(beta) 

    return gamma

def getGammaFromBeta(beta)
    """
    Obtain the gamma factor from a beta factor.

    @param beta The electron beta factor. Float or numpy array.

    @returns gamma The gamma factor. Float or numpy array.
    """

    gamma = 1 / np.sqrt(1 - beta**2)

    return gamma



