"""!
@file
File containing expressions for single pointing spectral distortions.
"""

import numpy as np

import MockSZ.ElectronDistributions as EDist
from MockSZ.Backgrounds import CMB

import matplotlib.pyplot as pt

def getIntensityMuRM(mu, Te, tau_e):
    """!
    Calculate the specific intensity of the CMBR distortion along a single line of sight.

    @param mu Range of frequencies over which to evaluate the intensity.
    @param Te Electron temperature of the cluster gas.
    @param tau_e Optical depth of cluster gas along line of sight. Note that this method assumes optically thin gases, i.e. tau_e << 1.
    """

    s_range = np.linspace(-2, 2)
