import numpy as np
from MockSZ.Constants import Constants as ct

class CMB(object):
    def __init__(self, Tcmb = ct.Tcmb):
        self.T = Tcmb

    def getSpecificIntensity(self, freqs):
        prefac = 2 * ct.h * freqs**3 / ct.c**2
        distri = (np.exp(ct.h * freqs / (ct.k * self.T)) - 1)**(-1)

        return prefac * distri
