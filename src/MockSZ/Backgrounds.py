import numpy as np
import MockSZ.Constants as ct

def getSpecificIntensityCMB(freqs):
    prefac = 2 * ct.h * freqs**3 / ct.c**2
    distri = (np.exp(ct.h * freqs / (ct.k * ct.Tcmb)) - 1)**(-1)

    return prefac * distri
