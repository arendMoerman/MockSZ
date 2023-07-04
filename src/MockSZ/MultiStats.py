import numpy as np

import MockSZ.SingleStats as MSingle
import MockSZ.Utils as MUtils
import MockSZ.ElectronDistributions as EDist

import matplotlib.pyplot as pt

def getP1_RM(s, Te, num_beta=100, num_mu=100):
    beta_lim = (np.exp(np.absolute(s)) - 1) / (np.exp(np.absolute(s)) + 1)

    dbeta = (1 - beta_lim) / num_beta

    integrand = np.zeros(s.shape)
    for i in range(num_beta):
        be = beta_lim + i*dbeta
        Psb = MSingle.getPsbThomson(s, be, num_mu, grid=False)
        pe = EDist.relativisticMaxwellian(be, Te)
        integrand += pe * Psb * dbeta
        
        #pt.plot(pe)
        #pt.show()
        
    return integrand
