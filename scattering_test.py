import numpy as np
import matplotlib.pyplot as pt

import MockSZ.Models as MModels

def get_zeros():
    num_arr = 3000  

    s = np.linspace(-1.2, 2.4, num=num_arr)                                                                                                                                     
    beta = np.array([0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.6])
    simObj = MModels.ScatteringKernels()
    ds = s[1] - s[0]

    fig, ax = pt.subplots(1,2, figsize=(10,5), gridspec_kw={"wspace":0.3})

    check = np.zeros(beta.size)
    zeros = []
    for i in range(beta.size):
        amplitudes = simObj.getSingleScattering(s, beta[i], num_mu=1000)
        
        ax[0].plot(s, amplitudes, label=r"$\beta$ = {}".format(beta[i]))

        check[i] = np.sum(amplitudes * ds, axis=0)

    ax[0].set_ylim(0, 12)
    ax[0].set_xlim(-1.2, 2*1.2)
        
    ax[0].set_ylabel(r"$P(s;\beta)$")
    ax[0].set_xlabel(r"$s$")
    ax[0].legend(frameon=False, prop={'size': 13},handlelength=1)
        
    ax[1].scatter(beta, np.log10(np.absolute(1 - check)), marker="X", color="black")
    ax[1].set_ylabel(r"$\log_{10}(1 - P_\mathrm{tot})$")
    ax[1].set_xlabel(r"$\beta$")

    pt.show()

if __name__ == "__main__":
    get_zeros()
