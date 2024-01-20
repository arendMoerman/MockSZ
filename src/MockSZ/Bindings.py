"""!
@file
Bindings for the ctypes interface for MockSZ. 
"""

import ctypes
import numpy as np
import os
import pathlib

import MockSZ.Threadmgr as TManager

def loadMockSZlib():
    """!
    Load the MockSZ shared library. Will detect the operating system and link the library accordingly.

    @returns lib The ctypes library containing the C/C++ functions.
    """

    path_cur = pathlib.Path(__file__).parent.resolve()
    try:
        lib = ctypes.CDLL(os.path.join(path_cur, "libmocksz.dll"))
    except:
        try:
            lib = ctypes.CDLL(os.path.join(path_cur, "libmocksz.so"))
        except:
            lib = ctypes.CDLL(os.path.join(path_cur, "libmocksz.dylib"))
    
    lib.MockSZ_getThomsonScatter.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    
    lib.MockSZ_getMaxwellJuttner.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double)]
    
    lib.MockSZ_getMultiScatteringMJ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    
    lib.MockSZ_getSignal_tSZ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_bool]
    
    lib.MockSZ_getSignal_kSZ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.c_int] 
    
    lib.MockSZ_getThomsonScatter.restype = None
    lib.MockSZ_getMaxwellJuttner.restype = None
    lib.MockSZ_getMultiScatteringMJ.restype = None
    lib.MockSZ_getSignal_tSZ.restype = None
    lib.MockSZ_getSignal_kSZ.restype = None


    return lib

def getThomsonScatter(s_arr, beta, num_mu):
    """!
    Binding for calculating single electron scattering kernel.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param beta Dimensionless electron velocity.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()

    cs_arr = (ctypes.c_double * s_arr.size)(*(s_arr.ravel().tolist()))
    cnum_s = ctypes.c_int(s_arr.size)
    cbeta = ctypes.c_double(beta)
    cnum_mu = ctypes.c_int(num_mu)

    coutput = (ctypes.c_double * s_arr.size)(*(np.zeros(s_arr.size).tolist()))
    
    args = [cs_arr, cnum_s, cbeta, coutput, cnum_mu]

    mgr.new_thread(target=lib.MockSZ_getThomsonScatter, args=args)

    output = np.ctypeslib.as_array(coutput, shape=s_arr.shape).astype(np.float64)

    return output

def getMaxwellJuttner(beta_arr, Te):
    """!
    Binding for calculating Maxwell-Juttner distribution over range of beta, given Te.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param beta Dimensionless electron velocity.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()

    cbeta_arr = (ctypes.c_double * beta_arr.size)(*(beta_arr.ravel().tolist()))
    cnum_beta = ctypes.c_int(beta_arr.size)
    cTe = ctypes.c_double(Te)

    coutput = (ctypes.c_double * beta_arr.size)(*(np.zeros(beta_arr.size).tolist()))
    
    args = [cbeta_arr, cnum_beta, cTe, coutput]

    mgr.new_thread(target=lib.MockSZ_getMaxwellJuttner, args=args)

    output = np.ctypeslib.as_array(coutput, shape=beta_arr.shape).astype(np.float64)

    return output

def getMultiScatteringMJ(s_arr, Te, n_beta=100):
    """!
    Binding for calculating relativistic powerlaw distribution over range of beta, given Te.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param beta Dimensionless electron velocity.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()

    cs_arr = (ctypes.c_double * s_arr.size)(*(s_arr.ravel().tolist()))
    cnum_s = ctypes.c_int(s_arr.size)
    cTe = ctypes.c_double(Te)
    cn_beta = ctypes.c_int(n_beta)

    coutput = (ctypes.c_double * s_arr.size)(*(np.zeros(s_arr.size).tolist()))
    
    args = [cs_arr, cnum_s, cTe, coutput, cn_beta]

    mgr.new_thread(target=lib.MockSZ_getMultiScatteringMJ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=s_arr.shape).astype(np.float64)

    return output

def getSinglePointing_tSZ(nu_arr, Te, n_s=500, n_beta=500, no_CMB=False):
    """!
    Binding for calculating relativistic powerlaw distribution over range of beta, given Te.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param beta Dimensionless electron velocity.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()
    
    cnu_arr = (ctypes.c_double * nu_arr.size)(*(nu_arr.ravel().tolist()))
    cnum_nu = ctypes.c_int(nu_arr.size)
    cTe = ctypes.c_double(Te)
    cn_s = ctypes.c_int(n_s)
    cn_beta = ctypes.c_int(n_beta)
    cno_CMB = ctypes.c_bool(no_CMB)

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, cTe, coutput, cn_s, cn_beta, cno_CMB]

    mgr.new_thread(target=lib.MockSZ_getSignal_tSZ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output

def getSinglePointing_kSZ(nu_arr, v_pec, n_mu=500):
    """!
    Binding for calculating relativistic powerlaw distribution over range of beta, given Te.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param v_pec Cluster peculiar velocity, in km/s.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()
    
    cnu_arr = (ctypes.c_double * nu_arr.size)(*(nu_arr.ravel().tolist()))
    cnum_nu = ctypes.c_int(nu_arr.size)
    cv_pec = ctypes.c_double(v_pec)
    cn_mu = ctypes.c_int(n_mu)

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, cv_pec, coutput, cn_mu]

    mgr.new_thread(target=lib.MockSZ_getSignal_kSZ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output
