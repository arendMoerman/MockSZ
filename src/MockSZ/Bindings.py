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
    
    lib.MockSZ_getPowerlaw.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double)]
    
    lib.MockSZ_getMultiScatteringMJ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    
    lib.MockSZ_getMultiScatteringPL.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    
    lib.MockSZ_getSignal_tSZ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_bool]
    
    lib.MockSZ_getSignal_ntSZ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_bool]
    
    lib.MockSZ_getSignal_kSZ.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), ctypes.c_int] 
    
    lib.MockSZ_getIsoBeta.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_int,
                                         ctypes.c_double, ctypes.c_double,
                                         ctypes.c_double, ctypes.c_double,
                                         ctypes.POINTER(ctypes.c_double), ctypes.c_bool] 
   
    lib.MockSZ_getCMB.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double)]

    lib.MockSZ_getThomsonScatter.restype = None
    lib.MockSZ_getMaxwellJuttner.restype = None
    lib.MockSZ_getPowerlaw.restype = None
    lib.MockSZ_getMultiScatteringMJ.restype = None
    lib.MockSZ_getMultiScatteringPL.restype = None
    lib.MockSZ_getSignal_tSZ.restype = None
    lib.MockSZ_getSignal_ntSZ.restype = None
    lib.MockSZ_getSignal_kSZ.restype = None
    lib.MockSZ_getIsoBeta.restype = None
    lib.MockSZ_getCMB.restype = None

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

def getPowerlaw(beta_arr, alpha):
    """!
    Binding for calculating powerlaw distribution over range of beta, given alpha.

    @param beta_arr Array containing dimensionless velocities over which to calculate powerlaw.
    @param alpha Slope of powerlaw.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()

    cbeta_arr = (ctypes.c_double * beta_arr.size)(*(beta_arr.ravel().tolist()))
    cnum_beta = ctypes.c_int(beta_arr.size)
    calpha = ctypes.c_double(alpha)

    coutput = (ctypes.c_double * beta_arr.size)(*(np.zeros(beta_arr.size).tolist()))
    
    args = [cbeta_arr, cnum_beta, calpha, coutput]

    mgr.new_thread(target=lib.MockSZ_getPowerlaw, args=args)

    output = np.ctypeslib.as_array(coutput, shape=beta_arr.shape).astype(np.float64)

    return output

def getMultiScatteringMJ(s_arr, Te, n_beta=500):
    """!
    Binding for calculating relativistic thermal distribution over range of beta, given Te.

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

def getMultiScatteringPL(s_arr, alpha, n_beta=500):
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
    calpha = ctypes.c_double(alpha)
    cn_beta = ctypes.c_int(n_beta)

    coutput = (ctypes.c_double * s_arr.size)(*(np.zeros(s_arr.size).tolist()))
    
    args = [cs_arr, cnum_s, calpha, coutput, cn_beta]

    mgr.new_thread(target=lib.MockSZ_getMultiScatteringPL, args=args)

    output = np.ctypeslib.as_array(coutput, shape=s_arr.shape).astype(np.float64)

    return output

def getSinglePointing_tSZ(nu_arr, Te, tau_e, n_s=500, n_beta=500, no_CMB=False):
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
    ctau_e = ctypes.c_double(tau_e)
    cn_s = ctypes.c_int(n_s)
    cn_beta = ctypes.c_int(n_beta)
    cno_CMB = ctypes.c_bool(no_CMB)

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, cTe, ctau_e, coutput, cn_s, cn_beta, cno_CMB]

    mgr.new_thread(target=lib.MockSZ_getSignal_tSZ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output

def getSinglePointing_ntSZ(nu_arr, alpha, tau_e=0.01, n_s=500, n_beta=500, no_CMB=False):
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
    calpha = ctypes.c_double(alpha)
    ctau_e = ctypes.c_double(tau_e)
    cn_s = ctypes.c_int(n_s)
    cn_beta = ctypes.c_int(n_beta)
    cno_CMB = ctypes.c_bool(no_CMB)

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, calpha, ctau_e, coutput, cn_s, cn_beta, cno_CMB]

    mgr.new_thread(target=lib.MockSZ_getSignal_ntSZ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output

def getSinglePointing_kSZ(nu_arr, v_pec, tau_e=0.01, n_mu=500):
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
    ctau_e = ctypes.c_double(tau_e)
    cn_mu = ctypes.c_int(n_mu)

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, cv_pec, ctau_e, coutput, cn_mu]

    mgr.new_thread(target=lib.MockSZ_getSignal_kSZ, args=args)

    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output

def getIsoBeta(Az, El, ibeta, ne0, thetac, Da, grid):
    """!
    Binding for calculating relativistic powerlaw distribution over range of beta, given Te.

    @param s_arr Array containing logarithmic frequency shifts over which to calculate scattering probability.
    @param v_pec Cluster peculiar velocity, in km/s.
    @param num_mu Number of direction cosines over which to calculate scattering probability.

    @returns output Array containing scattering probabilities, for each s in s_arr, given beta.
    """

    lib = loadMockSZlib()
    mgr = TManager.Manager()
    
    cAz = (ctypes.c_double * Az.size)(*(Az.ravel().tolist()))
    cEl = (ctypes.c_double * El.size)(*(El.ravel().tolist()))
    cnum_Az = ctypes.c_int(Az.size)
    cnum_El = ctypes.c_int(El.size)
    cibeta = ctypes.c_double(ibeta)
    cne0 = ctypes.c_double(ne0)
    cthetac = ctypes.c_double(thetac)
    cDa = ctypes.c_double(Da)
    cgrid = ctypes.c_bool(grid)

    if grid:
        n_out = Az.size * El.size
        out_shape = (Az.size, El.size)
    else:
        n_out = Az.size 
        out_shape = (Az.size,)

    coutput = (ctypes.c_double * n_out)(*(np.zeros(n_out).tolist()))
    
    args = [cAz, cEl, cnum_Az, cnum_El, cibeta, cne0, cthetac, cDa, coutput, cgrid]

    mgr.new_thread(target=lib.MockSZ_getIsoBeta, args=args)
    output = np.ctypeslib.as_array(coutput, shape=n_out).astype(np.float64)
    output = output.reshape(out_shape)

    return output

def getCMB(nu_arr):
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

    coutput = (ctypes.c_double * nu_arr.size)(*(np.zeros(nu_arr.size).tolist()))
    
    args = [cnu_arr, cnum_nu, coutput]

    mgr.new_thread(target=lib.MockSZ_getCMB, args=args)
    output = np.ctypeslib.as_array(coutput, shape=nu_arr.shape).astype(np.float64)

    return output
