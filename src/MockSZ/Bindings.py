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
    
    lib.MockSZ_getThomsonScatter.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                             ctypes.c_int, ctypes.c_double, 
                                             ctypes.POINTER(ctypes.c_double), 
                                             ctypes.c_double]
    
    lib.MockSZ_getMaxwellJuttner.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                             ctypes.c_int, ctypes.c_double, 
                                             ctypes.POINTER(ctypes.c_double), 
                                             ctypes.c_double]
    
    lib.MockSZ_getPowerlaw.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                       ctypes.c_int, ctypes.c_double, 
                                       ctypes.POINTER(ctypes.c_double), 
                                       ctypes.c_double]
    
    lib.MockSZ_getMultiScatteringMJ.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                                ctypes.c_int, ctypes.c_double, 
                                                ctypes.POINTER(ctypes.c_double), 
                                                ctypes.c_double]
    
    lib.MockSZ_getMultiScatteringPL.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                                ctypes.c_int, ctypes.c_double, 
                                                ctypes.POINTER(ctypes.c_double), 
                                                ctypes.c_double]
    
    lib.MockSZ_getSignal_tSZ.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_bool, ctypes.c_double]
    
    lib.MockSZ_getSignal_tSZ_beta2.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_bool, ctypes.c_double]
    
    lib.MockSZ_getSignal_ntSZ.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                          ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                          ctypes.POINTER(ctypes.c_double), 
                                          ctypes.c_bool, ctypes.c_double]
    
    lib.MockSZ_getSignal_kSZ.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_bool, ctypes.c_double] 
    
    lib.MockSZ_getSignal_kSZ_betatheta.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_bool, ctypes.c_double]
    
    lib.MockSZ_getSignal_kSZ_betat2heta.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_int, ctypes.c_double, ctypes.c_double, 
                                         ctypes.POINTER(ctypes.c_double), 
                                         ctypes.c_bool, ctypes.c_double]
    
    lib.MockSZ_getIsoBeta.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                      ctypes.POINTER(ctypes.c_double), 
                                      ctypes.c_int, ctypes.c_int,
                                      ctypes.c_double, ctypes.c_double,
                                      ctypes.c_double, ctypes.c_double,
                                      ctypes.POINTER(ctypes.c_double), ctypes.c_bool] 
   
    lib.MockSZ_getCMB.argtypes = [ctypes.POINTER(ctypes.c_double), 
                                  ctypes.c_int, ctypes.POINTER(ctypes.c_double)]

    lib.MockSZ_getThomsonScatter.restype = None
    lib.MockSZ_getMaxwellJuttner.restype = None
    lib.MockSZ_getPowerlaw.restype = None
    lib.MockSZ_getMultiScatteringMJ.restype = None
    lib.MockSZ_getMultiScatteringPL.restype = None
    lib.MockSZ_getSignal_tSZ.restype = None
    lib.MockSZ_getSignal_tSZ_beta2.restype = None
    lib.MockSZ_getSignal_ntSZ.restype = None
    lib.MockSZ_getSignal_kSZ.restype = None
    lib.MockSZ_getSignal_kSZ_betatheta.restype = None
    lib.MockSZ_getSignal_kSZ_betat2heta.restype = None
    lib.MockSZ_getIsoBeta.restype = None
    lib.MockSZ_getCMB.restype = None

    return lib

def getDistributionSingleParam(x_arr, param, acc, func):
    """!
    Binding for evaluating various single-parameter distributions used in MockSZ.

    @param x_arr Array of independent variables.
    @param param Parameter defining distribution.
    @param acc Accuracy of evaluation of distribution. 
        Note: if the distribution does not involve integration, acc is ignored.
    @param func Function from library.

    @returns output Array containing distribution.
    """
    
    mgr = TManager.Manager()
    
    cx_arr = (ctypes.c_double * x_arr.size)(*(x_arr.ravel().tolist()))
    cnum_x = ctypes.c_int(x_arr.size)
    cparam = ctypes.c_double(param)

    coutput = (ctypes.c_double * x_arr.size)(*(np.zeros(x_arr.size).tolist()))
    
    cacc = ctypes.c_double(acc)
    
    args = [cx_arr, cnum_x, cparam, coutput, cacc]

    mgr.new_thread(target=func, args=args)

    output = np.ctypeslib.as_array(coutput, shape=x_arr.shape).astype(np.float64)

    return output

def getDistributionTwoParam(x_arr, param1, param2, no_CMB, acc, func):
    """!
    Binding for evaluating various two-parameter distributions used in MockSZ.
    These include the actual tSZ, ntSZ and kSZ signal.

    @param x_arr Array of independent variables.
    @param param1 Parameter 1 defining distribution.
    @param param2 Parameter 2 defining distribution.
    @param acc Accuracy of evaluation of distribution. 
    @param func Function from library.

    @returns output Array containing distribution.
    """
    
    mgr = TManager.Manager()
    
    cx_arr = (ctypes.c_double * x_arr.size)(*(x_arr.ravel().tolist()))
    cnum_x = ctypes.c_int(x_arr.size)
    cparam1 = ctypes.c_double(param1)
    cparam2 = ctypes.c_double(param2)
    cno_CMB = ctypes.c_bool(no_CMB)
    cacc = ctypes.c_double(acc)

    coutput = (ctypes.c_double * x_arr.size)(*(np.zeros(x_arr.size).tolist()))
    
    args = [cx_arr, cnum_x, cparam1, cparam2, coutput, cno_CMB, cacc]
    
    mgr.new_thread(target=func, args=args)

    output = np.ctypeslib.as_array(coutput, shape=x_arr.shape).astype(np.float64)

    return output

def getIsoBeta(Az, El, ibeta, ne0, thetac, Da, grid):
    """!
    Binding for calculating an isothermal-beta optical depth screen. 

    @param Az Numpy array containing the range of Azimuth co-ordinates, in arcseconds.
    @param El Numpy array containing the range of Elevation co-ordinates, in arcseconds.
    @param ibeta Beta parameter for isothermal model.
    @param ne0 Central electron number density, in number / cm**3.
    @param thetac Angular cluster core radius, in arcseconds.
    @param Da Angular diameter distance to cluster, in Megaparsec.
    @param grid Whether or not to evaluate the model on a 2D grid spanned by Az and El, or on a 1D trace.
        If grid=True, the screen will be of size Az.size * El.size.
        If grid=False (default), it is required that Az.size ==  El.size, and the screen will be of size Az.size = El.size.

    @returns output The optical depth screen.
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
    Binding for calculating CMB blackbody.

    @param nu_arr Numpy array of frequencies for CMB, in Hz.

    @returns output Array containing CMB.
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
