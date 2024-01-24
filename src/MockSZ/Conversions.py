"""!
@file
Methods for unit conversions.
"""

import MockSZ.Constants as ct
import numpy as np

def keV_theta(Te):
    """!
    Get dimensionless electron temperature.

    @param Te Electron temperature in keV.
    
    @returns theta Dimensionless electron temperature.
    """

    theta = ct.k * keV_Temp(Te) / (ct.me * ct.c**2)
    return theta

def keV_Temp(energy_keV):
    """!
    Convert an energy in kilo electronvolt to temperature in Kelvin

    @param energy_keV Energy in kilo electronvolt.

    @returns T temperature in Kelvin.
    """

    T = energy_keV / ct.k * ct.eV * 1e3

    return T

def SI_JySr(I_nu):
    """
    Convert a specific brightness in SI units and convert to Jansky over steradian.

    @param I_nu Specific intensity in SI units.

    @returns JySr The specific intensity in Jansky / steradian
    """

    JySr = I_nu / 1e-26

    return JySr

def SI_Temp(I_nu, nu_arr):
    """!
    Take specific intensity in SI units.
    Convert to a brightness temperature in Kelvin, assuming Rayleigh-Jeans tail.

    @param I_nu Specific intensity in SI units.
    @param nu_arr Numpy array with frequencies of I_nu in Hz.

    @returns Tb Brightness temperature.
    """

    Tb = I_nu * ct.c**2 / (2 * ct.k * nu_arr**2)
    return Tb

def freq_x(nu_arr):
    """!
    Convert frequency in Hertz to dimensionless frequency using CMB temperature.

    @param nu_arr Numpy array with frequencies of I_nu in Hz.
    
    @returns x The dimensionless frequency.
    """

    x = ct.h * nu_arr / (ct.k * ct.Tcmb)

    return x

def x_freq(x):
    """!
    Convert dimensionless frequency to frequency in Hertz using CMB temperature.

    @param x Dimensionless frequency.
    
    @param nu_arr Numpy array with frequencies of I_nu in Hz.
    """

    nu_arr = x / ct.h * ct.k * ct.Tcmb

    return nu_arr

