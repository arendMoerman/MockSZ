"""!
@file
Methods for unit conversions.
"""

from typing import Union, Sequence

import numpy as np
import scipy.constants as const

Numbers = Union[Union[float, Sequence], Union[int, Sequence[int]]]


TCMB = 2.726 # K

def keV_theta(Te : Numbers) -> Numbers:
    """!
    Get dimensionless electron temperature.

    @param Te Electron temperature in keV.
    
    @returns theta Dimensionless electron temperature.
    """

    theta = const.k * keV_Temp(Te) / (const.m_e * const.c**2)
    return theta

def keV_Temp(energy_keV : Numbers) -> Numbers:
    """!
    Convert an energy in kilo electronvolt to temperature in Kelvin

    @param energy_keV Energy in kilo electronvolt.

    @returns T temperature in Kelvin.
    """

    T = energy_keV / const.k * const.eV * 1e3

    return T

def SI_JySr(I_nu : Numbers) -> Numbers:
    """
    Convert a specific brightness in SI units and convert to Jansky over steradian.

    @param I_nu Specific intensity in SI units.

    @returns JySr The specific intensity in Jansky / steradian
    """

    JySr = I_nu / 1e-26

    return JySr

def SI_Temp(I_nu   : Numbers, 
            nu_arr : Numbers) -> Numbers:
    """!
    Take specific intensity in SI units.
    Convert to a brightness temperature in Kelvin, assuming Rayleigh-Jeans tail.

    @param I_nu Specific intensity in SI units.
    @param nu_arr Numpy array with frequencies of I_nu in Hz.

    @returns Tb Brightness temperature.
    """

    Tb = I_nu * const.c**2 / (2 * const.k * nu_arr**2)
    return Tb

def freq_x(nu_arr : Numbers) -> Numbers:
    """!
    Convert frequency in Hertz to dimensionless frequency using CMB temperature.

    @param nu_arr Numpy array with frequencies of I_nu in Hz.
    
    @returns x The dimensionless frequency.
    """

    x = const.h * nu_arr / (const.k * TCMB)

    return x

def x_freq(x : Numbers) -> Numbers:
    """!
    Convert dimensionless frequency to frequency in Hertz using CMB temperature.

    @param x Dimensionless frequency.
    
    @param nu_arr Numpy array with frequencies of I_nu in Hz.
    """

    nu_arr = x / const.h * const.k * TCMB

    return nu_arr

