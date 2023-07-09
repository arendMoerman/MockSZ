"""!
@file
File containing methods for unit conversions.
The methods in this file are part of the public interface.
"""

from MockSZ.Constants import Constants
import numpy as np

def eV_Temp(energy_eV):
    """!
    Convert an energy in electronvolt to temperature in Kelvin
    
    @param energy_eV Energy in electronvolt.
    
    @returns T temperature in Kelvin.
    """

    ct = Constants()
    T = energy_eV / ct.k * ct.eV

    return T

def SI_JySr(I_freq):
    """
    Convert a specific brightness in SI units and convert to Jansky over steradian.

    @param I_freq Specific intensity in SI units.

    @returns JySr The specific intensity in Jansky / steradian
    """

    JySr = I_freq / 1e-26

    return JySr

def SI_Temp(I_freq, freqHz):
    """!
    Take specific intensity in SI units.
    Convert to a brightness temperature in Kelvin.

    @param I_freq Specific intensity in SI units.
    @param freqHz Frequencies of I_freq in Hz.

    @returns Tb Brightness temperature.
    """
    
    ct = Constants()

    Tb = ct.h * freqHz / ct.k / (np.log(1 + 2 * ct.h * freqHz**3 / (I_freq * ct.c**2)))
    return Tb

def pc_m(l_pc):
    """!
    Convert a length in parsecs to meters.

    @param l_pc Length in parsecs.
    
    @returns l_m Length in meters.
    """

    conv = 3.0857e16
    l_m = l_pc * conv

    return l_m

