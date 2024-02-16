/*! \file Conversions.h
    \brief Functions for converting quantities.
*/

#include "Constants.h"
#include "math.h"

#ifndef __Conversions_h
#define __Conversions_h

/**
 * Convert electron temperature in keV to Kelvin.
 *
 * @param energy_keV Electron temperature in keV.
 *
 * @returns Temperature in Kelvin.
 */
double keV_Temp(double energy_keV);

/**
 * Convert distance in parsec to meter.
 *
 * @param l_pc Distance in parsec.
 *
 * @returns Distance in meters.
 */
double pc_m(double l_pc);

/**
 * Convert temperature in kelvin to dimensionless temperature.
 *
 * The dimension is removed by dividing by electron rest energy.
 *
 * @param Te Temperature in Kelvin.
 *
 * @returns Dimensionless temperature.
 */
double Te_theta(double Te);

/**
 * Convert velocity to dimensionless beta.
 *
 * @param velocity Velocity in m / s.
 *
 * @returns Dimensionless velocity.
 */
double v_beta(double velocity);

/**
 * Convert dimensionless beta velocity to Lorentz (gamma) factor.
 *
 * @param beta Dimensionless velocity
 *
 * @returns Gamma factor.
 */
double beta_gamma(double beta);

/**
 * Convert velocity to Lorentz (gamma) factor.
 *
 * @param velocity Velocity in m / s.
 *
 * @returns Gamma factor.
 */
double v_gamma(double velocity);

/**
 * Convert frequency nu in Hertz to dimensionless frequency.
 * Note that the conversion is done w.r.t. the electron temperature, not the CMB temperature, which is the form for the conversion in the Python conversions.
 *
 * @param nu Frequency in Hertz.
 * @param T Elektron temperature in K.
 *
 * @returns x Dimensionless frequency.
 */
double nu_x(double nu, double T);

inline double keV_Temp(double energy_keV) {
    return energy_keV * 1e3 / KB * EV;
}

inline double pc_m(double l_pc) {
    return l_pc * 3.0857e16;
}

inline double Te_theta(double Te) {
    return KB * Te / (ME * CL * CL);
}

inline double v_beta(double velocity) {
    return velocity / CL;
}

inline double beta_gamma(double beta) {
    return 1 / sqrt(1 - beta*beta);
}

inline double v_gamma(double velocity) {
    return beta_gamma(v_beta(velocity)); 
}

inline double nu_x(double nu, double T) {
    return CH * nu / KB / T; 
}
#endif
