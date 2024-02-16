/*! \file Signal.h
    \brief Declarations of the single-pointing kSZ signal and utility functions for tSZ and ntSZ signals in MockSZ.
*/
#include <cstdio>

#include "Constants.h"
#include "Conversions.h"
#include "Stats.h"

#ifndef __Signal_h
#define __Signal_h

/**
 * Input structure for integral for kSZ effect.
 */
struct kSZ_params {
    double nu;      /*< Frequency at which to integrate kSZ signal.*/
    double v_pec;   /*< Velocity of gas in km / s.*/
    double tau_e;   /*< Optical depth of gas.*/
};

/**
 * Single-pointing signal assuming kinematic SZ effect.
 *
 * @param mu direction cosine between electron and photon (integration variable).
 * @param args Struct containing nu, v_pec and tau_e.
 *
 * @returns Integrated kSZ signal at frequency nu.
 */
double calcSignal_kSZ(double mu, void *args);

/**
 * Obtain CMB intensity at frequency nu.
 *
 * @param nu Frequency at which to calculate CMB.
 *
 * @returns CMB intensity at frequency nu.
 */
double get_CMB(double nu);

/**
 * Calculate number of s-points necessary for integrating scattering kernel.
 *
 * This function should be passed to the romberg integrator, which returns (number of) function evaluations.
 *
 * @param func Pointer to function (either getMultiScatteringMJ, or getMultiScatteringPL) that calculates scattering kernel.
 * @param s Logarithmic frequency shift.
 * @param arg Extra argument for func, either Te or alpha depending on scattering kernel.
 *
 * @returns integrated value (integrated over s) of scattering kernel.
 */
double get_n_eval(double (*func)(double, void*), double s, double arg);

/**
 * Wrapper routine for romberg integration of CMB-scattering kernel convolution.
 *
 * @param s Logarithmic frequency shift (integration variable).
 * @param args Extra arguments, contains nu and tau_e.
 *
 * @param returns Integrated (n)tSZ signal at frequency nu.
 */
double conv_CMB_scatt(double s, double *args);

#endif
