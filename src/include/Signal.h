/*! \file Signal.h
    \brief Declarations of the single-pointing SZ signals in MockSZ.
*/

#include "Constants.h"
#include "Conversions.h"
#include "Utils.h"
#include "Stats.h"

#ifndef __Signal_h
#define __Signal_h

/**
 * Single-pointing signal assuming thermal SZ effect.
 *
 * @param nu Array with frequencies at which to calculate tSZ signal, in Hz.
 * @param n_nu Number of frequencies in nu.
 * @param Te Electron temperature in keV.
 * @param tau_e Optical depth along sightline.
 * @param output Array for storing output.
 * @param n_s Number of logarithmic frequency shifts to include.
 * @param n_beta Number of dimensionless electron velocities to include.
 * @param no_CMB Whether to add CMB to tSZ signal or not.
 */
void calcSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);

/**
 * Single-pointing signal assuming non-thermal SZ effect.
 *
 * @param nu Array with frequencies at which to calculate ntSZ signal, in Hz.
 * @param n_nu Number of frequencies in nu.
 * @param alpha Slope of powerlaw.
 * @param tau_e Optical depth along sightline.
 * @param output Array for storing output.
 * @param n_s Number of logarithmic frequency shifts to include.
 * @param n_beta Number of dimensionless electron velocities to include.
 * @param no_CMB Whether to add CMB to ntSZ signal or not.
 */
void calcSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, int n_s, int n_beta, bool no_CMB);

/**
 * Single-pointing signal assuming kinematic SZ effect.
 *
 * @param nu Array with frequencies at which to calculate tSZ signal, in Hz.
 * @param n_nu Number of frequencies in nu.
 * @param v_pec Peculiar velocity in km /s.
 * @param tau_e Optical depth along sightline.
 * @param output Array for storing output.
 * @param n_mu Number of scattering direction cosines to include.
 */
void calcSignal_kSZ(double *nu, int n_nu, double v_pec, double tau_e, double *output, int n_mu);

#endif
