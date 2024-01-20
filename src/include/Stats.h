#include "Constants.h"
#include "Conversions.h"
#include "Utils.h"

#include <gsl/gsl_sf_bessel.h>

#ifndef __Stats_h
#define __Stats_h

/**
 * Generate probablity for a single electron at speed beta to generate a logarithmic frequency shift (given by s_arr).
 *
 * This function assumes an elastic scattering in electron rest-frame, and thus a Thomson cross section.
 *
 * @param s_arr Array of doubles containing s-values over which to calculate probability.
 * @param n_s Number of s-values in array.
 * @param beta Single double containing beta factor of electron.
 * @param output Array of doubles for storing results.
 * @param num_mu Number of direction cosines to evaluate for scattering. Defaults to 1000, which should be enough.
 */
void getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu=100);

/**
 * Generate probablity for a single electron at speed beta to generate a logarithmic frequency shift (given by s).
 *
 * This function assumes an elastic scattering in electron rest-frame, and thus a Thomson cross section.
 * Overloaded for single-argument s-values.
 *
 * @param s_arr Double containing s-value at which to calculate probability.
 * @param beta Single double containing beta factor of electron.
 * @param output Doubles for storing result.
 * @param num_mu Number of direction cosines to evaluate for scattering. Defaults to 1000, which should be enough.
 */
void getThomsonScatter(double s, double beta, double &output, int num_mu=100);

/**
 * Generate a Maxwell-Juttner (relativistic thermal) distribution.
 *
 * @param beta_arr Array of beta values over which to calculate distribution.
 * @param n_beta Number of beta values in array.
 * @param Te Electron temperature in KeV.
 * @param output Array for storing output values.
 */
void getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output);

/**
 * Generate a Maxwell-Juttner (relativistic thermal) distribution.
 *
 * @param beta Beta value at which to calculate distribution.
 * @param Te Electron temperature in KeV.
 * @param output Double for storing output values.
 */
void getMaxwellJuttner(double beta, double Te, double &output);

/**
 * Generate a multi-electron scattering kernel using a Maxwell-Juttner distribution.
 *
 * @param s_arr Array of s-values over which to calculate distribution.
 * @param n_s Number of s values in array.
 * @param n_beta Number of beta points to integrate over.
 * @param output Array for storing output values.
 */
// NOTE: Maybe make this one parallel? Parallelize over s_arr!
void getMultiScatteringMJ(double *s_arr, int n_s, int n_beta, double Te, double *output);

#endif
