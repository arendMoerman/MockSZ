/*! \file Stats.h
    \brief Declarations of several statistical distributions used in MockSZ.
*/

#include "Constants.h"
#include "Conversions.h"
#include "Utils.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

#ifndef __Stats_h
#define __Stats_h

/**
 * Generate probablity for a single electron at speed beta to generate a logarithmic frequency shift (given by s_arr).
 *
 * This function assumes Thomson scattering in electron rest-frame.
 *
 * @param s_arr Array of doubles containing s-values over which to calculate probability.
 * @param n_s Number of s-values in array.
 * @param beta Double containing beta factor of electron.
 * @param output Array of doubles for storing results.
 * @param num_mu Number of direction cosines to evaluate for scattering. Defaults to 500, which should be enough.
 */
void getThomsonScatter(double *s_arr, int n_s, double beta, double *output, int num_mu=500);

/**
 * Generate probablity for a single electron at speed beta to generate a logarithmic frequency shift (given by s).
 *
 * This function assumes Thomson scattering in electron rest-frame.
 * Overloaded for single-argument s-values.
 *
 * @param s Double containing s-value at which to calculate probability.
 * @param beta Double containing beta factor of electron.
 * @param output Double for storing result.
 * @param num_mu Number of direction cosines to evaluate for scattering. Defaults to 500, which should be enough.
 */
void getThomsonScatter(double s, double beta, double &output, int num_mu=500);

/**
 * Generate a Maxwell-Juttner (relativistic thermal) distribution.
 *
 * @param beta_arr Array of beta values over which to calculate distribution.
 * @param n_beta Number of beta values in array.
 * @param Te Electron temperature in keV.
 * @param output Array for storing output values.
 */
void getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output);

/**
 * Generate a Maxwell-Juttner (relativistic thermal) distribution.
 *
 * Overloaded for single-argument beta-values
 *
 * @param beta Beta value at which to calculate distribution.
 * @param Te Electron temperature in keV.
 * @param output Double for storing output values.
 */
void getMaxwellJuttner(double beta, double Te, double &output);

/**
 * Generate a powerlaw (relativistic nonthermal) distribution.
 *
 * @param beta_arr Array of beta values over which to calculate distribution.
 * @param n_beta Number of beta values in array.
 * @param alpha Slope of powerlaw.
 * @param output Array for storing output values.
 */
void getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output);

/**
 * Generate a powerlaw (relativistic nonthermal) distribution.
 *
 * Overloaded for single-argument beta-values
 *
 * @param beta Beta value at which to calculate distribution.
 * @param alpha Slope of powerlaw.
 * @param A Normalisation factor of powerlaw. Can be calculated outside hot section.
 * @param output Double for storing output values.
 */
void getPowerlaw(double beta, double alpha, double A, double &output);

/**
 * Generate a multi-electron scattering kernel using a Maxwell-Juttner distribution.
 *
 * @param s_arr Array of s-values over which to calculate distribution.
 * @param n_s Number of s values in array.
 * @param n_beta Number of beta points to integrate over.
 * @param Te Electron temperature in keV.
 * @param output Array for storing output values.
 */
void getMultiScatteringMJ(double *s_arr, int n_s, int n_beta, double Te, double *output);

/**
 * Generate a multi-electron scattering kernel using a powerlaw distribution.
 *
 * @param s_arr Array of s-values over which to calculate distribution.
 * @param n_s Number of s values in array.
 * @param n_beta Number of beta points to integrate over.
 * @param alpha Slope of powerlaw.
 * @param output Array for storing output values.
 */
void getMultiScatteringPL(double *s_arr, int n_s, int n_beta, double alpha, double *output);

/**
 * Generate an isothermal-beta model, from an azimuth and elevation array.
 *
 * Returns an array of shape azimuth * elevation, containing the optical depth for each pointing.
 *
 * @param Az Array containing azimuth points in arcsec.
 * @param El Array containing elevation points in arcsec.
 * @param n_Az Number of azimuth points.
 * @param n_El Number of elevation points.
 * @param ibeta Beta parameter of isothermal model.
 * @param ne0 Central electron number density, in electrons / cm**3.
 * @param thetac Core radius of cluster in arcsec.
 * @param Da Angular diameter distance in Megaparsec.
 * @param output Array for storing outputs.
 * @param grid Whether or not to evaluate on Az-El grid, or along Az-El trace.
 *      Note: if grid=false, n_Az must equal n_El, and output must equal either one.
 *      If grid=true, n_Az does not need to equal n_El, output should have size n_Az*n_El.
 */
void getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid);

#endif
