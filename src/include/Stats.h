/*! \file Stats.h
    \brief Declarations of several statistical distributions used in MockSZ.
*/

#include "Constants.h"
#include "Conversions.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cstdio>

#ifndef __Stats_h
#define __Stats_h

struct thom_params { double s; double beta; };
struct MS_params { double s; double param; };

/**
 * Calculate integration limits for integral over Thomson scattering cross section.
 *
 * @param s Logarithmic frequency shift.
 * @param beta Dimensionless electron velocity.
 * @param mu1 Double for storing lower integration limit.
 * @param mu2 Double for storing upper integration limit.
 */ 
void get_lims_mu(double s, double beta, double &mu1, double &mu2);

/**
 * Calculate probability for a scattering to give frequency shift s, given beta.
 *
 * The probability is calculated by integrating over all direction cosines.
 * This function assumes Thomson scattering in electron rest-frame.
 *
 * @param mu Direction cosine (integration variable).
 * @param args Pointer to struct containing s and beta.
 *
 * @returns Probability for a scattering to give a frequency shift s, for a given beta.
 */
double getThomsonScatter(double mu, void *args);

/**
 * Generate a Maxwell-Juttner (relativistic thermal) distribution.
 *
 * @param beta Beta value at which to calculate distribution.
 * @param Te Electron temperature in keV.
 *
 * @returns Probability for an electron to have velocity beta, given temperature.
 */
double getMaxwellJuttner(double beta, double Te);

/**
 * Generate a powerlaw (relativistic nonthermal) distribution.
 *
 * @param beta Beta value at which to calculate distribution.
 * @param alpha Slope of powerlaw.
 * @param A Normalisation factor of powerlaw. Can be calculated outside hot section.
 *
 * @returns Probability for an electron to have velocity beta, given alpha.
 */
double getPowerlaw(double beta, double alpha, double A);

/**
 * Generate a multi-electron scattering kernel using a Maxwell-Juttner distribution.
 *
 * @param beta Dimensionless electron velocity (integration variable).
 * @param args Pointer to struct containing s and Te.
 *
 * @returns Probability for frequency shift s, given an electron temperature.
 */
double getMultiScatteringMJ(double beta, void *args);

/**
 * Generate a multi-electron scattering kernel using a relativistic powerlaw distribution.
 *
 * @param beta Dimensionless electron velocity (integration variable).
 * @param args Pointer to struct containing s and alpha.
 *
 * @returns Probability for frequency shift s, given an alpha.
 */
double getMultiScatteringPL(double beta, void *args);

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

/**
 * Calculate normalisation constant for powerlaw distribution.
 *
 * @param gamma1 Lower limit on domain for powerlaw.
 * @param gamma2 Upper limit on domain for powerlaw. 
 * @param alpha Slope of powerlaw.
 * @param A Double for storing normalisation constant.
 */
void getNormPL(double &gamma1, double &gamma2, double &alpha, double &A);

#endif
