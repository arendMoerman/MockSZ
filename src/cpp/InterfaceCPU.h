/*! \file InterfaceCPU.h
    \brief Interface for ctypes bindings.
*/

#include "Stats.h"
#include "Signal.h"
#include "Romberg.h"

#ifdef _WIN32
#   define MOCKSZ_DLL __declspec(dllexport)
#else
#   define MOCKSZ_DLL
#endif

#ifndef __InterfaceCPU_h
#define __InterfaceCPU_h

extern "C"
{
    /**
     * Generate probablity for a single electron at speed beta to generate a logarithmic frequency shift (given by s_arr).
     *
     * This function assumes Thomson scattering in electron rest-frame.  
     * 
     * @param s_arr Array of doubles containing s-values over which to calculate probability.
     * @param n_s Number of s-values in array.
     * @param beta Double containing beta factor of electron.
     * @param output Array of doubles for storing results.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getThomsonScatter(double *s_arr, int n_s, double beta, double *output, double acc);

    /**
     * Generate a Maxwell-Juttner (relativistic thermal) distribution.
     *
     * @param beta_arr Array of beta values over which to calculate distribution.
     * @param n_beta Number of beta values in array.
     * @param Te Electron temperature in keV.
     * @param output Array for storing output values.
     * @param acc Accuracy of integrator. 
     *      Note that this is only passed for homogenity in the bindings and ignored in the actual function
     */
    MOCKSZ_DLL void MockSZ_getMaxwellJuttner(double *beta_arr, int n_beta, double Te, double *output, double acc); 
    
    /**
     * Generate a powerlaw (relativistic nonthermal) distribution.
     *
     * @param beta_arr Array of beta values over which to calculate distribution.
     * @param n_beta Number of beta values in array.
     * @param alpha Slope of powerlaw.
     * @param output Array for storing output values.
     *      Note that this is only passed for homogenity in the bindings and ignored in the actual function
     */
    MOCKSZ_DLL void MockSZ_getPowerlaw(double *beta_arr, int n_beta, double alpha, double *output, double acc);
    
    /**
     * Generate a multi-electron scattering kernel using a Maxwell-Juttner distribution.
     *
     * @param s_arr Array of s-values over which to calculate distribution.
     * @param n_s Number of s values in array.
     * @param Te Electron temperature in keV.
     * @param output Array for storing output values.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getMultiScatteringMJ(double *s_arr, int n_s, double Te, double *output, double acc);
    
    /**
     * Generate a multi-electron scattering kernel using a powerlaw distribution.
     *
     * @param s_arr Array of s-values over which to calculate distribution.
     * @param n_s Number of s values in array.
     * @param alpha Slope of powerlaw.
     * @param output Array for storing output values.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getMultiScatteringPL(double *s_arr, int n_s, double alpha, double *output, double acc);
    
    /**
     * Single-pointing signal assuming thermal SZ effect.
     *
     * @param nu Array with frequencies at which to calculate tSZ signal, in Hz.
     * @param n_nu Number of frequencies in nu.
     * @param Te Electron temperature in keV.
     * @param tau_e Optical depth along sightline.
     * @param output Array for storing output.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getSignal_tSZ(double *nu, int n_nu, double Te, double tau_e, double *output, double acc);
    
    /**
     * Single-pointing signal assuming non-thermal SZ effect.
     *
     * @param nu Array with frequencies at which to calculate ntSZ signal, in Hz.
     * @param n_nu Number of frequencies in nu.
     * @param alpha Slope of powerlaw.
     * @param tau_e Optical depth along sightline.
     * @param output Array for storing output.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getSignal_ntSZ(double *nu, int n_nu, double alpha, double tau_e, double *output, double acc);

    /**
     * Single-pointing signal assuming kinematic SZ effect.
     *
     * @param nu Array with frequencies at which to calculate tSZ signal, in Hz.
     * @param n_nu Number of frequencies in nu.
     * @param beta_pec_z Dimensionless peculiar velocity of cluster along sightline.
     * @param tau_e Optical depth along sightline.
     * @param output Array for storing output.
     * @param acc Accuracy of integrator.
     */
    MOCKSZ_DLL void MockSZ_getSignal_kSZ(double *nu, int n_nu, double beta_pec_z, double tau_e, double *output, double acc);
    
    /**
     * Correction (cross) terms up to second order in bulk velocity and electron temperature.
     *
     * @param nu Array with frequencies at which to calculate correction term, in Hz.
     * @param n_nu Number of frequencies in nu.
     * @param Te Electron temperature in keV.
     * @param beta_pec Dimensionless peculiar velocity of cluster.
     * @param output Array for storing output.
     * @param cosu Direction cosine between peculiar velocity and sightline.
     */
    MOCKSZ_DLL void MockSZ_getSignal_corrections(double *nu, int n_nu, double Te, double beta_pec, double *output, double cosu);

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
    MOCKSZ_DLL void MockSZ_getIsoBeta(double *Az, double *El, int n_Az, int n_El, double ibeta, double ne0, double thetac, double Da, double *output, bool grid);

    /**
     * Obtain value of CMB intensity at a range of frequencies.
     *
     * @param nu Array of frequencies in Hz.
     * @param n_nu Number of frequencies in nu.
     * @param output Array for storing outputs.
     *
     * @returns CMB intensities.
     */
    MOCKSZ_DLL void MockSZ_getCMB(double *nu, int n_nu, double *output);
}

#endif
