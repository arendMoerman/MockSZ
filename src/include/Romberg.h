/*! \file Romberg.h
    \brief Declarations of functions used for romberg integrations.
*/

#include <stdio.h>
#include <math.h>

#ifndef __Romberg_h
#define __Romberg_h

/**
 * Routine for calculating the integral of a function using Romberg integration.
 *
 * This specific version calculates the amount of function evaluations and stores these in an array.
 * In this way, we can perform the final integral in MockSZ alot more efficiently.
 *
 * @param f Pointer to function, which is (currently) always get_n_eval.
 * @param g Pointer to a scattering kernel function, either getMultiScatteringMJ or getMultiScatteringPL.
 * @param a Lower limit on integral.
 * @param b Upper limit on integral.
 * @param arg Extra argument to pass to f, either Te or alpha.
 * @param write_arr Array for storing the function evaluations.
 * @param max_steps Maximum number of steps before giving up.
 * @param acc Relative accuracy of romberg integrator.
 *
 * @returns Number of rows in Richardson table.
 */
int romberg_write(double (*f)(double (*g)(double, void*), double, double), double (*g)(double, void*), double a, double b, double arg, double *write_arr, size_t max_steps, double acc);

/**
 * Routine for calculating the integral of a function using Romberg integration.
 *
 * This specific version performs the Richardson extrapolation on a given function, for a given set of function evaluations and returns the final SZ signals.
 * To be more efficient, the scattering kernel is passed as an array, from which the values are read.
 * In this way, we avoid unnecessary re-evaluations of the scattering kernel.
 *
 * @param f Pointer to function, which is (currently) always conv_CMB_scatt.
 * @param a Lower limit on integral.
 * @param b Upper limit on integral.
 * @param args Extra arguments to pass to f, nu and tau_e.
 * @param read_arr Array for reading the scattering kernel evaluations.
 * @param n_eval Number of rows used for calculating scattering kernel.
 *
 * @returns The integrated SZ signal at frequency nu.
 */
double romberg_read(double (*f)(double, double*), double a, double b, double *args, double *read_arr, size_t n_eval);

#endif
