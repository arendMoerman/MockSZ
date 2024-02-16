/*! \file Constants.h
    \brief Often used physical constants.
    Also contains other constants used in MockSZ.
*/

#ifndef __Constants_h
#define __Constants_h

#define CH 6.62607015E-34           /* Planck constant */
#define KB 1.380649E-23             /* Boltzmann constant */
#define CL 2.99792458E8             /* Speed of light */
#define ME 9.1093837E-31            /* Electron mass, kg */
#define EV 1.602176634E-19          /* Electronvolt, in Joules */
#define ST 6.65245E-29              /* Thomson scattering cross section */

#define DBL_EPSILON 2.220446E-16    /* Double precision machine epsilon */
#define BETA1 0.99999999999         /* Limit on beta */
#define TCMB 2.726                  /* CMB temperature */

#define NW_INT 1000                 /* Number of subintervals for gsl integration routines. MIGHT WANT TO CUSTOMIZE*/

#define PI 3.14159265359            /* Pi. */

#define GQMODE GSL_INTEG_GAUSS31    /* Integrator type */
#define MEVALS 1000                  /* Maximum evaluations for romberg integrator.*/

#endif
