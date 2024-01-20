#include "Constants.h"
#include "math.h"

#ifndef __Conversions_h
#define __Conversions_h

double KeV_Temp(double energy_KeV);
double pc_m(double l_pc);
double Te_theta(double Te);
double v_beta(double velocity);
double beta_gamma(double beta);
double v_gamma(double velocity);


inline double KeV_Temp(double energy_KeV) {
    return energy_KeV * 1e3 / KB * EV;
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
#endif
